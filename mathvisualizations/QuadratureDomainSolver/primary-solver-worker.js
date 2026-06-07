// @ts-check
// =============================================================================
// primary-solver-worker.js -- Single warm Web Worker for the Inverse-tab
// primary solve path.
//
// Why: the full QD.solveInverseQD pipeline (multistart + continuation +
// deflation + diverse-seeds) runs for 50-500 ms on hard h's. On the main
// thread that freezes input + repaint, forcing a 250 ms debounce on
// solveAndRender. A single warm worker doing the same work keeps the UI
// responsive.
//
// API:
//   QD.PrimarySolverWorker.ensureReady() -> Promise<void>
//     Lazy-create the worker. Idempotent.
//
//   QD.PrimarySolverWorker.solve(hData, opts, { signal? }) -> Promise<result>
//     Post a solve request to the worker. Resolves with the same
//     { success, primary, alternates, attempts, error } shape that
//     QD.solveInverseQD returns on the main thread. `signal.aborted` (or
//     calling cancel()) terminates the in-flight job and rejects with
//     { aborted: true }.
//
//   QD.PrimarySolverWorker.cancel() -> void
//     Abort any in-flight solve. The worker is terminated and recreated on
//     the next solve() — the easiest way to interrupt deeply-nested Newton
//     iteration since the worker thread can't otherwise be preempted.
//
//   QD.PrimarySolverWorker.isBusy() -> bool
//
//   QD.PrimarySolverWorker.searchAlternates(hData, norm, known, opts)
//       -> Promise<solution[]>
//     Background alternate-solution search (QD.searchAlternates) on a SEPARATE
//     aux worker, so it never blocks or preempts an interactive primary solve
//     (A3). Posting a new request supersedes any prior in-flight aux job.
//
//   QD.PrimarySolverWorker.cancelAux() -> void
//     Abort any in-flight alternate search (terminates the aux worker).
//
// Fallback: if Worker / fetch fails (e.g. file:// origin), solve() runs
// on the main thread synchronously inside a Promise.resolve().then(...) so
// callers don't need to special-case anything. Logged once at console.warn.
//
// Reuses the same solver-source list as param-slice-pool.js, with its own
// worker-side handler.
// =============================================================================

(function (global) {
  'use strict';

  // Worker bundle source list is the single source of truth in
  // asset-manifest.js (P3.3 / A7), which index.html (and the node-test ctx)
  // load before this file. No inline fallback — a stale copy silently drops
  // newer solver families from the worker bundle (it omitted all 5 PQD files
  // for months). Fail loud instead.
  const SOLVER_SRC_FILES = global.QD_ASSET_MANIFEST && global.QD_ASSET_MANIFEST.WORKER_BUNDLE_FILES;
  if (!SOLVER_SRC_FILES) {
    throw new Error(
      'primary-solver-worker.js: QD_ASSET_MANIFEST.WORKER_BUNDLE_FILES unavailable — ' +
      'asset-manifest.js must load before this file.');
  }

  // Worker-side message handler. Self-contained string appended after the
  // bundled solver source. Listens for 'solve' messages; calls QD.solveInverseQD
  // with the supplied { hData, opts } payload; posts the result back keyed
  // by the same jobId.
  const WORKER_HANDLER = `
;(function () {
  'use strict';
  self.onmessage = function (e) {
    const msg = e.data;
    if (!msg) return;
    // 'solve' — full primary solve (QD.solveInverseQD). Used by the warm
    // primary worker. 'altSearch' — background alternate-solution search
    // (QD.searchAlternates). Used by the dedicated aux worker so it never
    // queues behind / preempts a primary solve. Both echo the same jobId.
    if (msg.kind === 'solve') {
      const { jobId, hData, opts } = msg;
      let result;
      try {
        result = self.QD.solveInverseQD(hData, opts || {});
      } catch (err) {
        self.postMessage({ kind: 'solve', jobId, error: String(err && err.stack || err) });
        return;
      }
      self.postMessage({ kind: 'solve', jobId, result });
    } else if (msg.kind === 'altSearch') {
      const { jobId, hData, norm, known, opts } = msg;
      let result;
      try {
        result = self.QD.searchAlternates(hData, norm, known || [], opts || {});
      } catch (err) {
        self.postMessage({ kind: 'altSearch', jobId, error: String(err && err.stack || err) });
        return;
      }
      self.postMessage({ kind: 'altSearch', jobId, result });
    } else if (msg.kind === 'liveSolve') {
      // 'liveSolve' — one warm Newton + reduced-sample checks (QD.liveSolveStep),
      // the per-drag-frame path. Runs on a dedicated live worker so a slow frame
      // never blocks the UI; superseded jobs are dropped on the main side.
      const { jobId, hData, initPhi, opts } = msg;
      let result;
      try {
        result = self.QD.liveSolveStep(hData, initPhi, opts || {});
      } catch (err) {
        self.postMessage({ kind: 'liveSolve', jobId, error: String(err && err.stack || err) });
        return;
      }
      self.postMessage({ kind: 'liveSolve', jobId, result });
    }
  };
})();
`;

  // Bundle is cached at module scope (not pool-level) — each call to
  // ensureReady() reuses the same bundle URL across worker restarts.
  let _bundlePromise = null;
  function getBundleURL() {
    if (_bundlePromise) return _bundlePromise;
    _bundlePromise = (async () => {
      const parts = [];
      // Worker scope: alias `window` -> `self` so solver files' `window.QD`
      // namespace registrations attach to the worker global. Same trick as
      // param-slice-pool.js.
      parts.push('var window = self;\n');
      // Cache-bust each source fetch with the release version so a Worker can
      // never run stale solver source from the browser HTTP cache after a
      // deploy (see asset-manifest.js CACHE_VERSION).
      const _ver = (global.QD_ASSET_MANIFEST && global.QD_ASSET_MANIFEST.CACHE_VERSION) || '0';
      for (const f of SOLVER_SRC_FILES) {
        const resp = await fetch(f + '?v=' + encodeURIComponent(_ver));
        if (!resp.ok) throw new Error('primary-solver-worker: failed to fetch ' + f + ' (' + resp.status + ')');
        parts.push('/*===== ' + f + ' =====*/\n');
        parts.push(await resp.text());
        parts.push('\n');
      }
      parts.push('/*===== worker handler =====*/\n');
      parts.push(WORKER_HANDLER);
      const blob = new Blob(parts, { type: 'application/javascript' });
      return URL.createObjectURL(blob);
    })();
    return _bundlePromise;
  }

  /**
   * @typedef {Object} InflightJob
   * @property {number} jobId
   * @property {(value: any) => void} resolve
   * @property {(reason: any) => void} reject
   * @property {(e: MessageEvent) => void} onMessage
   */

  /** @type {Worker|null} */
  let _worker = null;
  /** @type {Promise<void>|null} */
  let _readyPromise = null;
  let _nextJobId = 1;
  /** @type {InflightJob|null} */
  let _inflight = null;            // current outstanding solve, if any
  let _mainThreadFallback = false; // set true after the worker fails to load once.

  function _disposeWorker() {
    if (_worker) {
      try { _worker.terminate(); } catch (_) { /* ignore */ }
      _worker = null;
    }
    _readyPromise = null;
    if (_inflight) {
      _inflight.reject({ aborted: true });
      _inflight = null;
    }
  }

  async function ensureReady() {
    if (_mainThreadFallback) return;
    if (_worker) return;
    if (_readyPromise) { await _readyPromise; return; }
    _readyPromise = (async () => {
      const url = await getBundleURL();
      const w = new Worker(url);
      w.addEventListener('error', (ev) => {
        const detail = (ev.message || ev) + ' @ ' + (ev.filename || 'bundle') + ':' + (ev.lineno || '?');
        console.error('[primary-solver worker] error: ' + detail);
        // A worker-level error (bundle load/syntax error, crash, OOM) posts NO
        // {error} message, so without this the in-flight solve() promise would
        // never settle and the UI would spin "Solving…" forever. Reject it as a
        // REAL error (not an abort, so the pipeline surfaces it) and respawn next.
        if (_inflight) {
          const job = _inflight; _inflight = null;
          try { w.removeEventListener('message', job.onMessage); } catch (_) {}
          job.reject(new Error('solver worker crashed: ' + detail));
        }
        _disposeWorker();
      });
      w.addEventListener('messageerror', (ev) => {
        console.error('[primary-solver worker] messageerror (postMessage clone failed):', ev);
        if (_inflight) {
          const job = _inflight; _inflight = null;
          try { w.removeEventListener('message', job.onMessage); } catch (_) {}
          job.reject(new Error('solver worker message error (structured-clone failed)'));
        }
        _disposeWorker();
      });
      _worker = w;
    })().catch((err) => {
      console.warn('[primary-solver-worker] Worker unavailable (' + (err && err.message || err) +
        '). Falling back to main-thread solver. ' +
        'Serve via a local web server (e.g. `python -m http.server`) to enable.');
      _mainThreadFallback = true;
      _readyPromise = null;
    });
    await _readyPromise;
  }

  function solve(hData, opts, runOpts) {
    runOpts = runOpts || {};
    return ensureReady().then(() => {
      // Main-thread fallback path — used when the Worker bundle could not be
      // built. Yields one microtask so callers are still .then-able.
      if (_mainThreadFallback || !_worker) {
        return Promise.resolve().then(() => {
          // QD here is the main-thread QD namespace.
          return global.QD.solveInverseQD(hData, opts || {});
        });
      }

      // Cancel any prior in-flight job before posting a new one. The
      // previous caller's promise rejects with { aborted: true }.
      if (_inflight) {
        _inflight.reject({ aborted: true, superseded: true });
        try { _worker.removeEventListener('message', _inflight.onMessage); } catch (_) {}
        _inflight = null;
      }

      const jobId = _nextJobId++;
      return new Promise((resolve, reject) => {
        const onMessage = (e) => {
          const m = e.data;
          if (!m || m.kind !== 'solve' || m.jobId !== jobId) return;
          try { _worker.removeEventListener('message', onMessage); } catch (_) {}
          _inflight = null;
          if (m.error) reject(new Error(m.error));
          else resolve(m.result);
        };
        _inflight = { jobId, resolve, reject, onMessage };
        _worker.addEventListener('message', onMessage);

        // Forward AbortSignal -> cancel(). The signal can be supplied by the
        // caller (e.g. ui.js when a fresh edit supersedes the prior solve).
        const signal = runOpts.signal;
        if (signal) {
          if (signal.aborted) { cancel(); return; }
          const onAbort = () => { cancel(); };
          signal.addEventListener('abort', onAbort, { once: true });
        }

        _worker.postMessage({ kind: 'solve', jobId, hData, opts: opts || {} });
      });
    });
  }

  function cancel() {
    // Terminate-and-recreate is the cheap way to preempt deeply-nested Newton.
    // Subsequent solve() calls trigger a fresh ensureReady().
    _disposeWorker();
  }

  function isBusy() { return _inflight !== null; }

  // ---------------------------------------------------------------------------
  // Aux worker — background alternate-solution search (A3).
  //
  // A SEPARATE Worker instance from the same bundle, dedicated to
  // QD.searchAlternates. Keeping it separate from the primary `solve` worker
  // means a background alt-search can never queue behind (delaying) or preempt
  // an interactive primary solve, and it can be terminated freely on
  // supersession without forcing the primary worker to respawn. Previously the
  // alt-search ran synchronously on the main thread (setTimeout chunks), which
  // janked the 2D plot on every chunk.
  // ---------------------------------------------------------------------------
  /** @type {Worker|null} */
  let _auxWorker = null;
  /** @type {Promise<void>|null} */
  let _auxReady = null;
  let _auxNextJobId = 1;
  /** @type {InflightJob|null} */
  let _auxInflight = null;

  function _disposeAux() {
    if (_auxWorker) {
      try { _auxWorker.terminate(); } catch (_) { /* ignore */ }
      _auxWorker = null;
    }
    _auxReady = null;
    if (_auxInflight) {
      _auxInflight.reject({ aborted: true });
      _auxInflight = null;
    }
  }

  async function ensureAuxReady() {
    if (_mainThreadFallback) return;
    if (_auxWorker) return;
    if (_auxReady) { await _auxReady; return; }
    _auxReady = (async () => {
      const url = await getBundleURL();
      const w = new Worker(url);
      w.addEventListener('error', (ev) => {
        const detail = (ev.message || ev) + ' @ ' + (ev.filename || 'bundle') + ':' + (ev.lineno || '?');
        console.error('[primary-solver aux worker] error: ' + detail);
        if (_auxInflight) {
          const job = _auxInflight; _auxInflight = null;
          try { w.removeEventListener('message', job.onMessage); } catch (_) {}
          job.reject(new Error('alt-search worker crashed: ' + detail));
        }
        _disposeAux();
      });
      _auxWorker = w;
    })().catch((err) => {
      console.warn('[primary-solver-worker] Aux worker unavailable (' + (err && err.message || err) +
        '). Alternate search will run on the main thread.');
      _mainThreadFallback = true;
      _auxReady = null;
    });
    await _auxReady;
  }

  // Run one alternate-search pass off the main thread. Resolves with the array
  // of candidate solutions QD.searchAlternates returns. Posting a new request
  // supersedes any prior in-flight aux job (its promise rejects { aborted }).
  function searchAlternatesAsync(hData, norm, known, opts) {
    return ensureAuxReady().then(() => {
      if (_mainThreadFallback || !_auxWorker) {
        return Promise.resolve().then(() =>
          global.QD.searchAlternates(hData, norm, known || [], opts || {}));
      }
      if (_auxInflight) {
        _auxInflight.reject({ aborted: true, superseded: true });
        try { _auxWorker.removeEventListener('message', _auxInflight.onMessage); } catch (_) {}
        _auxInflight = null;
      }
      const jobId = _auxNextJobId++;
      return new Promise((resolve, reject) => {
        const onMessage = (e) => {
          const m = e.data;
          if (!m || m.kind !== 'altSearch' || m.jobId !== jobId) return;
          try { _auxWorker.removeEventListener('message', onMessage); } catch (_) {}
          _auxInflight = null;
          if (m.error) reject(new Error(m.error));
          else resolve(m.result);
        };
        _auxInflight = { jobId, resolve, reject, onMessage };
        _auxWorker.addEventListener('message', onMessage);
        _auxWorker.postMessage({ kind: 'altSearch', jobId, hData, norm, known: known || [], opts: opts || {} });
      });
    });
  }

  function cancelAux() { _disposeAux(); }
  function isAuxBusy() { return _auxInflight !== null; }

  // ---------------------------------------------------------------------------
  // Live worker — per-drag-frame warm-start solve (QD.liveSolveStep).
  //
  // A THIRD Worker instance from the same bundle, dedicated to the interactive
  // drag path. Kept separate from the primary `solve` worker so a burst of live
  // frames can't queue behind / preempt a debounced full solve (and vice
  // versa). Cancellation is cancel-and-replace by DROPPING THE LISTENER — not
  // terminate: each live job is a single bounded Newton (≈ms), so the worker
  // frees itself almost immediately and a respawn-per-frame would cost more than
  // it saves. The main side additionally guards with a token (ui-solve.js) so a
  // late result can't overwrite newer state.
  // ---------------------------------------------------------------------------
  /** @type {Worker|null} */
  let _liveWorker = null;
  /** @type {Promise<void>|null} */
  let _liveReady = null;
  let _liveNextJobId = 1;
  /** @type {InflightJob|null} */
  let _liveInflight = null;

  function _disposeLive() {
    if (_liveWorker) {
      try { _liveWorker.terminate(); } catch (_) { /* ignore */ }
      _liveWorker = null;
    }
    _liveReady = null;
    if (_liveInflight) {
      _liveInflight.reject({ aborted: true });
      _liveInflight = null;
    }
  }

  async function ensureLiveReady() {
    if (_mainThreadFallback) return;
    if (_liveWorker) return;
    if (_liveReady) { await _liveReady; return; }
    _liveReady = (async () => {
      const url = await getBundleURL();
      const w = new Worker(url);
      w.addEventListener('error', (ev) => {
        const detail = (ev.message || ev) + ' @ ' + (ev.filename || 'bundle') + ':' + (ev.lineno || '?');
        console.error('[primary-solver live worker] error: ' + detail);
        if (_liveInflight) {
          const job = _liveInflight; _liveInflight = null;
          try { w.removeEventListener('message', job.onMessage); } catch (_) {}
          job.reject(new Error('live-solve worker crashed: ' + detail));
        }
        _disposeLive();
      });
      _liveWorker = w;
    })().catch((err) => {
      console.warn('[primary-solver-worker] Live worker unavailable (' + (err && err.message || err) +
        '). Live drag solve will run on the main thread.');
      _mainThreadFallback = true;
      _liveReady = null;
    });
    await _liveReady;
  }

  // Run one live (drag-frame) solve off the main thread. `initPhi` is the
  // caller-chosen seed (warm clone of the previous phi, or a fresh init); `opts`
  // is forwarded to QD.liveSolveStep ({ newton, numSamples, wantOriginInside }).
  // Posting a new request supersedes any prior in-flight live job (its promise
  // rejects { aborted, superseded }). Falls back to a main-thread liveSolveStep
  // when the worker bundle is unavailable (e.g. file:// origin).
  function liveSolveAsync(hData, initPhi, opts) {
    return ensureLiveReady().then(() => {
      if (_mainThreadFallback || !_liveWorker) {
        return Promise.resolve().then(() =>
          global.QD.liveSolveStep(hData, initPhi, opts || {}));
      }
      if (_liveInflight) {
        _liveInflight.reject({ aborted: true, superseded: true });
        try { _liveWorker.removeEventListener('message', _liveInflight.onMessage); } catch (_) {}
        _liveInflight = null;
      }
      const jobId = _liveNextJobId++;
      return new Promise((resolve, reject) => {
        const onMessage = (e) => {
          const m = e.data;
          if (!m || m.kind !== 'liveSolve' || m.jobId !== jobId) return;
          try { _liveWorker.removeEventListener('message', onMessage); } catch (_) {}
          _liveInflight = null;
          if (m.error) reject(new Error(m.error));
          else resolve(m.result);
        };
        _liveInflight = { jobId, resolve, reject, onMessage };
        _liveWorker.addEventListener('message', onMessage);
        _liveWorker.postMessage({ kind: 'liveSolve', jobId, hData, initPhi, opts: opts || {} });
      });
    });
  }

  function cancelLive() { _disposeLive(); }
  function isLiveBusy() { return _liveInflight !== null; }

  // Expose under QD namespace if available, else stash on global.
  const ns = { ensureReady, solve, cancel, isBusy,
    // Background alternate-search (A3) — runs on a dedicated aux worker.
    searchAlternates: searchAlternatesAsync, cancelAux, isAuxBusy,
    // Live drag-frame solve — runs on a dedicated live worker (Tier-2 pole-drag).
    liveSolve: liveSolveAsync, cancelLive, isLiveBusy,
    // Diagnostics — used by tests / dev tools.
    _isMainThreadFallback() { return _mainThreadFallback; },
    _hasWorker()           { return _worker !== null; },
  };
  if (global.QD) global.QD.PrimarySolverWorker = ns;
  else if (global.module && global.module.exports) global.module.exports = ns;
  else global.QD_PrimarySolverWorker = ns;

})(typeof window !== 'undefined' ? window : (typeof self !== 'undefined' ? self : globalThis));
