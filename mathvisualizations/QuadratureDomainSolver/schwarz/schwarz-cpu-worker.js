// @ts-check
// =============================================================================
// schwarz-cpu-worker.js -- Dedicated Web Worker for the CPU Schwarz escape-time
// field (the fallback render path used when the GPU shader is unavailable or
// refuses a family, e.g. PQDs).
//
// Why: the CPU path builds an Int16 escape-count field by iterating σ per
// pixel (W×H, up to ~384²). Even time-sliced across requestAnimationFrame the
// work runs on the main thread and competes with input/paint. This worker
// rebuilds the Schwarz handle from the *serializable* φ + boundary samples
// (QD.Schwarz.buildSchwarzFromPhi takes plain data — no closures cross the
// wire) and computes the whole progressive pyramid off-thread, streaming one
// transferable field snapshot back per pass (4×4 → 2×2 → 1×1).
//
// The main thread still keeps its OWN sState.schwarz handle for interactivity
// (hover readout, orbit tracing, preimage tree, …) — only the heavy field
// compute is offloaded.
//
// API:
//   QD.SchwarzCpuWorker.isUsable() -> bool
//     Synchronous gate: false on file:// (workers can't fetch the bundle) or
//     when Worker is unavailable, or after a load failure. The Schwarz UI
//     checks this and falls back to its in-process pyramid when false.
//
//   QD.SchwarzCpuWorker.renderField(params, callbacks) -> { cancel() }
//     params:    { phi, boundaryPts, view:{cx,cy,scale,cssW,cssH}, W, H,
//                  maxIter, strides:[4,2,1] }
//     callbacks: { onPass(msg), onError(err), onUnavailable() }
//       onPass:  { stride, W, H, field:Int16Array, fieldKind:Uint8Array, done }
//                — fieldKind uses the same KIND+1 offset as the main thread
//                  (0 = unresolved; fundamental/escaped/interior/invalid = 1..4;
//                  outside = 5). One call per pyramid pass; `done` on the last.
//       onUnavailable: fired when no worker could be created — the caller
//                  should fall back to its in-process renderer.
//     A new renderField() preempts any in-flight render (terminate + recreate;
//     the bundle URL is cached, so respawn is cheap). The returned handle's
//     cancel() does the same.
//
//   QD.SchwarzCpuWorker.cancel() -> void   (terminate any in-flight render)
//
// Mirrors the bundle-building pattern of primary-solver-worker.js: fetch the
// shared WORKER_BUNDLE_FILES (cache-busted by CACHE_VERSION) plus
// schwarz-common.js, concatenate with `var window = self;` so the solver/
// schwarz files' `window.QD` registrations attach to the worker global, append
// a message handler, and Blob → object URL.
// =============================================================================

(function (global) {
  'use strict';

  // Worker bundle = the shared solver source (single source of truth in
  // asset-manifest.js) + the Schwarz math kernel. buildSchwarzFromPhi +
  // escapeTime live in schwarz-common.js and only need QD.Complex /
  // QD.LqdCommon / the family registry, all of which WORKER_BUNDLE_FILES
  // already provides. No inline fallback — fail loud if the manifest is
  // missing (index.html loads asset-manifest.js first).
  const _CORE = global.QD_ASSET_MANIFEST && global.QD_ASSET_MANIFEST.WORKER_BUNDLE_FILES;
  if (!_CORE) {
    throw new Error(
      'schwarz-cpu-worker.js: QD_ASSET_MANIFEST.WORKER_BUNDLE_FILES unavailable — ' +
      'asset-manifest.js must load before this file.');
  }
  const SRC_FILES = [..._CORE, 'schwarz/schwarz-common.js'];

  // Worker-side handler. Rebuilds the Schwarz handle from the posted φ +
  // boundary samples and runs the same progressive pyramid as the main-thread
  // fallback (schwarz-ui.js), posting a transferable snapshot after each pass.
  const WORKER_HANDLER = `
;(function () {
  'use strict';
  const S = self.QD && self.QD.Schwarz;
  self.onmessage = function (e) {
    const msg = e.data;
    if (!msg || msg.kind !== 'schwarzRender') return;
    const jobId = msg.jobId;
    if (!S || typeof S.buildSchwarzFromPhi !== 'function') {
      self.postMessage({ kind: 'schwarzError', jobId, error: 'QD.Schwarz unavailable in worker' });
      return;
    }
    let sw;
    try {
      // hData is unused by every adapter; pass null. φ + boundaryPts are plain data.
      sw = S.buildSchwarzFromPhi(msg.phi, null, msg.boundaryPts || []);
    } catch (err) {
      self.postMessage({ kind: 'schwarzError', jobId, error: String(err && err.stack || err) });
      return;
    }
    const W = msg.W, H = msg.H, maxIter = msg.maxIter;
    const v = msg.view;
    const domain = msg.domain || 'w';   // 'z' → sample 𝔻 and lift via w = φ(z)
    const cssW = v.cssW, cssH = v.cssH, cx = v.cx, cy = v.cy, scale = v.scale;
    const pxPerCellX = cssW / W, pxPerCellY = cssH / H;
    const field = new Int16Array(W * H);
    const kind  = new Uint8Array(W * H);     // KIND+1 offset; 0 = unresolved
    const strides = msg.strides || [4, 2, 1];
    for (let s = 0; s < strides.length; s++) {
      const stride = strides[s];
      for (let row = 0; row < H; row++) {
        // Per-row warm-start chain (matches the main-thread renderer): the
        // converged ψ-seed from the left neighbor seeds the next pixel's Newton.
        let leftSeed = null;
        for (let col = 0; col < W; col++) {
          if ((row % stride) !== 0 || (col % stride) !== 0) continue;
          const idx = row * W + col;
          if (kind[idx] && stride > 1) continue;        // already resolved coarser
          const px = (col + 0.5) * pxPerCellX;
          const py = (row + 0.5) * pxPerCellY;
          const aRe = cx + (px - cssW / 2) / scale;
          const aIm = cy - (py - cssH / 2) / scale;     // y flip (screen → world)
          let wpt;
          if (domain === 'z') {
            const r2 = aRe * aRe + aIm * aIm;
            if (sw.unbounded ? r2 <= 1 : r2 >= 1) { field[idx] = 0; kind[idx] = 4 + 1; continue; }
            wpt = sw.evalPhi({ re: aRe, im: aIm });
            if (!wpt || !isFinite(wpt.re) || !isFinite(wpt.im)) { field[idx] = 0; kind[idx] = 4 + 1; continue; }
          } else {
            wpt = { re: aRe, im: aIm };
          }
          if (!sw.isInOmega(wpt)) {
            field[idx] = 0;
            kind[idx]  = 4 + 1;                          // KIND_OUTSIDE + 1
          } else {
            const et = S.escapeTime(wpt, sw, { maxIter: maxIter, initialSeedHint: leftSeed });
            field[idx] = et.n;
            kind[idx] = ((et.kind === 'fundamental') ? 0 :
                         (et.kind === 'escaped')     ? 1 :
                         (et.kind === 'interior')    ? 2 : 3) + 1;
            if (et.firstZ) leftSeed = et.firstZ;
          }
        }
      }
      // Snapshot this pass and transfer copies (the worker keeps its own
      // field/kind arrays for the next, finer pass; transferring the live
      // buffers would detach them).
      const fCopy = field.slice();
      const kCopy = kind.slice();
      self.postMessage(
        { kind: 'schwarzPass', jobId, stride, W, H, field: fCopy, fieldKind: kCopy,
          done: (s === strides.length - 1) },
        [fCopy.buffer, kCopy.buffer]);
    }
  };
})();
`;

  // ---------------------------------------------------------------------------
  // Bundle builder (cached at module scope; reused across worker restarts).
  // ---------------------------------------------------------------------------
  let _bundlePromise = null;
  function getBundleURL() {
    if (_bundlePromise) return _bundlePromise;
    _bundlePromise = (async () => {
      const parts = ['var window = self;\n'];
      const _ver = (global.QD_ASSET_MANIFEST && global.QD_ASSET_MANIFEST.CACHE_VERSION) || '0';
      for (const f of SRC_FILES) {
        const resp = await fetch(f + '?v=' + encodeURIComponent(_ver));
        if (!resp.ok) throw new Error('schwarz-cpu-worker: failed to fetch ' + f + ' (' + resp.status + ')');
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

  /** @type {Worker|null} */
  let _worker = null;
  /** @type {Promise<void>|null} */
  let _readyPromise = null;
  let _nextJobId = 1;
  /** @type {{ jobId:number, onMessage:(e:MessageEvent)=>void }|null} */
  let _inflight = null;
  let _mainThreadFallback = false;

  function _disposeWorker() {
    if (_worker) {
      try { _worker.terminate(); } catch (_) { /* ignore */ }
      _worker = null;
    }
    _readyPromise = null;
    _inflight = null;
  }

  async function ensureReady() {
    if (_mainThreadFallback) return;
    if (_worker) return;
    if (_readyPromise) { await _readyPromise; return; }
    _readyPromise = (async () => {
      const url = await getBundleURL();
      const w = new Worker(url);
      w.addEventListener('error', (ev) => {
        console.error('[schwarz-cpu worker] error: '
          + (ev.message || ev) + ' @ ' + (ev.filename || 'bundle') + ':' + (ev.lineno || '?'));
      });
      _worker = w;
    })().catch((err) => {
      console.warn('[schwarz-cpu-worker] Worker unavailable (' + (err && err.message || err) +
        '). CPU Schwarz render will run on the main thread.');
      _mainThreadFallback = true;
      _readyPromise = null;
    });
    await _readyPromise;
  }

  // Synchronous viability gate for the caller's "use worker vs in-process"
  // decision. file:// can't fetch the bundle (same constraint as the other
  // workers), so report unusable there up front.
  function isUsable() {
    if (_mainThreadFallback) return false;
    if (typeof Worker === 'undefined') return false;
    if (typeof location !== 'undefined' && location.protocol === 'file:') return false;
    return true;
  }

  function renderField(params, cbs) {
    cbs = cbs || {};
    // Preempt any in-flight render. The worker can't be interrupted mid-pass,
    // so terminate + recreate (bundle URL is cached → cheap) to start the new
    // render immediately instead of queueing behind the stale one.
    if (_inflight) _disposeWorker();
    let cancelled = false;
    ensureReady().then(() => {
      if (cancelled) return;
      if (_mainThreadFallback || !_worker) { if (cbs.onUnavailable) cbs.onUnavailable(); return; }
      const jobId = _nextJobId++;
      const onMessage = (e) => {
        const m = e.data;
        if (!m || m.jobId !== jobId) return;
        if (m.kind === 'schwarzError') {
          _finish();
          if (cbs.onError) cbs.onError(m.error);
          return;
        }
        if (m.kind === 'schwarzPass') {
          if (cbs.onPass) cbs.onPass(m);
          if (m.done) _finish();
        }
      };
      function _finish() {
        try { _worker.removeEventListener('message', onMessage); } catch (_) { /* ignore */ }
        if (_inflight && _inflight.jobId === jobId) _inflight = null;
      }
      _inflight = { jobId, onMessage };
      _worker.addEventListener('message', onMessage);
      try {
        _worker.postMessage(Object.assign({ kind: 'schwarzRender', jobId }, params));
      } catch (err) {
        // Non-clonable payload (shouldn't happen — φ is plain data). Treat as
        // unavailable so the caller falls back to its in-process renderer.
        _finish();
        if (cbs.onUnavailable) cbs.onUnavailable();
      }
    }).catch(() => { if (cbs.onUnavailable) cbs.onUnavailable(); });
    return { cancel() { cancelled = true; _disposeWorker(); } };
  }

  function cancel() { _disposeWorker(); }
  function isBusy() { return _inflight !== null; }

  const ns = { ensureReady, isUsable, renderField, cancel, isBusy,
    _isMainThreadFallback() { return _mainThreadFallback; },
    _hasWorker() { return _worker !== null; },
  };
  if (global.QD) global.QD.SchwarzCpuWorker = ns;
  else if (global.module && global.module.exports) global.module.exports = ns;
  else global.QD_SchwarzCpuWorker = ns;

})(typeof window !== 'undefined' ? window : (typeof self !== 'undefined' ? self : globalThis));
