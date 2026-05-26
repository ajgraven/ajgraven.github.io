// =============================================================================
// param-slice-pool.js -- Web Worker pool for the Parameter-slice tab.
//
// Builds a self-contained Worker bundle by fetching the existing solver
// source files and concatenating them with param-slice-common.js plus a
// thin worker-side message handler.  No build step required — works
// directly from the static site.
//
// API:
//   const pool = await ParamSlicePool.create();
//   const handle = pool.runSweep({
//     scenario,         // { hData, norm, opts } (base, never mutated)
//     mode,             // mode string (used by classifier)
//     axes,             // 1-D: [{ ref, min, max, n }]
//                       // 2-D: [{ ref, min, max, n }, { ref, min, max, n }]
//     onTile,           // (tile) => void   // tile = { row, results: [{cls,...}, ...] }
//     onDone,           // ({ tilesDone, totalTiles, msTotal }) => void
//     onError,          // (err) => void  (optional)
//   });
//   // handle.cancel() stops further dispatch.
//
// One "tile" = one row of pixels (1-D mode: a single tile with all pixels).
// Workers warm-start within a row by chaining the previous pixel's valid
// φ into the next QD.newtonSolve call.
// =============================================================================

(function (global) {
  'use strict';

  // Same script set as index.html, in load order. Paths are relative to the
  // HTML document (== relative to /app/) so we can resolve them at runtime.
  const SOLVER_SRC_FILES = [
    'complex.js',
    'taylor.js',
    'solver.js',
    'solver-faber.js',
    'solver-qd.js',
    'solver-uqd.js',
    'solver-lqd-common.js',
    'solver-lqd.js',
    'solver-lqd-singular.js',
    'solver-uqd-lqd.js',
    'solver-uqd-lqd-singular.js',
    'parse-h.js',
    'param-slice/param-slice-common.js',
  ];

  // The worker-side message handler. Self-contained string so we can append
  // it after the bundled solver source.  Uses QD.* / ParamSlice globals
  // exposed by the bundled scripts on `self`.
  //
  // Message protocol (kind: 'tile'):
  //   { jobId, scenario, sweepPoints, warmHints? }
  // - scenario     : { hData, norm, opts, expectedFamilyTag }
  // - sweepPoints  : [[{ref,value}, ...], ...]   one entry per pixel
  // - warmHints    : optional array of φ objects (one per pixel) — when
  //                  provided, overrides the implicit per-chunk warm chain.
  const WORKER_HANDLER = `
;(function () {
  'use strict';
  const PS = self.ParamSlice;

  self.onmessage = function (e) {
    const msg = e.data;
    if (!msg || msg.kind !== 'tile') return;
    const { jobId, scenario, sweepPoints, warmHints } = msg;
    const expectedFamilyTag = scenario.expectedFamilyTag || undefined;
    const results = new Array(sweepPoints.length);
    // One scratch scenario per tile message — mutated in place between
    // pixels so we pay the cloneScenario cost only once per tile rather
    // than once per pixel. Invariant: all pixels in a tile target the
    // same ParamRefs (the sweep's axes), so successive applyParamInPlace
    // calls just overwrite each other's values cleanly.
    const scratch = PS.cloneScenario(scenario);
    let chainWarm = null;
    for (let i = 0; i < sweepPoints.length; i++) {
      const hint = (warmHints && warmHints[i]) || chainWarm;
      const r = PS.solveOnePointWithScratch(scratch, sweepPoints[i], hint, expectedFamilyTag);
      if (r.phiSerialized) chainWarm = r.phiSerialized;
      results[i] = r;
    }
    self.postMessage({ kind: 'tile', jobId, results });
  };
})();
`;

  // ---------------------------------------------------------------------------
  // Bundle builder
  // ---------------------------------------------------------------------------
  // Fetch + concat solver source files once per page load; cache the Blob URL.
  let _bundlePromise = null;
  function getBundleURL() {
    if (_bundlePromise) return _bundlePromise;
    _bundlePromise = (async () => {
      const parts = [];
      // Each solver file gates global registration on `typeof window !==
      // 'undefined'`. In a Worker scope `window` is undefined by default —
      // so QD never lands on `self`, and the family files fail their
      // `window.QD` lookup. Aliasing `window = self` makes the registrations
      // attach to the worker's global object, and lets the family files find
      // `window.QD` set by the preceding `solver.js`.
      parts.push('var window = self;\n');
      for (const f of SOLVER_SRC_FILES) {
        const resp = await fetch(f);
        if (!resp.ok) throw new Error('param-slice-pool: failed to fetch ' + f + ' (' + resp.status + ')');
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

  // ---------------------------------------------------------------------------
  // Pool
  // ---------------------------------------------------------------------------
  class Pool {
    constructor(workers, bundleURL) {
      this.workers = workers;
      this.bundleURL = bundleURL;
      this.idle = workers.slice();         // queue of available workers
      this.pending = [];                   // tile jobs awaiting a worker
      this.activeJobs = new Map();         // jobId → { worker, resolve }
      this._nextJobId = 1;
      this._cancelled = false;
    }

    _dispatch() {
      while (this.idle.length > 0 && this.pending.length > 0 && !this._cancelled) {
        const worker = this.idle.shift();
        const job = this.pending.shift();
        const jobId = this._nextJobId++;
        this.activeJobs.set(jobId, { worker, resolve: job.resolve });
        const onMessage = (e) => {
          const m = e.data;
          if (!m || m.jobId !== jobId) return;
          worker.removeEventListener('message', onMessage);
          this.activeJobs.delete(jobId);
          this.idle.push(worker);
          job.resolve(m.results);
          this._dispatch();
        };
        worker.addEventListener('message', onMessage);
        worker.postMessage({
          kind: 'tile',
          jobId,
          scenario: job.scenario,
          sweepPoints: job.sweepPoints,
          warmHints: job.warmHints || null,
        });
      }
    }

    submitTile(scenario, sweepPoints, warmHints) {
      return new Promise((resolve) => {
        this.pending.push({ scenario, sweepPoints, warmHints, resolve });
        this._dispatch();
      });
    }

    // Dispatch a batch of `points` (param-assignment arrays, one per pixel)
    // by splitting evenly across workers. Each worker gets a contiguous
    // chunk so it can chain its own implicit warm-start within the chunk.
    // `warmHints` is an optional parallel array of φ objects (or null) used
    // as the explicit per-pixel seed; nearby cells in the parent caller
    // produce these hints from already-evaluated valid pixels.
    async solveBatch(scenario, mode, points, warmHints) {
      if (!points || points.length === 0) return [];
      const scenarioWithTag = _attachFamilyTag(scenario, mode);
      const nChunks = Math.min(this.workers.length, points.length);
      const out = new Array(points.length);
      const chunkSize = Math.ceil(points.length / nChunks);
      const promises = [];
      for (let c = 0; c < nChunks; c++) {
        const start = c * chunkSize;
        const end = Math.min(points.length, start + chunkSize);
        if (start >= end) break;
        const chunkPoints = points.slice(start, end);
        const chunkHints  = warmHints ? warmHints.slice(start, end) : null;
        const p = this.submitTile(scenarioWithTag, chunkPoints, chunkHints).then((results) => {
          if (!results) return;
          for (let i = 0; i < results.length; i++) out[start + i] = results[i];
        });
        promises.push(p);
      }
      await Promise.all(promises);
      return out;
    }

    cancel() {
      this._cancelled = true;
      // Resolve outstanding pending jobs with empty results so the
      // orchestration code's awaits don't hang forever.
      for (const job of this.pending) job.resolve(null);
      this.pending.length = 0;
    }

    terminate() {
      this.cancel();
      for (const w of this.workers) {
        try { w.terminate(); } catch (e) { /* ignore */ }
      }
      this.workers.length = 0;
    }

    // Linear (non-adaptive) sweep — dispatch the full grid row-by-row.
    // Used by the 1-D path; the 2-D path uses solveBatch via the adaptive
    // renderer in param-slice-ui.js.
    runSweep({ scenario, mode, axes, onTile, onError }) {
      if (this._cancelled) this._cancelled = false;
      const n0 = axes[0].n;
      const has2 = axes.length === 2;
      const n1 = has2 ? axes[1].n : 1;
      const xs = sampleAxis(axes[0]);
      const ys = has2 ? sampleAxis(axes[1]) : [null];
      const scenarioWithTag = _attachFamilyTag(scenario, mode);

      const t0 = performance.now();
      let tilesDone = 0;
      const promises = [];
      for (let row = 0; row < n1; row++) {
        const yVal = ys[row];
        const sweepPoints = new Array(n0);
        for (let col = 0; col < n0; col++) {
          const pt = [{ ref: axes[0].ref, value: xs[col] }];
          if (has2) pt.push({ ref: axes[1].ref, value: yVal });
          sweepPoints[col] = pt;
        }
        const p = this.submitTile(scenarioWithTag, sweepPoints, null).then((results) => {
          if (results == null) return;
          tilesDone++;
          if (onTile) {
            try { onTile({ row, results, xs, yVal, tilesDone, totalTiles: n1 }); }
            catch (e) { if (onError) onError(e); }
          }
        }).catch((e) => { if (onError) onError(e); });
        promises.push(p);
      }
      const donePromise = Promise.all(promises).then(() => ({
        tilesDone, totalTiles: n1, msTotal: performance.now() - t0,
      }));
      return { cancel: () => this.cancel(), done: donePromise };
    }
  }

  // Attach the expected family tag to the scenario so workers can quickly
  // gate warm-start applicability (mismatched-family phis would crash).
  function _attachFamilyTag(scenario, mode) {
    const ps = global.ParamSlice;
    return Object.assign({}, scenario, {
      expectedFamilyTag: (ps && ps.MODE_FAMILY_TAG) ? ps.MODE_FAMILY_TAG[mode] : undefined,
    });
  }

  function sampleAxis(axis) {
    const { min, max, n } = axis;
    if (n === 1) return [(min + max) / 2];
    const out = new Array(n);
    const step = (max - min) / (n - 1);
    for (let i = 0; i < n; i++) out[i] = min + i * step;
    return out;
  }

  async function createPool(opts = {}) {
    const url = await getBundleURL();
    const n = Math.max(1, Math.min(opts.maxWorkers || (navigator.hardwareConcurrency || 4), 16));
    const workers = [];
    for (let i = 0; i < n; i++) {
      const w = new Worker(url);
      // Surface worker-level errors (script-parse, uncaught throw inside
      // the message handler) to the console so misconfigurations don't get
      // silently absorbed and mis-classified as solver capability refusals.
      w.addEventListener('error', (e) => {
        console.error('[param-slice worker] error: '
          + (e.message || e) + ' @ ' + (e.filename || 'bundle') + ':' + (e.lineno || '?'));
      });
      w.addEventListener('messageerror', (e) => {
        console.error('[param-slice worker] messageerror (postMessage clone failed):', e);
      });
      workers.push(w);
    }
    return new Pool(workers, url);
  }

  // ---------------------------------------------------------------------------
  // Main-thread fallback "pool"
  // ---------------------------------------------------------------------------
  // Used when the Worker bundle fails to build — most commonly when the page
  // is opened via file:// and the browser refuses fetch() of relative paths.
  // Exposes the same `.runSweep` shape as the real Pool; runs solves
  // synchronously on the main thread, yielding via setTimeout(0) between
  // rows so the UI stays responsive.
  class MainThreadPool {
    constructor() {
      this._cancelled = false;
      this.kind = 'main-thread';
    }
    cancel()    { this._cancelled = true; }
    terminate() { this.cancel(); }

    // solveBatch: process all points sequentially on the main thread, yielding
    // between chunks of `chunkYieldEvery` so the UI stays responsive. Per-point
    // explicit warmHints take precedence over the implicit chain.
    async solveBatch(scenario, mode, points, warmHints) {
      const ps = global.ParamSlice;
      if (!ps || !points || points.length === 0) return [];
      const scenarioWithTag = _attachFamilyTag(scenario, mode);
      const expectedFamilyTag = scenarioWithTag.expectedFamilyTag;
      const out = new Array(points.length);
      // Scratch-scenario reuse — one clone per batch, mutated in place between pixels.
      const scratch = ps.cloneScenario(scenarioWithTag);
      let chainWarm = null;
      const yieldEvery = 8;
      for (let i = 0; i < points.length; i++) {
        if (this._cancelled) break;
        const hint = (warmHints && warmHints[i]) || chainWarm;
        const r = ps.solveOnePointWithScratch(scratch, points[i], hint, expectedFamilyTag);
        if (r.phiSerialized) chainWarm = r.phiSerialized;
        out[i] = r;
        if ((i & (yieldEvery - 1)) === (yieldEvery - 1)) {
          await new Promise(res => setTimeout(res, 0));
        }
      }
      return out;
    }

    runSweep({ scenario, mode, axes, onTile, onError }) {
      const has2 = axes.length === 2;
      const n0 = axes[0].n;
      const n1 = has2 ? axes[1].n : 1;
      const xs = sampleAxis(axes[0]);
      const ys = has2 ? sampleAxis(axes[1]) : [null];
      const t0 = performance.now();
      let tilesDone = 0;
      this._cancelled = false;

      const donePromise = (async () => {
        for (let row = 0; row < n1; row++) {
          if (this._cancelled) break;
          const yVal = ys[row];
          const sweepPoints = new Array(n0);
          for (let col = 0; col < n0; col++) {
            const pt = [{ ref: axes[0].ref, value: xs[col] }];
            if (has2) pt.push({ ref: axes[1].ref, value: yVal });
            sweepPoints[col] = pt;
          }
          let results;
          try { results = await this.solveBatch(scenario, mode, sweepPoints, null); }
          catch (e) { if (onError) onError(e); continue; }
          tilesDone++;
          if (onTile) {
            try { onTile({ row, results, xs, yVal, tilesDone, totalTiles: n1 }); }
            catch (e) { if (onError) onError(e); }
          }
        }
        return { tilesDone, totalTiles: n1, msTotal: performance.now() - t0 };
      })();
      return { cancel: () => this.cancel(), done: donePromise };
    }
  }

  // Try the Worker pool first; on any failure (fetch blocked on file://,
  // CSP, etc.) fall back to the main-thread pool with a warning.
  async function createPoolWithFallback(opts = {}) {
    try {
      return await createPool(opts);
    } catch (e) {
      console.warn('[param-slice] Worker pool unavailable (' + (e.message || e) +
        '). Falling back to main-thread solver. Open the app via a local web ' +
        'server (e.g. `python -m http.server`) to enable the Worker pool.');
      return new MainThreadPool();
    }
  }

  global.ParamSlicePool = {
    create: createPoolWithFallback,
    createWorkerOnly: createPool,        // for tests / explicit opt-out of fallback
    MainThreadPool,                       // exposed so tests can drive it directly
    _getBundleURL: getBundleURL,
  };
})(typeof self !== 'undefined' ? self : window);
