// =============================================================================
// param-slice-render.js -- Adaptive 2-D parameter-slice render engine.
//
// Extracted from param-slice-ui.js by the Phase-3 UI modularization (item E).
// Exposes QD_UI.installParamSliceRender(psCtx); param-slice-ui.js captures
// runAdaptive2D into an IIFE-local binding so startRun calls it unchanged.
//
// runAdaptive2D drives the progressive quadtree sweep: a coarse stride pass then
// stride/2 refinement passes (subdividing only cells whose corners disagree),
// each dispatched to the worker pool and painted as it lands, finished by a
// nearest-neighbour coverage fill. It also owns the warm-hint spatial index
// (bucketed insertPhi/nearestPhi, published on sliceState for the hover preview)
// so refined pixels warm-start from a cached neighbour phi. All callbacks (pool,
// paintCellBlock, paintImage, logErrorSamples, cancelToken, onProgress) arrive
// via its single options argument; the body is a VERBATIM move. Deps via psCtx:
// sliceState + cancelLiveSolve. `PS` (= global.ParamSlice) + performance /
// console are read as globals.
// =============================================================================

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  global.QD_UI.installParamSliceRender = function installParamSliceRender(s) {
    const sliceState     = s.sliceState;
    const cancelLiveSolve = s.cancelLiveSolve;
    const PS             = global.ParamSlice;

  const REFINE_ITER_DELTA = 8;

  async function runAdaptive2D({ pool, baseScenario, mode, axes, xs, ys, n0, n1,
                                 paintCellBlock, paintImage, logErrorSamples,
                                 cancelToken, onProgress }) {
    // Persistent grid state — classification index + iteration count per cell.
    // Stored on sliceState (HANDOFF #33) so the hover-tooltip and the
    // Hovered-QD preview card can read them after the render completes.
    const classGrid = new Uint8Array(n0 * n1).fill(PS.UNKNOWN_CLASS);
    const iterGrid  = new Uint8Array(n0 * n1);
    sliceState.classGrid = classGrid;
    sliceState.iterGrid  = iterGrid;
    sliceState.gridDims  = { n0, n1 };
    // HANDOFF #35 race fix: invalidate the φ cache *immediately* when the
    // new grid dims are published, not 90 lines later when the new
    // nearestPhi closure is wired up. Otherwise the hover handler can briefly
    // read stale φs from the previous render via the still-live old closure.
    sliceState.nearestPhi = function () { return null; };
    sliceState.miniCache = { phi: null, pts: null };
    // HANDOFF #36: also cancel any pending live-solve.
    cancelLiveSolve();

    // Warm-hint spatial index: bucket the grid into ~16×16 buckets so
    // nearestPhi(c, r) is O(1) amortized instead of O(N).
    //
    // Each bucket holds an array of {c, r, phi} entries. Querying looks at
    // the 9 neighboring buckets (the home bucket + 8 neighbors); good
    // enough for nearest-neighbor since adjacent cells are guaranteed
    // within one bucket-width.
    const BUCKETS_PER_AXIS = 16;
    const bucketW = Math.max(1, Math.ceil(n0 / BUCKETS_PER_AXIS));
    const bucketH = Math.max(1, Math.ceil(n1 / BUCKETS_PER_AXIS));
    const bucketCols = Math.ceil(n0 / bucketW);
    const bucketRows = Math.ceil(n1 / bucketH);
    const phiBuckets = new Array(bucketCols * bucketRows);
    for (let i = 0; i < phiBuckets.length; i++) phiBuckets[i] = [];
    let phiCacheSize = 0;
    const PHI_CACHE_CAP = 4096;

    function bucketIdx(c, r) {
      const bc = Math.min(bucketCols - 1, Math.max(0, Math.floor(c / bucketW)));
      const br = Math.min(bucketRows - 1, Math.max(0, Math.floor(r / bucketH)));
      return br * bucketCols + bc;
    }
    function insertPhi(c, r, phi, iterCount) {
      phiBuckets[bucketIdx(c, r)].push({ c, r, phi, iterCount: iterCount | 0 });
      phiCacheSize++;
      // Evict half the cache uniformly when it overflows. Eviction is rare
      // and the spatial distribution roughly stays intact.
      if (phiCacheSize > PHI_CACHE_CAP) {
        for (let i = 0; i < phiBuckets.length; i++) {
          const b = phiBuckets[i];
          if (b.length > 4) b.splice(0, b.length >> 1);
        }
        phiCacheSize = 0;
        for (let i = 0; i < phiBuckets.length; i++) phiCacheSize += phiBuckets[i].length;
      }
    }

    // Choose the coarsest stride: a power of 2 ≤ min(n0,n1)/8, capped so we
    // don't degenerate to "sample every cell" on very small grids.
    let stride = 1;
    while ((stride << 1) <= Math.min(n0, n1) / 4) stride <<= 1;
    stride = Math.max(1, stride);
    const startStride = stride;

    // Estimate total work for progress reporting.
    const totalCellsAtFineGrid = n0 * n1;
    let cellsDone = 0;

    function nearestPhi(c, r) {
      // Scan the 9 buckets around (c, r). Each holds O(K) cached φs where
      // K ≈ PHI_CACHE_CAP / (BUCKETS_PER_AXIS^2) ≈ 16. Total ~144 distance
      // comparisons regardless of cache size — vs O(N) for a flat scan.
      //
      // Returns the nearest cached entry as { phi, iterCount } or null.
      // The iterCount lets the solver speculatively tighten its Newton
      // maxIter cap for the refined sub-pixel (see param-slice-common.js
      // `_solveScenarioBody`).
      if (phiCacheSize === 0) return null;
      const bc = Math.floor(c / bucketW);
      const br = Math.floor(r / bucketH);
      let best = null, bestD = Infinity;
      for (let dbr = -1; dbr <= 1; dbr++) {
        const brI = br + dbr;
        if (brI < 0 || brI >= bucketRows) continue;
        for (let dbc = -1; dbc <= 1; dbc++) {
          const bcI = bc + dbc;
          if (bcI < 0 || bcI >= bucketCols) continue;
          const b = phiBuckets[brI * bucketCols + bcI];
          for (let i = 0; i < b.length; i++) {
            const p = b[i];
            const d = (p.c - c) * (p.c - c) + (p.r - r) * (p.r - r);
            if (d < bestD) { bestD = d; best = p; }
          }
        }
      }
      // Fallback: if no neighbor bucket had anything (e.g. very sparse
      // valid region), do one full scan as a last resort. This is rare in
      // practice because the coarsest pass populates many cells uniformly.
      if (!best) {
        for (let i = 0; i < phiBuckets.length; i++) {
          const b = phiBuckets[i];
          for (let j = 0; j < b.length; j++) {
            const p = b[j];
            const d = (p.c - c) * (p.c - c) + (p.r - r) * (p.r - r);
            if (d < bestD) { bestD = d; best = p; }
          }
        }
      }
      return best ? { phi: best.phi, iterCount: best.iterCount } : null;
    }
    // Expose for the hover-tooltip + Hovered-QD preview card (HANDOFF #33).
    sliceState.nearestPhi = nearestPhi;

    function paintAtStride(s) {
      // For each cell of stride s whose corners agree (in class AND, for
      // VALID, in iter count to within REFINE_ITER_DELTA), fill the whole
      // pixel block with the top-left-corner color. Cells that fail this
      // test will be subdivided in the next refinement pass and their
      // sub-cells painted then; deferring keeps the coarse paint from
      // committing to a misleading top-left iter count for a whole block.
      for (let r = 0; r < n1; r += s) {
        for (let c = 0; c < n0; c += s) {
          const blockCols = Math.min(s, n0 - c);
          const blockRows = Math.min(s, n1 - r);
          const k = classGrid[r * n0 + c];
          if (k === PS.UNKNOWN_CLASS) continue;
          const homog = (s === 1) || PS.cellIsHomogeneous(
            classGrid, iterGrid, n0, n1, c, r, s,
            { iterDelta: REFINE_ITER_DELTA });
          if (homog) {
            const color = PS.colorFor({ cls: PS.IDX_TO_CLASS[k], iterations: iterGrid[r * n0 + c] });
            paintCellBlock(c, r, blockCols, blockRows, color);
          }
        }
      }
    }

    function storeResults(points, results) {
      logErrorSamples(results);
      for (let i = 0; i < results.length; i++) {
        const { c, r } = points[i];
        const idx = r * n0 + c;
        const cls = results[i].cls;
        const iters = Math.min(255, results[i].iterations || 0);
        classGrid[idx] = PS.CLASS_TO_IDX[cls];
        iterGrid[idx]  = iters;
        if (results[i].phiSerialized) {
          insertPhi(c, r, results[i].phiSerialized, iters);
        }
      }
    }

    function buildParams(points) {
      return points.map(({ c, r }) => [
        { ref: axes[0].ref, value: xs[c] },
        { ref: axes[1].ref, value: ys[r] },
      ]);
    }

    function dispatchPoints(points) {
      const params = buildParams(points);
      // Build per-point warm hints by looking up the nearest cached φ from
      // a previous pass. Wrap each hit with `_coarseIter` so the solver can
      // speculatively tighten its Newton maxIter cap. The wrapper is a
      // shallow copy of the φ — the underlying phi object is not mutated
      // (it gets cloned again inside `_solveScenarioBody` via QD.clonePhi).
      const hints = points.map(({ c, r }) => {
        const hit = nearestPhi(c, r);
        if (hit) return Object.assign({}, hit.phi, { _coarseIter: hit.iterCount });
        // No cached neighbour yet (the whole COARSE pass, and the first point of
        // each worker chunk): fall back to the precomputed seed φ so the pixel
        // warm-starts (~0.03 ms) instead of cold-solving (~5.6 ms for PQD). A
        // wrong-basin seed is safe — `_solveScenarioBody` retries cold if the
        // warm Newton fails. Seeding is what keeps PQD slices from paying ~200+
        // cold solves on the coarse pass alone.
        return seedPhi ? Object.assign({}, seedPhi) : null;
      });
      return pool.solveBatch(baseScenario, mode, params, hints);
    }

    // --- Coarse pass: sample every (stride * k, stride * k) corner.
    const coarsePoints = [];
    for (let r = 0; r < n1; r += startStride) {
      for (let c = 0; c < n0; c += startStride) {
        coarsePoints.push({ c, r });
      }
    }
    // Also ensure the right + bottom boundary corners are evaluated so the
    // cornersAgree check at later passes has well-defined neighbors.
    const includeEdge = (c, r) => {
      if (classGrid[r * n0 + c] === PS.UNKNOWN_CLASS &&
          !coarsePoints.some(p => p.c === c && p.r === r)) {
        coarsePoints.push({ c, r });
      }
    };
    for (let r = 0; r < n1; r += startStride) includeEdge(n0 - 1, r);
    for (let c = 0; c < n0; c += startStride) includeEdge(c, n1 - 1);
    includeEdge(n0 - 1, n1 - 1);

    if (cancelToken.cancelled) return;

    // Warm-start seed (PQD perf): solve the base scenario once on the main
    // thread to get ONE valid φ of the right family + pole structure, then use
    // it as the warm hint for every point that has no cached neighbour yet —
    // i.e. the entire coarse pass and each worker chunk's first point. Without
    // it, those ~200+ coarse points each cold-solve (~5.6 ms for PQD); with it
    // they warm-refine (~0.03 ms). The empty-point solve uses the scenario's
    // current (already-valid, user-visible) config; if it fails to solve we
    // simply skip seeding (seedPhi stays null → prior cold-start behaviour).
    let seedPhi = null;
    try {
      const familyTag = PS.MODE_FAMILY_TAG ? PS.MODE_FAMILY_TAG[mode] : undefined;
      const seedR = PS.solveOnePoint(baseScenario, [], null, familyTag);
      if (seedR && seedR.cls === PS.CLASS_VALID && seedR.phiSerialized) {
        seedPhi = seedR.phiSerialized;
      }
    } catch (e) { seedPhi = null; }

    const t0 = performance.now();
    const coarseResults = await dispatchPoints(coarsePoints);
    if (cancelToken.cancelled || !coarseResults) return;
    storeResults(coarsePoints, coarseResults);
    cellsDone += coarsePoints.length;
    paintAtStride(startStride);
    paintImage();
    onProgress(cellsDone, totalCellsAtFineGrid);
    console.log(`[param-slice] coarse pass: ${coarsePoints.length} samples in ${((performance.now()-t0)/1000).toFixed(2)}s`);

    // --- Refinement passes: stride/2 down to 1.
    while (stride > 1) {
      if (cancelToken.cancelled) return;
      const halfStride = stride >> 1;
      // Collect new sample points within cells whose corners disagree.
      const newPoints = [];
      const seen = new Set();
      for (let r = 0; r + stride < n1; r += stride) {
        for (let c = 0; c + stride < n0; c += stride) {
          if (PS.cellIsHomogeneous(classGrid, iterGrid, n0, n1, c, r, stride,
                                   { iterDelta: REFINE_ITER_DELTA })) continue;
          const subPts = PS.subdivisionPoints(c, r, stride, n0, n1);
          for (const p of subPts) {
            const key = p.r * n0 + p.c;
            if (classGrid[key] === PS.UNKNOWN_CLASS && !seen.has(key)) {
              newPoints.push(p);
              seen.add(key);
            }
          }
        }
      }
      if (newPoints.length === 0) {
        stride = halfStride;
        continue;
      }
      const tLevel = performance.now();
      const results = await dispatchPoints(newPoints);
      if (cancelToken.cancelled || !results) return;
      storeResults(newPoints, results);
      cellsDone += newPoints.length;
      paintAtStride(halfStride);
      paintImage();
      onProgress(cellsDone, totalCellsAtFineGrid);
      console.log(`[param-slice] stride=${halfStride}: ${newPoints.length} samples in ${((performance.now()-tLevel)/1000).toFixed(2)}s`);
      stride = halfStride;
    }

    // Final coverage-fill pass (HANDOFF #37). The adaptive mesh only
    // samples cell midpoints during refinement (never corners), so cells
    // whose top-left is at a non-coarse-stride position — e.g. (120, 0),
    // (124, 0), (126, 0) for n0=128, startStride=32 — can stay
    // UNKNOWN_CLASS after the cascade. paintAtStride(1) skips those
    // pixels and the dark canvas-init color (rgb(24,24,24)) shows
    // through, producing the regular black-dot grid the user reported.
    //
    // Plug the gaps by painting each remaining UNKNOWN pixel with its
    // nearest sampled neighbour's colour. Doesn't touch classGrid (so
    // hover-tooltip still honestly says "(no sample)" for these pixels)
    // — purely visual fill. One pass over a grid with isolated 1-pixel
    // gaps is sufficient.
    fillUnpaintedFromNeighbor();
    paintImage();

    function fillUnpaintedFromNeighbor() {
      const NEIGHBOR_OFFSETS = [
        [-1, 0], [1, 0], [0, -1], [0, 1],
        [-1, -1], [1, 1], [-1, 1], [1, -1],
      ];
      let filled = 0;
      for (let r = 0; r < n1; r++) {
        for (let c = 0; c < n0; c++) {
          if (classGrid[r * n0 + c] !== PS.UNKNOWN_CLASS) continue;
          for (const [dc, dr] of NEIGHBOR_OFFSETS) {
            const cc = c + dc, rr = r + dr;
            if (cc < 0 || cc >= n0 || rr < 0 || rr >= n1) continue;
            const kk = classGrid[rr * n0 + cc];
            if (kk === PS.UNKNOWN_CLASS) continue;
            const color = PS.colorFor({
              cls: PS.IDX_TO_CLASS[kk],
              iterations: iterGrid[rr * n0 + cc],
            });
            paintCellBlock(c, r, 1, 1, color);
            filled++;
            break;
          }
        }
      }
      if (filled > 0) {
        console.log(`[param-slice] coverage-fill: ${filled} orphan pixel${filled === 1 ? '' : 's'} painted from nearest neighbour.`);
      }
    }
  }

    return { runAdaptive2D };
  };
})(typeof window !== 'undefined' ? window : globalThis);
