// =============================================================================
// schwarz-render.js -- Progressive escape-time renderer for the Schwarz tab.
//
// Extracted from schwarz-ui.js by the Phase-3 UI modularization (item E).
// Exposes QD_UI.installSchwarzRender(sCtx); schwarz-ui.js captures the debounced
// entry point requestRecompute into an IIFE-local binding, so every call site is
// unchanged.
//
// This is the render dispatcher + CPU progressive pyramid: doRecompute picks the
// GPU path (one synchronous frame) or the CPU path (4x4 -> 2x2 -> 1x1, off-thread
// via QD.SchwarzCpuWorker when usable, else chunked across rAF). Bodies are
// VERBATIM moves. Deps via sCtx: sState, the paint fns (clearCanvas / paintAll /
// paintBoundaryOnTop / paintOrbit / paintZView / setProgress, from schwarz-paint.js), the
// GPU/geometry helpers (syncCanvasSize / activeRenderer / showGLLayer), and the
// KIND_* pixel-class constants. `QD` (incl. QD.SchwarzCpuWorker) is the global.
// =============================================================================

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  global.QD_UI.installSchwarzRender = function installSchwarzRender(s) {
    const sState             = s.sState;
    const clearCanvas        = s.clearCanvas;
    const paintAll           = s.paintAll;
    const paintBoundaryOnTop = s.paintBoundaryOnTop;
    const paintOrbit         = s.paintOrbit;
    const paintZView         = s.paintZView;
    const setProgress        = s.setProgress;
    const syncCanvasSize     = s.syncCanvasSize;
    const activeRenderer     = s.activeRenderer;
    const showGLLayer        = s.showGLLayer;
    const isSchwarzActive    = s.isSchwarzActive;
    const KIND_FUND          = s.KIND_FUND;
    const KIND_ESC           = s.KIND_ESC;
    const KIND_INT           = s.KIND_INT;
    const KIND_INV           = s.KIND_INV;
    const KIND_OUTSIDE       = s.KIND_OUTSIDE;

  let recomputeTimer = null;
  function requestRecompute() {
    if (recomputeTimer) clearTimeout(recomputeTimer);
    recomputeTimer = setTimeout(() => { recomputeTimer = null; doRecompute(); }, 80);
  }

  function doRecompute() {
    if (!sState.schwarz) { clearCanvas(); return; }
    // The debounced requestRecompute timer (~80 ms) can fire AFTER the user has
    // left the Schwarz tab. doRecompute re-bumps renderToken below, so the token
    // guard can't catch this; bail on the active-tab check instead, otherwise
    // the GPU branch would re-show the GL layer and the pyramid would blit over
    // whatever tab is now active. (The tab-leave handler cancels hover/click
    // timers but not this debounce.)
    if (!isSchwarzActive()) return;
    syncCanvasSize();

    // Invalidate any prior render up front (covers CPU→GPU: the GPU branch
    // below returns without touching renderToken, so a stale CPU-worker pass
    // could otherwise paint over the GPU image). Cancelling the worker also
    // frees it; a fresh render re-creates it on demand.
    const myToken = ++sState.renderToken;
    if (QD.SchwarzCpuWorker) QD.SchwarzCpuWorker.cancel();

    // GPU path: synchronous, complete in one frame. Used for both the plane and
    // z-disk views — the shader's u_viewMode branch lifts z → w = φ(z) for z-disk
    // (see schwarz-webgl.js). PQD families set capacityError in setPhi, so
    // activeRenderer() returns 'cpu' for them and they take the CPU path below.
    const inZ = sState.viewMode === 'z';
    if (activeRenderer() === 'gpu') {
      showGLLayer(true);
      const t0 = performance.now();
      try {
        sState.gpu.setColormap(sState.grid.colormap);
        sState.gpu.render(inZ ? sState.zView : sState.view, {
          maxIter:   sState.grid.maxIter,
          scaleMode: sState.grid.scaleMode,
          modK:      sState.grid.modK,
          viewMode:  inZ ? 'z' : 'w',
        });
      } catch (e) {
        // GPU render failed (e.g. context lost). Fall through to CPU path.
        sState.gpuMsg = 'GPU render failed; using CPU. ' + (e.message || e);
        // Continue below — CPU pyramid.
      }
      if (!sState.gpuMsg || sState.gpuMsg.indexOf('failed') === -1) {
        // Field/fieldKind aren't populated under GPU rendering — hover readout
        // will fall back to coordinates-only.
        sState.field = null; sState.fieldKind = null;
        // Overlays drawn on top of the GL field (no field clearing needed).
        // z-disk: paintZView in overlay-only mode (GL drew the field); plane:
        // boundary + orbit.
        if (inZ) {
          paintZView(true);
        } else {
          paintBoundaryOnTop();
          paintOrbit();
        }
        const ms = (performance.now() - t0).toFixed(0);
        setProgress('GPU render: ' + ms + ' ms' + (sState.gpuMsg ? '  (' + sState.gpuMsg + ')' : ''));
        return;
      }
    }

    // CPU progressive pyramid path. Hide the GL layer so a stale GPU image
    // doesn't peek through edge cases. (renderToken was already bumped at the
    // top of doRecompute; myToken is captured there.)
    showGLLayer(false);
    sState.rendering = true;
    setProgress('Pass 1/3 (coarse) ...');
    // Allocate field at target resolution.
    const res = sState.grid.resolution;
    const aspect = sState.view.cssW / sState.view.cssH;
    let W, H;
    if (aspect >= 1) { W = res; H = Math.max(1, Math.round(res / aspect)); }
    else             { H = res; W = Math.max(1, Math.round(res * aspect)); }
    sState.field     = new Int16Array(W * H);
    sState.fieldKind = new Uint8Array(W * H);
    sState.fieldW = W; sState.fieldH = H;

    // Prefer the dedicated CPU worker (A7) — it computes the whole pyramid
    // off-thread and streams a field snapshot per pass. Falls back to the
    // in-process pyramid (below) on file:// / no-Worker / clone failure.
    if (QD.SchwarzCpuWorker && QD.SchwarzCpuWorker.isUsable()) {
      _renderCpuViaWorker(myToken, W, H);
      return;
    }
    _renderCpuPyramid(myToken);
  }

  // In-process progressive pyramid (main-thread fallback for doRecompute).
  // Pass 1: every 4th pixel → 4×4 blocks; pass 2: 2×2; pass 3: per-pixel.
  function _renderCpuPyramid(myToken) {
    chainPass(myToken, 4, () =>
      chainPass(myToken, 2, () =>
        chainPass(myToken, 1, () => {
          if (myToken !== sState.renderToken) return;
          sState.rendering = false;
          setProgress('');
          paintAll();
        })));
  }

  // Off-thread CPU render via QD.SchwarzCpuWorker (A7). The worker rebuilds the
  // Schwarz handle from the serializable φ + boundary samples and posts one
  // transferable field snapshot per pyramid pass; we adopt each snapshot, fill
  // the coarse cells (reusing the same routine as the in-process path), and
  // repaint. Stale snapshots (token mismatch) are discarded. Any failure path
  // falls back to the in-process pyramid so the tab always renders.
  function _renderCpuViaWorker(myToken, W, H) {
    const inZ = sState.viewMode === 'z';
    const v = inZ ? sState.zView : sState.view;
    const params = {
      phi:         sState.phiSnapshot,
      boundaryPts: sState.boundarySnapshot || [],
      view:        { cx: v.cx, cy: v.cy, scale: v.scale, cssW: v.cssW, cssH: v.cssH },
      domain:      inZ ? 'z' : 'w',          // 'z' → sample 𝔻, lift via w = φ(z)
      W, H,
      maxIter:     sState.grid.maxIter,
      strides:     [4, 2, 1],
    };
    const fallback = () => {
      if (myToken !== sState.renderToken) return;
      _renderCpuPyramid(myToken);
    };
    sState._cpuWorkerHandle = QD.SchwarzCpuWorker.renderField(params, {
      onPass(m) {
        if (myToken !== sState.renderToken) return;   // superseded — discard
        sState.field = m.field; sState.fieldKind = m.fieldKind;
        sState.fieldW = m.W; sState.fieldH = m.H;
        if (m.stride > 1) fillFromCoarseSamples(m.stride);
        paintAll();
        if (m.done) { sState.rendering = false; setProgress(''); }
        else        { setProgress('Pass ' + ((4 / m.stride) | 0) + '/3 …'); }
      },
      onUnavailable: fallback,
      onError(e) { console.warn('[schwarz cpu worker]', e); fallback(); },
    });
  }

  function chainPass(token, stride, next) {
    if (token !== sState.renderToken) return;
    setProgress('Pass ' + (4 / stride | 0) + (stride === 1 ? '/3 (full)…' : '/3 (refining)…'));
    const W = sState.fieldW, H = sState.fieldH;
    const sw = sState.schwarz;
    const maxIter = sState.grid.maxIter;
    // Map field coords → world (plane) or z-disk (z-view).
    const inZ = sState.viewMode === 'z';
    const vv = inZ ? sState.zView : sState.view;
    const cssW = vv.cssW, cssH = vv.cssH;
    const cx = vv.cx, cy = vv.cy, scale = vv.scale;
    const pxPerCellX = cssW / W, pxPerCellY = cssH / H;

    let row = 0;
    // Per-row warm-start chain: the converged ψ-seed from the left neighbor
    // (same row, prior col) is reused as initialSeedHint for the current pixel.
    // Adjacent pixels in w-space land on adjacent z-values in 𝔻, so Newton
    // typically converges in 1–3 iters instead of 5–10. Reset at row start.
    let leftSeed = null;
    function chunk() {
      if (token !== sState.renderToken) return;
      const tStart = performance.now();
      while (row < H) {
        leftSeed = null;
        for (let col = 0; col < W; col++) {
          if ((row % stride) !== 0 || (col % stride) !== 0) continue;
          const idx = row * W + col;
          if (sState.fieldKind[idx] && stride > 1) continue;
          const px = (col + 0.5) * pxPerCellX;
          const py = (row + 0.5) * pxPerCellY;
          const aRe = cx + (px - cssW / 2) / scale;
          const aIm = cy - (py - cssH / 2) / scale;
          let wpt;
          if (inZ) {
            // z-disk sample: mask off the domain disk, else lift to w = φ(z).
            const r2 = aRe * aRe + aIm * aIm;
            if (sw.unbounded ? r2 <= 1 : r2 >= 1) {
              sState.field[idx] = 0; sState.fieldKind[idx] = KIND_OUTSIDE + 1; continue;
            }
            wpt = sw.evalPhi({ re: aRe, im: aIm });
            if (!wpt || !isFinite(wpt.re) || !isFinite(wpt.im)) {
              sState.field[idx] = 0; sState.fieldKind[idx] = KIND_OUTSIDE + 1; continue;
            }
          } else {
            wpt = { re: aRe, im: aIm };
          }
          if (!sw.isInOmega(wpt)) {
            sState.field[idx] = 0;
            sState.fieldKind[idx] = KIND_OUTSIDE + 1;
            // leftSeed stays — outside pixels don't update the chain.
          } else {
            const et = QD.Schwarz.escapeTime(wpt, sw, { maxIter, initialSeedHint: leftSeed });
            sState.field[idx] = et.n;
            sState.fieldKind[idx] =
              (et.kind === 'fundamental' ? KIND_FUND :
               et.kind === 'escaped'     ? KIND_ESC  :
               et.kind === 'interior'    ? KIND_INT  :
                                           KIND_INV) + 1;
            // Carry forward only if ψ converged to a usable seed.
            if (et.firstZ) leftSeed = et.firstZ;
          }
        }
        row++;
        if (performance.now() - tStart > 14) {
          requestAnimationFrame(chunk);
          paintAll();
          return;
        }
      }
      // After this pass: fill in any cells skipped by larger stride with the
      // nearest sampled value (for the coarse-display effect).
      fillFromCoarseSamples(stride);
      paintAll();
      next();
    }
    requestAnimationFrame(chunk);
  }

  // Fill un-resolved cells (kind === 0) with their nearest stride-aligned
  // neighbor's value, so the coarse pass shows blocky filled-in content.
  function fillFromCoarseSamples(stride) {
    const W = sState.fieldW, H = sState.fieldH;
    for (let row = 0; row < H; row++) {
      const rAnchor = row - (row % stride);
      for (let col = 0; col < W; col++) {
        const idx = row * W + col;
        if (sState.fieldKind[idx]) continue;
        const cAnchor = col - (col % stride);
        const aIdx = rAnchor * W + cAnchor;
        sState.field[idx]     = sState.field[aIdx];
        sState.fieldKind[idx] = sState.fieldKind[aIdx];
      }
    }
  }


    return { requestRecompute };
  };
})(typeof window !== 'undefined' ? window : globalThis);
