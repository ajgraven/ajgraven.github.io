// =============================================================================
// schwarz-features.js -- Feature-compute methods for the Schwarz tab.
//
// Extracted from schwarz-ui.js by the Phase-3 UI modularization (item E).
// Exposes QD_UI.installSchwarzFeatures(sCtx); schwarz-ui.js captures the
// returned functions into IIFE-local bindings (forward-declared near the top),
// so every card-builder event handler + the interaction-install dep list calls
// them by their original names, unchanged.
//
// These are the per-feature compute/recompute routines wired to the analysis,
// limit-set, and forward-dynamics cards (plus the PNG export action): the
// preimage-tree rebuild + stats, σ domain-coloring, the limit-set chaos game,
// σ level curves, critical orbits, the cycle finder, the orbit-family sweep,
// the z-panel ψ-pullback, and high-res PNG export. Bodies are VERBATIM moves.
//
// Deps via sCtx: sState + the paint fns (paintBoundaryOnTop / paintPreimageTree
// / paintLimitSet, from schwarz-paint.js) + the GPU/canvas helpers used by the
// exporter (activeRenderer / getCtx / getCanvas). `QD` is the global; `document`
// / `performance` / `URL` / `console` are browser globals.
// =============================================================================

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  global.QD_UI.installSchwarzFeatures = function installSchwarzFeatures(s) {
    const sState             = s.sState;
    const paintBoundaryOnTop = s.paintBoundaryOnTop;
    const activeRenderer     = s.activeRenderer;
    const getCtx             = s.getCtx;
    const getCanvas          = s.getCanvas;

  function _recomputeDomainColoring() {
    if (!sState.schwarz) { sState.domainColor = null; return; }
    const v = sState.view;
    const viewport = {
      reMin: v.cx - (v.cssW / 2) / v.scale,
      reMax: v.cx + (v.cssW / 2) / v.scale,
      imMin: v.cy - (v.cssH / 2) / v.scale,
      imMax: v.cy + (v.cssH / 2) / v.scale,
    };
    const W = 256, H = 256;
    try {
      const buf = QD.Schwarz.domainColoringField(sState.schwarz, viewport, { W, H });
      sState.domainColor = { buf, W, H, viewport };
    } catch (_) { sState.domainColor = null; }
  }

  // Re-build the tree from the existing seed (root of the previous tree)
  // when depth or budget changes. Cheap: σ⁻¹ runs are millisecond-scale.
  function _rebuildPreimageTreeIfActive() {
    if (sState.mode !== 'fractal') return;
    if (!sState.schwarz || !sState.preimageTree) return;
    const seed = sState.preimageTree.generations[0][0];
    sState.preimageTree = QD.Schwarz.buildPreimageTree(seed, sState.schwarz, {
      depth:        sState.preimageDepth,
      visualBudget: sState.preimageBudget,
    });
    // paintBoundaryOnTop draws the tree for the active view; no bare w-plane
    // paintPreimageTree() (it would mis-draw onto the z-disk in z-view).
    paintBoundaryOnTop();
    _refreshPreimageTreeStats();
  }

  function _refreshPreimageTreeStats() {
    const el = document.getElementById('schwarz-preimage-count');
    if (!el || !sState.preimageTree) { if (el) el.textContent = ''; return; }
    let total = 0;
    for (const g of sState.preimageTree.generations) total += g.length;
    const trunc = sState.preimageTree.truncatedByBudget ? ' (capped)' : '';
    el.textContent = total + ' pts' + trunc;
  }

  function _computeLimitSet() {
    if (!sState.schwarz) {
      const el = document.getElementById('schwarz-ls-status');
      if (el) el.textContent = 'No φ captured.';
      return;
    }
    const statusEl = document.getElementById('schwarz-ls-status');
    const dimEl    = document.getElementById('schwarz-ls-dim');
    if (statusEl) statusEl.textContent = 'Computing…';
    if (dimEl)    dimEl.textContent    = '';
    // Defer one frame so the "Computing…" text actually paints.
    setTimeout(() => {
      const t0 = performance.now();
      try {
        sState.limitSet = QD.Schwarz.sampleLimitSet(sState.schwarz, {
          n: sState.limitSetN,
          burnIn: 200,
        });
      } catch (err) {
        if (statusEl) statusEl.textContent = 'Error: ' + (err.message || err);
        return;
      }
      const t1 = performance.now();
      const n = sState.limitSet.length / 2;
      if (statusEl) statusEl.textContent = n + ' pts in ' + (t1 - t0).toFixed(0) + ' ms';
      // Dimension estimate.
      if (n >= 200) {
        try {
          const r = QD.Schwarz.boxCountingDimension(sState.limitSet);
          sState.limitSetDim = r.dim;
          if (dimEl) dimEl.textContent = 'dim ≈ ' + (isFinite(r.dim) ? r.dim.toFixed(3) : 'NaN');
        } catch (_) { /* ignore */ }
      }
      // paintBoundaryOnTop draws the limit set for the active view (paintZView
      // in z); no bare w-plane paintLimitSet() (it would mis-draw in the z-view).
      paintBoundaryOnTop();
    }, 30);
  }

  function _clearLimitSet() {
    sState.limitSet = null;
    sState.limitSetDim = null;
    const el = document.getElementById('schwarz-ls-status');
    const dimEl = document.getElementById('schwarz-ls-dim');
    if (el) el.textContent = '';
    if (dimEl) dimEl.textContent = '';
    paintBoundaryOnTop();
  }

  function _recomputeCriticalOrbits() {
    if (!sState.schwarz) { sState.criticalOrbits = null; return; }
    const seeds = QD.Schwarz.canonicalSeeds(sState.schwarz);
    const out = [];
    for (const s of seeds) {
      const orbit = QD.Schwarz.makeOrbit(s.w, sState.schwarz,
                                          { maxIter: sState.grid.maxIter });
      out.push({ label: s.label, orbit });
    }
    sState.criticalOrbits = out;
  }

  function _findCycles() {
    if (!sState.schwarz) return;
    const n = +(document.getElementById('schwarz-cycle-n').value || 2);
    const statusEl = document.getElementById('schwarz-cycle-count');
    if (statusEl) statusEl.textContent = '…';
    setTimeout(() => {
      const t0 = performance.now();
      let cycles = [];
      try {
        cycles = QD.Schwarz.findCycles(sState.schwarz, n, { gridSize: 18 });
      } catch (_) { /* ignore */ }
      const t1 = performance.now();
      sState.cycles = cycles;
      if (statusEl) statusEl.textContent =
        cycles.length + ' cycles in ' + (t1 - t0).toFixed(0) + ' ms';
      paintBoundaryOnTop();
    }, 30);
  }

  // S6 / F8: PNG export. Composites whatever is currently visible on the
  // Schwarz tab into a single PNG download. For GPU mode, that means the
  // fractal layer from #schwarz-gl-canvas + the overlay layer from #canvas.
  // For CPU / domain-coloring modes, only the 2D canvas is needed.
  //
  // High-res: when multiplier > 1, we briefly re-render the GPU canvas at
  // multiplier·display-size, re-paint the 2D overlay onto an off-screen
  // canvas of the same size, composite, export, restore.
  function _exportPng() {
    const mult = +(document.getElementById('schwarz-export-mult').value || 1);
    const view  = sState.view;
    const baseW = Math.round(view.cssW);
    const baseH = Math.round(view.cssH);
    const outW  = baseW * mult;
    const outH  = baseH * mult;

    // Off-screen composite canvas.
    const out    = document.createElement('canvas');
    out.width    = outW;
    out.height   = outH;
    const outCtx = out.getContext('2d');

    // --- 1) Fractal layer ---
    const glCanvas = document.getElementById('schwarz-gl-canvas');
    const onGpu    = activeRenderer() === 'gpu' && glCanvas && sState.gpu;
    if (onGpu && sState.mode === 'fractal') {
      if (mult > 1) {
        // Re-render at higher resolution. We construct a temporary view with
        // larger css dimensions so the renderer chooses larger drawing-buffer.
        const tmpView = Object.assign({}, view, { cssW: outW, cssH: outH });
        try {
          sState.gpu.setColormap(sState.grid.colormap);
          sState.gpu.render(tmpView, {
            maxIter:   sState.grid.maxIter,
            scaleMode: sState.grid.scaleMode,
            modK:      sState.grid.modK,
          });
          outCtx.drawImage(glCanvas, 0, 0, outW, outH);
        } catch (e) {
          console.warn('[export] high-res GPU render failed:', e);
          outCtx.drawImage(glCanvas, 0, 0, outW, outH);
        }
      } else {
        outCtx.drawImage(glCanvas, 0, 0, outW, outH);
      }
    } else {
      // CPU / domain-coloring: the 2D canvas already has the fractal layer
      // (or there isn't one). Nothing extra here.
      outCtx.fillStyle = '#fafafa';
      outCtx.fillRect(0, 0, outW, outH);
    }

    // --- 2) 2D overlay (boundary, orbits, markers, z-panel, etc.) ---
    const ctx2d = getCtx();
    if (ctx2d) {
      const mainCanvas = getCanvas();
      // Render the 2D layer at the target resolution. Easiest: scale the
      // existing canvas with imageSmoothingEnabled = false. Boundary lines
      // will be 1-px regardless of multiplier — acceptable for typical
      // print/share use; pure-vector boundaries would need a re-render.
      outCtx.imageSmoothingEnabled = false;
      outCtx.drawImage(mainCanvas, 0, 0, outW, outH);
    }

    // --- 3) Restore GPU canvas to its display size ---
    if (onGpu && mult > 1) {
      try {
        sState.gpu.render(view, {
          maxIter:   sState.grid.maxIter,
          scaleMode: sState.grid.scaleMode,
          modK:      sState.grid.modK,
        });
      } catch (_) { /* ignore */ }
    }

    // --- 4) Download ---
    out.toBlob((blob) => {
      if (!blob) return;
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      const ts = new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19);
      a.download = `qd-schwarz-${ts}-${outW}x${outH}.png`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    }, 'image/png');
  }

  // F4: ψ-pullback of the current w-orbit. Each w in sState.orbit is run
  // through sw.psi to get the corresponding z in 𝔻 (or 𝔻*); used to render
  // the z-history in the z-disk view (paintZOverlays).
  function _recomputeZPanelOrbit() {
    if (!sState.schwarz || !sState.orbit || sState.orbit.length === 0) {
      sState.zPanelOrbit = null;
      return;
    }
    const out = [];
    for (const w of sState.orbit) {
      let z;
      try { z = sState.schwarz.psi(w); } catch (_) { z = null; }
      if (z && isFinite(z.re) && isFinite(z.im)) out.push(z);
      else out.push(null);
    }
    sState.zPanelOrbit = out;
  }

  function _computeSweep() {
    if (!sState.schwarz) return;
    const N     = +(document.getElementById('schwarz-sweep-n').value     || 16);
    const depth = +(document.getElementById('schwarz-sweep-depth').value || 12);
    // Default sweep: horizontal line across the boundary bbox at y = centroid.
    const bdy = sState.schwarz._boundaryPts || [];
    let minRe = -1, maxRe = 1, cy = 0;
    if (bdy.length > 0) {
      minRe = Infinity; maxRe = -Infinity;
      let cyAcc = 0;
      for (const p of bdy) {
        if (p.re < minRe) minRe = p.re;
        if (p.re > maxRe) maxRe = p.re;
        cyAcc += p.im;
      }
      cy = cyAcc / bdy.length;
      const dx = maxRe - minRe;
      minRe += 0.1 * dx; maxRe -= 0.1 * dx;
    }
    const seeds = QD.Schwarz.sampleSweepSeeds('line',
      { from: { re: minRe, im: cy }, to: { re: maxRe, im: cy }, n: N });
    const out = [];
    for (const seed of seeds) {
      if (!sState.schwarz.isInOmega(seed)) { out.push([]); continue; }
      const orb = QD.Schwarz.makeOrbit(seed, sState.schwarz, { maxIter: depth });
      out.push(orb);
    }
    sState.sweepOrbits = out;
    paintBoundaryOnTop();
  }

  // Compute level curves on the current viewport. Triggered on toggle-on
  // and on pan/zoom (so contours follow the view).
  function _recomputeLevelCurves() {
    if (!sState.schwarz) { sState.levelCurves = null; return; }
    const v = sState.view;
    const viewport = {
      reMin: v.cx - (v.cssW / 2) / v.scale,
      reMax: v.cx + (v.cssW / 2) / v.scale,
      imMin: v.cy - (v.cssH / 2) / v.scale,
      imMax: v.cy + (v.cssH / 2) / v.scale,
    };
    try {
      sState.levelCurves = QD.Schwarz.computeSigmaLevelCurves(sState.schwarz,
        { gridSize: 96, viewport });
    } catch (_) { sState.levelCurves = null; }
  }

    return {
      _recomputeDomainColoring, _rebuildPreimageTreeIfActive,
      _refreshPreimageTreeStats, _computeLimitSet, _clearLimitSet,
      _recomputeCriticalOrbits, _findCycles, _exportPng,
      _recomputeZPanelOrbit, _computeSweep, _recomputeLevelCurves,
    };
  };
})(typeof window !== 'undefined' ? window : globalThis);
