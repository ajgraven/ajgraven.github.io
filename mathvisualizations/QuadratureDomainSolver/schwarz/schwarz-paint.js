// =============================================================================
// schwarz-paint.js -- Rendering-output (paint) layer for the Schwarz tab.
//
// Extracted from schwarz-ui.js by the Phase-3 UI modularization (item E).
// Exposes QD_UI.installSchwarzPaint(sCtx); schwarz-ui.js captures the public
// paint entry points (clearCanvas / paintAll / paintBoundaryOnTop /
// repaintField / paintOrbit / paintPreimageTree / paintLimitSet / setProgress /
// colormap) into IIFE-local bindings, so every call site is unchanged.
//
// This is the 2D-canvas output layer: the escape-time field, the boundary, the
// forward / preimage / limit-set / analysis overlays, and the colormaps. It
// reads sState + a few geometry/canvas helpers (getCtx / syncCanvasSize /
// worldToPixel) and the KIND_* pixel-class constants via sCtx; the colormap
// utilities (cpuComputeT / colormap / interpStops / CMAP) are pure and live
// here. The bodies are VERBATIM moves. `QD` is the solver global.
// =============================================================================

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  global.QD_UI.installSchwarzPaint = function installSchwarzPaint(s) {
    const sState          = s.sState;
    const getCtx          = s.getCtx;
    const syncCanvasSize  = s.syncCanvasSize;
    const worldToPixel    = s.worldToPixel;
    const zToPixel        = s.zToPixel;       // z-disk → pixel (z-view transform)
    const activeRenderer  = s.activeRenderer; // 'gpu' | 'cpu' (z field source)
    const KIND_FUND       = s.KIND_FUND;
    const KIND_ESC        = s.KIND_ESC;
    const KIND_INT        = s.KIND_INT;
    const KIND_INV        = s.KIND_INV;
    const KIND_OUTSIDE    = s.KIND_OUTSIDE;

  function clearCanvas() {
    const ctx = getCtx(); if (!ctx) return;
    syncCanvasSize();
    ctx.clearRect(0, 0, sState.view.cssW, sState.view.cssH);
    ctx.fillStyle = '#fafafa';
    ctx.fillRect(0, 0, sState.view.cssW, sState.view.cssH);
    ctx.fillStyle = '#777';
    ctx.font = '13px system-ui, sans-serif';
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.fillText('Solve on the Inverse tab and click “Use this φ” to begin.',
                 sState.view.cssW/2, sState.view.cssH/2);
  }

  function paintAll() {
    const ctx = getCtx(); if (!ctx) return;
    syncCanvasSize();
    // z-plane view is a self-contained render (tiling in 𝔻 + pulled-back
    // overlays); none of the w-side painters below run in z-mode.
    if (sState.viewMode === 'z') { paintZView(); return; }
    // Clear up front so fieldless overlay modes don't ghost. paintField()
    // clears + blits the escape field, but it returns early when there is no
    // field (e.g. domain-coloring, or a bare boundary/limit-set/level-curve
    // view). Without this, panning re-blits the cached domain-coloring image at
    // a shifted transform over the un-cleared previous frame → offset ghosts.
    // (paintBoundaryOnTop already clears for the GPU overlay path.)
    ctx.clearRect(0, 0, sState.view.cssW, sState.view.cssH);
    paintField();
    paintDomainColoring();               // S5 / F6: under the rest
    paintSigmaLevelCurves();             // S4 / F12: under the boundary
    paintBoundary();
    paintSweepOrbits();                  // S5 / H8
    paintCurveImage();                   // S5 / E11
    paintCriticalOrbits();               // S5 / H7
    paintOrbit();
    paintPreimageTree();
    paintCycles();                       // S5 / E10
    paintLimitSet();
    paintSigmaSingularities();           // S4 / F3: on top so markers are visible
  }
  function repaintField() { if (sState.field) paintAll(); }

  // Used after a GPU render: WebGL has already drawn the fractal to the
  // sibling #schwarz-gl-canvas. We clear the main 2D canvas to transparent
  // and draw only the boundary + orbit overlays on top. S1 adds the
  // preimage tree to the overlay chain.
  function paintBoundaryOnTop() {
    const ctx = getCtx(); if (!ctx) return;
    if (sState.viewMode === 'z') { paintZView(); return; }
    ctx.clearRect(0, 0, sState.view.cssW, sState.view.cssH);
    paintDomainColoring();               // S5 / F6
    paintSigmaLevelCurves();             // S4: contours under the boundary
    paintBoundary();
    paintSweepOrbits();                  // S5 / H8
    paintCurveImage();                   // S5 / E11
    paintCriticalOrbits();               // S5 / H7
    paintOrbit();                        // w-side orbit
    paintPreimageTree();
    paintCycles();                       // S5 / E10
    paintLimitSet();
    paintSigmaSingularities();           // S4: markers on top
  }

  // Cached off-screen canvas + ImageData buffer for CPU repaint. Re-created
  // only when (W, H) change — avoids allocating a few MB on every paint
  // during the progressive pyramid passes.
  let offC = null, offCtx = null, offImg = null;
  function ensureOffscreen(W, H) {
    if (offC && offC.width === W && offC.height === H) return;
    offC = document.createElement('canvas');
    offC.width = W; offC.height = H;
    offCtx = offC.getContext('2d');
    offImg = offCtx.createImageData(W, H);
  }

  function paintField() {
    const ctx = getCtx();
    const W = sState.fieldW, H = sState.fieldH;
    if (!sState.field || !W || !H) return;
    ensureOffscreen(W, H);
    const imgData = offImg;
    const maxIter = sState.grid.maxIter;
    const cmap = sState.grid.colormap;
    // In z-view, KIND_OUTSIDE marks pixels off the unit disk (not part of the
    // picture) → neutral backdrop; in plane view it's Ω^c, the fundamental tile.
    const zMode = sState.viewMode === 'z';
    for (let i = 0; i < W * H; i++) {
      const kind = sState.fieldKind[i];
      const n    = sState.field[i];
      let r = 0, g = 0, b = 0;
      if (kind === KIND_OUTSIDE + 1)        { if (zMode) { r = 224; g = 226; b = 232; } else { r = 245; g = 245; b = 248; } }   // z: off-disk / plane: fundamental tile
      else if (kind === KIND_INT + 1)       { r = 28;  g = 28;  b = 36;  }     // interior (tiling-set limit)
      else if (kind === KIND_ESC + 1)       { r = 80;  g = 80;  b = 90;  }     // escaping set
      else if (kind === KIND_INV + 1)       { r = 180; g = 90;  b = 90;  }     // bad pixel
      else if (kind === KIND_FUND + 1) {
        const t = cpuComputeT(n, maxIter, sState.grid.scaleMode, sState.grid.modK);
        const c = colormap(cmap, t);
        r = c[0]; g = c[1]; b = c[2];
      }
      const j = i * 4;
      imgData.data[j]   = r;
      imgData.data[j+1] = g;
      imgData.data[j+2] = b;
      imgData.data[j+3] = 255;
    }
    offCtx.putImageData(imgData, 0, 0);
    ctx.imageSmoothingEnabled = false;
    ctx.clearRect(0, 0, sState.view.cssW, sState.view.cssH);
    ctx.drawImage(offC, 0, 0, sState.view.cssW, sState.view.cssH);
  }

  function paintBoundary() {
    const pts = sState.boundarySnapshot;
    if (!pts || !pts.length) return;
    const ctx = getCtx();
    ctx.strokeStyle = '#1a3e7a';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    const p0 = worldToPixel(pts[0].re, pts[0].im);
    ctx.moveTo(p0.x, p0.y);
    for (let i = 1; i < pts.length; i++) {
      const p = worldToPixel(pts[i].re, pts[i].im);
      ctx.lineTo(p.x, p.y);
    }
    ctx.closePath();
    ctx.stroke();
  }

  // Draw one forward-orbit polyline + dots in the given style. `style.dash`
  // (optional) applies a dashed connecting line; it is reset afterwards so
  // other overlays aren't affected.
  function drawOrbitPolyline(pts, style) {
    if (!pts || pts.length === 0) return;
    const ctx = getCtx(); if (!ctx) return;
    ctx.strokeStyle = style.line;
    ctx.lineWidth   = style.lineWidth;
    if (ctx.setLineDash) ctx.setLineDash(style.dash || []);
    ctx.beginPath();
    for (let i = 0; i < pts.length; i++) {
      const p = worldToPixel(pts[i].re, pts[i].im);
      if (i === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
    }
    ctx.stroke();
    if (ctx.setLineDash) ctx.setLineDash([]);
    for (let i = 0; i < pts.length; i++) {
      const p = worldToPixel(pts[i].re, pts[i].im);
      ctx.beginPath();
      ctx.arc(p.x, p.y, i === 0 ? style.seedR : style.dotR, 0, 2 * Math.PI);
      ctx.fillStyle = i === 0 ? style.seedFill : style.dotFill;
      ctx.fill();
      if (i === 0 && style.seedHalo) {
        ctx.strokeStyle = '#fff';
        ctx.lineWidth = 1.2;
        ctx.stroke();
      }
    }
  }

  // Draw the transient hover orbit (light, dashed, underneath) then the
  // pinned orbit (solid green, on top). The hover orbit is the live preview;
  // the pinned orbit is the single-click-committed one (= sState.orbit).
  function paintOrbit() {
    drawOrbitPolyline(sState.hoverOrbit, {
      line: 'rgba(20, 160, 60, 0.45)', lineWidth: 1.0, dash: [3, 3],
      seedR: 3.2, dotR: 2.0, seedFill: 'rgba(16,138,64,0.7)',
      dotFill: 'rgba(20,160,60,0.45)', seedHalo: false,
    });
    drawOrbitPolyline(sState.pinnedOrbit, {
      line: 'rgba(20, 160, 60, 0.95)', lineWidth: 1.3,
      seedR: 4.5, dotR: 2.8, seedFill: '#108a40',
      dotFill: 'rgba(20, 160, 60, 0.85)', seedHalo: true,
    });
  }

  // ---------------------------------------------------------------------------
  // S1: paint the preimage tree (TODO #16 visualization).
  //
  // Layout: generations colored along a plasma-like ramp. Gen-0 (the seed)
  // is bright yellow; later generations interpolate toward magenta. Edges
  // are translucent thin lines parent→child. Dot radius shrinks slightly
  // with depth (5 → 1.5 px) so deeper generations don't dominate visually.
  // ---------------------------------------------------------------------------
  function paintPreimageTree() {
    const tree = sState.preimageTree;
    if (!tree || !tree.generations || tree.generations.length === 0) return;
    const ctx = getCtx(); if (!ctx) return;
    const N = tree.generations.length;

    // Plasma-like ramp (gen 0 yellow → gen N-1 deep purple), shared with the
    // z-panel mirror via the module-level preimageGenColor.
    const genColor = (g) => preimageGenColor(g, N);

    // Edges first (so dots sit on top).
    ctx.lineWidth = 0.9;
    for (const e of tree.edges) {
      const p0 = worldToPixel(tree.generations[e.fromGen][e.fromIdx].re,
                              tree.generations[e.fromGen][e.fromIdx].im);
      const p1 = worldToPixel(tree.generations[e.toGen  ][e.toIdx  ].re,
                              tree.generations[e.toGen  ][e.toIdx  ].im);
      ctx.strokeStyle = genColor(e.toGen);
      ctx.beginPath();
      ctx.moveTo(p0.x, p0.y);
      ctx.lineTo(p1.x, p1.y);
      ctx.stroke();
    }

    // Dots — gen 0 (seed) drawn with a white halo so it's easy to spot.
    for (let g = 0; g < N; g++) {
      const r = Math.max(1.5, 5 - g * 0.55);
      ctx.fillStyle = genColor(g);
      for (const w of tree.generations[g]) {
        const p = worldToPixel(w.re, w.im);
        ctx.beginPath();
        ctx.arc(p.x, p.y, r, 0, 2 * Math.PI);
        ctx.fill();
        if (g === 0) {
          ctx.strokeStyle = '#fff';
          ctx.lineWidth = 1.5;
          ctx.stroke();
        }
      }
    }
  }

  // S4 / F3: paint σ-poles + branch points overlay.
  function paintSigmaSingularities() {
    if (!sState.showSingularities || !sState.sigmaSingularities) return;
    const { poles, branchPoints } = sState.sigmaSingularities;
    const ctx = getCtx();
    // σ-poles: red filled circles with a thin label.
    ctx.font = '10px ui-monospace, Consolas, monospace';
    ctx.textBaseline = 'middle';
    ctx.textAlign    = 'left';
    for (const pl of poles) {
      const p = worldToPixel(pl.w.re, pl.w.im);
      if (!isFinite(p.x) || !isFinite(p.y)) continue;
      ctx.beginPath();
      ctx.arc(p.x, p.y, 6, 0, 2 * Math.PI);
      ctx.fillStyle   = '#d12d2d';
      ctx.fill();
      ctx.strokeStyle = '#fff';
      ctx.lineWidth   = 1.4;
      ctx.stroke();
      ctx.fillStyle = '#d12d2d';
      ctx.fillText(pl.label + ' (' + pl.kind + ')', p.x + 8, p.y);
    }
    // Branch points: blue triangle pointing up.
    for (const bp of branchPoints) {
      const p = worldToPixel(bp.w.re, bp.w.im);
      if (!isFinite(p.x) || !isFinite(p.y)) continue;
      ctx.beginPath();
      ctx.moveTo(p.x,     p.y - 6);
      ctx.lineTo(p.x + 6, p.y + 4);
      ctx.lineTo(p.x - 6, p.y + 4);
      ctx.closePath();
      ctx.fillStyle = bp.severity === 'critical' ? '#3a54ae' : '#1a3e7a';
      ctx.fill();
      ctx.strokeStyle = '#fff';
      ctx.lineWidth   = 1.2;
      ctx.stroke();
      ctx.fillStyle = '#1a3e7a';
      ctx.fillText(bp.label, p.x + 8, p.y + 4);
    }
  }

  // S4 / F12: paint |σ| (solid) + arg(σ) (dashed) level curves.
  function paintSigmaLevelCurves() {
    if (!sState.showLevelCurves || !sState.levelCurves) return;
    const ctx = getCtx();
    const lc = sState.levelCurves;
    // |σ| solid contours — muted teal.
    ctx.lineWidth   = 1.0;
    ctx.strokeStyle = 'rgba(20, 130, 130, 0.7)';
    if (ctx.setLineDash) ctx.setLineDash([]);
    for (const c of lc.abs) {
      ctx.beginPath();
      for (const s of c.segments) {
        const a = worldToPixel(s.x0, s.y0);
        const b = worldToPixel(s.x1, s.y1);
        ctx.moveTo(a.x, a.y);
        ctx.lineTo(b.x, b.y);
      }
      ctx.stroke();
    }
    // arg(σ) dashed contours — muted magenta.
    ctx.lineWidth   = 0.9;
    ctx.strokeStyle = 'rgba(170, 60, 130, 0.65)';
    if (ctx.setLineDash) ctx.setLineDash([4, 3]);
    for (const c of lc.arg) {
      ctx.beginPath();
      for (const s of c.segments) {
        const a = worldToPixel(s.x0, s.y0);
        const b = worldToPixel(s.x1, s.y1);
        ctx.moveTo(a.x, a.y);
        ctx.lineTo(b.x, b.y);
      }
      ctx.stroke();
    }
    if (ctx.setLineDash) ctx.setLineDash([]);
  }

  // Shared preimage-tree generation color ramp (plasma-like: gen 0 yellow →
  // last gen deep purple). Used by paintPreimageTree (w-plane) and the z-panel
  // mirror so both read identically.
  function preimageGenColor(g, N) {
    const t = N > 1 ? (g / (N - 1)) : 0;
    if (t < 0.5) {
      const u = t * 2;
      return `rgba(${Math.round(253 + (240 - 253) * u)},`
           + `${Math.round(231 + (132 - 231) * u)},`
           + `${Math.round(37  + (74  - 37 ) * u)},0.92)`;
    }
    const u = (t - 0.5) * 2;
    return `rgba(${Math.round(240 + (93  - 240) * u)},`
         + `${Math.round(132 + (1   - 132) * u)},`
         + `${Math.round(74  + (166 - 74 ) * u)},0.92)`;
  }

  // ---------------------------------------------------------------------------
  // z-view overlay mirroring (S6 / F4). Every Schwarz overlay is stored in the
  // w-plane; to show it in the z-disk view each w-point is pulled back through
  // ψ (sState.schwarz.psi — a per-point Newton solve) into 𝔻 (or 𝔻*). ψ is
  // VIEW-INDEPENDENT, so a pullback only changes when the source overlay (or the
  // captured φ) changes — never on pan/zoom. We therefore cache each pullback
  // keyed by (φ handle, source-array reference): the z-view repaints every frame
  // during a pan but reuses the cache, and a pullback re-runs only when its
  // overlay is recomputed/cleared. Points with no preimage in the target domain
  // (ψ → null) become gaps, exactly as the orbit already handles.
  // ---------------------------------------------------------------------------
  const LIMITSET_Z_CAP = 2000;      // subsample the cloud for the 180px inset
  const LEVELCURVE_Z_SEG_CAP = 800; // cap pulled-back segments per contour set
  const _zc = {
    phi: null, tree: null, limit: null, cycles: null,
    sweep: null, curve: null, critical: null, sing: null, lc: null,
  };
  function _zInvalidateIfPhiChanged() {
    if (_zc.phi !== sState.schwarz) {
      _zc.phi = sState.schwarz;
      for (const k in _zc) if (k !== 'phi') _zc[k] = null;
    }
  }
  function _psiOrNull(w) {
    if (!w || !sState.schwarz) return null;
    let z;
    try { z = sState.schwarz.psi(w); } catch (_) { return null; }
    return (z && isFinite(z.re) && isFinite(z.im)) ? z : null;
  }
  function _psiArr(pts) {
    const out = new Array(pts.length);
    for (let i = 0; i < pts.length; i++) out[i] = _psiOrNull(pts[i]);
    return out;
  }
  // Recompute compute(src) only when the source reference changed.
  function _zMemo(slot, src, compute) {
    const c = _zc[slot];
    if (c && c.src === src) return c.out;
    const out = src ? compute(src) : null;
    _zc[slot] = { src, out };
    return out;
  }
  const _zTree = () => _zMemo('tree', sState.preimageTree, (t) =>
    (t.generations ? { generations: t.generations.map(_psiArr), edges: t.edges } : null));
  const _zLimit = () => _zMemo('limit', sState.limitSet, (arr) => {
    const nPts = arr.length >> 1;
    const stride = Math.max(1, Math.ceil(nPts / LIMITSET_Z_CAP));
    const out = [];
    for (let i = 0; i < nPts; i += stride) {
      out.push(_psiOrNull({ re: arr[2 * i], im: arr[2 * i + 1] }));
    }
    return out;
  });
  const _zCycles = () => _zMemo('cycles', sState.cycles, (cs) =>
    cs.map(c => ({ period: c.period, z: _psiArr(c.points) })));
  const _zSweep = () => _zMemo('sweep', sState.sweepOrbits, (sw) => sw.map(_psiArr));
  const _zCurve = () => _zMemo('curve', sState.curveImage, (ci) => ci.map(_psiArr));
  const _zCritical = () => _zMemo('critical', sState.criticalOrbits, (cos) =>
    cos.map(o => ({ label: o.label, z: _psiArr(o.orbit) })));
  const _zSing = () => _zMemo('sing', sState.sigmaSingularities, (s) => ({
    poles:        (s.poles || []).map(p => _psiOrNull(p.w)),
    branchPoints: (s.branchPoints || []).map(b => _psiOrNull(b.w)),
  }));
  const _zLevelCurves = () => _zMemo('lc', sState.levelCurves, (lc) => {
    const pull = (contours) => {
      const segs = [];
      let budget = LEVELCURVE_Z_SEG_CAP;
      for (const c of contours) {
        for (const s of c.segments) {
          if (budget-- <= 0) return segs;
          const a = _psiOrNull({ re: s.x0, im: s.y0 });
          const b = _psiOrNull({ re: s.x1, im: s.y1 });
          if (a && b) segs.push([a, b]);
        }
      }
      return segs;
    };
    return { abs: pull(lc.abs || []), arg: pull(lc.arg || []) };
  });

  // Draw a pulled-back z-polyline / dot cloud. `zToPanel` is the z→pixel
  // transform supplied by paintZView (it takes a {re,im} object); null entries
  // break the line (gaps) and are skipped as dots.
  function _zDrawPolyline(ctx, zToPanel, pts, stroke, lw, dash) {
    if (!pts) return;
    ctx.strokeStyle = stroke; ctx.lineWidth = lw;
    if (ctx.setLineDash) ctx.setLineDash(dash || []);
    ctx.beginPath();
    let first = true;
    for (const z of pts) {
      if (!z) { first = true; continue; }
      const p = zToPanel(z);
      if (first) { ctx.moveTo(p.x, p.y); first = false; }
      else        ctx.lineTo(p.x, p.y);
    }
    ctx.stroke();
    if (ctx.setLineDash) ctx.setLineDash([]);
  }
  function _zDrawDots(ctx, zToPanel, pts, r, fill) {
    if (!pts) return;
    ctx.fillStyle = fill;
    for (const z of pts) {
      if (!z) continue;
      const p = zToPanel(z);
      ctx.beginPath(); ctx.arc(p.x, p.y, r, 0, 2 * Math.PI); ctx.fill();
    }
  }

  // ---------------------------------------------------------------------------
  // z-plane view (S6 / F4): the Schwarz tiling uniformized onto the unit disk 𝔻
  // (or its exterior 𝔻* for unbounded Ω), drawn as the MAIN plot. The escape-
  // time field — computed over the zView transform with w = φ(z) (see
  // schwarz-render.js) — comes from one of two sources:
  //   • CPU path: sState.field (blitted here via paintField), drawn on this 2D
  //     canvas — paint the backdrop + field.
  //   • GPU path: the WebGL layer below (u_viewMode=1 shader) already drew the
  //     field. We must skip the backdrop + paintField and leave this 2D canvas
  //     transparent so the GL field shows through, drawing only ∂𝔻 + axes +
  //     the ψ-pulled-back overlays on top.
  // The source is derived from activeRenderer() (NOT a caller-passed flag) so
  // that EVERY repaint path — doRecompute, renderImmediate, paintAll, hover-
  // orbit, resize — does the right thing. Otherwise a stray full paintZView()
  // (e.g. from a hover repaint) would fill the opaque backdrop over the GL
  // field, making the GPU render vanish except while actively dragging (where
  // renderImmediate keeps re-running the overlay-only paint). `overlayOnly`
  // remains as an explicit override for safety.
  // ---------------------------------------------------------------------------
  function paintZView(overlayOnly) {
    const ctx = getCtx(); if (!ctx) return;
    syncCanvasSize();
    const cssW = sState.zView.cssW, cssH = sState.zView.cssH;
    const glField = overlayOnly || (activeRenderer && activeRenderer() === 'gpu');
    ctx.clearRect(0, 0, cssW, cssH);
    if (!glField) {
      ctx.fillStyle = '#f3f4f7';               // backdrop shown until the field arrives
      ctx.fillRect(0, 0, cssW, cssH);
      paintField();                            // blits the z escape-time field if present
    }
    // Unit circle ∂𝔻 + faint axes through z = 0, in the live z-view transform.
    const o = zToPixel(0, 0);
    const Rpx = sState.zView.scale;            // 1 z-unit = scale px
    ctx.save();
    ctx.strokeStyle = '#d5d8e0'; ctx.lineWidth = 0.8;
    ctx.beginPath();
    ctx.moveTo(0, o.y); ctx.lineTo(cssW, o.y);
    ctx.moveTo(o.x, 0); ctx.lineTo(o.x, cssH);
    ctx.stroke();
    ctx.strokeStyle = '#1a3e7a'; ctx.lineWidth = 1.6;
    ctx.beginPath(); ctx.arc(o.x, o.y, Rpx, 0, 2 * Math.PI); ctx.stroke();
    // paintZOverlays calls its transform as fn(zObject) (the zToPanel(z)
    // contract); zToPixel takes (re, im), so adapt — passing zToPixel directly
    // made every overlay point map to (NaN, NaN) and silently vanish.
    paintZOverlays(ctx, (z) => zToPixel(z.re, z.im));
    ctx.restore();
  }

  // Draw the ψ-pulled-back overlays into the z-view via the supplied z→pixel
  // transform `zToPanel`. Each overlay is gated by the SAME condition that shows
  // it in the w-plane and uses its w-side palette; pullbacks are cached by
  // source identity (see _zc) so this stays cheap on repaint.
  function paintZOverlays(ctx, zToPanel) {
    _zInvalidateIfPhiChanged();

    // S4 / F12 level curves (under everything): |σ| teal solid, arg(σ) magenta dashed.
    if (sState.showLevelCurves && sState.levelCurves) {
      const zlc = _zLevelCurves();
      if (zlc) {
        ctx.lineWidth = 0.8; ctx.strokeStyle = 'rgba(20, 130, 130, 0.6)';
        if (ctx.setLineDash) ctx.setLineDash([]);
        ctx.beginPath();
        for (const [a, b] of zlc.abs) {
          const pa = zToPanel(a), pb = zToPanel(b);
          ctx.moveTo(pa.x, pa.y); ctx.lineTo(pb.x, pb.y);
        }
        ctx.stroke();
        ctx.strokeStyle = 'rgba(170, 60, 130, 0.55)';
        if (ctx.setLineDash) ctx.setLineDash([3, 2]);
        ctx.beginPath();
        for (const [a, b] of zlc.arg) {
          const pa = zToPanel(a), pb = zToPanel(b);
          ctx.moveTo(pa.x, pa.y); ctx.lineTo(pb.x, pb.y);
        }
        ctx.stroke();
        if (ctx.setLineDash) ctx.setLineDash([]);
      }
    }

    // S5 / H8 sweep family (hue ramp by seed position).
    if (sState.sweepOrbits && sState.sweepOrbits.length) {
      const zs = _zSweep();
      for (let i = 0; i < zs.length; i++) {
        const t = zs.length > 1 ? i / (zs.length - 1) : 0;
        _zDrawPolyline(ctx, zToPanel, zs[i], `hsla(${Math.round(360 * t)}, 65%, 45%, 0.6)`, 0.8);
      }
    }

    // S5 / E11 curve forward-image (blue → red ramp by iteration).
    if (sState.curveImage) {
      const zci = _zCurve();
      const K = zci.length;
      for (let it = 0; it < K; it++) {
        const t = K > 1 ? it / (K - 1) : 0;
        const r = Math.round(40 + (180 - 40) * t);
        const g = Math.round(80 + (50 - 80) * t);
        const b = Math.round(180 + (40 - 180) * t);
        _zDrawPolyline(ctx, zToPanel, zci[it], `rgba(${r},${g},${b},${0.85 - 0.5 * t})`, it === 0 ? 1.4 : 0.9);
      }
    }

    // S5 / H7 critical orbits.
    if (sState.showCriticalOrbits && sState.criticalOrbits) {
      const palette = ['#d4570b', '#7e2d8c', '#107a40', '#b8860b'];
      const zco = _zCritical();
      for (let i = 0; i < zco.length; i++) {
        const col = palette[i % palette.length];
        _zDrawPolyline(ctx, zToPanel, zco[i].z, col, 1.1);
        _zDrawDots(ctx, zToPanel, zco[i].z.slice(0, 1), 3, col);
      }
    }

    // S1 preimage tree (generation ramp; edges then dots).
    if (sState.preimageTree && sState.preimageTree.generations) {
      const zt = _zTree();
      if (zt) {
        const N = zt.generations.length;
        ctx.lineWidth = 0.7;
        for (const e of zt.edges) {
          const a = zt.generations[e.fromGen][e.fromIdx];
          const b = zt.generations[e.toGen][e.toIdx];
          if (!a || !b) continue;
          const pa = zToPanel(a), pb = zToPanel(b);
          ctx.strokeStyle = preimageGenColor(e.toGen, N);
          ctx.beginPath(); ctx.moveTo(pa.x, pa.y); ctx.lineTo(pb.x, pb.y); ctx.stroke();
        }
        for (let g = 0; g < N; g++) {
          _zDrawDots(ctx, zToPanel, zt.generations[g], Math.max(1, 3.2 - g * 0.4), preimageGenColor(g, N));
        }
      }
    }

    // S3 limit set (subsampled pale-yellow cloud).
    if (sState.limitSet && sState.limitSet.length) {
      const zl = _zLimit();
      if (zl) {
        ctx.fillStyle = 'rgba(255, 240, 80, 0.75)';
        for (const z of zl) {
          if (!z) continue;
          const p = zToPanel(z);
          ctx.fillRect(p.x - 0.5, p.y - 0.5, 1, 1);
        }
      }
    }

    // S5 / E10 cycles (dashed closed polyline + dots).
    if (sState.cycles && sState.cycles.length) {
      const palette = ['#108a40', '#b53030', '#5677a8', '#b8860b', '#7e2d8c'];
      const zcy = _zCycles();
      for (let i = 0; i < zcy.length; i++) {
        const col = palette[i % palette.length];
        _zDrawPolyline(ctx, zToPanel, zcy[i].z, col, 1.0, [3, 2]);
        _zDrawDots(ctx, zToPanel, zcy[i].z, 2.5, col);
      }
    }

    // z-orbit polyline + dots, matching the w-orbit color palette (green).
    const zOrbit = sState.zPanelOrbit;
    if (zOrbit && zOrbit.length > 0) {
      // Connecting line.
      ctx.strokeStyle = 'rgba(20, 160, 60, 0.95)';
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      let first = true;
      for (const z of zOrbit) {
        if (!z) { first = true; continue; }
        const p = zToPanel(z);
        if (first) { ctx.moveTo(p.x, p.y); first = false; }
        else        ctx.lineTo(p.x, p.y);
      }
      ctx.stroke();
      // Dots — seed (first) gets a halo.
      for (let i = 0; i < zOrbit.length; i++) {
        const z = zOrbit[i];
        if (!z) continue;
        const p = zToPanel(z);
        ctx.beginPath();
        ctx.arc(p.x, p.y, i === 0 ? 3.5 : 2.0, 0, 2 * Math.PI);
        ctx.fillStyle = i === 0 ? '#108a40' : 'rgba(20, 160, 60, 0.85)';
        ctx.fill();
        if (i === 0) {
          ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.0;
          ctx.stroke();
        }
      }
    }

    // S4 / F3 σ-singularities on top: poles (red) + branch points (blue).
    // Markers only (no labels, to keep the z-view uncluttered); w-side hues.
    if (sState.showSingularities && sState.sigmaSingularities) {
      const zsg = _zSing();
      if (zsg) {
        _zDrawDots(ctx, zToPanel, zsg.poles, 3, '#d12d2d');
        _zDrawDots(ctx, zToPanel, zsg.branchPoints, 3, '#1a3e7a');
      }
    }
  }

  // S5 / F6: paint domain-coloring field into the canvas. Caches an
  // offscreen W×H ImageData and stretches to the current viewport.
  let _dcOffCanvas = null, _dcOffCtx = null;
  function paintDomainColoring() {
    if (sState.mode !== 'domain-coloring' || !sState.domainColor) return;
    const dc = sState.domainColor;
    if (!_dcOffCanvas || _dcOffCanvas.width !== dc.W || _dcOffCanvas.height !== dc.H) {
      _dcOffCanvas = document.createElement('canvas');
      _dcOffCanvas.width  = dc.W;
      _dcOffCanvas.height = dc.H;
      _dcOffCtx = _dcOffCanvas.getContext('2d');
    }
    const img = _dcOffCtx.createImageData(dc.W, dc.H);
    img.data.set(dc.buf);
    _dcOffCtx.putImageData(img, 0, 0);
    const ctx = getCtx();
    // Map source viewport → screen pixel rect.
    const a = worldToPixel(dc.viewport.reMin, dc.viewport.imMax);
    const b = worldToPixel(dc.viewport.reMax, dc.viewport.imMin);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(_dcOffCanvas, a.x, a.y, b.x - a.x, b.y - a.y);
  }

  // S5 / H7: paint critical orbits (canonical-point σ-orbits).
  function paintCriticalOrbits() {
    if (!sState.showCriticalOrbits || !sState.criticalOrbits) return;
    const ctx = getCtx();
    const palette = ['#d4570b', '#7e2d8c', '#107a40', '#b8860b'];
    for (let i = 0; i < sState.criticalOrbits.length; i++) {
      const { label, orbit } = sState.criticalOrbits[i];
      if (!orbit || orbit.length === 0) continue;
      const col = palette[i % palette.length];
      ctx.strokeStyle = col; ctx.lineWidth = 1.4;
      ctx.beginPath();
      for (let k = 0; k < orbit.length; k++) {
        const p = worldToPixel(orbit[k].re, orbit[k].im);
        if (k === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
      }
      ctx.stroke();
      // Seed dot + label.
      const p0 = worldToPixel(orbit[0].re, orbit[0].im);
      ctx.beginPath(); ctx.arc(p0.x, p0.y, 4.5, 0, 2*Math.PI);
      ctx.fillStyle = col; ctx.fill();
      ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.2; ctx.stroke();
      ctx.fillStyle = col;
      ctx.font = '11px ui-monospace, Consolas, monospace';
      ctx.fillText(label, p0.x + 8, p0.y);
    }
  }

  // S5 / E11: paint user-drawn curve + its σ-iterates.
  function paintCurveImage() {
    const ctx = getCtx();
    // Draft (during shift-drag): light gray.
    if (sState.curveImageDraft && sState.curveImageDraft.length > 1) {
      ctx.strokeStyle = 'rgba(80, 80, 80, 0.6)';
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      for (let i = 0; i < sState.curveImageDraft.length; i++) {
        const p = worldToPixel(sState.curveImageDraft[i].re, sState.curveImageDraft[i].im);
        if (i === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
      }
      ctx.stroke();
    }
    if (!sState.curveImage) return;
    const K = sState.curveImage.length;
    for (let it = 0; it < K; it++) {
      const pts = sState.curveImage[it];
      if (!pts || pts.length < 2) continue;
      // Color ramp blue → red as iteration grows.
      const t = (K > 1) ? (it / (K - 1)) : 0;
      const r = Math.round(40  + (180 - 40)  * t);
      const g = Math.round(80  + (50  - 80)  * t);
      const b = Math.round(180 + (40  - 180) * t);
      ctx.strokeStyle = `rgba(${r},${g},${b},${0.85 - 0.5 * t})`;
      ctx.lineWidth = (it === 0) ? 1.8 : 1.1;
      ctx.beginPath();
      for (let i = 0; i < pts.length; i++) {
        const p = worldToPixel(pts[i].re, pts[i].im);
        if (i === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
      }
      ctx.stroke();
    }
  }

  // S5 / E10: paint cycles as filled markers + connecting polylines.
  function paintCycles() {
    if (!sState.cycles || sState.cycles.length === 0) return;
    const ctx = getCtx();
    const palette = ['#108a40', '#b53030', '#5677a8', '#b8860b', '#7e2d8c'];
    for (let i = 0; i < sState.cycles.length; i++) {
      const c = sState.cycles[i];
      const col = palette[i % palette.length];
      // Connecting cycle polyline (closed).
      if (c.points.length >= 2) {
        ctx.strokeStyle = col; ctx.lineWidth = 1.4;
        if (ctx.setLineDash) ctx.setLineDash([3, 2]);
        ctx.beginPath();
        for (let k = 0; k < c.points.length; k++) {
          const p = worldToPixel(c.points[k].re, c.points[k].im);
          if (k === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
        }
        ctx.closePath();
        ctx.stroke();
        if (ctx.setLineDash) ctx.setLineDash([]);
      }
      // Dots at each cycle point.
      for (const w of c.points) {
        const p = worldToPixel(w.re, w.im);
        ctx.beginPath(); ctx.arc(p.x, p.y, 5, 0, 2*Math.PI);
        ctx.fillStyle = col; ctx.fill();
        ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.3; ctx.stroke();
      }
      // Period label near the first point.
      const p0 = worldToPixel(c.points[0].re, c.points[0].im);
      ctx.fillStyle = col;
      ctx.font = '10px ui-monospace, Consolas, monospace';
      ctx.fillText('n=' + c.period, p0.x + 7, p0.y - 7);
    }
  }

  // S5 / H8: paint orbit-family sweep as a faint family of polylines.
  function paintSweepOrbits() {
    if (!sState.sweepOrbits || sState.sweepOrbits.length === 0) return;
    const ctx = getCtx();
    for (let i = 0; i < sState.sweepOrbits.length; i++) {
      const pts = sState.sweepOrbits[i];
      if (!pts || pts.length < 2) continue;
      const t = (sState.sweepOrbits.length > 1) ? (i / (sState.sweepOrbits.length - 1)) : 0;
      // Hue ramp around the wheel for seed position.
      const hue = Math.round(360 * t);
      ctx.strokeStyle = `hsla(${hue}, 65%, 45%, 0.65)`;
      ctx.lineWidth = 0.9;
      ctx.beginPath();
      for (let k = 0; k < pts.length; k++) {
        const p = worldToPixel(pts[k].re, pts[k].im);
        if (k === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
      }
      ctx.stroke();
    }
  }

  // S3: paint the limit-set point cloud as 1-px translucent dots.
  function paintLimitSet() {
    const pts = sState.limitSet;
    if (!pts || pts.length === 0) return;
    const ctx = getCtx();
    ctx.fillStyle = 'rgba(255, 240, 80, 0.7)';      // pale yellow
    for (let i = 0; i < pts.length; i += 2) {
      const p = worldToPixel(pts[i], pts[i + 1]);
      // Use a 1×1 px rect for speed (faster than arc on dense clouds).
      ctx.fillRect(p.x - 0.5, p.y - 0.5, 1, 1);
    }
  }

  function setProgress(msg) {
    const el = document.getElementById('schwarz-progress');
    if (el) el.textContent = msg;
  }

  // ---------------------------------------------------------------------------
  // Colormaps + scale modes. Tables match schwarz-webgl.js so CPU and GPU
  // outputs render the same colors for the same input.
  // ---------------------------------------------------------------------------
  function cpuComputeT(n, maxIter, scaleMode, modK) {
    if (scaleMode === 'log') {
      return Math.min(1, Math.max(0, Math.log(n + 1) / Math.log(maxIter + 1)));
    }
    if (scaleMode === 'sqrt') {
      return Math.min(1, Math.max(0, Math.sqrt(n / maxIter)));
    }
    if (scaleMode === 'modulo') {
      const k = Math.max(2, modK | 0);
      return ((n - 1) % k) / k;
    }
    let t = (n - 1) / Math.max(1, maxIter - 1);
    if (scaleMode === 'discrete') {
      t = (Math.floor(t * maxIter) + 0.5) / maxIter;
    }
    return Math.min(1, Math.max(0, t));
  }
  function colormap(name, t) {
    t = Math.max(0, Math.min(1, t));
    if (name === 'cyclic') {
      const tt = (t * 6) % 1;
      return interpStops(tt, CMAP.magma);
    }
    return interpStops(t, CMAP[name] || CMAP.magma);
  }
  function interpStops(t, stops) {
    const n = stops.length - 1;
    const f = t * n;
    const i = Math.min(n - 1, Math.floor(f));
    const u = f - i;
    const a = stops[i], b = stops[i + 1];
    return [
      Math.round(a[0] + (b[0] - a[0]) * u),
      Math.round(a[1] + (b[1] - a[1]) * u),
      Math.round(a[2] + (b[2] - a[2]) * u),
    ];
  }
  const CMAP = {
    magma:      [[0,0,4],[28,16,68],[79,18,123],[129,37,129],[181,54,122],[229,80,100],[251,135,97],[254,194,135],[252,253,191]],
    inferno:    [[0,0,4],[31,12,72],[85,15,109],[136,34,106],[186,54,85],[227,89,51],[249,140,10],[249,201,50],[252,255,164]],
    plasma:     [[13,8,135],[75,3,161],[125,3,168],[168,34,150],[203,70,121],[229,107,93],[248,148,65],[253,195,40],[240,249,33]],
    viridis:    [[68,1,84],[72,40,120],[62,73,137],[49,104,142],[38,130,142],[31,158,137],[53,183,121],[109,205,89],[180,222,44],[253,231,37]],
    cividis:    [[0,32,76],[0,52,110],[40,75,124],[80,100,128],[120,127,128],[161,156,124],[197,187,108],[230,219,84],[253,253,51]],
    turbo:      [[48,18,59],[71,118,238],[26,196,231],[26,231,153],[97,239,71],[202,231,33],[255,184,33],[255,113,33],[224,40,9],[122,4,2]],
    grayscale:  [[0,0,0],[64,64,64],[128,128,128],[192,192,192],[255,255,255]],
    rainbow:    [[148,0,211],[75,0,130],[0,0,255],[0,255,0],[255,255,0],[255,127,0],[255,0,0]],
    iceandfire: [[10,40,100],[60,120,200],[160,210,240],[245,245,245],[250,210,90],[235,120,40],[170,30,30]],
    twotone:    [[245,245,248],[120,130,200],[40,50,110],[20,30,70]],
  };


    return {
      clearCanvas, paintAll, repaintField, paintBoundaryOnTop, paintOrbit,
      paintPreimageTree, paintLimitSet, paintZView, setProgress,
    };
  };
})(typeof window !== 'undefined' ? window : globalThis);
