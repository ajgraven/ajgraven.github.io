// =============================================================================
// schwarz-ui.js -- Sidebar + canvas overlay for the "Schwarz dynamics" tab.
//
// On first tab activation, mounts a sidebar with:
//   • Source-of-φ card     (mirrors the Inverse tab's last successful solve)
//   • Render controls      (resolution / maxIter / colormap / scale / renderer)
//   • Click & hover info
//
// Two renderer paths share a dispatcher (activeRenderer()):
//   • GPU (default): WebGL 2 fragment shader from schwarz-webgl.js. One full
//     frame in 10-30 ms typical. Drag/zoom calls renderImmediate() per
//     mousemove for true interactive panning.
//   • CPU (fallback): progressive 4×4 → 2×2 → 1×1 pyramid. Preferred path is
//     the dedicated QD.SchwarzCpuWorker (A7), which computes the whole pyramid
//     off the main thread and streams one field snapshot per pass; when no
//     Worker is available (file:// / clone failure) it runs in-process,
//     chunked across requestAnimationFrame ticks (~14 ms per slice) so the UI
//     stays responsive (_renderCpuPyramid). Per-pixel warm seeds from the left
//     neighbor in raster order → Newton converges in 1-3 iterations on average.
//
// Layer stacking: the GPU canvas (#schwarz-gl-canvas) is inserted behind the
// main #canvas. The 2D canvas keeps its in-flow layout but gets z-index 1 +
// position:relative so it can host transparent overlays (boundary + orbit)
// on top of the GPU pixels. The GL layer hides on tab-out.
//
// Fractal-mode canvas interaction (plane view):
//   • hover        → live forward σ-orbit preview {w₀, σ(w₀), …} (toggleable;
//                    rAF-coalesced; only inside Ω where σ is defined)
//   • single click → PIN that orbit (deferred by CLICK_DELAY so a double-click
//                    can cancel it; dragMoved guard suppresses pins after pans)
//   • double click → seed a preimage tree at the click, gated to the tiling
//                    set (escapeTime kind 'fundamental'); drawn over the fractal
// The preimage tree is a fractal-mode overlay, not a separate mode.
//
// Hover readout: pixel coords + escape time (CPU mode keeps a per-pixel field;
// GPU mode does an ad-hoc per-cursor escapeTime — see onMouseMove).
// =============================================================================

(function () {
  'use strict';
  if (typeof QD === 'undefined') return;

  const sState = {
    schwarz: null,            // built QD.Schwarz handle
    phiSnapshot: null,        // captured phi (with shape info for label)
    hDataSnapshot: null,
    boundarySnapshot: null,
    mounted: false,

    view: {
      // world ↔ pixel transform, centered on Ω
      cx: 0, cy: 0, scale: 200,           // px per unit
      cssW: 600, cssH: 600,
    },

    // z-plane view transform (independent pan/zoom): the unit disk 𝔻 (bounded)
    // or its exterior 𝔻* (unbounded). Kept separate from `view` so the w-plane
    // and z-diskframings don't corrupt each other and so the w-side painters
    // (which read `view` via worldToPixel) can never draw through a z-transform.
    zView: {
      cx: 0, cy: 0, scale: 200,           // px per unit (z-coordinates)
      cssW: 600, cssH: 600,
    },

    grid: {
      resolution: 384,                     // active sample count along the canvas's shorter side
      maxIter: 24,
      colormap: 'magma',
      scaleMode: 'smooth',                 // 'smooth' | 'discrete' | 'log' | 'sqrt' | 'modulo'
      modK: 8,                             // modulo period (modulo mode only)
      renderer: 'auto',                    // 'auto' | 'gpu' | 'cpu'
    },

    // GPU renderer handle (null until first capture or if WebGL 2 missing).
    gpu:    null,
    gpuMsg: '',

    // Computed escape-time field for the current view.
    field: null,                           // Int16Array; length = gridW*gridH
    fieldW: 0, fieldH: 0,
    fieldKind: null,                       // Uint8Array of escape kinds (0=fund, 1=esc, 2=int, 3=invalid)
    rendering: false,
    renderToken: 0,

    // Forward σ-orbit overlay. `orbit` is the AUTHORITATIVE pinned orbit
    // (kept in sync with pinnedOrbit) because downstream consumers — the
    // z-panel (_recomputeZPanelOrbit), sphere view, sweep, and PNG export —
    // all read sState.orbit. The transient hover orbit lives separately in
    // hoverOrbit and must NOT be assigned to `orbit` (so the z-panel etc.
    // don't flicker as the cursor moves).
    orbit:            [],                  // pinned forward-orbit polyline (= pinnedOrbit)
    pinnedOrbit:      [],                  // single-click-pinned orbit
    hoverOrbit:       null,                // transient hover-preview orbit (null = none)
    hoverOrbitEnabled: true,               // "Live orbit on hover" toggle (default ON)
    _hoverRaf:        null,                // rAF id coalescing hover-orbit recompute
    _pendingHoverW:   null,                // latest hovered world point awaiting the rAF
    _clickTimer:      null,                // deferred single-click pin (cancelled by dblclick)

    // View toggle (HANDOFF #29): 'plane' = Schwarz dynamics on the w-plane
    // (the original Schwarz tab), 'z' = the same tiling uniformized onto the
    // unit disk via z = ψ(w) (CPU-rendered: w = φ(z) then escape-time),
    // 'sphere' = the iteration textured onto a Riemann sphere. The sphere
    // renderer is lazy-mounted via QD.SphereView on first switch to sphere mode.
    viewMode:   'plane',                   // 'plane' | 'z' | 'sphere'
    sphereView: null,                      // QD.SphereView handle

    // Render mode. The preimage tree is no longer its own mode — it is a
    // fractal-mode OVERLAY drawn on top of the GL fractal, seeded by
    // double-clicking anywhere in the tiling set (TODO #16). Plane-view only.
    mode:           'fractal',             // 'fractal' | 'domain-coloring'
    preimageTree:   null,                  // last-built tree (S1; fractal-mode overlay)
    preimageDepth:  4,                     // current depth setting
    preimageBudget: 4096,                  // current visual budget

    // Limit-set overlay (S3). Chaos game over σ⁻¹; rendered as a point
    // cloud on the 2D canvas overlay. dim_H estimate displayed in the
    // sidebar. Sized small by default to keep compute time < 1s.
    limitSet:       null,                  // Float64Array of [re,im,...] points
    limitSetDim:    null,                  // last-computed dim_H estimate
    limitSetN:      5000,

    // S4 Analysis overlays.
    showSingularities: false,              // F3: σ-poles + branch points on canvas
    sigmaSingularities: null,              // cached { poles, branchPoints }
    showLevelCurves:   false,              // F12: |σ| solid + arg(σ) dashed contours
    levelCurves:       null,               // cached { abs, arg } contour lists
    sigmaForm:         null,               // E13: { family, phiText, fText, sigmaText, ... }

    // S5 Forward-dynamics state.
    showCriticalOrbits: false,             // H7: σ-orbits of canonical points
    criticalOrbits:     null,              // cached [{ label, orbit: [{re,im}, …] }, …]

    // E11: user-drawn curve forward-image. shift-drag on canvas to draw.
    curveImageDraft:    null,              // [{re,im}, …] being drawn live
    curveImage:         null,              // [[{re,im},…], [{re,im},…], …] σ⁰…σᵏ
    curveImageDepth:    4,
    isDrawingCurve:     false,             // true while shift-drag is active

    // E10: cycles (period-n).
    cycles:             null,              // [{period, points}, …]
    cyclePeriodMax:     2,

    // H8: orbit-family sweep — pre-baked horizontal line through Ω center.
    sweepOrbits:        null,              // [[{re,im}, …], …]
    sweepN:             16,
    sweepDepth:         12,

    // Domain-coloring (F6) backing field. When mode='domain-coloring', this
    // overlays the GPU/CPU fractal. Cached per φ-viewport.
    domainColor:        null,              // { buf, W, H, viewport }

    // ψ-pullback of the pinned w-orbit, drawn in the z-plane view.
    zPanelOrbit:        null,              // [{re, im}, ...] z-history of sState.orbit
  };

  // Kinds enum
  const KIND_FUND = 0, KIND_ESC = 1, KIND_INT = 2, KIND_INV = 3, KIND_OUTSIDE = 4;

  // Interaction tuning. CLICK_DELAY defers the single-click orbit-pin long
  // enough for a double-click (which seeds the preimage tree) to cancel it;
  // the dblclick handler clears the pending timer. Exposed via the test hook.
  // The preimage-tree seed gate runs escapeTime with a generous iteration cap
  // so genuinely-fundamental (tiling-set) points that escape slowly aren't
  // mis-classified as non-escaping `interior`. Display/hover use the smaller
  // grid.maxIter — this larger cap is for the accept/reject decision only.
  function gateMaxIter() { return Math.max(256, (sState.grid.maxIter | 0) * 4); }

  // Forward bindings for the Phase-3 paint module (assigned by the install near
  // the end of this file; called by name throughout — see schwarz-paint.js).
  let clearCanvas, paintAll, repaintField, paintBoundaryOnTop, paintOrbit;
  let paintPreimageTree, paintLimitSet, paintZView, setProgress;
  let requestRecompute;
  let attachCanvasHandlers, onCanvasClick, onCanvasDblClick, onMouseMove;
  let runHoverOrbit, pinOrbitAt, _schwarzInter;
  // Forward bindings for the Phase-3 feature-compute module (assigned by the
  // install near the end of this file; the card-builder event handlers + the
  // interaction install call them by name — see schwarz-features.js).
  let _recomputeDomainColoring, _rebuildPreimageTreeIfActive, _refreshPreimageTreeStats;
  let _computeLimitSet, _clearLimitSet, _recomputeCriticalOrbits, _findCycles;
  let _exportPng, _recomputeZPanelOrbit, _computeSweep, _recomputeLevelCurves;

  // ---------------------------------------------------------------------------
  // Lazy mount
  // ---------------------------------------------------------------------------
  document.addEventListener('tab-changed', function (e) {
    if (!e.detail || e.detail.tab !== 'schwarz') {
      // Leaving the Schwarz tab — don't keep rendering, and hide BOTH GL
      // layers so they can't show through under another tab's drawing.
      sState.renderToken++;
      showGLLayer(false);
      // Cancel any pending interaction timers/frames so they can't fire
      // (and repaint into another tab) after we've left.
      if (sState._hoverRaf != null) { cancelAnimationFrame(sState._hoverRaf); sState._hoverRaf = null; }
      if (sState._clickTimer != null) { clearTimeout(sState._clickTimer); sState._clickTimer = null; }
      sState.hoverOrbit = null;
      // HANDOFF #34: also wipe our pixels from the shared 2D canvas so the
      // CPU-pyramid pixmap and orbit-polyline overlay don't briefly bleed
      // through into whichever tab takes over. The receiving tab's
      // tab-changed handler is responsible for repainting its own
      // background / axes (the QD tab now does this via plot.resize()).
      const ctx = getCtx();
      if (ctx) ctx.clearRect(0, 0, sState.view.cssW, sState.view.cssH);
      if (sState.sphereView) sState.sphereView.deactivate();
      return;
    }
    if (!sState.mounted) { mountSchwarzSidebar(); sState.mounted = true; }
    refreshSourceStatus();
    if (sState.viewMode === 'sphere') {
      // sphere mode: hide the Schwarz GL layer, activate sphere view.
      showGLLayer(false);
      _activateSphereView();
    } else {
      // 2D views (plane or z). Both can use the GPU (the shader's u_viewMode
      // branch handles z → w = φ(z)); doRecompute owns the final GL visibility,
      // this just avoids a flash. PQDs fall back to CPU via activeRenderer().
      if (sState.sphereView) sState.sphereView.deactivate();
      if (sState.schwarz) {
        const is2D = sState.viewMode === 'plane' || sState.viewMode === 'z';
        showGLLayer(is2D && activeRenderer() === 'gpu');
        requestRecompute();
      } else {
        showGLLayer(false);
        clearCanvas();
      }
    }
  });

  // Repaint the Schwarz view on window resize. The shared 2D #canvas backing
  // store is resized (and cleared) by the QD plot's global resize handler
  // (ui.js); without this the Schwarz overlay/field is left blank or stale at
  // the new size until the user pans/zooms. Only acts when the Schwarz tab is
  // active + in plane mode; rAF-coalesced. Sphere mode has its own resize
  // observer (sphere-ui.js). requestRecompute / _recomputeDomainColoring /
  // paintAll are forward-declared lets, assigned by the installs at the tail —
  // resolved by the time this fires.
  let _schwarzResizeRaf = null;
  window.addEventListener('resize', () => {
    if (_schwarzResizeRaf) return;
    _schwarzResizeRaf = requestAnimationFrame(() => {
      _schwarzResizeRaf = null;
      if (!isSchwarzActive() || !sState.schwarz || sState.viewMode === 'sphere') return;
      syncCanvasSize();
      if (sState.viewMode === 'plane' && sState.mode === 'domain-coloring') { _recomputeDomainColoring(); paintAll(); }
      else requestRecompute();
    });
  });

  function mountSchwarzSidebar() {
    const root = document.getElementById('controls-schwarz');
    if (!root) return;
    root.innerHTML = '';
    // Source φ first — it's the first tile a user touches when opening the tab
    // (capture a φ from the Inverse tab before anything else is meaningful).
    root.appendChild(makeSourceCard());
    root.appendChild(makeViewToggleCard());
    root.appendChild(makeModeCard());
    root.appendChild(makeOverlaysCard());
    root.appendChild(makeLimitSetCard());
    root.appendChild(makeAnalysisCard());
    root.appendChild(makeForwardCard());
    root.appendChild(makeRenderCard());
    root.appendChild(makeInfoCard());
    // Mount-time placeholder for SphereView's display + camera cards. The
    // cards are lazily appended into this container when the user first
    // toggles to sphere mode (via _activateSphereView → QD.SphereView.mount).
    const sphereSlot = document.createElement('div');
    sphereSlot.id = 'schwarz-sphere-slot';
    root.appendChild(sphereSlot);
    attachCanvasHandlers();
    _applyViewModeVisibility();
    attachSchwarzHelp();      // HANDOFF #33
  }

  function attachSchwarzHelp() {
    if (!window.QD || !window.QD.QoL || !window.QD.QoL.attachHelp) return;
    const H = window.QD.QoL.attachHelp;
    const root = document.getElementById('controls-schwarz');
    if (!root) return;
    // Attach by card id (not header index) so sidebar reordering can't scramble
    // which card a popover lands on. The Source φ tile is intentionally compact,
    // so its full description lives here in the "?" help.
    const help = (sel, html) => { const h = root.querySelector(sel + ' h2'); if (h) H(h, html); };
    help('#schwarz-source-card',
      `<b>Source φ.</b> The Schwarz reflection σ(w) = conj(F(ψ(w))) is built from a
       Riemann map φ produced by the inverse solver. Solve on the Inverse tab, then
       click <b>Use this φ</b> to snapshot it here; the tab iterates σ and colors each
       pixel by escape time. All six inverse families are supported (classical
       bounded / unbounded, and all four LQD variants).`);
    help('#schwarz-render-card',
      `<b>Render.</b> The CPU fallback draws a progressive 4×4→2×2→1×1 pyramid,
       computed in a background Web Worker when available (falling back to an
       in-page, animation-frame-sliced render) so the page stays responsive
       while it refines; GPU uses a WebGL 2 fragment shader for instant frames.
       <i>Colormap</i> + <i>scaleMode</i> change the escape-time → colour
       mapping; <i>maxIter</i> caps the σ-iteration before declaring a pixel
       interior; <i>mod k</i> emphasises orbit-period structure.`);
    help('#schwarz-view-card',
      `<b>View.</b> <b>plane</b> draws the σ-dynamics tiling on the w-plane (Ω).
       <b>z-disk</b> shows the SAME tiling uniformized onto the unit disk 𝔻 (or 𝔻*
       for unbounded Ω): each z is colored by the escape time of w = φ(z). Every
       overlay (orbit, preimage tree, limit set, …) is shown ψ-pulled-back into 𝔻.
       <b>sphere</b> textures the iteration onto the Riemann sphere. plane + z-disk
       have independent pan/zoom and full click/hover interaction.`);
    help('#schwarz-info-card',
      `<b>Click & hover.</b> Single-click in Ω to pin a σ-orbit; double-click in the
       tiling set to seed a preimage tree; hover reads the coordinate + escape time
       (and, in CPU mode, the pixel kind). Drag to pan, scroll to zoom. In the
       z-disk view every click/hover point z maps to w = φ(z) before iterating.`);
  }

  function makeViewToggleCard() {
    const card = document.createElement('section');
    card.className = 'card';
    card.id = 'schwarz-view-card';
    card.innerHTML = `
      <h2>View</h2>
      <div class="segmented" role="tablist" aria-label="View mode">
        <button class="seg-btn active" data-view="plane"  type="button">plane</button>
        <button class="seg-btn"        data-view="z"      type="button">z-disk</button>
        <button class="seg-btn"        data-view="sphere" type="button">sphere</button>
      </div>
    `;
    setTimeout(() => {
      card.querySelectorAll('.seg-btn').forEach(btn => {
        btn.addEventListener('click', () => setViewMode(btn.dataset.view));
      });
    }, 0);
    return card;
  }

  // ---------------------------------------------------------------------------
  // Mode card: fractal vs domain-coloring. Plane-view only. The preimage tree
  // is a fractal-mode overlay (seeded by double-click), not its own mode —
  // its depth/budget controls live in the fractal-mode options block below.
  // ---------------------------------------------------------------------------
  function makeModeCard() {
    const card = document.createElement('section');
    card.className = 'card view-2d';
    card.id = 'schwarz-mode-card';
    card.innerHTML = `
      <h2>Mode</h2>
      <div class="view-plane-only">
        <div class="segmented" role="tablist" aria-label="Schwarz mode">
          <button class="seg-btn active" data-mode="fractal"        type="button">fractal</button>
          <button class="seg-btn"        data-mode="domain-coloring" type="button">domain color</button>
        </div>
      </div>
      <div id="schwarz-mode-options-fractal" style="margin-top:8px;">
        <div style="font-size:12px; color:#555; margin-bottom:6px;">
          <b>Double-click</b> in the tiling set to seed a preimage tree (each
          generation applies σ⁻¹). <b>Single-click</b> pins the forward orbit;
          hover previews it.
        </div>
        <label style="display:block; font-size:12px; margin:4px 0;">
          <input type="checkbox" id="schwarz-hover-orbit" checked
                 style="vertical-align:middle; margin-right:4px;">
          Live orbit on hover
        </label>
        <label style="display:block; font-size:12px; margin:4px 0;">
          Tree depth:
          <input type="range" min="1" max="8" value="4" id="schwarz-preimage-depth"
                 style="vertical-align:middle; margin-left:6px; width:120px;">
          <span id="schwarz-preimage-depth-val" style="font-family:monospace;">4</span>
        </label>
        <label style="display:block; font-size:12px; margin:4px 0;">
          Visual budget:
          <select id="schwarz-preimage-budget" style="font-size:12px;">
            <option value="1024">1024</option>
            <option value="4096" selected>4096</option>
            <option value="16384">16384</option>
          </select>
          <span id="schwarz-preimage-count" style="font-family:monospace; margin-left:8px; color:#777;"></span>
        </label>
      </div>
    `;
    setTimeout(() => {
      card.querySelectorAll('.seg-btn').forEach(btn => {
        btn.addEventListener('click', () => setMode(btn.dataset.mode));
      });
      const hover = card.querySelector('#schwarz-hover-orbit');
      if (hover) {
        hover.addEventListener('change', () => {
          sState.hoverOrbitEnabled = hover.checked;
          if (!hover.checked) {
            // Turning the preview off — drop any live hover orbit immediately.
            if (sState._hoverRaf != null) { cancelAnimationFrame(sState._hoverRaf); sState._hoverRaf = null; }
            sState.hoverOrbit = null;
            paintBoundaryOnTop();
          }
        });
      }
      const depthSlider = card.querySelector('#schwarz-preimage-depth');
      const depthVal    = card.querySelector('#schwarz-preimage-depth-val');
      if (depthSlider) {
        depthSlider.addEventListener('input', () => {
          sState.preimageDepth = +depthSlider.value;
          if (depthVal) depthVal.textContent = String(sState.preimageDepth);
          _rebuildPreimageTreeIfActive();
        });
      }
      const budget = card.querySelector('#schwarz-preimage-budget');
      if (budget) {
        budget.addEventListener('change', () => {
          sState.preimageBudget = +budget.value;
          _rebuildPreimageTreeIfActive();
        });
      }
    }, 0);
    return card;
  }

  function setMode(mode) {
    if (mode !== 'fractal' && mode !== 'domain-coloring') return;
    if (mode === sState.mode) return;
    sState.mode = mode;
    document.querySelectorAll('#schwarz-mode-card .seg-btn').forEach(btn => {
      btn.classList.toggle('active', btn.dataset.mode === mode);
    });
    _applyModeOptionsVisibility();
    // All overlays (tree + orbits) are fractal-mode-only; clear them on any
    // mode transition so they don't bleed into domain-coloring or linger when
    // we re-enter fractal mode.
    sState.preimageTree = null;
    sState.orbit = []; sState.pinnedOrbit = []; sState.hoverOrbit = null;
    if (mode === 'domain-coloring') {
      // S5 / F6: render σ-domain-coloring into the 2D canvas; hide GL.
      showGLLayer(false);
      _recomputeDomainColoring();
      paintAll();
    } else {
      sState.domainColor = null;
      if (sState.schwarz) {
        showGLLayer(activeRenderer() === 'gpu');
        requestRecompute();
      }
    }
  }

  // Feature-compute (σ domain-coloring, preimage-tree rebuild + stats) ->
  // schwarz-features.js (Phase-3 item E). Installed near the end of this file;
  // the card builders + interaction install call the captured names.

  // ---------------------------------------------------------------------------
  // Overlays card: clear controls for the plotted overlays. Individual overlay
  // cards already carry per-feature "Clear" buttons (limit set, cycles, sweep,
  // curve), but the click-pinned/hover orbit and the double-click preimage tree
  // had no user-facing clear at all. This card adds a dedicated "Clear orbit"
  // and a master "Clear all overlays" that wipes every COMPUTED/DRAWN overlay
  // (orbit, preimage tree, limit set, cycles, sweep, curve image) in one click.
  // It deliberately leaves the checkbox DISPLAY TOGGLES (σ-poles, level curves,
  // canonical-point orbits, z-panel) as the user set them — those hide via
  // their own checkboxes, so wiping their data while the box stays checked would
  // be inconsistent. Plane-view only (overlays are a plane-view concept).
  // ---------------------------------------------------------------------------
  function makeOverlaysCard() {
    const card = document.createElement('section');
    card.className = 'card view-2d';
    card.id = 'schwarz-overlays-card';
    card.innerHTML = `
      <h2>Overlays</h2>
      <div style="font-size:12px; color:#555; margin-bottom:6px;">
        Clear plotted overlays. <b>Clear all</b> removes orbits, the preimage
        (inverse) tree, the limit set, cycles, sweeps and curve images; the
        display-toggle checkboxes below are left as set.
      </div>
      <button type="button" id="schwarz-clear-orbit" class="small">Clear orbit</button>
      <button type="button" id="schwarz-clear-overlays" class="small"
              style="margin-left:6px;">Clear all overlays</button>
    `;
    setTimeout(() => {
      card.querySelector('#schwarz-clear-orbit').addEventListener('click', clearOrbit);
      card.querySelector('#schwarz-clear-overlays').addEventListener('click', clearAllOverlays);
    }, 0);
    return card;
  }

  // Clear the forward σ-orbit (both the pinned polyline and any live hover
  // preview). The z-panel reads sState.orbit, so refresh its pullback when the
  // panel is open (else just drop it). paintBoundaryOnTop repaints overlays over
  // the GL/CPU field — the same repaint the per-feature clear buttons use.
  function clearOrbit() {
    if (sState._hoverRaf != null) { cancelAnimationFrame(sState._hoverRaf); sState._hoverRaf = null; }
    sState.orbit = [];
    sState.pinnedOrbit = [];
    sState.hoverOrbit = null;
    _recomputeZPanelOrbit();           // orbit now empty → clears the z-view pullback
    paintBoundaryOnTop();
  }

  // Master clear: every computed/drawn overlay. Leaves the checkbox display
  // toggles (showSingularities / showLevelCurves / showCriticalOrbits) and their
  // cached data untouched — those are dismissed via their own checkboxes
  // (clearing the data while the box stays checked would desync).
  function clearAllOverlays() {
    clearOrbit();                       // orbit (+ z-view pullback)
    sState.preimageTree = null;         // S1 inverse tree
    sState.limitSet = null;             // S3 limit set
    sState.limitSetDim = null;
    sState.cycles = null;               // E10 cycles
    sState.sweepOrbits = null;          // H8 sweep
    sState.curveImage = null;           // E11 curve forward-image
    sState.curveImageDraft = null;
    // Blank the status/count labels owned by the per-overlay cards so they
    // don't claim a now-cleared overlay still exists.
    ['schwarz-ls-status', 'schwarz-ls-dim', 'schwarz-cycle-count',
     'schwarz-preimage-count'].forEach(id => {
      const el = document.getElementById(id);
      if (el) el.textContent = '';
    });
    paintBoundaryOnTop();
  }

  // ---------------------------------------------------------------------------
  // Limit-set card (S3): chaos-game point cloud + dim_H readout.
  // Visible in plane view (it's an overlay on the fractal). Hidden in
  // sphere view via the view-plane-only class.
  // ---------------------------------------------------------------------------
  function makeLimitSetCard() {
    const card = document.createElement('section');
    card.className = 'card view-2d';
    card.id = 'schwarz-limit-set-card';
    card.innerHTML = `
      <h2>Limit set</h2>
      <div style="font-size:12px; color:#555; margin-bottom:6px;">
        Random walk through σ⁻¹ approximates the limit set (∂ tiling set).
      </div>
      <label style="display:block; font-size:12px; margin:4px 0;">
        Sample size:
        <select id="schwarz-ls-n">
          <option value="1000">1k</option>
          <option value="5000" selected>5k</option>
          <option value="20000">20k</option>
          <option value="100000">100k</option>
        </select>
      </label>
      <button type="button" id="schwarz-ls-compute"
              style="margin-top:4px; font-size:12px;">Compute limit set</button>
      <button type="button" id="schwarz-ls-clear"
              style="margin-top:4px; font-size:12px;">Clear</button>
      <div id="schwarz-ls-status" style="font-size:12px; margin-top:6px; color:#555;"></div>
      <div id="schwarz-ls-dim" style="font-size:12px; margin-top:4px; font-family:monospace;"></div>
    `;
    setTimeout(() => {
      const nSel = card.querySelector('#schwarz-ls-n');
      nSel.addEventListener('change', () => { sState.limitSetN = +nSel.value; });
      card.querySelector('#schwarz-ls-compute').addEventListener('click', _computeLimitSet);
      card.querySelector('#schwarz-ls-clear').addEventListener('click', _clearLimitSet);
    }, 0);
    return card;
  }

  // Limit-set compute/clear (_computeLimitSet / _clearLimitSet) ->
  // schwarz-features.js (Phase-3 item E).

  // ---------------------------------------------------------------------------
  // S4: Analysis card — explicit σ form (E13), singularities (F3), level
  // curves (F12). All overlay on the plane view; sphere mode hides this card.
  // ---------------------------------------------------------------------------
  function makeAnalysisCard() {
    const card = document.createElement('section');
    card.className = 'card view-2d';
    card.id = 'schwarz-analysis-card';
    card.innerHTML = `
      <h2>σ analysis</h2>
      <div style="font-size:12px; margin-bottom:6px;">
        <button type="button" id="schwarz-show-sigma-form"
                style="font-size:12px;">Show σ(w) form</button>
      </div>
      <div id="schwarz-sigma-form-out"
           style="display:none; font-size:12px; background:#f7f7f9; padding:6px;
                  border:1px solid #e0e0e8; border-radius:4px; max-height:200px;
                  overflow:auto; margin:4px 0;"></div>
      <label style="display:block; font-size:12px; margin:4px 0;">
        <input type="checkbox" id="schwarz-show-singularities"> Show σ-poles + branch points
        <span id="schwarz-sigma-sing-count"
              style="font-family:monospace; color:#777; margin-left:6px;"></span>
      </label>
      <label style="display:block; font-size:12px; margin:4px 0;">
        <input type="checkbox" id="schwarz-show-level-curves"> Show |σ| / arg(σ) level curves
      </label>
    `;
    setTimeout(() => {
      const btn = card.querySelector('#schwarz-show-sigma-form');
      const out = card.querySelector('#schwarz-sigma-form-out');
      btn.addEventListener('click', () => {
        if (!sState.schwarz) {
          out.style.display = '';
          out.textContent = '(no φ captured)';
          return;
        }
        const f = QD.Schwarz.explicitSigmaForm(sState.schwarz);
        sState.sigmaForm = f;
        out.style.display = '';
        // Typeset the LaTeX forms (already produced by explicitSigmaForm) with
        // KaTeX; fall back per-row to the ASCII text if KaTeX is unavailable or
        // a particular expression fails to parse.
        out.innerHTML =
          '<div class="ssf-fam"></div>' +
          '<div class="ssf-row ssf-phi"></div>' +
          '<div class="ssf-row ssf-f"></div>' +
          '<div class="ssf-row ssf-sigma"></div>';
        out.querySelector('.ssf-fam').textContent = 'family: ' + f.family;
        const rk = (sel, latex, text) => {
          const el = out.querySelector(sel);
          if (!el) return;
          if (window.katex && latex) {
            try { window.katex.render(latex, el, { displayMode: true, throwOnError: false }); return; }
            catch (e) { /* fall through to text */ }
          }
          el.textContent = text || latex || '';
        };
        rk('.ssf-phi',   f.phiLatex,   f.phiText);
        rk('.ssf-f',     f.fLatex,     f.fText);
        rk('.ssf-sigma', f.sigmaLatex, f.sigmaText);
      });
      card.querySelector('#schwarz-show-singularities').addEventListener('change', (e) => {
        sState.showSingularities = e.target.checked;
        if (e.target.checked) {
          if (sState.schwarz) {
            sState.sigmaSingularities = QD.Schwarz.findSigmaSingularities(sState.schwarz);
            const el = document.getElementById('schwarz-sigma-sing-count');
            const s = sState.sigmaSingularities;
            if (el && s) el.textContent = s.poles.length + ' poles, ' + s.branchPoints.length + ' bp';
          }
        } else {
          const el = document.getElementById('schwarz-sigma-sing-count');
          if (el) el.textContent = '';
        }
        paintBoundaryOnTop();
      });
      card.querySelector('#schwarz-show-level-curves').addEventListener('change', (e) => {
        sState.showLevelCurves = e.target.checked;
        if (e.target.checked && sState.schwarz) {
          _recomputeLevelCurves();
        }
        paintBoundaryOnTop();
      });
    }, 0);
    return card;
  }

  // ---------------------------------------------------------------------------
  // S5: Forward-dynamics card. Bundles H7 (critical orbits), E11 (curve
  // forward-image), E10 (cycle finder), H8 (orbit sweep).
  // ---------------------------------------------------------------------------
  function makeForwardCard() {
    const card = document.createElement('section');
    card.className = 'card view-2d';
    card.id = 'schwarz-forward-card';
    card.innerHTML = `
      <h2>Dynamics</h2>
      <label style="display:block; font-size:12px; margin:4px 0;">
        <input type="checkbox" id="schwarz-show-critical-orbits"> Show canonical-point orbits
      </label>
      <div style="font-size:12px; margin:8px 0 4px;"><b>Cycle finder:</b></div>
      <label style="display:block; font-size:12px; margin:4px 0;">
        Period: <select id="schwarz-cycle-n">
          <option value="1">1 (fixed)</option>
          <option value="2" selected>2</option>
          <option value="3">3</option>
          <option value="4">4</option>
          <option value="5">5</option>
          <option value="6">6</option>
        </select>
        <button type="button" id="schwarz-cycle-find"
                style="font-size:11px; margin-left:6px;">Find</button>
        <button type="button" id="schwarz-cycle-clear"
                style="font-size:11px; margin-left:4px;">Clear</button>
        <span id="schwarz-cycle-count"
              style="font-family:monospace; color:#777; margin-left:6px;"></span>
      </label>
      <details style="margin-top:10px;">
        <summary style="font-size:12px; cursor:pointer; color:#555;">Advanced</summary>
        <div style="font-size:12px; margin:6px 0 4px;">
          <b>Curve forward-image:</b> shift-drag in Ω to draw (plane view).
        </div>
        <label style="display:block; font-size:12px; margin:4px 0;">
          Iterations:
          <input type="range" min="1" max="10" value="4" id="schwarz-curve-depth"
                 style="vertical-align:middle; margin-left:6px; width:100px;">
          <span id="schwarz-curve-depth-val" style="font-family:monospace;">4</span>
          <button type="button" id="schwarz-curve-clear"
                  style="font-size:11px; margin-left:6px;">Clear</button>
        </label>
        <div style="font-size:12px; margin:8px 0 4px;"><b>Orbit-family sweep:</b></div>
        <label style="display:block; font-size:12px; margin:4px 0;">
          N: <input type="number" id="schwarz-sweep-n" value="16" min="2" max="64"
                    style="width:48px;">
          Depth: <input type="number" id="schwarz-sweep-depth" value="12" min="1" max="50"
                        style="width:48px;">
          <button type="button" id="schwarz-sweep-compute"
                  style="font-size:11px; margin-left:6px;">Sweep horizontal</button>
          <button type="button" id="schwarz-sweep-clear"
                  style="font-size:11px; margin-left:4px;">Clear</button>
        </label>
        <div style="font-size:12px; margin:8px 0 4px;"><b>Export PNG:</b></div>
        <label style="display:block; font-size:12px; margin:4px 0;">
          Multiplier:
          <select id="schwarz-export-mult">
            <option value="1">1× (display)</option>
            <option value="2" selected>2×</option>
            <option value="4">4×</option>
          </select>
          <button type="button" id="schwarz-export-png"
                  style="font-size:11px; margin-left:6px;">Export PNG</button>
        </label>
      </details>
    `;
    setTimeout(() => {
      // H7
      card.querySelector('#schwarz-show-critical-orbits').addEventListener('change', (e) => {
        sState.showCriticalOrbits = e.target.checked;
        if (e.target.checked && sState.schwarz) _recomputeCriticalOrbits();
        else sState.criticalOrbits = null;
        paintBoundaryOnTop();
      });
      // E11
      const depthSlider = card.querySelector('#schwarz-curve-depth');
      const depthVal    = card.querySelector('#schwarz-curve-depth-val');
      depthSlider.addEventListener('input', () => {
        sState.curveImageDepth = +depthSlider.value;
        depthVal.textContent = String(sState.curveImageDepth);
        // Re-iterate the latest captured curve if any.
        if (sState.curveImage && sState.curveImage[0]) {
          sState.curveImage = QD.Schwarz.iterateCurveForward(
            sState.curveImage[0], sState.schwarz, sState.curveImageDepth);
          paintBoundaryOnTop();
        }
      });
      card.querySelector('#schwarz-curve-clear').addEventListener('click', () => {
        sState.curveImage = null;
        sState.curveImageDraft = null;
        paintBoundaryOnTop();
      });
      // E10
      card.querySelector('#schwarz-cycle-find').addEventListener('click', _findCycles);
      card.querySelector('#schwarz-cycle-clear').addEventListener('click', () => {
        sState.cycles = null;
        const el = document.getElementById('schwarz-cycle-count');
        if (el) el.textContent = '';
        paintBoundaryOnTop();
      });
      // H8
      card.querySelector('#schwarz-sweep-compute').addEventListener('click', _computeSweep);
      card.querySelector('#schwarz-sweep-clear').addEventListener('click', () => {
        sState.sweepOrbits = null;
        paintBoundaryOnTop();
      });
      // F8: PNG export
      card.querySelector('#schwarz-export-png').addEventListener('click', _exportPng);
    }, 0);
    return card;
  }

  // Forward-dynamics feature-compute (_recomputeCriticalOrbits / _findCycles /
  // _computeSweep / _recomputeZPanelOrbit), σ level curves
  // (_recomputeLevelCurves), and high-res PNG export (_exportPng) ->
  // schwarz-features.js (Phase-3 item E).

  function setViewMode(mode) {
    if (mode !== 'plane' && mode !== 'z' && mode !== 'sphere') return;
    if (mode === sState.viewMode) return;
    sState.viewMode = mode;
    // Update the view-toggle highlight. Scope to the view card so the Mode
    // card's fractal/domain seg-btns (which carry data-mode, not data-view)
    // aren't deactivated as a side effect.
    document.querySelectorAll('#schwarz-view-card .seg-btn').forEach(btn => {
      btn.classList.toggle('active', btn.dataset.view === mode);
    });
    _applyViewModeVisibility();
    refreshSourceStatus();
    if (mode === 'sphere') {
      showGLLayer(false);
      _activateSphereView();
      return;
    }
    // 2D views (plane or z). Both can use the GPU — the shader's u_viewMode
    // branch lifts z → w = φ(z) for the z-disk view; doRecompute owns the final
    // GL visibility, this just avoids a flash on switch. (PQD families have no
    // GPU shader and fall back to CPU automatically via activeRenderer().)
    if (sState.sphereView) sState.sphereView.deactivate();
    if (!sState.schwarz) { showGLLayer(false); clearCanvas(); return; }
    showGLLayer((mode === 'plane' || mode === 'z') && activeRenderer() === 'gpu');
    requestRecompute();
  }

  function _applyViewModeVisibility() {
    const mode = sState.viewMode;
    const is2D = mode === 'plane' || mode === 'z';
    const set = (cls, show) => {
      document.querySelectorAll('#controls-schwarz ' + cls)
        .forEach(el => { el.style.display = show ? '' : 'none'; });
    };
    set('.view-plane-only',  mode === 'plane');   // plane-only (fractal/domain, GPU)
    set('.view-z-only',      mode === 'z');        // z-only (reserved)
    set('.view-2d',          is2D);                // plane AND z (shared overlays)
    set('.view-sphere-only', mode === 'sphere');
    _applyModeOptionsVisibility();
  }

  // The fractal-options block (hover toggle + preimage-tree controls + the
  // click/double-click hint) applies to BOTH the plane fractal view and the
  // z-disk view (z is tiling, i.e. fractal semantics). Show it whenever we're in
  // z or in plane-fractal mode; hide only in plane domain-coloring.
  function _applyModeOptionsVisibility() {
    const opts = document.getElementById('schwarz-mode-options-fractal');
    if (opts) opts.style.display = (sState.viewMode === 'z' || sState.mode === 'fractal') ? '' : 'none';
  }

  // Lazy-mount QD.SphereView the first time the user switches to sphere mode;
  // then push the current captured φ (if any) and broadcast the latest render
  // params. Subsequent invocations just activate the existing handle.
  function _activateSphereView() {
    if (!sState.sphereView) {
      if (!QD.SphereView || !QD.SphereView.mount) {
        console.warn('schwarz-ui: QD.SphereView unavailable; sphere view disabled.');
        return;
      }
      const sidebarRoot = document.getElementById('schwarz-sphere-slot')
                       || document.getElementById('controls-schwarz');
      sState.sphereView = QD.SphereView.mount({
        plotArea:   document.getElementById('plot-area'),
        mainCanvas: getCanvas(),
        sidebar:    sidebarRoot,
        isActive:   () => sState.viewMode === 'sphere',
      });
      if (!sState.sphereView || !sState.sphereView.isAvailable()) {
        console.warn('schwarz-ui: sphere view unavailable (WebGL 2 missing).');
        return;
      }
      // Re-apply visibility so the newly-built display/camera cards respect
      // the current view mode (.view-sphere-only).
      _applyViewModeVisibility();
      // Push current render params so colormap/maxIter/scale carry over.
      sState.sphereView.setRenderParams({
        maxIter:   sState.grid.maxIter,
        colormap:  sState.grid.colormap,
        scaleMode: sState.grid.scaleMode,
        modK:      sState.grid.modK,
      });
      // If we already have a captured φ, push it now.
      if (sState.phiSnapshot) {
        sState.sphereView.setPhi(sState.phiSnapshot,
                                  sState.hDataSnapshot,
                                  sState.boundarySnapshot);
      }
    }
    sState.sphereView.activate();
  }

  function makeSourceCard() {
    const card = document.createElement('section');
    card.className = 'card';
    card.id = 'schwarz-source-card';
    // Compact: the explanatory text lives in the "?" hover help (attachSchwarzHelp)
    // so this — the first/most-used tile — stays small (title + status + button).
    card.innerHTML = `
      <h2>Source φ</h2>
      <div id="schwarz-src-status" class="hint" style="color:#333; margin-top:2px;">
        (no φ captured)
      </div>
      <div id="schwarz-bounded-warning" class="hint"
           style="display:none; color:#b8860b; margin-top:6px; background:#fffbe6;
                  border:1px solid #e8c840; border-radius:4px; padding:4px 8px;">
        ⚠ Bounded Ω: the sphere view is well-defined but visually uninformative
        (K fills most of the southern hemisphere). The render is still shown.
      </div>
      <button id="schwarz-capture" class="small" style="margin-top:8px;">Use this φ</button>
    `;
    setTimeout(() => {
      document.getElementById('schwarz-capture').addEventListener('click', captureFromInverseTab);
    }, 0);
    return card;
  }

  function makeRenderCard() {
    const card = document.createElement('section');
    card.className = 'card';
    card.id = 'schwarz-render-card';
    card.innerHTML = `
      <h2>Render</h2>
      <div class="row view-2d">
        <label>Resolution:
          <select id="schwarz-resolution">
            <option value="192">192</option>
            <option value="256">256</option>
            <option value="384" selected>384</option>
            <option value="512">512</option>
            <option value="768">768</option>
          </select>
        </label>
      </div>
      <div class="row" style="margin-top:8px;">
        <label>Max iterations:
          <input id="schwarz-maxiter" type="number" min="1" max="200" value="24" style="width:72px;">
        </label>
      </div>
      <div class="row" style="margin-top:8px;">
        <label>Colormap:
          <select id="schwarz-colormap">
            <option value="magma" selected>magma</option>
            <option value="inferno">inferno</option>
            <option value="plasma">plasma</option>
            <option value="viridis">viridis</option>
            <option value="cividis">cividis</option>
            <option value="turbo">turbo</option>
            <option value="grayscale">grayscale</option>
            <option value="rainbow">rainbow</option>
            <option value="iceandfire">ice & fire</option>
            <option value="twotone">two-tone</option>
            <option value="cyclic">cyclic (magma)</option>
          </select>
        </label>
      </div>
      <div class="row" style="margin-top:8px;">
        <label>Scale:
          <select id="schwarz-scalemode" title="How escape time n maps to colormap position">
            <option value="smooth" selected>smooth</option>
            <option value="discrete">discrete (per-n bands)</option>
            <option value="log">log</option>
            <option value="sqrt">sqrt</option>
            <option value="modulo">modulo (cyclic)</option>
          </select>
        </label>
        <label style="margin-left:8px;" id="schwarz-modk-wrap">
          K:
          <input id="schwarz-modk" type="number" min="2" max="64" value="8" style="width:52px;">
        </label>
      </div>
      <div class="row view-2d" style="margin-top:8px;">
        <label>Renderer:
          <select id="schwarz-renderer">
            <option value="auto" selected>auto (GPU if available)</option>
            <option value="gpu">GPU (WebGL 2)</option>
            <option value="cpu">CPU (fallback)</option>
          </select>
        </label>
      </div>
      <div class="row view-2d" style="margin-top:10px;">
        <button id="schwarz-recompute" class="small">Recompute</button>
        <button id="schwarz-fit" class="small" style="margin-left:6px;">Fit</button>
      </div>
      <div id="schwarz-progress" class="hint view-2d" style="margin-top:8px; min-height:1.2em;"></div>
    `;
    setTimeout(() => {
      document.getElementById('schwarz-renderer').addEventListener('change', e => {
        sState.grid.renderer = e.target.value;
        showGLLayer(activeRenderer() === 'gpu');
        requestRecompute();
      });
      document.getElementById('schwarz-resolution').addEventListener('change', e => {
        sState.grid.resolution = +e.target.value;
        requestRecompute();
      });
      document.getElementById('schwarz-maxiter').addEventListener('change', e => {
        sState.grid.maxIter = Math.max(1, Math.min(200, +e.target.value || 24));
        // Broadcast to sphere view too — same shared sliders for both renderers.
        if (sState.sphereView) sState.sphereView.setRenderParams({ maxIter: sState.grid.maxIter });
        if (sState.viewMode !== 'sphere') requestRecompute();   // plane + z re-iterate
      });
      // Recolor the current field without recomputing: GPU re-renders (plane or
      // z-disk — renderImmediate handles both), CPU repaints the existing field.
      const recolor = () => {
        const is2D = sState.viewMode === 'plane' || sState.viewMode === 'z';
        if (is2D && activeRenderer() === 'gpu') renderImmediate();
        else if (sState.viewMode !== 'sphere') repaintField();
      };
      document.getElementById('schwarz-colormap').addEventListener('change', e => {
        sState.grid.colormap = e.target.value;
        if (sState.sphereView) sState.sphereView.setRenderParams({ colormap: sState.grid.colormap });
        recolor();
      });
      document.getElementById('schwarz-scalemode').addEventListener('change', e => {
        sState.grid.scaleMode = e.target.value;
        updateModKVisibility();
        if (sState.sphereView) sState.sphereView.setRenderParams({ scaleMode: sState.grid.scaleMode });
        recolor();
      });
      document.getElementById('schwarz-modk').addEventListener('change', e => {
        sState.grid.modK = Math.max(2, Math.min(64, +e.target.value || 8));
        if (sState.sphereView) sState.sphereView.setRenderParams({ modK: sState.grid.modK });
        if (sState.grid.scaleMode === 'modulo') recolor();
      });
      updateModKVisibility();
      document.getElementById('schwarz-recompute').addEventListener('click', requestRecompute);
      document.getElementById('schwarz-fit').addEventListener('click', () => {
        if (sState.viewMode === 'z') fitToDisk(); else fitToOmega();
      });
    }, 0);
    return card;
  }

  function updateModKVisibility() {
    const w = document.getElementById('schwarz-modk-wrap');
    if (w) w.style.display = (sState.grid.scaleMode === 'modulo') ? '' : 'none';
  }

  function makeInfoCard() {
    const card = document.createElement('section');
    card.className = 'card view-2d';
    card.id = 'schwarz-info-card';
    card.innerHTML = `
      <h2>Click & hover</h2>
      <div class="hint">
        <b>Single-click</b> in Ω → pin the forward σ-orbit.<br>
        <b>Double-click</b> in the tiling set → seed a preimage tree.<br>
        <b>Hover</b> → coordinates + escape time; <b>drag</b> pans, <b>wheel</b> zooms.<br>
        In the <b>z-disk</b> view, clicks map through z ↦ φ(z).
      </div>
      <div id="schwarz-readout" class="hint" style="font-family:ui-monospace,Consolas,monospace; margin-top:8px; min-height:1.4em;">
        —
      </div>
    `;
    return card;
  }

  // ---------------------------------------------------------------------------
  // Source-φ capture from the Inverse tab's state.current.primary.
  // ---------------------------------------------------------------------------
  function refreshSourceStatus() {
    const el = document.getElementById('schwarz-src-status');
    if (!el) return;
    if (sState.phiSnapshot) {
      const phi = sState.phiSnapshot;
      let famLabel;
      switch (phi.family) {
        case 'boundedLQD':            famLabel = 'bounded LQD';            break;
        case 'boundedLQD_singular':   famLabel = 'bounded singular LQD';   break;
        case 'unboundedLQD':          famLabel = 'unbounded LQD';          break;
        case 'unboundedLQD_singular': famLabel = 'unbounded singular LQD'; break;
        case 'powerQD':               famLabel = `bounded PQD (α=${phi.alpha || '?'})`; break;
        case 'powerQD_singular':      famLabel = `bounded singular PQD (α=${phi.alpha || '?'})`; break;
        case 'unboundedPQD':          famLabel = `unbounded PQD (α=${phi.alpha || '?'})`; break;
        case 'unboundedPQD_singular': famLabel = `unbounded singular PQD (α=${phi.alpha || '?'})`; break;
        default: famLabel = phi.unbounded ? 'unbounded QD' : 'bounded QD';
      }
      const polyLen = (phi.polyA || phi.F || []).length;
      const branchLen = (phi.branches || []).reduce((a, b) => a + (b.A || []).length, 0);
      const bits = [];
      bits.push('branch terms=' + branchLen);
      if (polyLen) bits.push('Laurent m=' + polyLen);
      if (phi.c != null && phi.unbounded) bits.push('c=' + phi.c);
      if (phi.z0) bits.push('z₀=' + QD.Complex.format(phi.z0));
      el.innerHTML = `<b>Captured:</b> ${famLabel}, ${bits.join(', ')}`;
      el.style.color = '#1a3e7a';
    } else {
      el.textContent = '(no φ captured)';
      el.style.color = '#777';
    }
    // Bounded warning: shown only in sphere view AND when φ is bounded
    // (because the sphere view is visually uninformative there).
    const warnEl = document.getElementById('schwarz-bounded-warning');
    if (warnEl) {
      const showWarn = sState.viewMode === 'sphere'
                    && !!sState.phiSnapshot
                    && !sState.phiSnapshot.unbounded;
      warnEl.style.display = showWarn ? '' : 'none';
    }
  }

  function captureFromInverseTab() {
    // Prefer the QD.PrimarySolution envelope (HANDOFF — P0 refactor); fall
    // back to the legacy state.current read for any unanticipated load order.
    const envelope = (QD.PrimarySolution && QD.PrimarySolution.get())
                  || (typeof state !== 'undefined' ? state.current : null);
    if (!envelope || !envelope.success) {
      alert('No successful Inverse-tab solution yet. Solve on the Inverse tab first.');
      return;
    }
    const sol = envelope.primary;
    if (!sol || !sol.phi) { alert('Inverse-tab primary solution missing φ.'); return; }
    // All six families are supported: classical bounded/unbounded (which
    // don't set phi.family — see HANDOFF gotcha #1) plus all four LQD
    // families. No allowlist to enforce.
    sState.phiSnapshot = clonePhi(sol.phi);
    // hData lives on the envelope, NOT on the primary sol — keep them
    // separate so we can read either path without surprise.
    const hData = envelope.hData;
    sState.hDataSnapshot = hData ? cloneHData(hData) : null;
    // Boundary samples: re-derive from φ via the adaptive sampler. 512
    // samples is the larger of plane (was 384) and sphere (512) defaults —
    // both views share the snapshot now (HANDOFF #29).
    try {
      sState.boundarySnapshot = QD.sampleBoundaryAdaptive
        ? QD.sampleBoundaryAdaptive(sState.phiSnapshot, 512, 800).map(p => ({re:p.w.re, im:p.w.im}))
        : QD.sampleBoundary(sState.phiSnapshot, 512);
    } catch (e) {
      alert('Failed to sample ∂Ω: ' + (e.message || e));
      return;
    }
    sState.schwarz = QD.Schwarz.buildSchwarzFromPhi(
      sState.phiSnapshot, sState.hDataSnapshot, sState.boundarySnapshot);
    sState.orbit = []; sState.pinnedOrbit = []; sState.hoverOrbit = null;
    sState.preimageTree = null;

    // Try to bring up the GPU renderer (idempotent: only created once).
    ensureGPU();
    if (sState.gpu) {
      const okGpu = sState.gpu.setPhi(sState.phiSnapshot, {
        boundaryPts: sState.boundarySnapshot,
        escapeR:     sState.schwarz.escapeR,
      });
      if (!okGpu) {
        sState.gpuMsg = sState.gpu.capacityError() || 'GPU rejected this φ.';
      } else {
        sState.gpu.setColormap(sState.grid.colormap);
        sState.gpuMsg = '';
      }
    }
    // Push to sphere view too (only if it's already been lazy-mounted).
    // If the user hasn't toggled to sphere yet, the snapshot will be pushed
    // on first activation (see _activateSphereView).
    if (sState.sphereView) {
      sState.sphereView.setPhi(sState.phiSnapshot,
                                sState.hDataSnapshot,
                                sState.boundarySnapshot);
    }

    refreshSourceStatus();
    frameDisk();                       // pre-frame the z-view so a later switch is ready
    if (sState.viewMode === 'plane') {
      fitToOmega();
    } else if (sState.viewMode === 'z') {
      requestRecompute();
    } else if (sState.sphereView) {
      sState.sphereView.requestRender();
    }
  }

  function ensureGPU() {
    if (sState.gpu) return;
    if (!QD.Schwarz.createGPURenderer) return;
    // We need our own canvas — the existing #canvas is already in 2D mode
    // for the Inverse tab and a canvas can only hold one context type.
    // We add a sibling, positioned under #canvas, that the WebGL renderer
    // paints into. #canvas keeps its 2D context and is used for overlays
    // (boundary + orbit) layered on top.
    const plotArea = document.getElementById('plot-area');
    const mainC    = getCanvas();
    if (!plotArea || !mainC) return;
    let glC = document.getElementById('schwarz-gl-canvas');
    if (!glC) {
      glC = document.createElement('canvas');
      glC.id = 'schwarz-gl-canvas';
      // Stack the GL canvas behind #canvas. #plot-area already has
      // position:relative (style.css). We give the GL canvas absolute
      // positioning + z-index 0. The 2D canvas keeps its in-flow
      // layout but receives a z-index 1 via position:relative — this
      // way z-index works WITHOUT removing it from normal flow (which
      // would break sizing and was the cause of the drag-flicker
      // "I see the QD" symptom).
      glC.style.cssText =
        'position:absolute; left:0; top:0; width:100%; height:100%; '
        + 'pointer-events:none; z-index:0;';
      mainC.style.position = 'relative';   // keep in flow but enable z-index
      mainC.style.zIndex   = '1';
      mainC.style.background = 'transparent';
      plotArea.insertBefore(glC, mainC);
      // WebGL context-loss recovery (attached once, with the canvas). On loss
      // every GL object dies, so drop the renderer → activeRenderer() falls back
      // to CPU until the context returns; on restore, recreate the renderer and
      // re-render. (The renderer's own listener calls preventDefault, which is
      // what lets the browser fire 'restored' at all.)
      glC.addEventListener('webglcontextlost', () => {
        sState.gpu = null;
        sState.gpuMsg = 'GPU context lost; using CPU until it is restored.';
      }, false);
      glC.addEventListener('webglcontextrestored', () => {
        sState.gpu = null; sState.gpuMsg = '';
        ensureGPU();
        if (sState.schwarz && isSchwarzActive() && sState.viewMode === 'plane') {
          requestRecompute();
        }
      }, false);
    }
    try {
      sState.gpu = QD.Schwarz.createGPURenderer(glC);
      if (!sState.gpu) sState.gpuMsg = 'WebGL 2 unavailable; using CPU renderer.';
    } catch (e) {
      sState.gpu = null;
      sState.gpuMsg = 'GPU init failed: ' + (e.message || e);
    }
  }
  function showGLLayer(show) {
    const glC = document.getElementById('schwarz-gl-canvas');
    if (glC) glC.style.display = show ? '' : 'none';
  }

  function activeRenderer() {
    // 'auto': prefer GPU when present & no capacity-error; else CPU.
    if (sState.grid.renderer === 'cpu') return 'cpu';
    if (sState.grid.renderer === 'gpu') return sState.gpu ? 'gpu' : 'cpu';
    return (sState.gpu && !sState.gpu.capacityError()) ? 'gpu' : 'cpu';
  }

  function clonePhi(phi) {
    // Delegate to the canonical QD.clonePhi so EVERY family-specific field
    // (alpha for PQDs, lqdBeta/lqdGamma for LQDs, z0/gamma/q, polyA, branches)
    // is carried through a single source of truth. Layer on only the UI-only
    // `F` field (Laurent companion the Schwarz inverse path attaches) that the
    // canonical clone does not know about. Maintaining a second hand-rolled
    // clone here is what silently dropped lqdBeta (HANDOFF #26) and then alpha
    // (PQD captures → NaN powers); the delegation removes that drift class.
    const out = QD.clonePhi(phi);
    out.F = phi.F ? phi.F.map(c => ({ re: c.re, im: c.im })) : undefined;
    return out;
  }
  function cloneHData(h) {
    return {
      poles: (h.poles || []).map(p => ({
        a: { re: p.a.re, im: p.a.im },
        principal: p.principal.map(c => ({ re: c.re, im: c.im })),
      })),
      polyPart: (h.polyPart || []).map(c => ({ re: c.re, im: c.im })),
    };
  }

  // ---------------------------------------------------------------------------
  // Canvas plumbing (we own the shared canvas while this tab is active).
  // ---------------------------------------------------------------------------
  function getCanvas() { return document.getElementById('canvas'); }
  function getCtx()    { const c = getCanvas(); return c ? c.getContext('2d') : null; }

  // Canvas interaction -> schwarz-interaction.js (Phase-3 item E). The drag
  // state + attachCanvasHandlers + the handlers below are installed near the
  // end of this file; renderImmediate (used by render too) stays here.
  // Synchronous GPU re-render. Used during drag/zoom in GPU mode for both the
  // plane and z-disk views (the shader's u_viewMode branch lifts z → w = φ(z)).
  function renderImmediate() {
    if (!sState.schwarz || !sState.gpu || activeRenderer() !== 'gpu') return;
    // Belt-and-suspenders: never re-show the GL layer / paint when the Schwarz
    // tab isn't active (e.g. a late control event after a tab switch).
    if (!isSchwarzActive()) return;
    const inZ = sState.viewMode === 'z';
    showGLLayer(true);
    try {
      sState.gpu.setColormap(sState.grid.colormap);
      sState.gpu.render(inZ ? sState.zView : sState.view, {
        maxIter:   sState.grid.maxIter,
        scaleMode: sState.grid.scaleMode,
        modK:      sState.grid.modK,
        viewMode:  inZ ? 'z' : 'w',
      });
      if (inZ) {
        paintZView(true);
      } else {
        paintBoundaryOnTop();
        paintOrbit();
      }
    } catch (e) {
      // Fall through silently; the next debounced recompute will surface
      // any persistent error.
    }
  }

  // Canvas interaction (clearOverlay / isSchwarzActive / wheel / mousemove /
  // hover / click / dblclick / pin) -> schwarz-interaction.js (Phase-3 item E).

  // ---------------------------------------------------------------------------
  // Coordinate transforms.
  // ---------------------------------------------------------------------------
  function pixelToWorld(sx, sy) {
    const { cx, cy, scale, cssW, cssH } = sState.view;
    return {
      re: cx + (sx - cssW / 2) / scale,
      im: cy - (sy - cssH / 2) / scale,
    };
  }
  function worldToPixel(re, im) {
    const { cx, cy, scale, cssW, cssH } = sState.view;
    return {
      x: cssW / 2 + (re - cx) * scale,
      y: cssH / 2 - (im - cy) * scale,
    };
  }
  // z-plane (unit-disk) transforms — parallel to pixelToWorld/worldToPixel but
  // reading sState.zView, so the z-view has its own independent pan/zoom.
  function pixelToZ(sx, sy) {
    const { cx, cy, scale, cssW, cssH } = sState.zView;
    return {
      re: cx + (sx - cssW / 2) / scale,
      im: cy - (sy - cssH / 2) / scale,
    };
  }
  function zToPixel(re, im) {
    const { cx, cy, scale, cssW, cssH } = sState.zView;
    return {
      x: cssW / 2 + (re - cx) * scale,
      y: cssH / 2 - (im - cy) * scale,
    };
  }
  function fitToOmega() {
    if (!sState.boundarySnapshot || !sState.boundarySnapshot.length) return;
    const b = QD.Schwarz.polygonBounds(sState.boundarySnapshot);
    sState.view.cx = b.center.re;
    sState.view.cy = b.center.im;
    syncCanvasSize();
    const margin = sState.phiSnapshot && sState.phiSnapshot.unbounded ? 2.2 : 1.25;
    sState.view.scale = Math.min(sState.view.cssW, sState.view.cssH) / (2 * b.radius * margin);
    requestRecompute();
  }
  // Frame the z-view on the unit disk (bounded: 𝔻 fills ~82% of the shorter
  // half; unbounded: show |z| up to ~2.5 so 𝔻* is visible). Sets the transform
  // only; fitToDisk also kicks a recompute.
  function frameDisk() {
    syncCanvasSize();
    sState.zView.cx = 0;
    sState.zView.cy = 0;
    const unb  = !!(sState.phiSnapshot && sState.phiSnapshot.unbounded);
    const half = Math.min(sState.zView.cssW, sState.zView.cssH) / 2;
    sState.zView.scale = unb ? (half * 0.9) / 2.5 : half * 0.82;
  }
  function fitToDisk() { frameDisk(); requestRecompute(); }
  function syncCanvasSize() {
    const c = getCanvas();
    if (!c) return;
    const rect = c.getBoundingClientRect();
    sState.view.cssW = Math.max(50, rect.width);
    sState.view.cssH = Math.max(50, rect.height);
    // The z-view shares the same physical canvas, so its css dims track view's.
    sState.zView.cssW = sState.view.cssW;
    sState.zView.cssH = sState.view.cssH;
  }

  // True when the Schwarz tab is the active tab (its controls panel is shown).
  // Single predicate for the render-entry / resize / context-restore guards, so
  // a late timer / global event can't paint into another tab. (Hoisted, so the
  // resize listener above and the installSchwarzRender injection below can both
  // use it regardless of source order.)
  function isSchwarzActive() {
    const p = document.getElementById('controls-schwarz');
    return !!p && !p.hidden;
  }

  // ---------------------------------------------------------------------------
  // Progressive renderer.
  // ---------------------------------------------------------------------------
  // Progressive renderer -> schwarz-render.js (Phase-3 item E). requestRecompute
  // (the debounced entry) + doRecompute + the CPU pyramid (_renderCpuPyramid /
  // _renderCpuViaWorker / chainPass / fillFromCoarseSamples) are installed after
  // the paint install below; the rest of this file calls requestRecompute by the
  // forward-declared name. The module reads sState + the paint fns + a few
  // GPU/geometry helpers via sCtx.

  // ---------------------------------------------------------------------------
  // Painting.
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Painting + colormaps -> schwarz-paint.js (Phase-3 item E). The 2D-canvas
  // rendering-output layer is installed here; the rest of this file calls the
  // captured names (forward-declared near the top). The module reads sState +
  // getCtx / syncCanvasSize / worldToPixel + the KIND_* constants via sCtx.
  // ---------------------------------------------------------------------------
  ({
    clearCanvas, paintAll, repaintField, paintBoundaryOnTop, paintOrbit,
    paintPreimageTree, paintLimitSet, paintZView, setProgress,
  } = window.QD_UI.installSchwarzPaint({
    sState, getCtx, syncCanvasSize, worldToPixel, zToPixel, activeRenderer,
    KIND_FUND, KIND_ESC, KIND_INT, KIND_INV, KIND_OUTSIDE,
  }));

  // Progressive renderer (installed after paint so its paint deps are on sCtx).
  ({ requestRecompute } = window.QD_UI.installSchwarzRender({
    sState, clearCanvas, paintAll, paintBoundaryOnTop, paintOrbit, paintZView,
    setProgress, syncCanvasSize, activeRenderer, showGLLayer, isSchwarzActive,
    KIND_FUND, KIND_ESC, KIND_INT, KIND_INV, KIND_OUTSIDE,
  }));

  // Feature-compute methods (installed after paint+render so their paint deps
  // are on sCtx; before interaction, which destructures the recompute hooks).
  ({
    _recomputeDomainColoring, _rebuildPreimageTreeIfActive, _refreshPreimageTreeStats,
    _computeLimitSet, _clearLimitSet, _recomputeCriticalOrbits, _findCycles,
    _exportPng, _recomputeZPanelOrbit, _computeSweep, _recomputeLevelCurves,
  } = window.QD_UI.installSchwarzFeatures({
    sState, paintBoundaryOnTop, paintPreimageTree, paintLimitSet,
    activeRenderer, getCtx, getCanvas,
  }));

  // Canvas interaction (installed after paint+render so its renderer/paint deps
  // are available; reads the per-feature recompute hooks + geometry via sCtx).
  _schwarzInter = window.QD_UI.installSchwarzInteraction({
    sState, getCanvas, getCtx, pixelToWorld, worldToPixel, pixelToZ, zToPixel,
    syncCanvasSize, activeRenderer, renderImmediate, requestRecompute,
    paintBoundaryOnTop, paintOrbit, paintAll, paintPreimageTree,
    gateMaxIter, _recomputeLevelCurves, _recomputeDomainColoring,
    _recomputeZPanelOrbit, _refreshPreimageTreeStats,
    KIND_FUND, KIND_ESC, KIND_INT, KIND_INV, KIND_OUTSIDE,
  });
  ({ attachCanvasHandlers, onCanvasClick, onCanvasDblClick, onMouseMove,
     runHoverOrbit, pinOrbitAt } = _schwarzInter);

  // ---------------------------------------------------------------------------
  // Test-only hook (see node-test.js). Opt-in via a window sentinel so a normal
  // browser load NEVER attaches it. Exposes the fractal-mode interaction
  // handlers + state so the click/dblclick disambiguation and the tiling-set
  // seed gate can be unit-tested without mounting the sidebar.
  // ---------------------------------------------------------------------------
  if (typeof window !== 'undefined' && window.__SCHWARZ_UI_TEST_HOOK__) {
    window.__schwarzUiTest = {
      sState, setMode, onCanvasClick, onCanvasDblClick, onMouseMove,
      runHoverOrbit, pinOrbitAt,
      get CLICK_DELAY() { return _schwarzInter.getClickDelay(); },
      set CLICK_DELAY(v) { _schwarzInter.setClickDelay(v); },
    };
  }

})();
