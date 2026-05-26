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
//   • CPU (fallback): progressive 4×4 → 2×2 → 1×1 pyramid on the main thread,
//     chunked across requestAnimationFrame ticks (~14 ms per slice) so the
//     UI stays responsive. Per-pixel warm seeds from the left neighbor in
//     raster order → Newton converges in 1-3 iterations on average.
//
// Layer stacking: the GPU canvas (#schwarz-gl-canvas) is inserted behind the
// main #canvas. The 2D canvas keeps its in-flow layout but gets z-index 1 +
// position:relative so it can host transparent overlays (boundary + orbit)
// on top of the GPU pixels. The GL layer hides on tab-out.
//
// Double-click-to-orbit: plots {w₀, σ(w₀), σ²(w₀), …} as a polyline of small
// dots. Single click is reserved for the start of a drag (dragMoved guard
// suppresses orbit drops after pans).
//
// Hover readout: pixel coords + escape time (CPU mode only — GPU doesn't
// keep a per-pixel field array; hover shows coords only).
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

    orbit: [],                             // last-clicked orbit polyline

    // View toggle (HANDOFF #29): 'plane' = Schwarz dynamics on the w-plane
    // (the original Schwarz tab), 'sphere' = same iteration textured onto a
    // Riemann sphere. The sphere renderer is lazy-mounted via QD.SphereView
    // on first switch to sphere mode.
    viewMode:   'plane',                   // 'plane' | 'sphere'
    sphereView: null,                      // QD.SphereView handle
  };

  // Kinds enum
  const KIND_FUND = 0, KIND_ESC = 1, KIND_INT = 2, KIND_INV = 3, KIND_OUTSIDE = 4;

  // ---------------------------------------------------------------------------
  // Lazy mount
  // ---------------------------------------------------------------------------
  document.addEventListener('tab-changed', function (e) {
    if (!e.detail || e.detail.tab !== 'schwarz') {
      // Leaving the Schwarz tab — don't keep rendering, and hide BOTH GL
      // layers so they can't show through under another tab's drawing.
      sState.renderToken++;
      showGLLayer(false);
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
    if (sState.viewMode === 'plane') {
      if (sState.sphereView) sState.sphereView.deactivate();
      if (sState.schwarz) {
        showGLLayer(activeRenderer() === 'gpu');
        requestRecompute();
      } else {
        showGLLayer(false);
        clearCanvas();
      }
    } else {
      // sphere mode: hide the Schwarz GL layer, activate sphere view.
      showGLLayer(false);
      _activateSphereView();
    }
  });

  function mountSchwarzSidebar() {
    const root = document.getElementById('controls-schwarz');
    if (!root) return;
    root.innerHTML = '';
    root.appendChild(makeViewToggleCard());
    root.appendChild(makeSourceCard());
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
    // The view-toggle card has no h2; the others do.
    const headers = root.querySelectorAll('section.card h2');
    if (headers[0]) H(headers[0],
      `<b>Source φ.</b> The Schwarz dynamics tab iterates σ(w) = φ(1/φ⁻¹(w)),
       the Schwarz reflection associated with Ω. Capture a φ from the QD tab
       (after solving). Each pixel is colored by escape time of σ-iteration.`);
    if (headers[1]) H(headers[1],
      `<b>Render.</b> CPU pyramid uses a Worker that progressively refines
       resolution; GPU uses a WebGL 2 fragment shader for instant frames.
       <i>Colormap</i> + <i>scaleMode</i> change the escape-time → colour
       mapping; <i>maxIter</i> caps the σ-iteration before declaring a pixel
       interior; <i>mod k</i> emphasises orbit-period structure.`);
    if (headers[2]) H(headers[2],
      `<b>Click & hover.</b> Click pixels to trace and overlay individual σ
       orbits. Hover to read the w-plane coordinate + escape time (and, in
       CPU mode, the pixel kind). In plane view, drag to pan, scroll to zoom.`);
  }

  function makeViewToggleCard() {
    const card = document.createElement('section');
    card.className = 'card';
    card.innerHTML = `
      <div class="segmented" role="tablist" aria-label="View mode">
        <button class="seg-btn active" data-view="plane" type="button">plane</button>
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

  function setViewMode(mode) {
    if (mode !== 'plane' && mode !== 'sphere') return;
    if (mode === sState.viewMode) return;
    sState.viewMode = mode;
    // Update segmented-control highlight.
    document.querySelectorAll('#controls-schwarz .seg-btn').forEach(btn => {
      btn.classList.toggle('active', btn.dataset.view === mode);
    });
    _applyViewModeVisibility();
    if (mode === 'plane') {
      if (sState.sphereView) sState.sphereView.deactivate();
      if (sState.schwarz) {
        showGLLayer(activeRenderer() === 'gpu');
        requestRecompute();
      } else {
        showGLLayer(false);
        clearCanvas();
      }
    } else {
      showGLLayer(false);
      _activateSphereView();
    }
    refreshSourceStatus();
  }

  function _applyViewModeVisibility() {
    const planeShow  = sState.viewMode === 'plane';
    const sphereShow = sState.viewMode === 'sphere';
    document.querySelectorAll('#controls-schwarz .view-plane-only')
      .forEach(el => { el.style.display = planeShow ? '' : 'none'; });
    document.querySelectorAll('#controls-schwarz .view-sphere-only')
      .forEach(el => { el.style.display = sphereShow ? '' : 'none'; });
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
    card.innerHTML = `
      <h2>Source φ (from Inverse tab)</h2>
      <div class="hint">
        The Schwarz reflection σ(w) = conj(F(ψ(w))) is built from a Riemann
        map φ produced by the inverse solver. Solve on the Inverse tab,
        then click <b>Use this φ</b> to snapshot it here.
        <br>All six inverse families supported (classical bounded / unbounded,
        and all four LQD variants).
      </div>
      <div id="schwarz-src-status" class="hint" style="color:#333; margin-top:8px;">
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
    card.innerHTML = `
      <h2>Render</h2>
      <div class="row view-plane-only">
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
      <div class="row view-plane-only" style="margin-top:8px;">
        <label>Renderer:
          <select id="schwarz-renderer">
            <option value="auto" selected>auto (GPU if available)</option>
            <option value="gpu">GPU (WebGL 2)</option>
            <option value="cpu">CPU (fallback)</option>
          </select>
        </label>
      </div>
      <div class="row view-plane-only" style="margin-top:10px;">
        <button id="schwarz-recompute" class="small">Recompute</button>
        <button id="schwarz-fit" class="small" style="margin-left:6px;">Fit to Ω</button>
      </div>
      <div id="schwarz-progress" class="hint view-plane-only" style="margin-top:8px; min-height:1.2em;"></div>
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
        if (sState.viewMode === 'plane') requestRecompute();
      });
      document.getElementById('schwarz-colormap').addEventListener('change', e => {
        sState.grid.colormap = e.target.value;
        if (sState.sphereView) sState.sphereView.setRenderParams({ colormap: sState.grid.colormap });
        if (sState.viewMode === 'plane') {
          // For GPU, just re-render (very fast). For CPU, repaint existing field.
          if (activeRenderer() === 'gpu') renderImmediate();
          else repaintField();
        }
      });
      document.getElementById('schwarz-scalemode').addEventListener('change', e => {
        sState.grid.scaleMode = e.target.value;
        updateModKVisibility();
        if (sState.sphereView) sState.sphereView.setRenderParams({ scaleMode: sState.grid.scaleMode });
        if (sState.viewMode === 'plane') {
          if (activeRenderer() === 'gpu') renderImmediate();
          else repaintField();
        }
      });
      document.getElementById('schwarz-modk').addEventListener('change', e => {
        sState.grid.modK = Math.max(2, Math.min(64, +e.target.value || 8));
        if (sState.sphereView) sState.sphereView.setRenderParams({ modK: sState.grid.modK });
        if (sState.viewMode === 'plane' && sState.grid.scaleMode === 'modulo') {
          if (activeRenderer() === 'gpu') renderImmediate();
          else repaintField();
        }
      });
      updateModKVisibility();
      document.getElementById('schwarz-recompute').addEventListener('click', requestRecompute);
      document.getElementById('schwarz-fit').addEventListener('click', fitToOmega);
    }, 0);
    return card;
  }

  function updateModKVisibility() {
    const w = document.getElementById('schwarz-modk-wrap');
    if (w) w.style.display = (sState.grid.scaleMode === 'modulo') ? '' : 'none';
  }

  function makeInfoCard() {
    const card = document.createElement('section');
    card.className = 'card view-plane-only';
    card.innerHTML = `
      <h2>Click & hover</h2>
      <div class="hint">
        <b>Double-click</b> on Ω → plot the orbit {w₀, σ(w₀), σ²(w₀), …}<br>
        <b>Hover</b> → show pixel coords + escape time.<br>
        <b>Drag</b> to pan; <b>wheel</b> to zoom.
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
    if (typeof state === 'undefined' || !state.current || !state.current.success) {
      alert('No successful Inverse-tab solution yet. Solve on the Inverse tab first.');
      return;
    }
    const sol = state.current.primary;
    if (!sol || !sol.phi) { alert('Inverse-tab primary solution missing φ.'); return; }
    // All six families are supported: classical bounded/unbounded (which
    // don't set phi.family — see HANDOFF gotcha #1) plus all four LQD
    // families. No allowlist to enforce.
    sState.phiSnapshot = clonePhi(sol.phi);
    // hData lives on state.current, NOT on the primary sol — keep them
    // separate so we can read either path without surprise.
    const hData = state.current.hData;
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
    sState.orbit = [];

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
    if (sState.viewMode === 'plane') {
      fitToOmega();
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
    // Deep-enough copy for our purposes. The inverse solver uses `polyA` for
    // the Laurent part on unbounded; some other paths use `F`. Carry both.
    // `lqdBeta` (polynomial-h β-correction, HANDOFF #22) and `lqdGamma`
    // (higher-order pole at origin, HANDOFF #24) must also be preserved so
    // the Schwarz adapters can evaluate the full φ.
    const out = {
      family:    phi.family,
      unbounded: phi.unbounded,
      w0:        phi.w0 ? { re: phi.w0.re, im: phi.w0.im } : undefined,
      c:         phi.c,
      q:         phi.q     ? { re: phi.q.re,     im: phi.q.im }     : undefined,
      gamma:     phi.gamma ? { re: phi.gamma.re, im: phi.gamma.im } : undefined,
      z0:        phi.z0    ? { re: phi.z0.re,    im: phi.z0.im }    : undefined,
      polyA:     phi.polyA ? phi.polyA.map(c => ({ re: c.re, im: c.im })) : undefined,
      F:         phi.F     ? phi.F    .map(c => ({ re: c.re, im: c.im })) : undefined,
      branches:  phi.branches ? phi.branches.map(b => ({
        z: { re: b.z.re, im: b.z.im },
        A: b.A.map(a => ({ re: a.re, im: a.im })),
      })) : [],
      lqdBeta:   phi.lqdBeta  ? phi.lqdBeta .map(c => ({ re: c.re, im: c.im })) : [],
      lqdGamma:  phi.lqdGamma ? phi.lqdGamma.map(c => ({ re: c.re, im: c.im })) : [],
    };
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

  let dragging = false, dragMoved = false, lastX = 0, lastY = 0;
  function attachCanvasHandlers() {
    const c = getCanvas();
    if (!c) return;
    c.addEventListener('mousemove', onMouseMove);
    c.addEventListener('mouseleave', () => {
      const r = document.getElementById('schwarz-readout');
      if (r) r.textContent = '—';
    });
    // Orbit on DOUBLE click only. Single click is reserved for click-and-drag
    // panning (and for not adding an unwanted orbit).
    c.addEventListener('dblclick', onDoubleClickOrbit);
    c.addEventListener('wheel', onWheel, { passive: false });
    c.addEventListener('mousedown', e => {
      if (e.button !== 0) return;
      dragging = true; dragMoved = false;
      lastX = e.clientX; lastY = e.clientY;
      c.style.cursor = 'grabbing';
    });
    window.addEventListener('mousemove', e => {
      if (!dragging || !isSchwarzActive()) return;
      const dx = e.clientX - lastX, dy = e.clientY - lastY;
      if (dx !== 0 || dy !== 0) dragMoved = true;
      lastX = e.clientX; lastY = e.clientY;
      sState.view.cx -= dx / sState.view.scale;
      sState.view.cy += dy / sState.view.scale;          // screen y is flipped
      // GPU is fast enough (10-30 ms typical) to render every mousemove
      // without debounce. CPU mode debounces because the pyramid is slow.
      if (activeRenderer() === 'gpu') renderImmediate();
      else { clearOverlay(); requestRecompute(); }
    });
    window.addEventListener('mouseup', () => {
      if (!dragging) return;
      dragging = false;
      c.style.cursor = '';
      // After a real drag (mouse moved), trigger a fresh render so the
      // final position is sharp even in CPU mode.
      if (dragMoved && activeRenderer() !== 'gpu') requestRecompute();
    });
  }

  // Synchronous GPU re-render. Used during drag/zoom in GPU mode.
  function renderImmediate() {
    if (!sState.schwarz || !sState.gpu || activeRenderer() !== 'gpu') return;
    showGLLayer(true);
    try {
      sState.gpu.setColormap(sState.grid.colormap);
      sState.gpu.render(sState.view, {
        maxIter:   sState.grid.maxIter,
        scaleMode: sState.grid.scaleMode,
        modK:      sState.grid.modK,
      });
      paintBoundaryOnTop();
      paintOrbit();
    } catch (e) {
      // Fall through silently; the next debounced recompute will surface
      // any persistent error.
    }
  }

  function clearOverlay() {
    const ctx = getCtx(); if (!ctx) return;
    syncCanvasSize();
    ctx.clearRect(0, 0, sState.view.cssW, sState.view.cssH);
  }
  function isSchwarzActive() {
    const panel = document.getElementById('controls-schwarz');
    return panel && !panel.hidden;
  }
  function onWheel(e) {
    if (!isSchwarzActive() || !sState.schwarz) return;
    e.preventDefault();
    const c = getCanvas();
    const rect = c.getBoundingClientRect();
    const sx = e.clientX - rect.left, sy = e.clientY - rect.top;
    const w = pixelToWorld(sx, sy);
    const k = (e.deltaY > 0) ? 0.85 : 1.18;
    sState.view.scale *= k;
    // Keep the world point under the cursor pinned in screen space.
    const after = pixelToWorld(sx, sy);
    sState.view.cx += w.re - after.re;
    sState.view.cy += w.im - after.im;
    if (activeRenderer() === 'gpu') renderImmediate();
    else { clearOverlay(); requestRecompute(); }
  }
  function onMouseMove(e) {
    if (!isSchwarzActive() || !sState.schwarz) return;
    const c = getCanvas();
    const rect = c.getBoundingClientRect();
    const sx = e.clientX - rect.left, sy = e.clientY - rect.top;
    const w = pixelToWorld(sx, sy);
    let info = `w = (${w.re.toFixed(3)}, ${w.im.toFixed(3)})`;
    if (sState.field && sState.fieldW > 0) {
      const gx = Math.floor(sx / sState.view.cssW  * sState.fieldW);
      const gy = Math.floor(sy / sState.view.cssH  * sState.fieldH);
      if (gx >= 0 && gx < sState.fieldW && gy >= 0 && gy < sState.fieldH) {
        const idx = gy * sState.fieldW + gx;
        const n = sState.field[idx];
        const kind = sState.fieldKind ? sState.fieldKind[idx] : KIND_OUTSIDE;
        info += '  ' + describeKind(kind, n);
      }
    } else if (activeRenderer() === 'gpu' && QD.Schwarz && QD.Schwarz.escapeTime) {
      // GPU-mode parity (HANDOFF #33): the field array isn't populated in
      // GPU mode, so do an ad-hoc per-cursor CPU iteration. Cheap — at most
      // `maxIter` σ-evals (μs scale on typical scenarios).
      try {
        const et = QD.Schwarz.escapeTime(w, sState.schwarz, { maxIter: sState.grid.maxIter });
        // escapeTime returns kind as a string; map to the KIND_* enum
        // used by describeKind. Note: a pixel that was already outside Ω
        // returns kind='fundamental' n=0 — display it as KIND_OUTSIDE.
        const kindMap = { fundamental: KIND_FUND, escaped: KIND_ESC,
                          interior: KIND_INT, invalid: KIND_INV };
        let kindI = (et && kindMap[et.kind] != null) ? kindMap[et.kind] : KIND_OUTSIDE;
        if (kindI === KIND_FUND && (et.n | 0) === 0) kindI = KIND_OUTSIDE;
        info += '  ' + describeKind(kindI, et ? (et.n | 0) : 0);
      } catch (err) { /* swallow; coordinate readout still shown */ }
    }
    const r = document.getElementById('schwarz-readout');
    if (r) r.textContent = info;
  }
  function describeKind(kind, n) {
    switch (kind) {
      case KIND_FUND:    return 'escape time n=' + n;
      case KIND_ESC:     return 'in escaping set';
      case KIND_INT:     return 'still in Ω after maxIter (tiling-set interior)';
      case KIND_INV:     return 'Newton diverged';
      case KIND_OUTSIDE: return 'in Ω^c (fundamental tile)';
      default:           return '';
    }
  }
  function onDoubleClickOrbit(e) {
    if (!isSchwarzActive() || !sState.schwarz) return;
    // dblclick can fire after the user double-clicks at the end of a drag —
    // guard against that with the dragMoved sentinel.
    if (dragMoved) return;
    const c = getCanvas();
    const rect = c.getBoundingClientRect();
    const sx = e.clientX - rect.left, sy = e.clientY - rect.top;
    const w = pixelToWorld(sx, sy);
    if (!sState.schwarz.isInOmega(w)) {
      sState.orbit = [];
    } else {
      sState.orbit = QD.Schwarz.makeOrbit(w, sState.schwarz, { maxIter: sState.grid.maxIter });
    }
    // Just redraw the overlay; the GPU fractal layer doesn't need re-render.
    if (activeRenderer() === 'gpu') {
      paintBoundaryOnTop();
      paintOrbit();
    } else {
      paintAll();
    }
  }

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
  function syncCanvasSize() {
    const c = getCanvas();
    if (!c) return;
    const rect = c.getBoundingClientRect();
    sState.view.cssW = Math.max(50, rect.width);
    sState.view.cssH = Math.max(50, rect.height);
  }

  // ---------------------------------------------------------------------------
  // Progressive renderer.
  // ---------------------------------------------------------------------------
  let recomputeTimer = null;
  function requestRecompute() {
    if (recomputeTimer) clearTimeout(recomputeTimer);
    recomputeTimer = setTimeout(() => { recomputeTimer = null; doRecompute(); }, 80);
  }

  function doRecompute() {
    if (!sState.schwarz) { clearCanvas(); return; }
    syncCanvasSize();

    // GPU path: synchronous, complete in one frame.
    if (activeRenderer() === 'gpu') {
      showGLLayer(true);
      const t0 = performance.now();
      try {
        sState.gpu.setColormap(sState.grid.colormap);
        sState.gpu.render(sState.view, {
          maxIter:   sState.grid.maxIter,
          scaleMode: sState.grid.scaleMode,
          modK:      sState.grid.modK,
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
        // Boundary + orbit overlays drawn on top (no field clearing needed).
        paintBoundaryOnTop();
        paintOrbit();
        const ms = (performance.now() - t0).toFixed(0);
        setProgress('GPU render: ' + ms + ' ms' + (sState.gpuMsg ? '  (' + sState.gpuMsg + ')' : ''));
        return;
      }
    }

    // CPU progressive pyramid path. Hide the GL layer so a stale GPU image
    // doesn't peek through edge cases.
    showGLLayer(false);
    const myToken = ++sState.renderToken;
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

    // Pass 1: every 4th pixel → fill 4×4 blocks.
    // Pass 2: every 2nd pixel → fill 2×2 blocks (only un-resolved cells).
    // Pass 3: per-pixel.
    chainPass(myToken, 4, () =>
      chainPass(myToken, 2, () =>
        chainPass(myToken, 1, () => {
          if (myToken !== sState.renderToken) return;
          sState.rendering = false;
          setProgress('');
          paintAll();
        })));
  }

  function chainPass(token, stride, next) {
    if (token !== sState.renderToken) return;
    setProgress('Pass ' + (4 / stride | 0) + (stride === 1 ? '/3 (full)…' : '/3 (refining)…'));
    const W = sState.fieldW, H = sState.fieldH;
    const sw = sState.schwarz;
    const maxIter = sState.grid.maxIter;
    // Map field coords → world.
    const cssW = sState.view.cssW, cssH = sState.view.cssH;
    const cx = sState.view.cx, cy = sState.view.cy, scale = sState.view.scale;
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
          const wRe = cx + (px - cssW / 2) / scale;
          const wIm = cy - (py - cssH / 2) / scale;
          const wpt = { re: wRe, im: wIm };
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

  // ---------------------------------------------------------------------------
  // Painting.
  // ---------------------------------------------------------------------------
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
    paintField();
    paintBoundary();
    paintOrbit();
  }
  function repaintField() { if (sState.field) paintAll(); }

  // Used after a GPU render: WebGL has already drawn the fractal to the
  // sibling #schwarz-gl-canvas. We clear the main 2D canvas to transparent
  // and draw only the boundary + orbit overlays on top.
  function paintBoundaryOnTop() {
    const ctx = getCtx(); if (!ctx) return;
    ctx.clearRect(0, 0, sState.view.cssW, sState.view.cssH);
    paintBoundary();
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
    for (let i = 0; i < W * H; i++) {
      const kind = sState.fieldKind[i];
      const n    = sState.field[i];
      let r = 0, g = 0, b = 0;
      if (kind === KIND_OUTSIDE + 1)        { r = 245; g = 245; b = 248; }     // fundamental tile
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

  function paintOrbit() {
    const pts = sState.orbit;
    if (!pts || pts.length === 0) return;
    const ctx = getCtx();
    // Connecting line.
    ctx.strokeStyle = 'rgba(20, 160, 60, 0.95)';
    ctx.lineWidth = 1.3;
    ctx.beginPath();
    for (let i = 0; i < pts.length; i++) {
      const p = worldToPixel(pts[i].re, pts[i].im);
      if (i === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
    }
    ctx.stroke();
    // Dots.
    for (let i = 0; i < pts.length; i++) {
      const p = worldToPixel(pts[i].re, pts[i].im);
      ctx.beginPath();
      ctx.arc(p.x, p.y, i === 0 ? 4.5 : 2.8, 0, 2 * Math.PI);
      ctx.fillStyle = i === 0 ? '#108a40' : 'rgba(20, 160, 60, 0.85)';
      ctx.fill();
      if (i === 0) {
        ctx.strokeStyle = '#fff';
        ctx.lineWidth = 1.2;
        ctx.stroke();
      }
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

})();
