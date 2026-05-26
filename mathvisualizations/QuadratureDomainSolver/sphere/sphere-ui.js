// =============================================================================
// sphere-ui.js  — Sphere VIEW ADAPTER (mounted by the Schwarz tab).
//
// As of HANDOFF #29 the Riemann-sphere tab no longer exists; the sphere is a
// view-mode toggle inside the Schwarz dynamics tab. This module exposes
// `QD.SphereView.mount(opts)` — the Schwarz tab calls it once when the user
// first switches to sphere view, and uses the returned handle to drive the
// sphere renderer (activate/deactivate/setPhi/setRenderParams/...).
//
// Layer stacking (unchanged from the standalone sphere tab):
//   The sphere GL canvas (#sphere-gl-canvas) is at z-index 2, ON TOP of the
//   main #canvas (z-index 1), with pointer-events:auto. The sphere "owns" all
//   pointer input while sphere-view is active. On deactivate the GL canvas is
//   hidden (display:none) so the underlying Schwarz GL canvas (z-index 0) and
//   the 2D overlay canvas (z-index 1) become interactive again.
//
// Camera: orbit (left-drag), zoom (wheel), reset (double-click / Reset btn).
// Default view: azimuth 45°, elevation 30°, distance 2.5 — looks from the
// "northeast-above" direction so both the north pole (∞, top) and the
// compact set K (bottom) are visible simultaneously.
//
// Hover readout: raycasts against the unit sphere, converts the hit point via
// inverse stereographic to w ∈ C, and displays (x,y,z) + w in a floating
// tooltip (#sphere-hover-tooltip) plus a sidebar mirror (#sp-hover).
// =============================================================================

(function () {
  'use strict';
  if (typeof QD === 'undefined') return;

  const DEF_CAMERA = { azimuth: Math.PI / 4, elevation: Math.PI / 6, distance: 2.5 };

  // Mount a sphere view into the given host. Returns null if WebGL 2 / the
  // sphere renderer is unavailable on this browser.
  function mount(opts) {
    opts = opts || {};
    const plotArea   = opts.plotArea   || document.getElementById('plot-area');
    const mainCanvas = opts.mainCanvas || document.getElementById('canvas');
    const sidebar    = opts.sidebar;            // where the display+camera cards go
    const isActive   = opts.isActive  || (() => true);
    if (!plotArea || !mainCanvas) return null;

    // ---- View state (private to this mount instance) -----------------------
    const state = {
      phiSnapshot:      null,
      boundarySnapshot: null,
      polesSnapshot:    [],
      renderer:         null,
      glCanvas:         null,
      frameRequested:   false,
      camera: {
        azimuth:   DEF_CAMERA.azimuth,
        elevation: DEF_CAMERA.elevation,
        distance:  DEF_CAMERA.distance,
      },
      params: {
        maxIter:   64,
        colormap:  'magma',
        scaleMode: 'smooth',
        modK:      8,
        texSize:   1024,
      },
      display: {
        rimDarken:     0.30,
        showBoundary:  true,
        showPoles:     true,
        showNorthPole: true,
      },
      drag: { active: false, lastX: 0, lastY: 0 },
    };

    // ---- Build display + camera sidebar cards ------------------------------
    // The host owns the source-φ + render-params cards; this adapter only
    // contributes the sphere-specific display toggles + camera-reset card.
    if (sidebar) {
      sidebar.appendChild(_makeDisplayCard(state));
      sidebar.appendChild(_makeCameraCard(state));
    }

    // ---- Create GL canvas + renderer ---------------------------------------
    const glC = _ensureGLCanvas(plotArea, mainCanvas);
    state.glCanvas = glC;
    _ensureHoverTooltip(plotArea);
    _attachGLCanvasEvents(glC, state, isActive);
    _attachWindowHandlers(state, isActive);
    _ensureResizeObserver(glC, state, isActive);

    if (!QD.Sphere || !QD.Sphere.createRenderer) {
      console.warn('sphere-ui: QD.Sphere.createRenderer unavailable.');
      return _disabledHandle();
    }
    try {
      state.renderer = QD.Sphere.createRenderer(glC);
    } catch (e) {
      console.error('sphere-ui: renderer init failed:', e);
      state.renderer = null;
    }
    if (!state.renderer) return _disabledHandle();

    // Apply the initial render+display params so the placeholder sphere shows
    // up correctly even before any φ is captured.
    state.renderer.setRenderParams(state.params);
    state.renderer.setDisplayParams(state.display);

    // ---- Handle exposed to the host ----------------------------------------
    return {
      isAvailable: () => true,
      activate() {
        _showGLLayer(glC, true);
        _requestRender(state);
      },
      deactivate() {
        _showGLLayer(glC, false);
        _hideHoverTooltip();
        state.frameRequested = false;   // any RAF in flight will no-op
      },
      setPhi(phi, hData, boundaryPts) {
        state.phiSnapshot      = phi;
        state.boundarySnapshot = boundaryPts;
        state.polesSnapshot    = _extractPoles(hData);

        // escape radius shared with Schwarz CPU/GPU path
        let escapeR;
        try {
          const sw = QD.Schwarz && QD.Schwarz.buildSchwarzFromPhi
            ? QD.Schwarz.buildSchwarzFromPhi(phi, hData, boundaryPts)
            : null;
          escapeR = sw ? sw.escapeR : undefined;
        } catch (e) { escapeR = undefined; }

        const ok = state.renderer.setPhi(phi, {
          boundaryPts: boundaryPts || [],
          escapeR,
        });
        if (!ok) {
          console.warn('sphere-ui: setPhi rejected (GPU capacity exceeded).');
          return false;
        }
        if (state.polesSnapshot.length && state.renderer.setPolePts) {
          state.renderer.setPolePts(state.polesSnapshot);
        }
        state.renderer.setRenderParams(state.params);
        state.renderer.setDisplayParams(state.display);
        _requestRender(state);
        return true;
      },
      setRenderParams(p) {
        Object.assign(state.params, p);
        if (state.renderer) {
          state.renderer.setRenderParams(state.params);
          _requestRender(state);
        }
      },
      setDisplayParams(d) {
        Object.assign(state.display, d);
        if (state.renderer) {
          state.renderer.setDisplayParams(state.display);
          _requestRender(state);
        }
      },
      resetCamera() { _resetCamera(state); },
      requestRender() { _requestRender(state); },
      markFractalDirty() {
        if (state.renderer && state.renderer.markFractalDirty) {
          state.renderer.markFractalDirty();
          _requestRender(state);
        }
      },
      destroy() {
        _showGLLayer(glC, false);
        _hideHoverTooltip();
        if (state.renderer && state.renderer.destroy) state.renderer.destroy();
        state.renderer = null;
      },
    };
  }

  // Returned when WebGL 2 / the renderer is unavailable. The host can detect
  // this via isAvailable() === false and surface a graceful error to the user.
  function _disabledHandle() {
    const noop = () => {};
    return {
      isAvailable: () => false,
      activate: noop, deactivate: noop,
      setPhi: () => false,
      setRenderParams: noop, setDisplayParams: noop,
      resetCamera: noop, requestRender: noop,
      markFractalDirty: noop, destroy: noop,
    };
  }

  // ===========================================================================
  // GL canvas + tooltip DOM management
  // ===========================================================================
  function _ensureGLCanvas(plotArea, mainC) {
    let glC = document.getElementById('sphere-gl-canvas');
    if (glC) return glC;
    glC = document.createElement('canvas');
    glC.id = 'sphere-gl-canvas';
    // z-index 2 puts the sphere ON TOP of the main #canvas (z-index 1).
    // pointer-events:auto means GL canvas captures all pointer input while
    // shown. The Schwarz GL canvas (z-index 0) and the 2D overlay canvas
    // (z-index 1) are unreachable as long as this canvas is displayed. On
    // deactivate the canvas is hidden via display:none so events fall
    // through to the lower layers again.
    glC.style.cssText =
      'position:absolute; left:0; top:0; width:100%; height:100%; '
      + 'pointer-events:auto; z-index:2; cursor:grab; display:none;';
    mainC.style.position   = 'relative';
    mainC.style.zIndex     = '1';
    mainC.style.background = 'transparent';
    plotArea.insertBefore(glC, mainC.nextSibling);   // sibling after main
    return glC;
  }
  function _showGLLayer(glC, show) {
    if (glC) glC.style.display = show ? '' : 'none';
  }
  function _ensureHoverTooltip(plotArea) {
    if (document.getElementById('sphere-hover-tooltip')) return;
    const tt = document.createElement('div');
    tt.id = 'sphere-hover-tooltip';
    tt.style.display = 'none';
    if (plotArea) plotArea.appendChild(tt);
  }
  function _showHoverTooltip(clientX, clientY, text) {
    const tt = document.getElementById('sphere-hover-tooltip');
    if (!tt) return;
    const plotArea = document.getElementById('plot-area');
    if (!plotArea) return;
    const rect = plotArea.getBoundingClientRect();
    let lx = clientX - rect.left + 14;
    let ly = clientY - rect.top  + 14;
    tt.textContent = text;
    tt.style.display = 'block';
    const ttW = tt.offsetWidth, ttH = tt.offsetHeight;
    if (lx + ttW > rect.width  - 4) lx = clientX - rect.left - ttW - 12;
    if (ly + ttH > rect.height - 4) ly = clientY - rect.top  - ttH - 12;
    if (lx < 4) lx = 4;
    if (ly < 4) ly = 4;
    tt.style.left = lx + 'px';
    tt.style.top  = ly + 'px';
  }
  function _hideHoverTooltip() {
    const tt = document.getElementById('sphere-hover-tooltip');
    if (tt) tt.style.display = 'none';
  }

  // ===========================================================================
  // Sidebar cards (display + camera)
  // ===========================================================================
  function _makeDisplayCard(state) {
    const card = document.createElement('section');
    card.className = 'card view-sphere-only';
    card.innerHTML = `
      <h2>Sphere display</h2>
      <div class="row">
        <label><input type="checkbox" id="sp-show-boundary" checked> Show boundary ∂K</label>
      </div>
      <div class="row">
        <label><input type="checkbox" id="sp-show-poles" checked> Show poles &amp; ∞ marker</label>
      </div>
      <div class="row" style="margin-top:6px;">
        <label>Rim shading:
          <input id="sp-rim" type="range" min="0" max="0.5" step="0.01" value="0.3"
                 style="width:100px; margin-left:6px;">
          <span id="sp-rim-val">0.30</span>
        </label>
      </div>
      <div class="row" style="margin-top:6px;">
        <label>Texture resolution:
          <select id="sp-texsize">
            <option value="512">512²</option>
            <option value="1024" selected>1024²</option>
            <option value="2048">2048²</option>
          </select>
        </label>
      </div>
    `;
    setTimeout(() => {
      document.getElementById('sp-show-boundary').addEventListener('change', function () {
        state.display.showBoundary = this.checked;
        if (state.renderer) { state.renderer.setDisplayParams(state.display); _requestRender(state); }
      });
      document.getElementById('sp-show-poles').addEventListener('change', function () {
        state.display.showPoles = state.display.showNorthPole = this.checked;
        if (state.renderer) { state.renderer.setDisplayParams(state.display); _requestRender(state); }
      });
      const rimSlider = document.getElementById('sp-rim');
      const rimVal    = document.getElementById('sp-rim-val');
      rimSlider.addEventListener('input', function () {
        state.display.rimDarken = parseFloat(this.value);
        rimVal.textContent = state.display.rimDarken.toFixed(2);
        if (state.renderer) { state.renderer.setDisplayParams(state.display); _requestRender(state); }
      });
      document.getElementById('sp-texsize').addEventListener('change', function () {
        state.params.texSize = parseInt(this.value) || 1024;
        if (state.renderer) { state.renderer.setRenderParams(state.params); _requestRender(state); }
      });
    }, 0);
    return card;
  }
  function _makeCameraCard(state) {
    const card = document.createElement('section');
    card.className = 'card view-sphere-only';
    card.innerHTML = `
      <h2>Camera &amp; readout</h2>
      <div class="hint">Left-drag to orbit · Scroll to zoom · Double-click to reset</div>
      <div class="row" style="margin-top:6px;">
        <button id="sp-reset-cam" class="small">Reset view</button>
        <button id="sp-recompute-sphere" class="small" style="margin-left:6px;">Recompute</button>
      </div>
      <div id="sp-hover" class="hint" style="margin-top:8px; font-family: ui-monospace, Consolas, monospace; font-size:11px; min-height:2em;">
        Hover over the sphere for coordinates.
      </div>
    `;
    setTimeout(() => {
      document.getElementById('sp-reset-cam').addEventListener('click', () => _resetCamera(state));
      document.getElementById('sp-recompute-sphere').addEventListener('click', () => {
        if (state.renderer && state.renderer.markFractalDirty) {
          state.renderer.markFractalDirty();
          _requestRender(state);
        }
      });
    }, 0);
    return card;
  }

  // ===========================================================================
  // Render loop
  // ===========================================================================
  function _requestRender(state) {
    if (state.frameRequested) return;
    state.frameRequested = true;
    requestAnimationFrame(() => {
      state.frameRequested = false;
      if (!state.renderer || !state.glCanvas) return;
      const dpr = window.devicePixelRatio || 1;
      const W   = Math.max(1, Math.floor(state.glCanvas.clientWidth  * dpr));
      const H   = Math.max(1, Math.floor(state.glCanvas.clientHeight * dpr));
      state.renderer.render(state.camera, { W, H });
    });
  }

  // ===========================================================================
  // Camera + event handlers
  // ===========================================================================
  function _resetCamera(state) {
    state.camera.azimuth   = DEF_CAMERA.azimuth;
    state.camera.elevation = DEF_CAMERA.elevation;
    state.camera.distance  = DEF_CAMERA.distance;
    _requestRender(state);
  }

  function _attachGLCanvasEvents(glC, state, isActive) {
    glC.addEventListener('mousedown', e => {
      if (!isActive()) return;
      if (e.button !== 0) return;
      state.drag.active = true;
      state.drag.lastX  = e.clientX;
      state.drag.lastY  = e.clientY;
      glC.style.cursor = 'grabbing';
      e.preventDefault();
    });
    glC.addEventListener('mousemove', e => _onHover(e, state, isActive));
    glC.addEventListener('mouseleave', _hideHoverTooltip);
    glC.addEventListener('wheel', e => {
      if (!isActive()) return;
      e.preventDefault();   // {passive:false} below; stops page scroll
      state.camera.distance *= Math.pow(1.1, e.deltaY * 0.03);
      state.camera.distance  = Math.max(1.1, Math.min(10, state.camera.distance));
      _requestRender(state);
    }, { passive: false });
    glC.addEventListener('dblclick', () => _resetCamera(state));
    glC.addEventListener('contextmenu', e => e.preventDefault());
  }
  function _attachWindowHandlers(state, isActive) {
    window.addEventListener('mousemove', e => {
      if (!state.drag.active || !isActive()) return;
      const dx = e.clientX - state.drag.lastX;
      const dy = e.clientY - state.drag.lastY;
      state.drag.lastX = e.clientX;
      state.drag.lastY = e.clientY;
      state.camera.azimuth   -= dx * 0.008;
      state.camera.elevation  = Math.max(
        -Math.PI / 2 + 0.04,
        Math.min(Math.PI / 2 - 0.04, state.camera.elevation + dy * 0.008));
      _requestRender(state);
    });
    window.addEventListener('mouseup', () => {
      if (!state.drag.active) return;
      state.drag.active = false;
      if (state.glCanvas) state.glCanvas.style.cursor = 'grab';
    });
    document.addEventListener('keydown', e => {
      if ((e.key === 'r' || e.key === 'R') && isActive()) _resetCamera(state);
    });
  }
  function _ensureResizeObserver(glC, state, isActive) {
    if (typeof ResizeObserver === 'undefined') {
      window.addEventListener('resize', () => {
        if (isActive() && state.renderer) _requestRender(state);
      });
      return;
    }
    new ResizeObserver(() => {
      if (isActive() && state.renderer) _requestRender(state);
    }).observe(glC);
  }

  // ===========================================================================
  // Hover readout — raycasts against unit sphere
  // ===========================================================================
  function _onHover(e, state, isActive) {
    if (!isActive()) return;
    const readoutEl = document.getElementById('sp-hover');
    const glC = state.glCanvas;
    if (!glC) return;
    const rect = glC.getBoundingClientRect();
    const mx   = e.clientX - rect.left;
    const my   = e.clientY - rect.top;
    const W    = rect.width, H = rect.height;
    if (W < 1 || H < 1) return;

    const cam  = state.camera;
    const cosEl = Math.cos(cam.elevation), sinEl = Math.sin(cam.elevation);
    const cosAz = Math.cos(cam.azimuth),   sinAz = Math.sin(cam.azimuth);
    const d    = cam.distance;
    const eye  = [d*cosEl*cosAz, d*cosEl*sinAz, d*sinEl];

    const fwd = _normalize(_neg(eye));
    const rgt = _normalize(_cross(fwd, [0,0,1]));
    const up  = _cross(rgt, fwd);

    const fovY  = Math.PI / 3;
    const tanH  = Math.tan(fovY * 0.5);
    const ndcX  =  2 * mx / W - 1;
    const ndcY  = -(2 * my / H - 1);
    const aspect = W / H;

    const dir = _normalize([
      fwd[0] + rgt[0] * ndcX * tanH * aspect + up[0] * ndcY * tanH,
      fwd[1] + rgt[1] * ndcX * tanH * aspect + up[1] * ndcY * tanH,
      fwd[2] + rgt[2] * ndcX * tanH * aspect + up[2] * ndcY * tanH,
    ]);

    const b = 2 * _dot(eye, dir);
    const c = _dot(eye, eye) - 1;
    const D = b * b - 4 * c;
    if (D < 0) { if (readoutEl) readoutEl.textContent = '—'; _hideHoverTooltip(); return; }
    const t = (-b - Math.sqrt(D)) / 2;
    if (t < 0) { if (readoutEl) readoutEl.textContent = '—'; _hideHoverTooltip(); return; }

    const hit = { x: eye[0]+t*dir[0], y: eye[1]+t*dir[1], z: eye[2]+t*dir[2] };

    let wStr;
    const SC = typeof SphereCommon !== 'undefined' ? SphereCommon : null;
    if (SC) {
      const w = SC.unprojectFromSphere(hit, 1e-3);
      wStr = w ? ('w = ' + _fmtC(w.re, w.im, 3)) : 'w = ∞';
    } else { wStr = ''; }

    const txt = `(x,y,z) = (${hit.x.toFixed(3)}, ${hit.y.toFixed(3)}, ${hit.z.toFixed(3)})\n${wStr}`;
    if (readoutEl) readoutEl.textContent = txt;
    _showHoverTooltip(e.clientX, e.clientY, txt);
  }

  // ===========================================================================
  // Helpers
  // ===========================================================================
  function _extractPoles(hData) {
    const out = [];
    if (!hData || !hData.poles) return out;
    for (const p of hData.poles) {
      if (p.a && typeof p.a.re === 'number') {
        out.push({ re: p.a.re, im: p.a.im });
      }
    }
    return out;
  }

  function _fmtC(re, im, d) {
    const r = re.toFixed(d);
    if (Math.abs(im) < 5e-4) return r;
    const iStr = (im < 0 ? ' − ' : ' + ') + Math.abs(im).toFixed(d) + 'i';
    return r + iStr;
  }
  function _normalize(a) {
    const len = Math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    if (len < 1e-300) return [0,0,1];
    return [a[0]/len, a[1]/len, a[2]/len];
  }
  function _neg(a) { return [-a[0], -a[1], -a[2]]; }
  function _dot(a, b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
  function _cross(a, b) {
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
  }

  // ---------------------------------------------------------------------------
  // Public export
  // ---------------------------------------------------------------------------
  QD.SphereView = { mount };

})();
