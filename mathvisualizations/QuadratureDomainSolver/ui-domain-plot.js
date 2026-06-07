// =============================================================================
// ui-domain-plot.js -- DomainPlot canvas class (Inverse-tab plot)
//
// Extracted from ui.js by the P0 split (§1.1). Pulled out as a factory so it
// can receive its ui.js closures (state, modeDescriptor, formatTick, sub) via
// dependency injection rather than relying on shared script-scope variables
// (which don't cross <script> tags).
//
// ui.js installs the class via:
//   const DomainPlot = window.QD_UI.installDomainPlot({ state, modeDescriptor,
//                                                       formatTick, sub });
//
// Inside the class, the four injected names are simple closure variables —
// the class body itself is identical to what previously lived in ui.js.
// =============================================================================

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  global.QD_UI.installDomainPlot = function installDomainPlot(deps) {
    const state          = deps.state;
    const modeDescriptor = deps.modeDescriptor;
    const formatTick     = deps.formatTick;
    const sub            = deps.sub;

    // Hover hit-test radius for the pole-proximity annotation in the readout
    // (HANDOFF #33 / #35). Larger than the click hit-radius (9 px) so the
    // cursor doesn't need to be pixel-perfect over the pole dot.
    const POLE_HOVER_HIT_RADIUS_PX = 12;

class DomainPlot {
  constructor(canvas, readout) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.readout = readout;
    this.dpr = window.devicePixelRatio || 1;

    this.view = { cx: 0, cy: 0, scale: 100 };  // pixels per unit
    this.data = null;

    // Vector-field perf: the field re-samples h(w) over the whole visible grid
    // on every repaint (no cache), which is the dominant cost while dragging /
    // panning / zooming. We defer it during an active gesture and recompute once
    // when the gesture settles (~150 ms idle, or immediately on mouseup).
    this._vfInteracting = false;
    this._vfSettleTimer = null;

    // Callbacks for click-drag on quadrature-node dots. Set by ui.js.
    this.onPoleDrag    = null;   // (idx, worldPoint) -> void
    this.onPoleDragEnd = null;   // (idx)            -> void
    // Double-click on empty plot space -> add a new simple pole there. Set by
    // ui.js (which owns state.poles + the solve pipeline). (worldPoint) -> void
    this.onAddPole     = null;

    this.attachEvents();
    this.resize();
  }

  // Returns the index of the pole dot under (x, y) in CSS pixels, or -1 if
  // none is within the hit-test radius.
  _hitTestPole(x, y, radius = 9) {
    if (!this.data || !this.data.poles) return -1;
    let bestI = -1, bestD2 = radius * radius;
    for (let i = 0; i < this.data.poles.length; i++) {
      const sp = this.toScreen(this.data.poles[i].re, this.data.poles[i].im);
      const dx = sp.x - x, dy = sp.y - y;
      const d2 = dx * dx + dy * dy;
      if (d2 <= bestD2) { bestD2 = d2; bestI = i; }
    }
    return bestI;
  }

  resize() {
    const rect = this.canvas.getBoundingClientRect();
    this.cssW = rect.width;
    this.cssH = rect.height;
    this.canvas.width  = Math.round(this.cssW * this.dpr);
    this.canvas.height = Math.round(this.cssH * this.dpr);
    this.render();
  }

  toScreen(re, im) {
    return {
      x: this.cssW / 2 + (re - this.view.cx) * this.view.scale,
      y: this.cssH / 2 - (im - this.view.cy) * this.view.scale,
    };
  }
  toWorld(x, y) {
    return {
      re: this.view.cx + (x - this.cssW / 2) / this.view.scale,
      im: this.view.cy - (y - this.cssH / 2) / this.view.scale,
    };
  }

  attachEvents() {
    let panning = false, lastX = 0, lastY = 0;
    let draggingPole = -1;          // index, or -1

    // Belt-and-suspenders gate: the QD/LQD inverse tab is the only one that
    // owns the main #canvas as a 2D drawing surface; the Schwarz and Riemann-
    // sphere tabs overlay their own GL canvases.  These handlers must early-
    // return when those other tabs are active, otherwise a drag/wheel in
    // those tabs would trigger this.render() and repaint the QD plot over
    // the other tab's content.  (The sphere tab additionally puts its GL
    // canvas on top with pointer-events:auto so most events never reach
    // this listener — this check is the second line of defense, and the
    // only line of defense for the Schwarz tab whose GL canvas sits below.)
    //
    // Exception: an in-progress pan or pole-drag started on the QD tab is
    // allowed to complete even if the user mid-drag switched tabs (avoids
    // a stuck "panning=true" state).
    function _qdTabActive() {
      const btn = document.querySelector('.tab-btn.active');
      return btn && btn.dataset.tab === 'qd';
    }

    this.canvas.addEventListener('mousedown', e => {
      if (!_qdTabActive()) return;
      const rect = this.canvas.getBoundingClientRect();
      const x = e.clientX - rect.left, y = e.clientY - rect.top;
      // Hit-test a pole first; if hit, take ownership of this drag for the
      // pole rather than starting a pan.
      const hit = this._hitTestPole(x, y);
      if (hit >= 0) {
        draggingPole = hit;
        this.canvas.style.cursor = 'grabbing';
        e.preventDefault();
        return;
      }
      panning = true; lastX = e.clientX; lastY = e.clientY;
      this.canvas.style.cursor = 'grabbing';
    });

    window.addEventListener('mousemove', e => {
      // Allow an already-in-progress drag/pole-drag to finish even after a
      // tab switch.  Brand-new hover/readout updates require the QD tab.
      const inProgressDrag = panning || draggingPole >= 0;
      if (!inProgressDrag && !_qdTabActive()) return;

      const rect = this.canvas.getBoundingClientRect();
      const x = e.clientX - rect.left, y = e.clientY - rect.top;

      if (draggingPole >= 0) {
        this._markVfInteracting();
        const w = this.toWorld(x, y);
        if (this.onPoleDrag) this.onPoleDrag(draggingPole, w);
        return;
      }
      if (panning) {
        this._markVfInteracting();
        const dx = e.clientX - lastX, dy = e.clientY - lastY;
        this.view.cx -= dx / this.view.scale;
        this.view.cy += dy / this.view.scale;
        lastX = e.clientX; lastY = e.clientY;
        this.render();
      }

      // Mouse-coordinate readout + hover cursor (pointer over a pole, grab
      // elsewhere) — only when the cursor is over the canvas.
      if (x >= 0 && x <= this.cssW && y >= 0 && y <= this.cssH) {
        const w = this.toWorld(x, y);
        let text = `w = ${w.re.toFixed(4)} ${w.im >= 0 ? '+' : '-'} ${Math.abs(w.im).toFixed(4)}i`;
        // Append nearby-pole annotation when the cursor is within the
        // hit-test radius (HANDOFF #33). Useful for quickly identifying
        // which pole a residue belongs to as the user scans.
        const hitIdx = this._hitTestPole(x, y, POLE_HOVER_HIT_RADIUS_PX);
        if (hitIdx >= 0) {
          const a = this.data.poles[hitIdx];
          const SUBS = ['₀','₁','₂','₃','₄','₅','₆','₇','₈','₉'];
          const sub = String(hitIdx + 1).split('').map(d => SUBS[+d] || d).join('');
          text += `  ·  near pole a${sub} = ${a.re.toFixed(3)}${a.im >= 0 ? '+' : '-'}${Math.abs(a.im).toFixed(3)}i`;
        }
        this.readout.textContent = text;
        if (!panning && draggingPole < 0) {
          this.canvas.style.cursor = hitIdx >= 0 ? 'pointer' : 'grab';
        }
      }
    });

    window.addEventListener('mouseup', () => {
      // Gesture over → redraw the deferred vector field now (no-op if no gesture
      // was active, so plain clicks don't force an extra repaint).
      this._settleVectorField();
      if (draggingPole >= 0) {
        const idx = draggingPole;
        draggingPole = -1;
        this.canvas.style.cursor = 'grab';
        if (this.onPoleDragEnd) this.onPoleDragEnd(idx);
        return;
      }
      panning = false;
      this.canvas.style.cursor = 'grab';
    });

    this.canvas.addEventListener('wheel', e => {
      if (!_qdTabActive()) return;
      e.preventDefault();
      const rect = this.canvas.getBoundingClientRect();
      const x = e.clientX - rect.left, y = e.clientY - rect.top;
      const wBefore = this.toWorld(x, y);
      const factor = Math.exp(-e.deltaY * 0.001);
      this.view.scale = Math.max(1e-3, Math.min(1e7, this.view.scale * factor));
      const wAfter = this.toWorld(x, y);
      this.view.cx += wBefore.re - wAfter.re;
      this.view.cy += wBefore.im - wAfter.im;
      this._markVfInteracting();   // wheel has no mouseup; the settle timer redraws
      this.render();
    }, { passive: false });

    // Double-click on empty plot space drops a new simple pole at that w
    // (coefficient 1) — the inverse of the click-drag-to-move gesture. A
    // double-click that lands on an existing pole dot is ignored so we never
    // stack a duplicate node on one the user was aiming to grab. ui.js owns
    // state.poles + the solve, so we just hand it the world coordinate.
    this.canvas.addEventListener('dblclick', e => {
      if (!_qdTabActive()) return;
      const rect = this.canvas.getBoundingClientRect();
      const x = e.clientX - rect.left, y = e.clientY - rect.top;
      if (this._hitTestPole(x, y) >= 0) return;   // on an existing pole — ignore
      e.preventDefault();
      if (this.onAddPole) this.onAddPole(this.toWorld(x, y));
    });
  }

  setData(d) {
    this.data = d;
    this.render();
  }
  clear() {
    this.data = null;
    this.render();
  }

  reset() {
    this.view = { cx: 0, cy: 0, scale: 100 };
    this.render();
  }

  fit() {
    if (!this.data || !this.data.boundaryPts || this.data.boundaryPts.length === 0) return;
    let minRe = Infinity, maxRe = -Infinity, minIm = Infinity, maxIm = -Infinity;
    for (const p of this.data.boundaryPts) {
      if (p.re < minRe) minRe = p.re; if (p.re > maxRe) maxRe = p.re;
      if (p.im < minIm) minIm = p.im; if (p.im > maxIm) maxIm = p.im;
    }
    for (const p of this.data.poles) {
      if (p.re < minRe) minRe = p.re; if (p.re > maxRe) maxRe = p.re;
      if (p.im < minIm) minIm = p.im; if (p.im > maxIm) maxIm = p.im;
    }
    const dx = Math.max(1e-6, maxRe - minRe);
    const dy = Math.max(1e-6, maxIm - minIm);
    const pad = 0.15;
    const sx = this.cssW / (dx * (1 + 2*pad));
    const sy = this.cssH / (dy * (1 + 2*pad));
    this.view.scale = Math.min(sx, sy);
    this.view.cx = (minRe + maxRe) / 2;
    this.view.cy = (minIm + maxIm) / 2;
    this.render();
  }

  // Live-update one pole marker's position WITHOUT waiting for a re-solve.
  // Used during a pole drag so the dot tracks the cursor 1:1 even when the
  // solver (which redraws the boundary + every marker via setData/showSolution)
  // lags a frame or two behind. The next solve's showSolution() rewrites
  // this.data.poles from the canonical hData, so the transient write here is
  // harmless. render() is rAF-coalesced, so calling this every mousemove is
  // cheap. No-op if there's no data yet or the index is out of range.
  setLivePole(idx, w) {
    if (!this.data || !this.data.poles) return;
    if (idx < 0 || idx >= this.data.poles.length) return;
    this.data.poles[idx] = { re: w.re, im: w.im };
    this.render();
  }

  // Mark the start/continuation of an interactive gesture (pan / pole-drag /
  // wheel-zoom): suppress the expensive vector field for the duration and arm a
  // settle timer that redraws it once the gesture goes idle (covers wheel-zoom,
  // which has no mouseup). render() itself is unchanged — the boundary, markers
  // and grid still repaint every frame; only drawVectorField() is gated.
  _markVfInteracting() {
    this._vfInteracting = true;
    if (this._vfSettleTimer) clearTimeout(this._vfSettleTimer);
    this._vfSettleTimer = setTimeout(() => {
      this._vfSettleTimer = null;
      this._vfInteracting = false;
      this.render();
    }, 150);
  }

  // Settle immediately (e.g. on mouseup). No-op when no gesture was in progress
  // (a plain click) so we don't force an extra repaint.
  _settleVectorField() {
    if (this._vfSettleTimer) { clearTimeout(this._vfSettleTimer); this._vfSettleTimer = null; }
    if (this._vfInteracting) { this._vfInteracting = false; this.render(); }
  }

  // Public repaint entry point. Coalesces bursts of render() calls (pan,
  // wheel-zoom, resize, setData) into a single paint per animation frame so a
  // fast drag can't queue dozens of full redraws in one frame (A8). The actual
  // drawing is _renderNow(); it always reads the latest this.data / this.view,
  // so coalescing never paints stale state.
  render() {
    if (this._renderScheduled) return;
    this._renderScheduled = true;
    const raf = (typeof requestAnimationFrame === 'function')
      ? requestAnimationFrame
      : (fn) => setTimeout(() => fn(), 16);
    raf(() => { this._renderScheduled = false; this._renderNow(); });
  }

  _renderNow() {
    const c = this.ctx;
    c.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    c.clearRect(0, 0, this.cssW, this.cssH);
    c.fillStyle = '#fafafa';
    c.fillRect(0, 0, this.cssW, this.cssH);

    this.drawGrid();
    this.drawAxes();

    // Vector field underlays the domain so the boundary remains crisp. Skipped
    // mid-gesture (drag / pan / zoom) — it re-samples h(w) over the whole grid
    // with no cache, the dominant cost while interacting; _settleVectorField /
    // the settle timer redraw it once the gesture ends.
    if (state.vectorFieldMode !== 'off' && this.data && this.data.hData
        && !this._vfInteracting) {
      this.drawVectorField();
    }

    if (this.data && this.data.boundaryPts && this.data.boundaryPts.length > 0) {
      this.drawBoundary();
    }
    // Optional dashed overlay (e.g. for Direct-tab round-trip diagnostics).
    if (this.data && this.data.overlayBoundary && this.data.overlayBoundary.length > 0) {
      this.drawOverlayBoundary();
    }
    if (this.data && this.data.poles)  this.drawPoles();
    if (this.data && this.data.w0)     this.drawW0();
    // Critical-set image overlay (zeros of φ', mapped to w-plane).  Drawn
    // last so the markers sit on top of everything; lazy-computed and
    // cached on state.current.criticalSet to avoid recomputing per render.
    // Curvature heat-strip overlay: color ∂Ω by |κ| (brightest at a cusp).
    // Lazy-computed + cached on state.current.observables (shared with the
    // "Geometry & accuracy" card). Drawn over the boundary but under the
    // critical-set / cusp markers.
    if (state.showCurvature && this.data && this.data.phi) {
      this.drawCurvature();
    }
    if (state.showCriticalSet && this.data && this.data.phi) {
      this.drawCriticalSet();
    }
    // Boundary cusp markers (zeros of φ′ on ∂𝔻), from the async classifyCusps
    // result stashed on state.current.cuspProps. Drawn last, on top.
    if (this.data && this.data.phi) this.drawCusps();
    // Annotated-phenomena overlay (#9): harmonic-measure / curvature peaks +
    // symmetry axes. Opt-in; reuses the cached observables + symmetry results.
    if (state.showPhenomena && this.data && this.data.phi) this.drawPhenomenaAnnotations();
    // Faber-polynomial roots overlay (UQD only): roots of the Faber polynomials
    // F_n of the bounded complement K, plotted in the ζ = image plane (where ∂Ω
    // lives). They cluster inside K, the "hole" of the unbounded domain. The
    // payload is pushed by ui-faber.js onto state.faberRoots. Drawn last (top).
    if (state.showFaberRoots && state.faberRoots && this.data && this.data.unbounded) {
      this.drawFaberRoots();
    }

    // Empty-state hint (item 3/4): when there's no solved boundary and no poles
    // to drag, the canvas would be an unexplained blank — tell the user how to
    // begin instead.
    const noBoundary = !(this.data && this.data.boundaryPts && this.data.boundaryPts.length > 0);
    const noPoles = !(this.data && this.data.poles && this.data.poles.length > 0);
    if (noBoundary && noPoles) this.drawEmptyState();
  }

  drawEmptyState() {
    const c = this.ctx;
    c.save();
    c.textAlign = 'center';
    c.fillStyle = '#9aa3b2';
    c.font = '600 15px system-ui, -apple-system, sans-serif';
    c.fillText('No domain yet', this.cssW / 2, this.cssH / 2 - 10);
    c.font = '13px system-ui, -apple-system, sans-serif';
    c.fillText('Pick a Preset or Thesis example, or double-click on the plot to add a pole.',
               this.cssW / 2, this.cssH / 2 + 14);
    c.restore();
  }

  drawOverlayBoundary() {
    const c = this.ctx;
    const pts = this.data.overlayBoundary;
    c.save();
    c.beginPath();
    const p0 = this.toScreen(pts[0].re, pts[0].im);
    c.moveTo(p0.x, p0.y);
    for (let i = 1; i < pts.length; i++) {
      const p = this.toScreen(pts[i].re, pts[i].im);
      c.lineTo(p.x, p.y);
    }
    c.closePath();
    c.strokeStyle = '#d4a017';                            // gold
    c.lineWidth = 1.8;
    if (c.setLineDash) c.setLineDash([6, 4]);
    c.stroke();
    c.restore();
  }

  // ----- Vector-field overlay: V(w) = conj(h(w)) ---------------------------
  // Samples h on a grid that is anchored to WORLD coordinates (multiples of
  // a "nice" step in the w-plane), so that as the user pans/zooms each arrow
  // stays glued to a specific point in the domain. The step is chosen so
  // that adjacent grid points are roughly 26 px apart on screen at the
  // current zoom, snapped to a 1 / 2 / 5 / 10 sequence times a power of 10.
  //
  // Each arrow's length and opacity scale by tanh(log10(1 + |h|)) so a few
  // large values near poles don't dominate the visual.
  //
  // Screen-direction note: conj(h) = Re(h) - i Im(h) is the vector
  // (Re h, -Im h) in the math plane. Screen y is flipped relative to the
  // imaginary axis, so the screen-space arrow direction is just (Re h, Im h).
  drawVectorField() {
    const c = this.ctx;
    const hData = this.data.hData;
    if (!hData) return;
    // h(w) = poly part + sum of principal parts. Bail only if BOTH are empty;
    // a poly-only unbounded h still has a vector field worth drawing.
    const polyLen = (hData.polyPart && hData.polyPart.length) || 0;
    if (hData.poles.length === 0 && polyLen === 0) return;

    // Pick a nice world-coordinate step targeting ~26 px on screen.
    const targetPx = 26;
    const targetWorld = targetPx / this.view.scale;
    const exp10 = Math.floor(Math.log10(targetWorld));
    const frac = targetWorld / Math.pow(10, exp10);
    let mult;
    if      (frac < 1.5) mult = 1;
    else if (frac < 3)   mult = 2;
    else if (frac < 7)   mult = 5;
    else                 mult = 10;
    const worldStep = mult * Math.pow(10, exp10);
    const stepPx = worldStep * this.view.scale;
    const arrowMaxLen = stepPx * 0.72;

    // Visible world bounds (screen y is flipped vs imaginary axis).
    const wTL = this.toWorld(0, 0);
    const wBR = this.toWorld(this.cssW, this.cssH);
    const minRe = Math.min(wTL.re, wBR.re);
    const maxRe = Math.max(wTL.re, wBR.re);
    const minIm = Math.min(wTL.im, wBR.im);
    const maxIm = Math.max(wTL.im, wBR.im);
    const iMin = Math.floor(minRe / worldStep);
    const iMax = Math.ceil (maxRe / worldStep);
    const jMin = Math.floor(minIm / worldStep);
    const jMax = Math.ceil (maxIm / worldStep);

    c.lineWidth = 1;
    c.lineCap = 'round';

    // Pole screen positions are constant across grid samples — compute once so
    // the per-sample proximity test is a cheap squared-distance check rather
    // than O(poles) toScreen() + hypot() calls at every grid point.
    const poleScreens = hData.poles.map(p => this.toScreen(p.a.re, p.a.im));

    for (let i = iMin; i <= iMax; i++) {
      for (let j = jMin; j <= jMax; j++) {
        const wRe = i * worldStep;
        const wIm = j * worldStep;
        const screen = this.toScreen(wRe, wIm);
        const sx = screen.x;
        const sy = screen.y;
        // Skip if the grid point is well off-canvas.
        if (sx < -arrowMaxLen || sx > this.cssW + arrowMaxLen ||
            sy < -arrowMaxLen || sy > this.cssH + arrowMaxLen) continue;
        const w = { re: wRe, im: wIm };

        // Skip arrows that fall too close to any pole in screen-space; the
        // field diverges and the visual is meaningless there.
        let skip = false;
        for (let k = 0; k < poleScreens.length; k++) {
          const ddx = poleScreens[k].x - sx, ddy = poleScreens[k].y - sy;
          if (ddx * ddx + ddy * ddy < 81) { skip = true; break; }   // <9px, squared
        }
        if (skip) continue;

        // Evaluate h(w) = Σ_{l=0..m_∞} C_{∞,l} w^l + Σ_j Σ_s C_{j,s}/(w − a_j)^s
        let hr = 0, hi = 0;
        let bad = false;
        // Polynomial part of h (unbounded only; otherwise polyPart is empty).
        const polyPart = hData.polyPart || [];
        if (polyPart.length > 0) {
          let wPowRe = 1, wPowIm = 0;                       // w^0
          for (let l = 0; l < polyPart.length; l++) {
            const cR = polyPart[l].re, cI = polyPart[l].im;
            hr += cR * wPowRe - cI * wPowIm;
            hi += cR * wPowIm + cI * wPowRe;
            if (l + 1 < polyPart.length) {
              const nR = wPowRe * w.re - wPowIm * w.im;
              const nI = wPowRe * w.im + wPowIm * w.re;
              wPowRe = nR; wPowIm = nI;
            }
          }
        }
        // Finite-pole part.
        for (const pole of hData.poles) {
          const dx = w.re - pole.a.re, dy = w.im - pole.a.im;
          const d2 = dx*dx + dy*dy;
          if (d2 < 1e-14) { bad = true; break; }
          // dPowRe + dPowIm·i  =  (w - a)^{s+1}  for s = 0, 1, ...
          let dPowRe = dx, dPowIm = dy;
          for (let s = 0; s < pole.principal.length; s++) {
            const cR = pole.principal[s].re, cI = pole.principal[s].im;
            const den = dPowRe*dPowRe + dPowIm*dPowIm;
            if (den < 1e-30) { bad = true; break; }
            // C / (dPowRe + dPowIm·i)  =  C · (dPowRe − dPowIm·i) / den
            hr += (cR * dPowRe + cI * dPowIm) / den;
            hi += (cI * dPowRe - cR * dPowIm) / den;
            // Advance (w-a)^{s+1} -> (w-a)^{s+2}
            const nRe = dPowRe * dx - dPowIm * dy;
            const nIm = dPowRe * dy + dPowIm * dx;
            dPowRe = nRe; dPowIm = nIm;
          }
          if (bad) break;
        }
        if (bad) continue;

        // Per-family extra contributions to h (e.g. q/w pole at origin for
        // singular LQDs). Default is null = no extra terms.
        const desc = modeDescriptor();
        if (desc.extraHContrib) {
          const extra = desc.extraHContrib(w, hData, this.data && this.data.phi, state);
          hr += extra.re;
          hi += extra.im;
        }

        // Pólya field is V = conj(h). External-potential field depends on
        // the QD family:
        //   classical:       V = w − conj(h)        (∇ of |w|² − 2 Re H)
        //   PQD (weight α):  V = (1/α)·w·|w|^{2(α−1)} − conj(h)
        //                                            (∇ of (1/α)|w|^{2α}/2 − 2 Re H)
        //   LQD:             V = ln|w|²/conj(w) − conj(h)
        //                                            (∇ of (1/2)ln²|w|² − 2 Re H)
        // Math-plane vectors (before screen y-flip):
        //   conj(h)      = (Re h, −Im h)
        //   w − conj(h)  = (Re w − Re h, Im w + Im h)
        //   ln|w|²/conj(w) = ln|w|² · (Re w, Im w) / |w|²
        //                  = (Re w · ln|w|²/|w|², Im w · ln|w|²/|w|²)
        // Screen y is flipped vs the imaginary axis, so the screen-direction
        // negates the math y-component.
        let fieldX, fieldY;
        if (state.vectorFieldMode === 'external') {
          if (desc.externalFieldKind === 'lqd') {
            // V = ln|w|² / conj(w) − conj(h). Clip near origin per descriptor
            // (singular LQDs need a larger floor since 0 ∈ Ω).
            const absW2 = w.re * w.re + w.im * w.im;
            if (absW2 < desc.vectorFieldOriginAbs2Floor) continue;
            const logScale = Math.log(absW2) / absW2;
            fieldX = w.re * logScale - hr;
            fieldY = -(w.im * logScale) - hi;
            // ln|w|²/conj(w) = ln|w|² · w/|w|², so its math-Im = (Im w) · ln|w|²/|w|²;
            // screen-y = −(math-Im); and conj(h)'s math-Im = −Im h, screen-y = +Im h.
            // V = (math)  (Re w · L − Re h, Im w · L + Im h)  with L = ln|w|²/|w|²
            // screen      (Re w · L − Re h, −(Im w · L + Im h))
            //           = (Re w · L − Re h, −Im w · L − Im h)
            // (Re-deriving here matches the assignment above.)
          } else if (desc.externalFieldKind === 'pqd') {
            // V = (1/α)·w·|w|^{2(α−1)} − conj(h). The |w|^{2(α−1)} weight makes
            // the external field family-specific (analogous to the LQD ln|w|²
            // factor). Reduces to the classical w − conj(h) at α = 1.
            const absW2 = w.re * w.re + w.im * w.im;
            if (absW2 < desc.vectorFieldOriginAbs2Floor) continue;
            const alpha = (this.data && this.data.phi && this.data.phi.alpha) || +state.alpha || 1;
            const s = Math.pow(absW2, alpha - 1) / alpha;   // (1/α)·|w|^{2(α−1)}
            fieldX =  w.re * s - hr;
            fieldY = -(w.im * s) - hi;
          } else {
            fieldX =  w.re - hr;
            fieldY = -w.im - hi;
          }
        } else {
          // Pólya field V = conj(h): screen vector is (Re h, Im h).
          fieldX = hr;
          fieldY = hi;
        }
        const mag = Math.hypot(fieldX, fieldY);
        if (!isFinite(mag) || mag === 0) continue;

        // Length + opacity from tanh(log10(1 + |V|)) — short/dim for small
        // values, saturating for very large ones.
        const sat = Math.tanh(Math.log10(1 + mag));
        const len = arrowMaxLen * (0.25 + 0.75 * sat);
        const alpha = 0.18 + 0.55 * sat;

        const dirX = fieldX / mag;
        const dirY = fieldY / mag;

        // Center the arrow on the grid point.
        const baseX = sx - 0.35 * len * dirX;
        const baseY = sy - 0.35 * len * dirY;
        const tipX  = sx + 0.65 * len * dirX;
        const tipY  = sy + 0.65 * len * dirY;

        c.strokeStyle = `rgba(58, 84, 124, ${alpha.toFixed(3)})`;
        c.beginPath();
        c.moveTo(baseX, baseY);
        c.lineTo(tipX, tipY);
        c.stroke();

        // Arrowhead
        const ahLen = Math.max(3, len * 0.32);
        const ahAng = 0.45;
        const ang = Math.atan2(dirY, dirX);
        c.beginPath();
        c.moveTo(tipX, tipY);
        c.lineTo(tipX - ahLen * Math.cos(ang - ahAng), tipY - ahLen * Math.sin(ang - ahAng));
        c.moveTo(tipX, tipY);
        c.lineTo(tipX - ahLen * Math.cos(ang + ahAng), tipY - ahLen * Math.sin(ang + ahAng));
        c.stroke();
      }
    }
  }

  // Pick a "nice" tick spacing for the current scale
  niceStep() {
    const target = 80 / this.view.scale;            // ~80 px per major grid line
    const exp = Math.floor(Math.log10(target));
    const frac = target / Math.pow(10, exp);
    let step;
    if      (frac < 1.5) step = 1;
    else if (frac < 3)   step = 2;
    else if (frac < 7)   step = 5;
    else                 step = 10;
    return step * Math.pow(10, exp);
  }

  drawGrid() {
    const c = this.ctx;
    const step = this.niceStep();
    const tl = this.toWorld(0, 0);
    const br = this.toWorld(this.cssW, this.cssH);
    const minRe = Math.floor(tl.re / step) * step;
    const maxRe = Math.ceil(br.re / step) * step;
    const minIm = Math.floor(br.im / step) * step;
    const maxIm = Math.ceil(tl.im / step) * step;

    c.strokeStyle = '#e8eaef';
    c.lineWidth = 1;
    c.beginPath();
    for (let r = minRe; r <= maxRe + 1e-9; r += step) {
      const x = this.toScreen(r, 0).x;
      c.moveTo(x, 0); c.lineTo(x, this.cssH);
    }
    for (let i = minIm; i <= maxIm + 1e-9; i += step) {
      const y = this.toScreen(0, i).y;
      c.moveTo(0, y); c.lineTo(this.cssW, y);
    }
    c.stroke();

    // tick labels
    c.fillStyle = '#777';
    c.font = '10px ui-monospace, "SF Mono", Consolas, monospace';
    c.textAlign = 'left';
    c.textBaseline = 'top';
    const y0 = Math.max(2, Math.min(this.cssH - 12, this.toScreen(0, 0).y + 2));
    for (let r = minRe; r <= maxRe + 1e-9; r += step) {
      if (Math.abs(r) < step * 1e-6) continue;
      c.fillText(formatTick(r, step), this.toScreen(r, 0).x + 2, y0);
    }
    c.textAlign = 'left';
    const x0 = Math.max(2, Math.min(this.cssW - 30, this.toScreen(0, 0).x + 2));
    for (let i = minIm; i <= maxIm + 1e-9; i += step) {
      if (Math.abs(i) < step * 1e-6) continue;
      c.fillText(formatTick(i, step) + 'i', x0, this.toScreen(0, i).y + 2);
    }
  }

  drawAxes() {
    const c = this.ctx;
    c.strokeStyle = '#bbb';
    c.lineWidth = 1;
    c.beginPath();
    const yAxisX = this.toScreen(0, 0).x;
    const xAxisY = this.toScreen(0, 0).y;
    c.moveTo(0, xAxisY); c.lineTo(this.cssW, xAxisY);
    c.moveTo(yAxisX, 0); c.lineTo(yAxisX, this.cssH);
    c.stroke();
  }

  drawBoundary() {
    const c = this.ctx;
    const pts = this.data.boundaryPts;
    c.beginPath();
    const p0 = this.toScreen(pts[0].re, pts[0].im);
    c.moveTo(p0.x, p0.y);
    for (let i = 1; i < pts.length; i++) {
      const p = this.toScreen(pts[i].re, pts[i].im);
      c.lineTo(p.x, p.y);
    }
    c.closePath();

    const ok = this.data.univalent;
    if (this.data.unbounded) {
      // Unbounded: shade the bounded complement K (= inside of the boundary
      // curve) in a contrasting muted color and outline ∂Ω.
      c.fillStyle = ok ? 'rgba(180, 195, 220, 0.45)' : 'rgba(220, 180, 180, 0.45)';
      c.fill('evenodd');
      c.strokeStyle = ok ? '#1a3e7a' : '#b53030';
      c.lineWidth = 1.8;
      c.stroke();
    } else {
      // Bounded: shade Ω (= inside of the curve) in the standard tint.
      c.fillStyle   = ok ? 'rgba(86, 119, 168, 0.16)' : 'rgba(181, 48, 48, 0.14)';
      c.fill('evenodd');
      c.strokeStyle = ok ? '#1a3e7a' : '#b53030';
      c.lineWidth = 1.6;
      c.stroke();
    }
  }

  drawPoles() {
    const c = this.ctx;
    c.font = '11px system-ui, sans-serif';
    c.textBaseline = 'middle';
    for (let i = 0; i < this.data.poles.length; i++) {
      const p = this.data.poles[i];
      const s = this.toScreen(p.re, p.im);
      c.beginPath();
      c.arc(s.x, s.y, 5.5, 0, 2*Math.PI);
      c.fillStyle = '#b53030';
      c.fill();
      c.strokeStyle = '#fff'; c.lineWidth = 1.5;
      c.stroke();
      c.fillStyle = '#b53030';
      c.textAlign = 'left';
      c.fillText('a' + sub(i+1), s.x + 7, s.y);
    }
  }

  drawW0() {
    const c = this.ctx;
    const s = this.toScreen(this.data.w0.re, this.data.w0.im);
    c.strokeStyle = '#1a3e7a';
    c.lineWidth = 1.6;
    c.beginPath();
    c.moveTo(s.x - 5, s.y); c.lineTo(s.x + 5, s.y);
    c.moveTo(s.x, s.y - 5); c.lineTo(s.x, s.y + 5);
    c.stroke();
    c.fillStyle = '#1a3e7a';
    c.font = '11px system-ui, sans-serif';
    c.textAlign = 'left';
    c.textBaseline = 'top';
    c.fillText('φ(0)', s.x + 6, s.y + 4);
  }

  // -------------------------------------------------------------------------
  // Curvature heat-strip: color ∂Ω by signed curvature |κ| from QD.boundary-
  // Observables. The per-node curvature[] is aligned with the boundary points
  // w[], so each segment is stroked with a color ramped by |κ|/maxκ (cool→hot);
  // a cusp (κ → ∞) shows as the hottest stretch. Result is cached on
  // state.current.observables (keyed by phi identity), shared with the card.
  // -------------------------------------------------------------------------
  drawCurvature() {
    const phi = this.data.phi;
    if (!phi || typeof QD === 'undefined' || !QD.boundaryObservables) return;
    if (!state.current) return;

    const cur = state.current.observables;
    let obs = cur && cur.obs;
    if (!(obs && obs._phiRef === phi)) {
      try { obs = QD.boundaryObservables(phi, { samples: 720 }); } catch (e) { return; }
      obs._phiRef = phi;
      state.current.observables = Object.assign({}, cur, { obs });
    }
    if (!obs.w || obs.w.length < 2 || !(obs.maxCurvature > 0)) return;

    const c = this.ctx;
    const w = obs.w, kappa = obs.curvature, N = w.length;
    // Robust color reference: at a cusp |κ| → ∞, so normalizing by the TRUE max
    // crushes every other segment to ≈0 and the whole boundary reads uniform
    // blue. Use a high percentile of |κ| as the reference instead — the cusp
    // region then saturates to red while the rest of ∂Ω gets a real gradient.
    const absK = [];
    for (let i = 0; i < N; i++) { const v = Math.abs(kappa[i]); if (isFinite(v)) absK.push(v); }
    absK.sort((p, q) => p - q);
    let denom = absK.length
      ? absK[Math.min(absK.length - 1, Math.floor(absK.length * 0.90))]
      : obs.maxCurvature;
    if (!(denom > 0)) denom = (obs.maxCurvature > 0 ? obs.maxCurvature : 1);
    c.save();
    c.lineWidth = 3;
    c.lineCap = 'round';
    for (let i = 0; i < N; i++) {
      const a = this.toScreen(w[i].re, w[i].im);
      const b = this.toScreen(w[(i + 1) % N].re, w[(i + 1) % N].im);
      // |κ| at this node, normalized to the robust reference (cusp clamps to 1);
      // sqrt eases the dynamic range so the ramp toward the cusp stays visible.
      let t = Math.abs(kappa[i]) / denom;
      if (!isFinite(t)) t = 1;
      t = Math.max(0, Math.min(1, Math.sqrt(t)));
      const hue = 210 - 210 * t;            // 210° (cool blue) → 0° (hot red)
      c.strokeStyle = `hsl(${hue.toFixed(0)}, 85%, 45%)`;
      c.beginPath();
      c.moveTo(a.x, a.y);
      c.lineTo(b.x, b.y);
      c.stroke();
    }
    c.restore();
  }

  // -------------------------------------------------------------------------
  // Critical-set overlay: w-plane images of {z : φ'(z) = 0}.
  //
  // Computed lazily on first request (and on toggle-on after a fresh solve),
  // cached on state.current.criticalSet so subsequent pans/zooms don't pay
  // the Newton cost.  The cache is keyed by reference identity of
  // this.data.phi — if a new solve produces a new phi object, the cache
  // is recomputed automatically.
  //
  // Visual encoding:
  //   severity 'critical' (zero of φ' strictly inside the relevant disk):
  //     red filled disk, 6 px radius.  This is the bad case — φ is
  //     non-univalent.
  //   severity 'near'     (|z| within 0.05 of the unit circle):
  //     orange filled disk, 5 px radius.  Imminent-degeneracy warning.
  //   severity 'safe'     (zero of φ' outside the relevant disk):
  //     small gray hollow circle, 3.5 px radius.  Background info.
  //
  // Each marker carries a 1-letter tag showing severity.
  // -------------------------------------------------------------------------
  drawCriticalSet() {
    const phi = this.data.phi;
    if (!phi || typeof QD === 'undefined' || !QD.findCriticalPoints) return;

    // Look up / refresh cache.
    if (!state.current) return;
    const cached = state.current.criticalSet;
    let cs;
    if (cached && cached._phiRef === phi) {
      cs = cached;
    } else {
      try {
        cs = QD.findCriticalPoints(phi);
      } catch (e) {
        return;       // silent on solver error — overlay is purely diagnostic
      }
      cs._phiRef = phi;
      state.current.criticalSet = cs;
    }
    if (!cs.points || cs.points.length === 0) return;

    const c = this.ctx;
    c.save();
    c.font = '10px ui-monospace, Consolas, monospace';
    c.textBaseline = 'middle';
    c.textAlign    = 'left';

    for (const p of cs.points) {
      const s = this.toScreen(p.w.re, p.w.im);
      // Skip if off-screen by a large margin (saves draw work; markers very
      // far outside the visible region clutter the corner-clip).
      if (s.x < -40 || s.x > this.cssW + 40) continue;
      if (s.y < -40 || s.y > this.cssH + 40) continue;

      let fill, stroke, r, tag;
      switch (p.severity) {
        case 'critical':
          fill   = '#d12d2d';
          stroke = '#ffffff';
          r      = 6.0;
          tag    = '!';
          break;
        case 'near':
          fill   = '#d97706';   // orange
          stroke = '#ffffff';
          r      = 5.0;
          tag    = '~';
          break;
        case 'safe':
        default:
          fill   = null;        // hollow
          stroke = '#888888';
          r      = 3.5;
          tag    = '';
          break;
      }
      c.beginPath();
      c.arc(s.x, s.y, r, 0, 2 * Math.PI);
      if (fill) {
        c.fillStyle = fill;
        c.fill();
      }
      c.strokeStyle = stroke;
      c.lineWidth   = (p.severity === 'safe') ? 1.2 : 1.6;
      c.stroke();

      if (tag) {
        c.fillStyle = (p.severity === 'critical') ? '#d12d2d' : '#a85706';
        c.fillText(tag + ' |z|=' + p.absZ.toFixed(3), s.x + r + 3, s.y);
      }
    }
    c.restore();
  }

  // -------------------------------------------------------------------------
  // Faber-polynomial roots overlay (ζ = image plane). state.faberRoots:
  //   { mode:'all'|'single', N, n, sets:[{ n, roots:[Complex], converged }] }
  // 'all'    → union of roots of F_1..F_N as teal hollow circles, alpha fading
  //            with order so the low orders read clearly (no labels — too many).
  // 'single' → roots of one F_n as violet filled diamonds (distinct from the
  //            critical-set disks / cusp triangles / phenomena diamonds), labelled.
  // -------------------------------------------------------------------------
  drawFaberRoots() {
    const fr = state.faberRoots;
    if (!fr || !fr.sets || !fr.sets.length) return;
    const c = this.ctx;
    c.save();
    c.font = '10px ui-monospace, Consolas, monospace';
    c.textBaseline = 'middle';
    c.textAlign    = 'left';

    const onScreen = (s) => s.x >= -40 && s.x <= this.cssW + 40 && s.y >= -40 && s.y <= this.cssH + 40;

    if (fr.mode === 'single') {
      const set = fr.sets[0];
      if (!set || !set.roots) { c.restore(); return; }
      const r = 4.5;
      let labelled = false;
      for (const root of set.roots) {
        const s = this.toScreen(root.re, root.im);
        if (!onScreen(s)) continue;
        c.beginPath();
        c.moveTo(s.x, s.y - r); c.lineTo(s.x + r, s.y);
        c.lineTo(s.x, s.y + r); c.lineTo(s.x - r, s.y);
        c.closePath();
        c.fillStyle   = '#7c3aed';        // violet
        c.fill();
        c.strokeStyle = '#ffffff';
        c.lineWidth   = 1.2;
        c.stroke();
        if (!labelled) {
          c.fillStyle = '#5b21b6';
          c.fillText('F' + this._sub(set.n) + ' roots', s.x + r + 3, s.y);
          labelled = true;
        }
      }
    } else {
      // 'all' — union over orders, teal hollow circles fading with order.
      const N = fr.sets.length;
      for (let i = 0; i < N; i++) {
        const set = fr.sets[i];
        if (!set || !set.roots) continue;
        const alpha = 0.35 + 0.55 * (set.n / Math.max(1, fr.N));   // higher order → brighter
        c.strokeStyle = 'rgba(13,148,136,' + alpha.toFixed(3) + ')';   // teal
        c.lineWidth   = 1.1;
        for (const root of set.roots) {
          const s = this.toScreen(root.re, root.im);
          if (!onScreen(s)) continue;
          c.beginPath();
          c.arc(s.x, s.y, 3, 0, 2 * Math.PI);
          c.stroke();
        }
      }
    }
    c.restore();
  }

  // Unicode subscript for small non-negative integers (for "Fₙ" labels).
  _sub(n) {
    const map = '₀₁₂₃₄₅₆₇₈₉';
    return String(n).split('').map(d => map[+d] || d).join('');
  }

  // -------------------------------------------------------------------------
  // Boundary cusp markers — the w-plane locations of (incipient) cusps from
  // QD.classifyCusps, stashed async on state.current.cuspProps. A filled
  // magenta triangle marks an ACTUAL boundary cusp (φ′-zero on ∂𝔻); a hollow
  // gray triangle marks an INCIPIENT one (φ′-zero near but not on ∂𝔻). Each is
  // labelled with its (p,q) type. Only drawn for the primary solution (the
  // cuspProps are computed for it); alternates suppress the markers to avoid
  // mismatched overlays.
  // -------------------------------------------------------------------------
  drawCusps() {
    if (!state.current || !state.current.cuspProps) return;
    if ((state.selectedSolutionIdx || 0) !== 0) return;     // primary only
    const cusps = state.current.cuspProps.cusps || [];
    if (!cusps.length) return;

    const c = this.ctx;
    c.save();
    c.font = '10px ui-monospace, Consolas, monospace';
    c.textBaseline = 'middle';
    c.textAlign    = 'left';

    for (const cu of cusps) {
      const s = this.toScreen(cu.w.re, cu.w.im);
      if (s.x < -40 || s.x > this.cssW + 40) continue;
      if (s.y < -40 || s.y > this.cssH + 40) continue;
      const r = 6;
      // Upward triangle centered on the cusp tip.
      c.beginPath();
      c.moveTo(s.x, s.y - r);
      c.lineTo(s.x - r * 0.9, s.y + r * 0.7);
      c.lineTo(s.x + r * 0.9, s.y + r * 0.7);
      c.closePath();
      if (cu.isCusp) {
        c.fillStyle = '#b5179e';   // magenta — an actual cusp
        c.fill();
        c.strokeStyle = '#ffffff';
        c.lineWidth = 1.4;
        c.stroke();
      } else {
        c.strokeStyle = '#888888'; // hollow — incipient
        c.lineWidth = 1.2;
        c.stroke();
      }
      c.fillStyle = cu.isCusp ? '#b5179e' : '#888888';
      c.fillText('(' + cu.type[0] + ',' + cu.type[1] + ')', s.x + r + 3, s.y);
    }
    c.restore();
  }

  // -------------------------------------------------------------------------
  // Annotated-phenomena overlay (#9) — labels the phenomena the cusp/critical
  // overlays don't: the harmonic-measure hot spot (the tip, where ρ = 1/(2π|φ′|)
  // peaks ⇔ |φ′| is smallest), the max-curvature point on ∂Ω, and the domain's
  // symmetry axes. Reads the already-cached observables sweep + the
  // detectSymmetry result (state.current.symmetry); adds no solve cost. Primary
  // solution only (the caches are computed for it).
  // -------------------------------------------------------------------------
  drawPhenomenaAnnotations() {
    if (!state.current || !this.data || !this.data.phi) return;
    if ((state.selectedSolutionIdx || 0) !== 0) return;     // primary only
    const phi = this.data.phi;
    const c = this.ctx;
    const sub = (n) => String(n).replace(/[0-9]/g, d => '₀₁₂₃₄₅₆₇₈₉'[+d]);

    // ---- symmetry axes (dashed) + a D_n / Z_n / circle badge ----------------
    const sym = state.current.symmetry;
    if (sym && (sym.reflectionAxes.length || sym.rotationalOrder > 1 || sym.continuous)) {
      const ctr = sym.center || { re: 0, im: 0 };
      const sCtr = this.toScreen(ctr.re, ctr.im);
      const L = (this.cssW + this.cssH) / this.view.scale;   // world length covering the canvas
      c.save();
      c.setLineDash([6, 5]);
      c.strokeStyle = 'rgba(124, 58, 237, 0.55)';            // faint violet
      c.lineWidth = 1.2;
      for (const a of sym.reflectionAxes) {
        const dx = Math.cos(a) * L, dy = Math.sin(a) * L;
        const p1 = this.toScreen(ctr.re - dx, ctr.im - dy);
        const p2 = this.toScreen(ctr.re + dx, ctr.im + dy);
        c.beginPath(); c.moveTo(p1.x, p1.y); c.lineTo(p2.x, p2.y); c.stroke();
      }
      c.setLineDash([]);
      const nAx = sym.reflectionAxes.length;
      let badge = '';
      if (sym.continuous) badge = '○ rotational';
      else if (sym.rotationalOrder > 1) badge = (nAx ? 'D' : 'Z') + sub(sym.rotationalOrder);
      else if (nAx) badge = 'mirror';
      if (badge) {
        c.font = '11px ui-monospace, Consolas, monospace';
        c.fillStyle = 'rgba(124, 58, 237, 0.95)';
        c.textBaseline = 'middle'; c.textAlign = 'left';
        c.fillText(badge, sCtr.x + 6, sCtr.y - 7);
      }
      c.restore();
    }

    // ---- harmonic-measure peak + max-curvature point ------------------------
    const obs = state.current.observables && state.current.observables.obs;
    if (obs && obs._phiRef === phi && typeof QD !== 'undefined' && QD.evalPhi) {
      c.save();
      c.font = '10px ui-monospace, Consolas, monospace';
      c.textBaseline = 'middle'; c.textAlign = 'left';
      const mark = (theta, color, label) => {
        let wpt; try { wpt = QD.evalPhi({ re: Math.cos(theta), im: Math.sin(theta) }, phi); }
        catch (e) { return; }
        if (!wpt || !isFinite(wpt.re) || !isFinite(wpt.im)) return;
        const s = this.toScreen(wpt.re, wpt.im);
        if (s.x < -40 || s.x > this.cssW + 40 || s.y < -40 || s.y > this.cssH + 40) return;
        c.beginPath();                                       // diamond marker
        c.moveTo(s.x, s.y - 5); c.lineTo(s.x + 5, s.y);
        c.lineTo(s.x, s.y + 5); c.lineTo(s.x - 5, s.y); c.closePath();
        c.fillStyle = color; c.fill();
        c.strokeStyle = '#ffffff'; c.lineWidth = 1.2; c.stroke();
        // Label sits BELOW the marker (with a short leader). A boundary cusp at
        // this same point already labels its (p,q) type to the RIGHT at marker
        // height, so dropping the phenomena label avoids overlapping it.
        const lx = s.x + 8, ly = s.y + 15;
        c.strokeStyle = color; c.lineWidth = 1;
        c.beginPath(); c.moveTo(s.x + 3, s.y + 4); c.lineTo(lx - 2, ly - 3); c.stroke();
        c.fillStyle = color; c.fillText(label, lx, ly);
      };
      const thHM = obs.minAbsPhiPrimeTheta, thK = obs.argMaxCurvatureTheta;
      const maxK = obs.maxCurvature;
      const angClose = (a, b) => {
        let d = Math.abs(((a - b) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI));
        return Math.min(d, 2 * Math.PI - d) < 0.09;          // ~5°
      };
      if (thHM != null && thK != null && maxK > 0 && angClose(thHM, thK)) {
        mark(thHM, '#7c3aed', 'tip: ρ-peak · max|κ|=' + maxK.toFixed(2));
      } else {
        if (thHM != null) mark(thHM, '#0d9488', 'harmonic-measure peak');
        if (thK != null && maxK > 0) mark(thK, '#dc2626', 'max |κ|=' + maxK.toFixed(2));
      }
      c.restore();
    }
  }
}

    return DomainPlot;
  };

})(typeof window !== 'undefined' ? window
                                 : (typeof globalThis !== 'undefined' ? globalThis : this));
