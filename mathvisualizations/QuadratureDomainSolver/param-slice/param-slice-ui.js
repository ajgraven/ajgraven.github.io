// =============================================================================
// param-slice-ui.js -- Sidebar + canvas overlay for the "Parameter slice" tab.
//
// On first tab activation, mounts a sidebar with:
//   • Source-of-scenario card     (mirrors the Inverse tab's current state)
//   • Axis pickers (X, optional Y)
//   • Range/resolution controls for each active axis
//   • Render + cancel buttons
//   • Legend + progress
//
// Owns the shared #canvas while active. Draws pixel data via putImageData
// (categorical classification color, see param-slice-common.js).
//
// Click a pixel → re-solve at that parameter value, push φ into the QD tab
// via window.QD_UI.loadPhiIntoQdTab, switch to QD tab.
// =============================================================================

(function (global) {
  'use strict';
  if (typeof QD === 'undefined') return;
  const PS = global.ParamSlice;
  if (!PS) return;

  const sliceState = {
    mounted: false,
    canvas: null,
    pool: null,                // Worker pool (created on first render)
    poolPromise: null,
    activeJob: null,           // { cancel, done } from runSweep
    lastImageData: null,
    lastTiles: null,           // [{ row, results }, ...]  for click-load
    lastAxes: null,            // [{ ref, min, max, n }, optional 2nd]
    lastScenario: null,        // { hData, norm, opts, mode }
    progress: { done: 0, total: 1, t0: 0 },
    // Adaptive-mesh φ cache + lookup (hoisted from runAdaptive2D in
    // HANDOFF #33 so the hover handler can read it after the render
    // completes). nearestPhi is bound to the most recent runAdaptive2D's
    // grid dimensions; classGrid / iterGrid let tooltip and hover-card
    // read the per-pixel result without re-running the solver.
    classGrid: null,           // Uint8Array(n0*n1) of CLASS_TO_IDX values
    iterGrid:  null,           // Uint8Array(n0*n1) of iteration counts
    nearestPhi: function () { return null; },  // (col, row) → {phi, iterCount} | null
    gridDims: { n0: 0, n1: 0 },
    // Mini-canvas memoisation (HANDOFF #35): identity-keyed cache of the
    // last hovered φ's boundary samples, so scrubbing across many pixels
    // pointing at the same nearest φ doesn't re-sample 128 points per frame.
    miniCache: { phi: null, pts: null },
  };

  // ---------------------------------------------------------------------------
  // Lazy mount
  // ---------------------------------------------------------------------------
  document.addEventListener('tab-changed', function (e) {
    if (!e.detail || e.detail.tab !== 'param-slice') {
      detachCanvasHandler();
      return;
    }
    if (!sliceState.mounted) { mountSidebar(); sliceState.mounted = true; }
    attachCanvasHandler();
    refreshAxisOptions();
    refreshScenarioStatus();
    refreshQuadratureDataCard();
    // If we have a cached image, re-draw it (canvas might have been
    // overwritten by the QD tab in the meantime).
    if (sliceState.lastImageData) repaint();
  });

  // ---------------------------------------------------------------------------
  // Sidebar construction
  // ---------------------------------------------------------------------------
  function mountSidebar() {
    const root = document.getElementById('controls-param-slice');
    if (!root) return;
    root.innerHTML = '';
    root.appendChild(makeIntroCard());
    root.appendChild(makeScenarioCard());
    root.appendChild(makeQuadratureDataCard());   // HANDOFF #33
    root.appendChild(makeAxisCard());
    root.appendChild(makeRunCard());
    root.appendChild(makeHoveredQDCard());        // HANDOFF #33
    root.appendChild(makeLegendCard());
    attachHelpButtons();                          // HANDOFF #33
  }

  function makeIntroCard() {
    const card = document.createElement('section');
    card.className = 'card';
    card.id = 'ps-intro-card';
    card.innerHTML = `
      <h2>Parameter-slice cartography</h2>
      <div class="hint">
        Sweep one or two parameters of <code>h(w)</code>; at each sample,
        attempt the inverse-solve; color each pixel by whether the resulting
        Ω is a valid QD/LQD.  Click any green pixel to load that φ into the
        Inverse tab.
        <br><br>
        The base scenario is whatever you have currently set in the
        <b>QD / LQD</b> tab (poles, residues, c, q, w₀, mode).  Re-open this
        tab after editing to refresh the axis options.
      </div>
    `;
    return card;
  }

  function makeScenarioCard() {
    const card = document.createElement('section');
    card.className = 'card';
    card.id = 'ps-scenario-card';
    card.innerHTML = `
      <h2>Base scenario</h2>
      <div id="ps-scenario-status" class="hint" style="color:#333;">(no scenario)</div>
    `;
    return card;
  }

  function makeAxisCard() {
    const card = document.createElement('section');
    card.className = 'card';
    card.id = 'ps-axis-card';
    card.innerHTML = `
      <h2>Axes</h2>
      <div class="row" style="margin-bottom:6px;">
        <label>X axis:
          <select id="ps-axis-x" style="min-width:160px;"></select>
        </label>
      </div>
      <div class="row" id="ps-x-range-row" style="margin-bottom:8px;">
        <label>min <input id="ps-x-min" type="number" step="any" value="-1" style="width:64px;"></label>
        <label style="margin-left:6px;">max <input id="ps-x-max" type="number" step="any" value="1" style="width:64px;"></label>
        <label style="margin-left:6px;">n <input id="ps-x-n" type="number" min="2" max="512" value="128" style="width:60px;"></label>
      </div>
      <div class="row" style="margin-bottom:6px;">
        <label>Y axis:
          <select id="ps-axis-y" style="min-width:160px;">
            <option value="">— (1-D sweep)</option>
          </select>
        </label>
      </div>
      <div class="row" id="ps-y-range-row" style="margin-bottom:8px; display:none;">
        <label>min <input id="ps-y-min" type="number" step="any" value="-1" style="width:64px;"></label>
        <label style="margin-left:6px;">max <input id="ps-y-max" type="number" step="any" value="1" style="width:64px;"></label>
        <label style="margin-left:6px;">n <input id="ps-y-n" type="number" min="2" max="512" value="96" style="width:60px;"></label>
      </div>
    `;
    setTimeout(() => {
      document.getElementById('ps-axis-x').addEventListener('change', onAxisChange);
      document.getElementById('ps-axis-y').addEventListener('change', onAxisChange);
    }, 0);
    return card;
  }

  // Quality presets — per-pixel identity-check rigor. Higher samples /
  // tighter tolerance = fewer false-positive `identity-fail` (yellow)
  // pixels at the cost of per-pixel solve time. See HANDOFF #32.
  const QUALITY_PRESETS = {
    fast:     { univalenceSamples: 32,  identityTol: 1e-5 },
    standard: { univalenceSamples: 128, identityTol: 1e-6 },
    rigorous: { univalenceSamples: 512, identityTol: 1e-7 },
  };
  function readQualityPreset() {
    const el = document.getElementById('ps-quality');
    const key = (el && el.value) || 'standard';
    return QUALITY_PRESETS[key] || QUALITY_PRESETS.standard;
  }

  function makeRunCard() {
    const card = document.createElement('section');
    card.className = 'card';
    card.id = 'ps-run-card';
    card.innerHTML = `
      <h2>Run</h2>
      <div class="row" style="margin-bottom:8px;">
        <label>Quality:
          <select id="ps-quality" style="margin-left:4px;">
            <option value="fast">Fast (N=32, tol=1e-5)</option>
            <option value="standard" selected>Standard (N=128, tol=1e-6)</option>
            <option value="rigorous">Rigorous (N=512, tol=1e-7)</option>
          </select>
        </label>
      </div>
      <div class="row">
        <button id="ps-run" class="small">Render slice</button>
        <button id="ps-cancel" class="small" style="margin-left:6px;" disabled>Cancel</button>
      </div>
      <div id="ps-progress" class="hint" style="margin-top:8px; min-height:1.2em;"></div>
    `;
    setTimeout(() => {
      document.getElementById('ps-run').addEventListener('click', startRun);
      document.getElementById('ps-cancel').addEventListener('click', cancelRun);
    }, 0);
    return card;
  }

  function makeLegendCard() {
    const card = document.createElement('section');
    card.className = 'card';
    card.id = 'ps-legend-card';
    const rows = PS.CLASS_ORDER.map(cls => {
      const [r,g,b] = PS.CLASS_COLORS[cls];
      const label = PS.CLASS_LABELS[cls];
      return `<div style="display:flex; align-items:center; gap:6px; margin:3px 0;">
        <span style="display:inline-block; width:14px; height:14px;
                     background:rgb(${r},${g},${b}); border:1px solid #888;"></span>
        <span style="font-size:12px;">${label}</span>
      </div>`;
    }).join('');
    card.innerHTML = `<h2>Legend</h2>${rows}
      <div class="hint" style="margin-top:6px;">
        Valid-QD pixels are dimmed with iteration count (bright = fast convergence).
        Click any pixel to load that φ into the Inverse tab.
      </div>`;
    return card;
  }

  // -------------------------------------------------------------------------
  // Quadrature-data card — shows the h(w) form of the scenario being plotted
  // plus the axis assignments. HANDOFF #33.
  // -------------------------------------------------------------------------
  function makeQuadratureDataCard() {
    const card = document.createElement('section');
    card.className = 'card';
    card.id = 'ps-quad-card';
    card.innerHTML = `
      <h2>Quadrature data</h2>
      <div class="hint">The h(w) being sliced, and which parameter each axis varies.</div>
      <div id="ps-quad-h-row" style="display:flex; align-items:flex-start; gap:6px; margin-top:6px;">
        <span style="font-size:12px; color:#5677a8; font-family:ui-monospace,monospace;">h(w) =</span>
        <code id="ps-quad-h" style="font-size:11px; flex:1; white-space:pre-wrap; word-break:break-all;">(open the QD tab once)</code>
      </div>
      <div id="ps-quad-axes" style="margin-top:8px; font-size:12px;"></div>
    `;
    setTimeout(() => {
      // Add a copy button next to the h(w) text.
      const row = document.getElementById('ps-quad-h-row');
      if (row && global.QD && global.QD.QoL && global.QD.QoL.copyButton) {
        row.appendChild(global.QD.QoL.copyButton(() => {
          const el = document.getElementById('ps-quad-h');
          return el ? el.textContent : '';
        }, { title: 'Copy h(w) expression' }));
      }
    }, 0);
    return card;
  }

  function refreshQuadratureDataCard() {
    const hEl = document.getElementById('ps-quad-h');
    if (!hEl) return;
    const snap = snapshotScenario();
    if (!snap) { hEl.textContent = '(open the QD tab once)'; return; }
    // QD.formatH expects the parse-h.js residue-shape; convert from the
    // solver's principal-shape used by hData.
    let text = '0';
    try {
      const polesForFormat = (snap.hData.poles || []).map(p => ({
        a: p.a, residues: p.principal || [],
      }));
      text = QD.formatH({ poles: polesForFormat, polyCoeffs: snap.hData.polyPart || [] });
    } catch (e) { text = '(formatting failed)'; }
    hEl.textContent = text;
    // Axis assignments — show whatever's currently selected.
    const axesEl = document.getElementById('ps-quad-axes');
    if (!axesEl) return;
    const xSel = document.getElementById('ps-axis-x');
    const ySel = document.getElementById('ps-axis-y');
    if (!xSel || !ySel) { axesEl.innerHTML = ''; return; }
    const params = PS.listAvailableParams({ hData: snap.hData, norm: snap.norm }, snap.mode);
    const xLabel = (xSel.value !== '' && params[+xSel.value]) ? params[+xSel.value].label : '—';
    const yLabel = (ySel.value !== '' && params[+ySel.value]) ? params[+ySel.value].label : null;
    const rows = [`<div><b>X axis:</b> <code style="font-size:11px;">${escapeHTML(xLabel)}</code></div>`];
    if (yLabel) rows.push(`<div><b>Y axis:</b> <code style="font-size:11px;">${escapeHTML(yLabel)}</code></div>`);
    else        rows.push(`<div style="color:#888;"><i>1-D sweep (no Y axis)</i></div>`);
    axesEl.innerHTML = rows.join('');
  }

  // -------------------------------------------------------------------------
  // Hovered-QD card — mini-canvas preview + per-pixel readout. Updates on
  // mousemove via refreshHoveredQDCard. HANDOFF #33.
  // -------------------------------------------------------------------------
  function makeHoveredQDCard() {
    const card = document.createElement('section');
    card.className = 'card';
    card.id = 'ps-hover-card';
    card.innerHTML = `
      <h2>Hovered QD</h2>
      <div class="hint">Hover any pixel above to preview the QD at those parameters.</div>
      <div style="display:flex; gap:10px; align-items:flex-start; margin-top:6px;">
        <canvas id="ps-hover-canvas" width="160" height="160"
                style="width:160px; height:160px; background:#0d111a;
                       border:1px solid #c5c9d2; border-radius:4px;"></canvas>
        <div style="flex:1; min-width:0;">
          <div id="ps-hover-text" style="font-size:11px; line-height:1.5;
               color:#444; font-family:ui-monospace,monospace; word-break:break-all;">
            (hover the slice)
          </div>
          <button id="ps-hover-send" class="small" style="margin-top:8px;" disabled>Send to inverse</button>
        </div>
      </div>
    `;
    setTimeout(() => {
      const btn = document.getElementById('ps-hover-send');
      if (btn) btn.addEventListener('click', () => {
        if (_lastHovered) loadPixel(_lastHovered.col, _lastHovered.row,
                                    _lastHovered.xVal, _lastHovered.yVal);
      });
    }, 0);
    return card;
  }

  // Tracks the last cell shown in the Hovered-QD card so the "Send" button
  // knows what to load.
  let _lastHovered = null;

  function refreshHoveredQDCard(cellInfo) {
    const txtEl = document.getElementById('ps-hover-text');
    const btn   = document.getElementById('ps-hover-send');
    const canvas = document.getElementById('ps-hover-canvas');
    if (!txtEl || !canvas) return;
    if (!cellInfo) {
      _lastHovered = null;
      txtEl.innerHTML = '<i style="color:#888;">(hover the slice)</i>';
      if (btn) btn.disabled = true;
      _clearMiniCanvas(canvas);
      return;
    }
    _lastHovered = cellInfo;
    if (btn) btn.disabled = false;
    const lbl = PS.CLASS_LABELS[cellInfo.cls] || cellInfo.cls;
    const lines = [];
    lines.push(`<b>${escapeHTML(cellInfo.xAxis.label || 'x')}</b> = ${fmtNum(cellInfo.xVal)}`);
    if (cellInfo.yAxis) {
      lines.push(`<b>${escapeHTML(cellInfo.yAxis.label || 'y')}</b> = ${fmtNum(cellInfo.yVal)}`);
    }
    lines.push(`<span style="color:#5677a8;">${escapeHTML(lbl)}</span>` +
               (cellInfo.cls === PS.CLASS_VALID ? ` · ${cellInfo.iter} iters` : ''));
    txtEl.innerHTML = lines.map(l => `<div>${l}</div>`).join('');
    // Render the boundary if we have a φ for (or near) this pixel.
    const hit = sliceState.nearestPhi(cellInfo.col, cellInfo.row);
    _drawMiniBoundary(canvas, hit && hit.phi);
  }

  function _clearMiniCanvas(canvas) {
    const ctx = canvas.getContext('2d');
    ctx.fillStyle = '#0d111a';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = '#666';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText('—', canvas.width/2, canvas.height/2 + 4);
  }
  function _drawMiniBoundary(canvas, phi) {
    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height;
    ctx.fillStyle = '#0d111a';
    ctx.fillRect(0, 0, W, H);
    if (!phi) {
      ctx.fillStyle = '#666';
      ctx.font = '10px sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText('no cached φ', W/2, H/2 + 4);
      return;
    }
    // Identity-keyed memoisation (HANDOFF #35): scrubbing the cursor across
    // adjacent pixels often hits the SAME nearest-cached φ for many frames.
    // Without this, we re-sample 128 boundary points per hover frame. The
    // cache holds one slot; nearestPhi returns the same φ object until the
    // next render replaces the bucket, so === identity is the right key.
    let pts = null;
    if (sliceState.miniCache && sliceState.miniCache.phi === phi) {
      pts = sliceState.miniCache.pts;
    } else {
      try {
        pts = (QD.sampleBoundary || function () { return null; })(phi, MINI_BOUNDARY_SAMPLES);
      } catch (e) { pts = null; }
      if (pts && sliceState.miniCache) {
        sliceState.miniCache.phi = phi;
        sliceState.miniCache.pts = pts;
      }
    }
    if (!pts || pts.length === 0) {
      ctx.fillStyle = '#a44';
      ctx.font = '10px sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText('render failed', W/2, H/2 + 4);
      return;
    }
    // Compute bounding box.
    let xmin = +Infinity, xmax = -Infinity, ymin = +Infinity, ymax = -Infinity;
    for (const p of pts) {
      if (p.re < xmin) xmin = p.re; if (p.re > xmax) xmax = p.re;
      if (p.im < ymin) ymin = p.im; if (p.im > ymax) ymax = p.im;
    }
    const pad = 0.1 * Math.max(xmax - xmin, ymax - ymin, 0.01);
    xmin -= pad; xmax += pad; ymin -= pad; ymax += pad;
    const sx = W / (xmax - xmin), sy = H / (ymax - ymin);
    const s = Math.min(sx, sy);
    const ox = (W - s * (xmax - xmin)) / 2;
    const oy = (H - s * (ymax - ymin)) / 2;
    const toX = (re) => ox + (re - xmin) * s;
    const toY = (im) => H - (oy + (im - ymin) * s);    // flip Y so +imag is up
    // Filled interior.
    ctx.beginPath();
    ctx.moveTo(toX(pts[0].re), toY(pts[0].im));
    for (let i = 1; i < pts.length; i++) ctx.lineTo(toX(pts[i].re), toY(pts[i].im));
    ctx.closePath();
    ctx.fillStyle = 'rgba(86, 170, 120, 0.32)';
    ctx.fill();
    ctx.strokeStyle = '#5acc8e';
    ctx.lineWidth = 1.2;
    ctx.stroke();
    // Poles (small dots) if present on phi.poles.
    if (phi.poles) {
      ctx.fillStyle = '#ddd';
      for (const p of phi.poles) {
        const a = (p && p.a) || (p && p.pole) || null;
        if (!a) continue;
        const x = toX(a.re), y = toY(a.im);
        if (x < 0 || x > W || y < 0 || y > H) continue;
        ctx.beginPath(); ctx.arc(x, y, 2, 0, Math.PI * 2); ctx.fill();
      }
    }
  }

  // -------------------------------------------------------------------------
  // Help-button content (HANDOFF #33). Wired by attachHelpButtons after
  // all cards are mounted; finds each by anchor selector.
  // -------------------------------------------------------------------------
  function attachHelpButtons() {
    if (!global.QD || !global.QD.QoL || !global.QD.QoL.attachHelp) return;
    const H = global.QD.QoL.attachHelp;
    const hOf = (cardId) => {
      const card = document.getElementById(cardId);
      return card ? card.querySelector('h2') : null;
    };
    H(hOf('ps-intro-card'),
      `<b>Parameter-slice cartography.</b> Pick one or two parameters of the
       current QD scenario and sweep them across a range. Each pixel runs an
       independent inverse-solve and colors by classification. Click any pixel
       to send that φ into the QD tab.`);
    H(hOf('ps-scenario-card'),
      `<b>Scenario snapshot</b>: the h(w) currently set in the QD tab. The
       parameter slice runs against this base scenario, varying only the
       parameters you pick on the Axes card.`);
    H(hOf('ps-quad-card'),
      `Shows the analytic form of h(w) (a sum of rational and polynomial
       terms) and which axis varies which parameter. Copy the formula via the
       button next to it.`);
    H(hOf('ps-axis-card'),
      `<b>Axes.</b> Pick one parameter for the X axis and optionally a second
       for the Y axis. <code>n</code> is the resolution along each axis. The
       coarse pass renders at ~16×16 and refines adaptively until <code>n</code>
       is reached.`);
    H(hOf('ps-run-card'),
      `<b>Run.</b> <i>Quality</i> sets how rigorously the identity check
       verifies each pixel: <code>Fast</code> (N=32, tol=1e-5) renders quickly
       but may yield false-positive <i>identity-fail</i> pixels;
       <code>Standard</code> (N=128, tol=1e-6) matches the inverse tab's
       tolerance; <code>Rigorous</code> (N=512, tol=1e-7) is publication-grade
       but ~3-4× slower.`);
    H(hOf('ps-hover-card'),
      `Live preview of the QD at the pixel you're hovering over the slice.
       Drawn from the cached φ closest to the cursor; if no cached φ exists
       (e.g. for an <i>invalid</i> pixel) the preview is empty. Click
       "Send to inverse" to push the scenario into the QD tab.`);
    H(hOf('ps-legend-card'),
      `<b>Legend.</b> Each color is a per-pixel classification outcome.
       Brightness inside <i>Valid QD</i> scales with iteration count
       (brighter = fewer Newton iterations).`);
  }

  // ---------------------------------------------------------------------------
  // Scenario capture: pull from QD-tab ui.js state (window.QD_UI.snapshotScenario)
  // ---------------------------------------------------------------------------
  function snapshotScenario() {
    if (typeof global.QD_UI !== 'object' || !global.QD_UI.snapshotScenario) return null;
    try { return global.QD_UI.snapshotScenario(); }
    catch (e) { return null; }
  }

  function refreshScenarioStatus() {
    const el = document.getElementById('ps-scenario-status');
    if (!el) return;
    const snap = snapshotScenario();
    if (!snap) { el.textContent = '(scenario unavailable — open the QD tab once first)'; return; }
    const { hData, norm, mode } = snap;
    const nPoles = hData.poles.length;
    const polyDeg = (hData.polyPart && hData.polyPart.length) ? hData.polyPart.length - 1 : -1;
    const parts = [`mode: <b>${mode}</b>`, `poles: ${nPoles}`, `poly degree: ${polyDeg}`];
    if (norm.c != null)  parts.push(`c=${(+norm.c).toFixed(3)}`);
    if (norm.w0)         parts.push(`w₀=(${norm.w0.re.toFixed(2)},${norm.w0.im.toFixed(2)})`);
    if (norm.q)          parts.push(`q=(${norm.q.re.toFixed(2)},${norm.q.im.toFixed(2)})`);
    el.innerHTML = parts.join('  •  ');
  }

  // ---------------------------------------------------------------------------
  // Axis picker
  // ---------------------------------------------------------------------------
  function refreshAxisOptions() {
    const snap = snapshotScenario();
    if (!snap) return;
    const params = PS.listAvailableParams({ hData: snap.hData, norm: snap.norm }, snap.mode);
    const xSel = document.getElementById('ps-axis-x');
    const ySel = document.getElementById('ps-axis-y');
    if (!xSel || !ySel) return;
    const prevX = xSel.value, prevY = ySel.value;
    xSel.innerHTML = params.map((p, i) =>
      `<option value="${i}">${p.label}</option>`).join('');
    ySel.innerHTML = `<option value="">— (1-D sweep)</option>` + params.map((p, i) =>
      `<option value="${i}">${p.label}</option>`).join('');
    // Preserve previous selection by label match if possible.
    const restore = (sel, prev) => {
      if (!prev) return;
      for (const opt of sel.options) {
        if (opt.textContent === prev || opt.value === prev) { sel.value = opt.value; return; }
      }
    };
    restore(xSel, prevX);
    restore(ySel, prevY);
    onAxisChange();
  }

  function onAxisChange() {
    const ySel = document.getElementById('ps-axis-y');
    document.getElementById('ps-y-range-row').style.display = ySel.value ? '' : 'none';
    // Default range to (cur ± 1) when the axis changes and the current
    // user-entered range hasn't been touched in this session.
    const snap = snapshotScenario();
    if (!snap) return;
    const params = PS.listAvailableParams({ hData: snap.hData, norm: snap.norm }, snap.mode);
    const xIdx = +document.getElementById('ps-axis-x').value;
    const xParam = params[xIdx];
    if (xParam) {
      const cur = PS.readParam(snap, xParam.ref);
      if (Math.abs(+document.getElementById('ps-x-min').value - (cur - 1)) > 1e9 ||
          isNaN(+document.getElementById('ps-x-min').value)) {
        // leave the field alone — user owns it after first interaction
      }
      // Only auto-set defaults when fields are still at the initial -1/1.
      const xMin = document.getElementById('ps-x-min');
      const xMax = document.getElementById('ps-x-max');
      if (xMin.dataset.userTouched !== '1') xMin.value = (cur - 1).toFixed(3);
      if (xMax.dataset.userTouched !== '1') xMax.value = (cur + 1).toFixed(3);
    }
    if (ySel.value) {
      const yParam = params[+ySel.value];
      if (yParam) {
        const cur = PS.readParam(snap, yParam.ref);
        const yMin = document.getElementById('ps-y-min');
        const yMax = document.getElementById('ps-y-max');
        if (yMin.dataset.userTouched !== '1') yMin.value = (cur - 1).toFixed(3);
        if (yMax.dataset.userTouched !== '1') yMax.value = (cur + 1).toFixed(3);
      }
    }
    refreshQuadratureDataCard();          // HANDOFF #33
  }
  // Mark range fields as user-touched on first edit so we stop overwriting them.
  document.addEventListener('input', e => {
    if (e.target && /^ps-[xy]-(min|max|n)$/.test(e.target.id)) {
      e.target.dataset.userTouched = '1';
    }
  });
  // QoL (HANDOFF #33): Enter in any axis-range input triggers a render
  // so the user can iterate axis params without reaching for the button.
  document.addEventListener('keydown', e => {
    if (e.key !== 'Enter' || !e.target) return;
    if (!/^ps-[xy]-(min|max|n)$/.test(e.target.id)) return;
    e.preventDefault();
    const runBtn = document.getElementById('ps-run');
    if (runBtn && !runBtn.disabled) runBtn.click();
  });

  // ---------------------------------------------------------------------------
  // Canvas handling
  // ---------------------------------------------------------------------------
  function getCanvas() { return document.getElementById('canvas'); }

  function paintImage(imageData) {
    const canvas = getCanvas();
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    // Center the slice image in the canvas. Image is pre-scaled to the
    // canvas's CSS-pixel square; we just blit.
    ctx.fillStyle = '#1a1a1a';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    const dx = Math.max(0, Math.floor((canvas.width  - imageData.width)  / 2));
    const dy = Math.max(0, Math.floor((canvas.height - imageData.height) / 2));
    ctx.putImageData(imageData, dx, dy);
    sliceState.canvasDestX = dx;
    sliceState.canvasDestY = dy;
  }

  function repaint() {
    if (sliceState.lastImageData) paintImage(sliceState.lastImageData);
  }

  let _clickAttached = false;
  let _hoverDetach   = null;
  function attachCanvasHandler() {
    if (_clickAttached) return;
    const canvas = getCanvas();
    if (!canvas) return;
    canvas.addEventListener('click', onCanvasClick);
    _clickAttached = true;
    // Hover tooltip + Hovered-QD preview wiring (HANDOFF #33).
    if (global.QD && global.QD.QoL && global.QD.QoL.attachHoverTooltip) {
      _hoverDetach = global.QD.QoL.attachHoverTooltip(canvas, hoverFormatter);
    }
  }
  function detachCanvasHandler() {
    if (!_clickAttached) return;
    const canvas = getCanvas();
    if (canvas) canvas.removeEventListener('click', onCanvasClick);
    _clickAttached = false;
    if (_hoverDetach) { _hoverDetach(); _hoverDetach = null; }
    refreshHoveredQDCard(null);
  }

  // Decode a CSS-pixel position on the canvas into a slice cell, returning
  // { col, row, xVal, yVal, xAxis, yAxis } or null if outside the slice
  // image. Used by both the click handler and the hover formatter so they
  // stay in lockstep.
  function pixelToCell(cssX, cssY) {
    if (!sliceState.lastImageData || !sliceState.lastAxes) return null;
    const dpr = window.devicePixelRatio || 1;
    const cx = cssX * dpr;
    const cy = cssY * dpr;
    const px = cx - (sliceState.canvasDestX || 0);
    const py = cy - (sliceState.canvasDestY || 0);
    const W = sliceState.lastImageData.width;
    const H = sliceState.lastImageData.height;
    if (px < 0 || py < 0 || px >= W || py >= H) return null;
    const axes = sliceState.lastAxes;
    const xAxis = axes[0];
    const yAxis = axes.length === 2 ? axes[1] : null;
    const colN = xAxis.n;
    const rowN = yAxis ? yAxis.n : 1;
    const col = Math.min(colN - 1, Math.max(0, Math.floor((px / W) * colN)));
    const row = yAxis ? Math.min(rowN - 1, Math.max(0, Math.floor((py / H) * rowN))) : 0;
    const xVal = lerp(xAxis.min, xAxis.max, col / Math.max(1, colN - 1));
    const yVal = yAxis ? lerp(yAxis.min, yAxis.max, row / Math.max(1, rowN - 1)) : null;
    return { col, row, xVal, yVal, xAxis, yAxis };
  }

  function onCanvasClick(e) {
    const canvas = getCanvas();
    if (!canvas) return;
    const rect = canvas.getBoundingClientRect();
    const cell = pixelToCell(e.clientX - rect.left, e.clientY - rect.top);
    if (!cell) return;
    loadPixel(cell.col, cell.row, cell.xVal, cell.yVal);
  }

  // Build the tooltip HTML for the cell at (cssX, cssY) on the canvas, and
  // sync the Hovered-QD preview card to the same cell. Returns null when
  // the cursor is off the slice (so the tooltip hides).
  function hoverFormatter(cssX, cssY) {
    const cell = pixelToCell(cssX, cssY);
    if (!cell) { refreshHoveredQDCard(null); return null; }
    const cellInfo = cellSummary(cell);
    refreshHoveredQDCard(cellInfo);
    const parts = [];
    parts.push(`<span class="tt-row">${escapeHTML(cell.xAxis.label || 'x')} = ${fmtNum(cell.xVal)}</span>`);
    if (cell.yAxis) {
      parts.push(`<span class="tt-row">${escapeHTML(cell.yAxis.label || 'y')} = ${fmtNum(cell.yVal)}</span>`);
    }
    if (cellInfo) {
      const lbl = PS.CLASS_LABELS[cellInfo.cls] || cellInfo.cls;
      const iterPart = (cellInfo.cls === PS.CLASS_VALID) ? ` · ${cellInfo.iter} iters` : '';
      parts.push(`<span class="tt-class">${escapeHTML(lbl)}${iterPart}</span>`);
    } else {
      parts.push(`<span class="tt-class">(no sample)</span>`);
    }
    return parts.join('');
  }

  function cellSummary(cell) {
    const dims = sliceState.gridDims;
    if (!dims || !sliceState.classGrid) return null;
    // The class/iter grids are sized to the adaptive-mesh resolution
    // (n0 × n1), which equals the slice image width × height. cell.col /
    // cell.row are in axis-space; map them into grid-space the same way.
    const colN = cell.xAxis.n;
    const rowN = cell.yAxis ? cell.yAxis.n : 1;
    if (dims.n0 !== colN || dims.n1 !== rowN) return null;
    const idx = cell.row * dims.n0 + cell.col;
    const k = sliceState.classGrid[idx];
    if (k === PS.UNKNOWN_CLASS) return null;
    const cls = PS.IDX_TO_CLASS[k];
    const iter = sliceState.iterGrid[idx];
    return { col: cell.col, row: cell.row, xVal: cell.xVal, yVal: cell.yVal,
             xAxis: cell.xAxis, yAxis: cell.yAxis, cls, iter };
  }

  function fmtNum(v) {
    if (v == null || !isFinite(v)) return String(v);
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-3 || a >= 1e5)) return v.toExponential(3);
    return v.toFixed(4);
  }
  // escapeHTML: delegates to QD.QoL.escapeHTML (HANDOFF #35 consolidation).
  // Falls back to a local impl if qol.js failed to load (defensive).
  function escapeHTML(s) {
    return (global.QD && global.QD.QoL && global.QD.QoL.escapeHTML)
      ? global.QD.QoL.escapeHTML(s)
      : String(s).replace(/[&<>"']/g, c => (
          { '&':'&amp;', '<':'&lt;', '>':'&gt;', '"':'&quot;', "'":'&#39;' }[c]));
  }

  function loadPixel(col, row, xVal, yVal) {
    const scenario = sliceState.lastScenario;
    if (!scenario) return;
    // Mutate a scenario clone with the pixel's parameter values.
    const s = PS.cloneScenario(scenario);
    PS.applyParamInPlace(s, sliceState.lastAxes[0].ref, xVal);
    if (sliceState.lastAxes[1]) {
      PS.applyParamInPlace(s, sliceState.lastAxes[1].ref, yVal);
    }
    // Hand to ui.js. It will switch tabs + re-solve.
    if (global.QD_UI && typeof global.QD_UI.loadScenarioIntoQdTab === 'function') {
      global.QD_UI.loadScenarioIntoQdTab(s, scenario.mode);
    }
  }

  function lerp(a, b, t) { return a + (b - a) * t; }

  // ---------------------------------------------------------------------------
  // Render orchestration
  // ---------------------------------------------------------------------------
  async function startRun() {
    const snap = snapshotScenario();
    if (!snap) { setProgress('No scenario.'); return; }
    const xSel = document.getElementById('ps-axis-x');
    const ySel = document.getElementById('ps-axis-y');
    if (!xSel.value && xSel.value !== '0') { setProgress('Pick an X axis.'); return; }

    const params = PS.listAvailableParams({ hData: snap.hData, norm: snap.norm }, snap.mode);
    const xRef = params[+xSel.value].ref;
    const yRef = ySel.value ? params[+ySel.value].ref : null;

    const xAxis = readAxis('x');
    const yAxis = yRef ? readAxis('y') : null;
    if (!xAxis || (yRef && !yAxis)) { setProgress('Invalid axis range.'); return; }
    xAxis.ref = xRef;
    xAxis.label = PS.formatParamLabel(xRef);
    if (yAxis) { yAxis.ref = yRef; yAxis.label = PS.formatParamLabel(yRef); }
    const axes = yAxis ? [xAxis, yAxis] : [xAxis];

    // Build base scenario for workers. Per-pixel opts are still cheaper
    // than the full Inverse-tab solver (no continuation, no diverse seeds),
    // but the identity-check rigor is now user-tunable via the Run-card
    // Quality dropdown — see QUALITY_PRESETS / readQualityPreset above and
    // HANDOFF #32 for the false-positive-fix rationale.
    const quality = readQualityPreset();
    const baseScenario = {
      hData: snap.hData,
      norm:  snap.norm,
      opts: {
        numRestarts: 1,                              // warm-start chain replaces multistart
        univalenceSamples: quality.univalenceSamples,
        identityTol:       quality.identityTol,
        findAlternates: false,
        newton: { maxIter: 40, tolerance: 1e-9 },
        usePhases: {
          direct: true,
          continuation: false,            // warm-start IS the continuation between pixels
          multistart: true,               // kept for cold pixels with no warm hint
          diverse: false,
          deflation: false,
        },
      },
    };

    // Allocate the image buffer (canvas-square sized to the smaller of
    // canvas-min-dim and a sensible cap so 256² doesn't blow up memory).
    const canvas = getCanvas();
    const W = Math.min(canvas.width, canvas.height, 1024);
    const H = W;
    const img = new ImageData(W, H);
    // Fill with dark background so unsampled regions are visible while rendering.
    for (let i = 0; i < img.data.length; i += 4) {
      img.data[i] = 24; img.data[i+1] = 24; img.data[i+2] = 24; img.data[i+3] = 255;
    }
    sliceState.lastImageData = img;
    sliceState.lastAxes      = axes;
    sliceState.lastScenario  = { ...snap, mode: snap.mode };
    sliceState.lastTiles     = new Array(axes[axes.length === 2 ? 1 : 0].n).fill(null);
    paintImage(img);

    // Ensure pool ready. createPool falls back to a main-thread pool when
    // the Worker bundle can't be built (e.g. fetch blocked on file://).
    let pool;
    try { pool = await ensurePool(); }
    catch (e) { setProgress('Pool init failed: ' + e.message); return; }
    if (pool.kind === 'main-thread') {
      setProgress('Note: running on main thread (Worker bundle unavailable — see console).');
    }

    setRunningUI(true);
    sliceState._errSamples = new Map();   // reset per-run error log
    sliceState.progress.t0 = performance.now();
    sliceState.progress.done = 0;
    sliceState.progress.total = axes.length === 2 ? axes[1].n : 1;
    updateProgressText();

    const colN = xAxis.n;
    const rowN = axes.length === 2 ? axes[1].n : 1;
    const xs = sampleAxis(xAxis);
    const ys = axes.length === 2 ? sampleAxis(axes[1]) : [null];

    // Shared paint helpers — paint a (col, row) grid block at the given
    // stride into the image buffer.
    const cellPxW = (colN > 0) ? (W / colN) : W;
    const cellPxH = (rowN > 0) ? (H / rowN) : H;
    function paintCellBlock(col, row, blockCols, blockRows, color) {
      const x0 = Math.floor(col * cellPxW);
      const y0 = Math.floor(row * cellPxH);
      const x1 = Math.min(W, Math.floor((col + blockCols) * cellPxW));
      const y1 = Math.min(H, Math.floor((row + blockRows) * cellPxH));
      for (let y = y0; y < y1; y++) {
        let idx = (y * W + x0) * 4;
        for (let x = x0; x < x1; x++) {
          img.data[idx] = color[0];
          img.data[idx+1] = color[1];
          img.data[idx+2] = color[2];
          img.data[idx+3] = 255;
          idx += 4;
        }
      }
    }
    function logErrorSamples(results) {
      for (const r of results) {
        if (r && r.errSample && sliceState._errSamples.size < 5
            && !sliceState._errSamples.has(r.errSample)) {
          sliceState._errSamples.set(r.errSample, 1);
          console.warn('[param-slice] first occurrence of error: ' + JSON.stringify(r.errSample));
        }
      }
    }

    // Cancellation handle the run loop checks between dispatches.
    const cancelToken = { cancelled: false };
    sliceState.activeJob = {
      cancel: () => { cancelToken.cancelled = true; if (sliceState.pool && sliceState.pool.cancel) sliceState.pool.cancel(); },
    };

    const tRunStart = performance.now();
    try {
      if (rowN === 1) {
        // ---- 1-D path: single linear sweep, no adaptive subdivision ----
        const points = xs.map(v => [{ ref: xAxis.ref, value: v }]);
        sliceState.progress.total = 1;
        const results = await pool.solveBatch(baseScenario, snap.mode, points, null);
        if (!cancelToken.cancelled && results) {
          logErrorSamples(results);
          for (let col = 0; col < colN; col++) {
            paintCellBlock(col, 0, 1, 1, PS.colorFor(results[col]));
          }
          paintImage(img);
          sliceState.progress.done = 1;
          updateProgressText();
        }
      } else {
        // ---- 2-D path: adaptive quadtree mesh refinement ----
        await runAdaptive2D({
          pool, baseScenario, mode: snap.mode, axes, xs, ys,
          n0: colN, n1: rowN,
          paintCellBlock, paintImage: () => paintImage(img),
          logErrorSamples, cancelToken,
          onProgress: (done, total) => {
            sliceState.progress.done = done;
            sliceState.progress.total = total;
            updateProgressText();
          },
        });
      }
      const dt = (performance.now() - tRunStart) / 1000;
      let extra = '';
      if (sliceState._errSamples.size > 0) {
        const top = Array.from(sliceState._errSamples.keys())[0];
        extra = ` — first error: ${top.slice(0, 120)} (see console for details)`;
      }
      setProgress(`Done in ${dt.toFixed(1)} s${extra}`);
    } catch (e) {
      console.warn('[param-slice] render error', e);
      setProgress('Render error: ' + (e.message || e));
    } finally {
      setRunningUI(false);
      sliceState.activeJob = null;
    }
  }

  // ---------------------------------------------------------------------------
  // Adaptive 2-D renderer
  // ---------------------------------------------------------------------------
  // Quadtree subdivision: render coarsely first (≈16×16 samples), then refine
  // only cells whose 4 corners disagree in classification OR whose iter-count
  // spread inside a uniformly-VALID block exceeds REFINE_ITER_DELTA (the
  // iter-gradient trigger — see param-slice-common.js's `cellIsHomogeneous`).
  // Each new sample gets an explicit warm-start hint from its nearest
  // previously-evaluated valid neighbour (cross-cell warm-start chaining)
  // and the hint's `_coarseIter` lets the solver speculatively tighten
  // its Newton maxIter cap.
  //
  // `n0 × n1` is the target sample resolution. The image buffer is laid out
  // independently; paintCellBlock maps cell coords to pixel coords.
  //
  // Refinement-trigger knob — change here to tune visible blockiness vs
  // sample count. 8 iters ≈ a ~10–15% perceived brightness shift in the
  // VALID-class colormap.
  const REFINE_ITER_DELTA = 8;
  // Mini-canvas boundary-sample count (HANDOFF #35). 128 is indistinguishable
  // from the inverse tab's 512+ render at 160×160 px; bump if jagged edges
  // are visible on high-DPR displays.
  const MINI_BOUNDARY_SAMPLES = 128;

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
        if (!hit) return null;
        return Object.assign({}, hit.phi, { _coarseIter: hit.iterCount });
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
  }

  // Sample an axis range into n evenly-spaced values.
  function sampleAxis(axis) {
    const { min, max, n } = axis;
    if (n === 1) return [(min + max) / 2];
    const out = new Array(n);
    const step = (max - min) / (n - 1);
    for (let i = 0; i < n; i++) out[i] = min + i * step;
    return out;
  }

  function cancelRun() {
    if (sliceState.activeJob) sliceState.activeJob.cancel();
    // HANDOFF #35: drop the stale φ cache so the hover-preview
    // doesn't show a QD from the aborted render.
    sliceState.nearestPhi = function () { return null; };
    sliceState.miniCache = { phi: null, pts: null };
    setProgress('Cancelled.');
  }

  function readAxis(which) {
    const min = +document.getElementById(`ps-${which}-min`).value;
    const max = +document.getElementById(`ps-${which}-max`).value;
    const n   = Math.max(2, Math.min(512, Math.floor(+document.getElementById(`ps-${which}-n`).value || 64)));
    if (!isFinite(min) || !isFinite(max) || min === max) return null;
    return { min: Math.min(min, max), max: Math.max(min, max), n };
  }

  function setRunningUI(running) {
    const runBtn = document.getElementById('ps-run');
    const cancelBtn = document.getElementById('ps-cancel');
    if (runBtn)    runBtn.disabled    = running;
    if (cancelBtn) cancelBtn.disabled = !running;
  }

  function setProgress(msg) {
    const el = document.getElementById('ps-progress');
    if (el) el.textContent = msg;
  }

  function updateProgressText() {
    const { done, total, t0 } = sliceState.progress;
    const dt = (performance.now() - t0) / 1000;
    const rate = done / Math.max(0.001, dt);
    const eta = rate > 0 ? Math.max(0, (total - done) / rate) : 0;
    setProgress(`row ${done}/${total}  •  ${dt.toFixed(1)}s elapsed  •  ${eta.toFixed(1)}s ETA`);
  }

  async function ensurePool() {
    if (sliceState.pool) return sliceState.pool;
    if (!sliceState.poolPromise) {
      sliceState.poolPromise = global.ParamSlicePool.create();
    }
    sliceState.pool = await sliceState.poolPromise;
    return sliceState.pool;
  }
})(typeof window !== 'undefined' ? window : globalThis);
