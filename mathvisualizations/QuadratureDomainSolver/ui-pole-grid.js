// =============================================================================
// ui-pole-grid.js -- Pole-control + polynomial-coefficient grid renderers.
//
// Extracted from ui.js by the Phase-3 UI modularization (item E). Exposes a
// QD_UI.installPoleGrid(uiCtx) factory; ui.js captures:
//   ({ renderPolesList, renderPolyCoefList } =
//        window.QD_UI.installPoleGrid(uiCtx));
//
// These two functions only BUILD DOM from state; all the slider-bookkeeping
// helpers they use (residueKey / magMaxFor / fmtArg / the magSliderMax cache),
// plus escapeHTML / escapeAttr / sub, stay in ui.js (they're shared with the
// pole input handlers + other cards) and are read off `ui`. The h-text mirror
// helpers (modeAllowsPoly / refreshHText) live in ui-h-text.js and are read as
// `ui.*` AT CALL TIME (the two modules call each other; this avoids any
// install-order dependency). `Complex` / `QD` are solver globals.
// =============================================================================

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  global.QD_UI.installPoleGrid = function installPoleGrid(ui) {
    const state       = ui.state;
    const $           = ui.$;
    const sub         = ui.sub;
    const escapeHTML  = ui.escapeHTML;
    const escapeAttr  = ui.escapeAttr;
    const residueKey  = ui.residueKey;
    const magMaxFor   = ui.magMaxFor;
    const fmtArg      = ui.fmtArg;

    // ---------- Render the pole controls ------------------------------------
    function renderPolesList() {
      const list = $('#poles-list');
      list.innerHTML = '';

      state.poles.forEach((pole, idx) => {
        // Each pole is a collapsible <details>, collapsed by default. The summary
        // shows the pole index + its location so collapsed poles stay identifiable.
        // Event delegation (#poles-list) targets `.pole` via closest(), so a
        // <details class="pole"> works unchanged.
        const div = document.createElement('details');
        div.className = 'pole';
        div.dataset.idx = idx;
        div.innerHTML = `
          <summary class="pole-header">
            <span class="pole-num">Pole ${idx + 1}</span>
            <span class="pole-loc">a = ${escapeHTML(pole.a)}</span>
            <button type="button" class="small danger" data-action="remove" title="Remove this pole">×</button>
          </summary>
          <div class="row">
            <label>a${sub(idx+1)} =
              <input type="text" class="cnum" data-field="a" value="${escapeAttr(pole.a)}"
                     aria-label="Pole ${idx + 1} location (complex)">
            </label>
          </div>
          <div class="row">
            <label>Order:
              <input type="number" min="1" max="6" value="${pole.order}" data-field="order" style="width: 56px;"
                     aria-label="Pole ${idx + 1} order">
            </label>
          </div>
          <div class="residues"></div>
        `;
        const residuesEl = $('.residues', div);
        for (let s = 0; s < pole.order; s++) {
          const cval = Complex.parse(pole.residues[s] || '0') || { re: 0, im: 0 };
          const mag = Math.hypot(cval.re, cval.im);
          const arg = Math.atan2(cval.im, cval.re);
          const key = residueKey(idx, s);
          const magMax = magMaxFor(key, mag);

          const block = document.createElement('div');
          block.className = 'residue-block';
          block.dataset.s = s;
          block.innerHTML = `
            <div class="residue-row">
              <span class="label-fixed">C${sub(idx+1)}${sub(s+1)}</span>
              =
              <input type="text" class="cnum residue" data-field="residue" data-s="${s}" value="${escapeAttr(pole.residues[s] || '')}"
                     aria-label="Pole ${idx + 1} residue C${idx + 1},${s + 1} (complex)">
            </div>
            <div class="slider1d-row">
              <label>|C|</label>
              <input type="range" class="slider1d slider1d-mag" data-s="${s}"
                     min="0" max="${magMax}" step="any" value="${mag}"
                     aria-label="Pole ${idx + 1} residue ${s + 1} magnitude">
              <span class="slider1d-val mag-val">${mag.toFixed(3)}</span>
            </div>
            <div class="slider1d-row">
              <label>arg</label>
              <input type="range" class="slider1d slider1d-arg" data-s="${s}"
                     min="${-Math.PI}" max="${Math.PI}" step="any" value="${arg}"
                     aria-label="Pole ${idx + 1} residue ${s + 1} argument">
              <span class="slider1d-val arg-val">${fmtArg(arg)}</span>
            </div>
          `;
          residuesEl.appendChild(block);
        }
        list.appendChild(div);
      });
      if (typeof ui.refreshHText === 'function') ui.refreshHText();
    }

    // Render the polynomial-part coefficient list. One block per C_{∞,l} for
    // l = 0..polyDegree, with magnitude/argument sliders matching the residue
    // rows. Visible in any mode where polynomial-h is meaningful (classical
    // unbounded + both unbounded-LQD variants — see modeAllowsPoly).
    function renderPolyCoefList() {
      const list = $('#poly-coefs-list');
      if (!list) return;
      list.innerHTML = '';
      const deg = state.polyDegree;
      if (!ui.modeAllowsPoly(state.mode) || deg < 0) return;

      // Ensure polyCoeffs has at least deg+1 entries (pad with '0').
      while (state.polyCoeffs.length < deg + 1) state.polyCoeffs.push('0');
      state.polyCoeffs.length = deg + 1;          // truncate any extras

      for (let l = 0; l <= deg; l++) {
        const cval = QD.Complex.parse(state.polyCoeffs[l] || '0') || { re: 0, im: 0 };
        const mag = Math.hypot(cval.re, cval.im);
        const arg = Math.atan2(cval.im, cval.re);
        const key = `poly-coef-${l}`;
        const magMax = magMaxFor(key, mag);
        const block = document.createElement('div');
        block.className = 'residue-block';
        block.dataset.polyL = l;
        block.innerHTML = `
          <div class="residue-row">
            <span class="label-fixed">C<sub>∞,${l}</sub></span>
            =
            <input type="text" class="cnum poly-coef" data-poly-l="${l}" value="${escapeAttr(state.polyCoeffs[l] || '')}"
                   aria-label="Polynomial-part coefficient C∞,${l} (complex)">
          </div>
          <div class="slider1d-row">
            <label>|C|</label>
            <input type="range" class="slider1d slider1d-poly-mag" data-poly-l="${l}"
                   min="0" max="${magMax}" step="any" value="${mag}"
                   aria-label="Polynomial coefficient ${l} magnitude">
            <span class="slider1d-val poly-mag-val">${mag.toFixed(3)}</span>
          </div>
          <div class="slider1d-row">
            <label>arg</label>
            <input type="range" class="slider1d slider1d-poly-arg" data-poly-l="${l}"
                   min="${-Math.PI}" max="${Math.PI}" step="any" value="${arg}"
                   aria-label="Polynomial coefficient ${l} argument">
            <span class="slider1d-val poly-arg-val">${fmtArg(arg)}</span>
          </div>
        `;
        list.appendChild(block);
      }
      if (typeof ui.refreshHText === 'function') ui.refreshHText();
    }

    return { renderPolesList, renderPolyCoefList };
  };
})(typeof window !== 'undefined' ? window : globalThis);
