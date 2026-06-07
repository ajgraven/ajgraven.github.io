// =============================================================================
// ui-h-text.js -- Custom h(w) text input ⇄ structured-grid mirror.
//
// Extracted from ui.js by the Phase-3 UI modularization (item E). Exposes a
// QD_UI.installHText(uiCtx) factory; ui.js captures:
//   ({ modeAllowsPoly, refreshHText, setHTextMsg, parseAndApplyHText } =
//        window.QD_UI.installHText(uiCtx));
//
// The #h-text input is a two-way-coupled mirror of the structured pole grid and
// polynomial-part coefficient list. refreshHText() rebuilds the text from state;
// parseAndApplyHText() goes the other direction via QD.parseH.
//
// Dependencies on ui.js / sibling modules are read as `ui.*` AT CALL TIME:
// renderPolesList / renderPolyCoefList (ui-pole-grid.js — these two modules call
// each other), and syncPolyDegreeInput / markAsCustom / scheduleSolve / the
// magSliderMax cache (all ui.js). `QD` is the solver namespace and `math` is the
// math.js CDN global (parser-only).
// =============================================================================

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  // Caps mirror the structured grid: pole order ∈ [1,6] (ui.js order input) and
  // polynomial degree ∈ [-1,6] (ui.js #poly-degree). The text path must enforce
  // the same limits so a pasted/typed `1/(w-2)^50` or huge polynomial can't slip
  // a high order past the grid into the solver (Taylor truncation + an O((2m)²)
  // dense Newton step per pole).
  const POLE_ORDER_MAX = 6;
  const POLY_DEGREE_MAX = 6;

  global.QD_UI.installHText = function installHText(ui) {
    const state = ui.state;

    // Polynomial part of h is meaningful exactly in the three unbounded family
    // panels. Keep this predicate centralized so refreshHText / parseAndApplyHText
    // agree with what the mode descriptors expose (cards.poly).
    function modeAllowsPoly(mode) {
      return mode === 'unbounded' ||
             mode === 'pqd-unbounded' ||
             mode === 'pqd-unbounded-singular' ||
             mode === 'lqd-unbounded' ||
             mode === 'lqd-unbounded-singular';
    }

    function refreshHText() {
      const inp = document.getElementById('h-text');
      if (!inp) return;
      try {
        const poles = state.poles.map(po => {
          const a = QD.Complex.parse(po.a) || { re: 0, im: 0 };
          const residues = po.residues.slice(0, po.order).map(r =>
            QD.Complex.parse(r) || { re: 0, im: 0 });
          return { a, order: po.order, residues };
        });
        let polyCoeffs = [];
        if (modeAllowsPoly(state.mode) && state.polyDegree >= 0) {
          polyCoeffs = state.polyCoeffs.slice(0, state.polyDegree + 1).map(s =>
            QD.Complex.parse(s) || { re: 0, im: 0 });
        }
        inp.value = QD.formatH({ poles, polyCoeffs });
        setHTextMsg('');
      } catch (e) {
        // Defensive: never let formatter errors break the panel.
      }
    }

    function setHTextMsg(msg, kind) {
      const el = document.getElementById('h-text-msg');
      if (!el) return;
      el.textContent = msg || '';
      el.style.color = (kind === 'warn') ? '#9a6a00' : '#b53030';
    }

    function parseAndApplyHText() {
      const inp = document.getElementById('h-text');
      if (!inp) return;
      const expr = inp.value.trim();
      if (!expr) { setHTextMsg('Enter an expression in w.'); return; }
      let parsed;
      try {
        parsed = QD.parseH(expr, math, { mode: state.mode });
      } catch (e) {
        setHTextMsg(e.message || String(e));
        return;
      }

      let clamped = false;   // true if any order / degree exceeded the grid cap

      // Convert parsed.poles (Complex-typed) back to the state's string form.
      if (parsed.poles.length === 0) {
        // Need at least one row in the grid so the user can extend it.
        state.poles = [{ a: '0', order: 1, residues: ['0'] }];
      } else {
        state.poles = parsed.poles.map(p => {
          let order = p.order;
          let residues = p.residues;
          if (order > POLE_ORDER_MAX) {
            clamped = true;
            order = POLE_ORDER_MAX;
            residues = residues.slice(0, POLE_ORDER_MAX);
          }
          return {
            a: QD.Complex.format(p.a),
            order,
            residues: residues.map(c => QD.Complex.format(c)),
          };
        });
      }

      if (modeAllowsPoly(state.mode)) {
        if (parsed.polyCoeffs.length > 0) {
          let coeffs = parsed.polyCoeffs;
          if (coeffs.length - 1 > POLY_DEGREE_MAX) {
            clamped = true;
            coeffs = coeffs.slice(0, POLY_DEGREE_MAX + 1);
          }
          state.polyCoeffs = coeffs.map(c => QD.Complex.format(c));
          state.polyDegree = coeffs.length - 1;
        } else {
          state.polyDegree = -1;
          state.polyCoeffs = [];
        }
        ui.syncPolyDegreeInput();
      }

      for (const k of Object.keys(ui.magSliderMax)) delete ui.magSliderMax[k];
      ui.renderPolesList();
      ui.renderPolyCoefList();
      ui.markAsCustom();
      if (clamped) {
        setHTextMsg('Truncated to the supported limit: pole order ≤ ' + POLE_ORDER_MAX +
          ', polynomial degree ≤ ' + POLY_DEGREE_MAX + '.', 'warn');
      } else if (parsed.warnings && parsed.warnings.length) {
        setHTextMsg('Parsed with warning: ' + parsed.warnings[0], 'warn');
      } else {
        setHTextMsg('');
      }
      ui.scheduleSolve();
    }

    return { modeAllowsPoly, refreshHText, setHTextMsg, parseAndApplyHText };
  };
})(typeof window !== 'undefined' ? window : globalThis);
