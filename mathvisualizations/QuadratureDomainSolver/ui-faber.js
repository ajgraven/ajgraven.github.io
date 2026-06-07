// =============================================================================
// ui-faber.js  —  Faber-polynomials analysis card (UQD).
//
// Installs via QD_UI.installFaber(uiCtx). For a classical UNBOUNDED QD
// (family 'unboundedQD') it computes the Faber polynomials F_n of the bounded
// complement K = ℂ∖Ω from φ's Laurent expansion at ∞ (QD.FaberAnalysis), shows
// them (formula list + expandable coefficient table) with capacity / degree /
// leading-coeff context and per-order convergence flags, and on demand pushes
// their roots to the main canvas via uiCtx.setFaberRoots (roots cluster in K).
//
// The card is hidden for any non-UQD solve (bounded, PQD, LQD) — those carry a
// power/Blaschke weight, so φ's Laurent is not the plain exterior-map expansion
// and the clean Faber-of-K identity does not apply. UI-only; no solver math.
// =============================================================================

(function () {
  'use strict';

  function esc(s) {
    return String(s).replace(/[&<>"]/g, c =>
      ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  }
  // Subscripts index an object (Fₙ, c₀); superscripts are exponents (ζ², 1/cⁿ).
  function subDigits(n) {
    const map = '₀₁₂₃₄₅₆₇₈₉';
    return String(n).split('').map(d => map[+d] || d).join('');
  }
  function supDigits(n) {
    const map = '⁰¹²³⁴⁵⁶⁷⁸⁹';
    return String(n).split('').map(d => map[+d] || d).join('');
  }
  // Format a complex value with the true minus sign (U+2212), matching the
  // formula list (formatFaberPoly) so the card doesn't mix '-' and '−'.
  function fmtC(z, digits) {
    return (window.QD.Complex.format(z, { digits }) || '').replace(/-/g, '−');
  }
  function num(x, p) {
    if (typeof x !== 'number' || !isFinite(x)) return String(x);
    return (Math.abs(x) >= 1e4 || (x !== 0 && Math.abs(x) < 1e-3))
      ? x.toExponential(3) : String(+x.toPrecision(p || 5));
  }

  function installFaber(ctx) {
    const QD = window.QD;
    const $ = ctx.$;
    const FA = QD && QD.FaberAnalysis;

    const card    = $('#faber-card');
    const content = $('#faber-content');
    if (!card || !content || !FA) return {};

    const orderInp  = $('#faber-order');
    const singleInp = $('#faber-single-n');
    const modeAll   = $('#faber-mode-all');
    const modeSingle= $('#faber-mode-single');
    const showRoots = $('#faber-show-roots');
    const showPolys = $('#faber-show-polys');
    const capSpan   = $('#faber-cap');

    let activePhi = null;                 // current UQD φ (null ⇒ card hidden / stale)
    let cache = null;                      // { phi, N, c, c0, coeffs, orders }
    let debounceTimer = null;

    // Help popover on the header. Prose lives in QD.Strings.faber (ui-strings.js).
    const STR = (QD.Strings && QD.Strings.faber) || {};
    if (QD.QoL && QD.QoL.attachHelp && STR.help) {
      const h = card.querySelector('h2');
      if (h) QD.QoL.attachHelp(h, STR.help);
    }

    function clampInt(v, lo, hi, dflt) {
      let n = parseInt(v, 10);
      if (!isFinite(n)) n = dflt;
      return Math.max(lo, Math.min(hi, n));
    }

    // Heavy step: build F_1..F_N and root-find each order. Cached on (phi, N).
    function compute(phi, N) {
      const { c, c0, coeffs } = FA.faberPolynomials(phi, N);
      const orders = [];
      for (let n = 1; n <= N; n++) {
        const r = FA.polynomialRoots(coeffs[n]);
        orders.push({ n, roots: r.roots, converged: r.converged });
      }
      return { phi, N, c, c0, coeffs, orders };
    }

    function pushRoots(N, mode, sn) {
      if (!showRoots.checked) { ctx.setFaberRoots(null); return; }
      if (mode === 'single') {
        const o = cache.orders[sn - 1];
        ctx.setFaberRoots({ mode: 'single', N, n: sn,
          sets: o ? [{ n: o.n, roots: o.roots, converged: o.converged }] : [] });
      } else {
        ctx.setFaberRoots({ mode: 'all', N,
          sets: cache.orders.map(o => ({ n: o.n, roots: o.roots, converged: o.converged })) });
      }
    }

    function leadingCoeffStr(c, n) {
      const lc = 1 / Math.pow(c, n);
      return '1/c' + supDigits(n) + ' ≈ ' + num(lc, 5);   // exponent → superscript
    }

    function render(N, mode, sn) {
      const c = cache.c, c0 = cache.c0, coeffs = cache.coeffs;
      if (capSpan) capSpan.textContent = 'cap(K) = c ≈ ' + num(c, 6);

      const parts = [];
      parts.push('<div class="geom-row"><span class="key">cap(K) = c:</span> ' + num(c, 6) + '</div>');
      parts.push('<div class="geom-row"><span class="key">c₀ (Laurent const):</span> ' +
        esc(fmtC(c0, 5)) + '</div>');

      // Polynomial display (formula list + expandable coefficient table for the focus order).
      if (showPolys.checked) {
        const focus = (mode === 'single') ? sn : N;
        const listMax = (mode === 'single') ? sn : N;
        const listMin = (mode === 'single') ? sn : 1;
        let formulas = '';
        for (let n = listMin; n <= listMax; n++) {
          formulas += '<div class="faber-formula">F' + subDigits(n) + '(ζ) = ' +
            esc(FA.formatFaberPoly(coeffs[n])) + '</div>';
        }
        parts.push('<div style="margin-top:6px;">' + formulas + '</div>');

        // Coefficient table for the focus polynomial.
        const Ff = coeffs[focus];
        let rows = '';
        for (let k = Ff.length - 1; k >= 0; k--) {
          rows += '<tr><td>ζ' + supDigits(k) + '</td><td>' +    // ζ^k → superscript exponent
            esc(fmtC(Ff[k], 6)) + '</td></tr>';
        }
        parts.push(
          '<details style="margin-top:6px;"><summary class="key" style="cursor:pointer;">' +
          'Coefficients of F' + subDigits(focus) + '</summary>' +
          '<table class="faber-coef"><tr><th>term</th><th>coefficient</th></tr>' + rows + '</table>' +
          '<div class="hint">degree ' + focus + ', leading coeff = ' + leadingCoeffStr(c, focus) + '</div>' +
          '</details>');
      }

      // Convergence flags (the one explicitly-requested property list).
      const bad = cache.orders.filter(o => (mode === 'single' ? o.n === sn : true) && !o.converged);
      if (bad.length === 0) {
        parts.push('<div class="geom-row"><span class="ok">✓</span> ' +
          '<span class="key">root-finder converged</span> for all orders shown</div>');
      } else {
        parts.push('<div class="geom-row"><span class="warn">⚠</span> ' +
          '<span class="key">non-convergence</span> at F' +
          bad.map(o => subDigits(o.n)).join(', F') +
          ' <span class="hint">(high-degree conditioning — roots approximate)</span></div>');
      }

      content.innerHTML = parts.join('');
    }

    function refresh() {
      if (!activePhi) return;
      const N = clampInt(orderInp.value, 1, 30, 6);
      if (String(N) !== orderInp.value) orderInp.value = String(N);
      const mode = (modeSingle && modeSingle.checked) ? 'single' : 'all';
      let sn = clampInt(singleInp.value, 1, N, Math.min(6, N));
      if (sn > N) sn = N;

      try {
        if (!cache || cache.phi !== activePhi || cache.N !== N) cache = compute(activePhi, N);
        render(N, mode, sn);
        pushRoots(N, mode, sn);
      } catch (e) {
        content.innerHTML = '<div class="warn">' + esc(STR.unavailablePrefix || 'Faber analysis unavailable: ') +
          esc((e && e.message) || String(e)) + '</div>';
        ctx.setFaberRoots(null);
      }
    }

    function refreshDebounced() {
      if (debounceTimer) clearTimeout(debounceTimer);
      debounceTimer = setTimeout(refresh, 150);
    }

    // Controls.
    if (orderInp)   orderInp.addEventListener('input', refreshDebounced);
    if (singleInp)  singleInp.addEventListener('input', refreshDebounced);
    if (modeAll)    modeAll.addEventListener('change', refresh);
    if (modeSingle) modeSingle.addEventListener('change', refresh);
    if (showPolys)  showPolys.addEventListener('change', refresh);
    if (showRoots)  showRoots.addEventListener('change', () => {
      // Keep the Layers toggle in sync, then recompute/push.
      const t = $('#faber-roots-toggle');
      if (t && t.checked !== showRoots.checked) { t.checked = showRoots.checked; }
      ctx.state.showFaberRoots = showRoots.checked;
      refresh();
    });

    // A manual edit (pole drag / preset) invalidates the current solution: clear
    // the on-plot roots and grey the card until the next solve lands. Keep the
    // card visible if the domain is still a UQD (the re-solve will refresh it).
    document.addEventListener('qd-customized', () => {
      cache = null;
      ctx.setFaberRoots(null);
      if (activePhi) content.innerHTML = '<div class="hint">' + esc(STR.pending || 'solving… Faber analysis pending') + '</div>';
    });

    // Refresh on every fresh primary solution. Show only for classical UQD.
    if (QD.PrimarySolution && QD.PrimarySolution.subscribe) {
      QD.PrimarySolution.subscribe((env) => {
        const phi = env && env.success && env.primary && env.primary.phi;
        // Classical unbounded QD only. The solver leaves phi.family UNSET for the
        // classical UQD; the power/log-weighted families (PQD/LQD) tag themselves
        // ('unboundedPQD'/'unboundedLQD'…) and carry a weight marker (alpha /
        // lqdBeta / z0 / gamma / q). Exclude those — their φ-Laurent is not the
        // plain exterior-map expansion, so the Faber-of-K identity doesn't hold.
        const isUQD = !!(phi && phi.unbounded
          && (!phi.family || phi.family === 'unboundedQD')
          && phi.alpha == null && phi.lqdBeta == null
          && phi.z0 == null && phi.gamma == null && phi.q == null);
        if (!isUQD) {
          activePhi = null; cache = null;
          card.classList.add('hidden');
          ctx.setFaberRoots(null);
          return;
        }
        activePhi = phi;
        card.classList.remove('hidden');
        refresh();
      });
    }

    return {};
  }

  window.QD_UI = window.QD_UI || {};
  window.QD_UI.installFaber = installFaber;
})();
