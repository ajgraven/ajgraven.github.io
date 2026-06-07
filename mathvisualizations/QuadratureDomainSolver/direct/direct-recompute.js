// =============================================================================
// direct-recompute.js -- Recompute-and-render pipeline for the Direct tab.
//
// Extracted from direct-ui.js by the Phase-3 UI modularization (item E).
// Exposes QD_UI.installDirectRecompute(dCtx); direct-ui.js captures
// recomputeAndRender into an IIFE-local binding so every card-builder handler +
// _activate call it unchanged.
//
// recomputeAndRender dispatches on directState.mode to recomputeBounded /
// recomputeUnbounded / recomputeNumerical, each of which builds h via the
// QD.Direct.* forward kernels, renders it (displayH), and pushes the sampled
// boundary to the shared plot (pushBoundaryToPlot). sampleBoundedPhi samples the
// weight-honoring boundary; the complex->string/katex helpers + formatNumLocal
// are display-only. Bodies are VERBATIM moves. Deps via dCtx: directState +
// parseComplex + isMounted (the host's mount-guard accessor). `QD` / `window` /
// `document` / `math` / `katex` are globals.
// =============================================================================
'use strict';

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  global.QD_UI.installDirectRecompute = function installDirectRecompute(d) {
    const directState  = d.directState;
    const parseComplex = d.parseComplex;
    const isMounted    = d.isMounted;

  // Thin wrappers around QD.Complex.format (moved with displayH, their only user).
  function complexToString(c)  { return QD.Complex.format(c); }
  function complexToKatex(c)   { return QD.Complex.format(c); }

  function sampleBoundedPhi(P, Q, N) {
    const evalRat = (z) => {
      const pv = QD.Direct.evalPolyAscending(P, z), qv = QD.Direct.evalPolyAscending(Q, z);
      const d2 = qv.re * qv.re + qv.im * qv.im;
      return { re: (pv.re * qv.re + pv.im * qv.im) / d2, im: (pv.im * qv.re - pv.re * qv.im) / d2 };
    };
    const w = directState.weight;
    const pts = new Array(N);
    // Singular: φ also carries the Blaschke factor b_{z₀}(z) (0 ∈ Ω).
    const sing = directState.singular && w !== 'classical';
    let z0 = null, absz0 = 1;
    if (sing) { try { z0 = parseComplex(directState.z0); absz0 = QD.Complex.abs(z0); } catch (e) { z0 = null; } }
    const cmul = (a, b) => ({ re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re });
    if (w === 'power') {
      const alpha = parseFloat(directState.alpha);
      const R0 = evalRat({ re: 0, im: 0 });
      const anchorArg0 = Math.atan2(R0.im, R0.re);
      const evalRHashRaw = (z) => evalRat(z);
      for (let n = 0; n < N; n++) {
        const t = 2 * Math.PI * n / N, z = { re: Math.cos(t), im: Math.sin(t) };
        const argZ = QD.PqdCommon.argContAt(null, z, evalRHashRaw, anchorArg0, { re: 0, im: 0 });
        const R = evalRat(z), mag = Math.pow(R.re * R.re + R.im * R.im, 0.5 / alpha);
        let pt = { re: mag * Math.cos(argZ / alpha), im: mag * Math.sin(argZ / alpha) };
        if (sing && z0) pt = cmul(QD.LqdCommon.blaschkeEval(z, z0), pt);   // ×b_{z₀}
        pts[n] = pt;
      }
      return pts;
    }
    if (w === 'log') {
      let w0; try { w0 = parseComplex(directState.logW0); } catch (e) { w0 = { re: 1, im: 0 }; }
      for (let n = 0; n < N; n++) {
        const t = 2 * Math.PI * n / N, z = { re: Math.cos(t), im: Math.sin(t) };
        const r = evalRat(z), e = Math.exp(r.re);
        const er = { re: e * Math.cos(r.im), im: e * Math.sin(r.im) };
        let pt = { re: w0.re * er.re - w0.im * er.im, im: w0.re * er.im + w0.im * er.re };
        // Singular: φ = γ·b·exp(r#) = (b/|z₀|)·(w₀·exp(r#)), since γ = w₀/|z₀|.
        if (sing && z0) { const bf = QD.LqdCommon.blaschkeEval(z, z0); pt = cmul({ re: bf.re / absz0, im: bf.im / absz0 }, pt); }
        pts[n] = pt;
      }
      return pts;
    }
    for (let n = 0; n < N; n++) {
      const t = 2 * Math.PI * n / N;
      pts[n] = evalRat({ re: Math.cos(t), im: Math.sin(t) });
    }
    return pts;
  }

  function recomputeAndRender() {
    if (!isMounted()) return;
    const root = document.getElementById('controls-direct');
    if (!root) return;

    const hDisp   = root.querySelector('.dir-h-display');
    const hKatex  = root.querySelector('.dir-h-katex');
    const errBox  = root.querySelector('.dir-error');
    if (hDisp) hDisp.textContent = '';
    if (hKatex) hKatex.innerHTML = '';
    if (errBox) errBox.textContent = '';

    if (directState.mode === 'bounded') {
      recomputeBounded(root, hDisp, hKatex, errBox);
    } else if (directState.mode === 'unbounded') {
      recomputeUnbounded(root, hDisp, hKatex, errBox);
    } else {
      recomputeNumerical(root, hDisp, hKatex, errBox);
    }

    // Empty-state guidance: with no h produced (blank input / nothing computed
    // yet) and no error shown, replace the blank output box with a one-line hint.
    if (hDisp && !hDisp.textContent.trim()
        && (!hKatex || !hKatex.innerHTML.trim())
        && (!errBox || !errBox.textContent.trim())) {
      hDisp.innerHTML = '<span style="color:var(--c-muted)">Enter or paste φ(z) above — the resulting h(w) appears here.</span>';
    }
  }

  function recomputeNumerical(root, hDisp, hKatex, errBox) {
    const card    = root.querySelector('.dir-phi-card-numerical');
    const status  = card && card.querySelector('.dir-phi-num-status');
    const msg     = card && card.querySelector('.dir-phi-num-msg');
    const diag    = card && card.querySelector('.dir-phi-num-diag');
    const warnBox = card && card.querySelector('.dir-phi-num-warnings');
    if (status)  status.textContent = '';
    if (msg)     msg.textContent = '';
    if (diag)    diag.textContent = '';
    if (warnBox) warnBox.textContent = '';

    const exprStr = directState.numExpr.trim();
    if (!exprStr) {
      if (msg) { msg.style.color = '#888'; msg.textContent = 'enter an expression'; }
      return;
    }
    if (typeof window === 'undefined' || !window.math || !window.math.parse) {
      if (errBox) errBox.textContent = 'math.js not loaded';
      return;
    }

    // Parse expression once; build a numeric phiFn.
    let node;
    try { node = window.math.parse(exprStr); }
    catch (err) {
      if (status) { status.textContent = '✗'; status.style.color = '#b53030'; }
      if (msg) { msg.style.color = '#b53030'; msg.textContent = 'parse: ' + (err.message || err); }
      return;
    }
    let compiled;
    try { compiled = node.compile(); }
    catch (err) {
      if (status) { status.textContent = '✗'; status.style.color = '#b53030'; }
      if (msg) { msg.style.color = '#b53030'; msg.textContent = 'compile: ' + (err.message || err); }
      return;
    }
    const phiFn = z => {
      const v = compiled.evaluate({ z: window.math.complex(z.re, z.im) });
      if (typeof v === 'number') return { re: v, im: 0 };
      if (v && typeof v.re === 'number' && typeof v.im === 'number') return { re: v.re, im: v.im };
      throw new Error('expression did not evaluate to a complex/number');
    };

    let result;
    try {
      result = QD.Direct.numericalBoundedQD(phiFn, {
        numSamples: 256,
        maxOrder: directState.numMaxOrder,
      });
    } catch (err) {
      if (status) { status.textContent = '✗'; status.style.color = '#b53030'; }
      if (msg) { msg.style.color = '#b53030'; msg.textContent = err.message || String(err); }
      return;
    }

    directState.lastH = result.hData;
    directState.lastC = 0;                          // bounded mode (numerical reduces to bounded)

    if (status) { status.textContent = '✓'; status.style.color = '#2a8f2a'; }
    if (msg) {
      msg.style.color = '#2a8f2a';
      msg.textContent = 'truncated at degree ' + result.truncationOrder
                      + ' (analyticity score = ' + result.analyticityScore.toExponential(2) + ')';
    }
    if (diag) {
      const lines = ['Recovered Taylor coefficients of φ at z=0:'];
      for (let k = 0; k <= Math.min(result.truncationOrder, 6); k++) {
        const c = result.taylorCoeffs[k];
        lines.push('  c_' + k + ' = ' + formatNumLocal(c.re) + (c.im >= 0 ? '+' : '') + formatNumLocal(c.im) + 'i');
      }
      if (result.truncationOrder > 6) lines.push('  …');
      diag.textContent = lines.join('\n');
    }
    if (warnBox && result.warnings.length) warnBox.textContent = '⚠ ' + result.warnings.join('; ');

    displayH(hDisp, hKatex, result.hData, /*isUnbounded=*/false);

    // Live ∂Ω preview: sample using the user's phiFn directly.
    try {
      const N = 400;
      const pts = new Array(N);
      for (let n = 0; n < N; n++) {
        const theta = 2 * Math.PI * n / N;
        pts[n] = phiFn({ re: Math.cos(theta), im: Math.sin(theta) });
      }
      pushBoundaryToPlot(pts, false);
    } catch (e) { /* preview is best-effort */ }
  }

  function formatNumLocal(x) {
    if (!isFinite(x)) return String(x);
    if (Math.abs(x) < 1e-12) return '0';
    if (Math.abs(x - Math.round(x)) < 1e-10) return String(Math.round(x));
    return Number(x.toPrecision(6)).toString();
  }

  function recomputeBounded(root, hDisp, hKatex, errBox) {
    const warnBox = root.querySelector('.dir-phi-warnings');
    if (warnBox) warnBox.textContent = '';
    directState.lastQ = null;   // reset singular origin residues (set below if applicable)
    directState.lastOriginRes = null;
    directState.lastWeight = directState.weight;   // family of this solve (Send/Verify dispatch)

    const weighted = directState.weight !== 'classical';
    if (directState.coeffsKind === 'rational' || weighted) {
      let P, Q;
      try {
        if (directState.coeffsKind === 'rational') {
          P = directState.coeffsNum.map(parseComplex);
          Q = directState.coeffsDen.map(parseComplex);
        } else {
          // Weighted mode with a polynomial-kind entry: kernel = polynomial / 1.
          P = directState.coeffs.map(parseComplex);
          Q = [{ re: 1, im: 0 }];
        }
      } catch (err) {
        if (errBox) errBox.textContent = 'Coefficient parse error: ' + err.message;
        return;
      }
      // Weighted (PQD / LQD) forward kernels treat P/Q as the rational KERNEL
      // R#/r#; classical treats it as φ = P/Q.
      let result;
      // Singular (0 ∈ Ω) forward is supported for both Power and Log weights.
      const sing = directState.singular && directState.weight !== 'classical';
      let z0 = null;
      if (sing) {
        try { z0 = parseComplex(directState.z0); } catch (e) { if (errBox) errBox.textContent = 'z₀ parse error: ' + e.message; return; }
      }
      if (directState.weight === 'power') {
        const alpha = parseFloat(directState.alpha);
        if (!(alpha > 0) || Math.abs(alpha - 1) < 1e-9) {
          if (errBox) errBox.textContent = 'α must be a real number > 0 and ≠ 1.';
          directState.lastH = null; return;
        }
        try { result = sing ? QD.Direct.boundedPowerQDSingular({ num: P, den: Q }, alpha, z0)
                            : QD.Direct.boundedPowerQD({ num: P, den: Q }, alpha); }
        catch (err) { if (errBox) errBox.textContent = err.message; directState.lastH = null; return; }
      } else if (directState.weight === 'log') {
        let w0;
        try { w0 = parseComplex(directState.logW0); } catch (e) { if (errBox) errBox.textContent = 'w₀ parse error: ' + e.message; return; }
        try { result = sing ? QD.Direct.boundedLogQDSingular({ num: P, den: Q }, w0, z0)
                            : QD.Direct.boundedLogQD({ num: P, den: Q }, w0); }
        catch (err) { if (errBox) errBox.textContent = err.message; directState.lastH = null; return; }
      } else {
        try { result = QD.Direct.boundedQDRational(P, Q); }
        catch (err) { if (errBox) errBox.textContent = err.message; directState.lastH = null; return; }
      }
      // Singular PQD realizability ⟺ univalence (Thm 4.3.3); report when φ fails it.
      if (result.univalent === false) {
        if (errBox) errBox.textContent = (result.warnings && result.warnings[0]) || 'Not realizable (φ not univalent).';
        directState.lastH = null; return;
      }
      directState.lastH = result.hData;
      directState.lastC = 0;
      // Stash the singular origin term (LQD q, or PQD r₀ = t − ΣC_{j,1}) + the built
      // φ (for the family-verifier Verify on the FULL h).
      directState.lastQ = (result.q && (result.q.re || result.q.im)) ? result.q : null;
      directState.lastOriginRes = (result.originResidue && (result.originResidue.re || result.originResidue.im)) ? result.originResidue : null;
      directState.lastPhi = result.phi || null;
      directState.lastSingular = sing;
      directState.lastWeight = directState.weight;     // for the Send/Verify dispatch
      const warns = (result.warnings || []).slice();
      const fmtc = (z) => z.re.toFixed(4) + (z.im >= 0 ? '+' : '') + z.im.toFixed(4) + 'i';
      if (directState.lastQ) warns.push('origin term: h has q/w with q = ' + fmtc(directState.lastQ));
      if (directState.lastOriginRes) warns.push('origin term: h has r₀/w with r₀ = ' + fmtc(directState.lastOriginRes) + ' (= ∫|w|^{2(α−1)}dA − ΣC)');
      if (warnBox) warnBox.textContent = warns.length ? '⚠ ' + warns.join('; ') : '';

      displayH(hDisp, hKatex, result.hData);
      pushBoundaryToPlot(sampleBoundedPhi(P, Q, 400), false);
      return;
    }

    // Polynomial path.
    let coeffs;
    try { coeffs = directState.coeffs.map(parseComplex); }
    catch (err) {
      if (errBox) errBox.textContent = 'Coefficient parse error: ' + err.message;
      return;
    }
    let result;
    try { result = QD.Direct.boundedQD(coeffs); }
    catch (err) {
      if (errBox) errBox.textContent = err.message;
      directState.lastH = null;
      return;
    }
    directState.lastH = result.hData;
    directState.lastC = 0;
    if (warnBox && result.warnings.length) warnBox.textContent = '⚠ ' + result.warnings.join('; ');

    displayH(hDisp, hKatex, result.hData);
    pushBoundaryToPlot(QD.Direct.sampleBoundaryPolynomial(coeffs, 400), false);
  }

  function recomputeUnbounded(root, hDisp, hKatex, errBox) {
    const warnBox = root.querySelector('.dir-phi-uns-warnings');
    if (warnBox) warnBox.textContent = '';
    directState.lastQ = null; directState.lastOriginRes = null;
    directState.lastPhi = null; directState.lastSingular = false;

    // ---- Weighted (PQD / LQD) unbounded forward kernels (Thm 4.3.7). The unified
    // weight/singular axes live in the Domain-type card. ----
    if (directState.weight === 'power' || directState.weight === 'log') {
      let P, Q;
      try { P = directState.unsKernelNum.map(parseComplex); Q = directState.unsKernelDen.map(parseComplex); }
      catch (err) { if (errBox) errBox.textContent = 'Kernel parse error: ' + err.message; return; }
      const rH = { num: P, den: Q };
      const sing = directState.singular;
      let result;
      try {
        if (directState.weight === 'power') {
          const alpha = parseFloat(directState.unsAlpha);
          if (!(alpha > 0) || Math.abs(alpha - 1) < 1e-9) { if (errBox) errBox.textContent = 'α must be > 0 and ≠ 1.'; directState.lastH = null; return; }
          if (sing) {
            let z0hint = null; try { z0hint = parseComplex(directState.unsZ0); } catch (e) { z0hint = null; }
            result = QD.Direct.unboundedPowerQDSingular(rH, alpha, z0hint ? { z0: z0hint } : {});
          } else {
            result = QD.Direct.unboundedPowerQD(rH, alpha);
          }
        } else { // log
          let c; try { const pc = parseComplex(directState.cValue); if (!(pc.re > 0) || Math.abs(pc.im) > 1e-12) throw new Error('c must be a positive real'); c = pc.re; }
          catch (e) { if (errBox) errBox.textContent = 'c parse error: ' + e.message; return; }
          if (sing) {
            let z0; try { z0 = parseComplex(directState.unsZ0); } catch (e) { if (errBox) errBox.textContent = 'z₀ parse error: ' + e.message; return; }
            result = QD.Direct.unboundedLogQDSingular(rH, c, z0);
          } else {
            result = QD.Direct.unboundedLogQD(rH, c);
          }
        }
      } catch (err) { if (errBox) errBox.textContent = err.message; directState.lastH = null; return; }

      if (result.univalent === false) {
        if (errBox) errBox.textContent = (result.warnings && result.warnings[0]) || 'Not realizable (φ not univalent).';
        directState.lastH = null; return;
      }
      if (result.originInside) {
        if (errBox) errBox.textContent = (result.warnings && result.warnings[0]) || '0 ∈ Ω — use the singular kernel.';
        directState.lastH = null; return;
      }
      directState.lastH = result.hData;
      directState.lastC = result.c;
      directState.lastPhi = result.phi || null;
      directState.lastSingular = sing;
      directState.lastWeight = directState.weight;
      directState.lastQ = (result.q && (result.q.re || result.q.im)) ? result.q : null;
      const warns = (result.warnings || []).slice();
      const fmtc = (z) => z.re.toFixed(4) + (z.im >= 0 ? '+' : '') + z.im.toFixed(4) + 'i';
      if (typeof result.c === 'number') warns.push('c = φ′(∞) = ' + result.c.toFixed(5));
      if (result.z0) warns.push('z₀ (origin preimage) = ' + fmtc(result.z0));
      if (directState.lastQ) warns.push('origin term: h has q/w with q = ' + fmtc(directState.lastQ));
      if (warnBox) warnBox.textContent = warns.length ? '⚠ ' + warns.join('; ') : '';

      displayH(hDisp, hKatex, result.hData, /*isUnbounded=*/true, result.c);
      // Boundary plot: sample φ on |z|=1 via the family evaluator.
      try {
        const fam = QD.Family[result.phi.family];
        const N = 400, pts = new Array(N);
        for (let n = 0; n < N; n++) { const th = 2 * Math.PI * n / N; pts[n] = fam.evalPhi({ re: Math.cos(th), im: Math.sin(th) }, result.phi); }
        pushBoundaryToPlot(pts, /*unbounded=*/true);
      } catch (e) { /* plot is best-effort */ }
      return;
    }

    // ---- Classical unbounded QD (Laurent c + F). ----
    directState.lastWeight = 'classical';
    let c;
    try {
      const parsed = parseComplex(directState.cValue);
      if (Math.abs(parsed.im) > 1e-12 || parsed.re <= 0 || !isFinite(parsed.re)) {
        throw new Error("c must be a positive real number");
      }
      c = parsed.re;
    } catch (err) {
      if (errBox) errBox.textContent = 'c parse error: ' + err.message;
      return;
    }

    let F;
    try { F = directState.Fcoeffs.map(parseComplex); }
    catch (err) {
      if (errBox) errBox.textContent = 'F coefficient parse error: ' + err.message;
      return;
    }

    let result;
    try { result = QD.Direct.unboundedQD(c, F); }
    catch (err) {
      if (errBox) errBox.textContent = err.message;
      directState.lastH = null;
      return;
    }
    directState.lastH = result.hData;
    directState.lastC = c;
    if (warnBox && result.warnings.length) warnBox.textContent = '⚠ ' + result.warnings.join('; ');

    displayH(hDisp, hKatex, result.hData, /*isUnbounded=*/true, c);
    pushBoundaryToPlot(QD.Direct.sampleBoundaryLaurent(c, F, 400), /*unbounded=*/true);
  }

  function displayH(hDisp, hKatex, hData, isUnbounded, cValue) {
    const lines = ['h(w) = '];
    const polyPart = hData.polyPart || [];
    if (polyPart.length > 0) {
      const polyTerms = [];
      for (let l = 0; l < polyPart.length; l++) {
        const c = polyPart[l];
        if (Math.abs(c.re) < 1e-14 && Math.abs(c.im) < 1e-14) continue;
        polyTerms.push('  ' + complexToString(c) + (l === 0 ? '' : ' · w' + (l === 1 ? '' : '^' + l)));
      }
      if (polyTerms.length) lines.push.apply(lines, polyTerms);
    }
    for (const pole of hData.poles) {
      for (let k = 0; k < pole.principal.length; k++) {
        const C = pole.principal[k];
        const denPow = (k === 0) ? '' : '^' + (k + 1);
        lines.push('  ' + complexToString(C) + ' / (w − ' + complexToString(pole.a) + ')' + denPow);
      }
    }
    if (lines.length === 1) lines.push('  0');
    if (hDisp) hDisp.textContent = lines.join('\n');

    if (hKatex && window.katex) {
      try {
        let body = '';
        let first = true;
        for (let l = 0; l < polyPart.length; l++) {
          const c = polyPart[l];
          if (Math.abs(c.re) < 1e-14 && Math.abs(c.im) < 1e-14) continue;
          const cstr = complexToKatex(c);
          body += (first ? '' : ' + ') + (l === 0 ? cstr
                  : (l === 1 ? cstr + '\\,w' : cstr + '\\,w^{' + l + '}'));
          first = false;
        }
        for (const pole of hData.poles) {
          for (let k = 0; k < pole.principal.length; k++) {
            const C = pole.principal[k];
            const cstr = complexToKatex(C);
            const den = (k === 0)
              ? `(w - (${complexToKatex(pole.a)}))`
              : `(w - (${complexToKatex(pole.a)}))^{${k + 1}}`;
            body += (first ? '' : ' + ') + '\\frac{' + cstr + '}{' + den + '}';
            first = false;
          }
        }
        if (first) body = '0';
        // Display mode: full-size fractions, left-aligned via the .rm-sym CSS
        // (was inline/script style, which rendered the fractions tiny).
        window.katex.render('h(w) = ' + body, hKatex, { displayMode: true, throwOnError: false });
      } catch (e) { /* leave text fallback */ }
    }
  }

  function pushBoundaryToPlot(pts, unbounded) {
    const setBdy = window.QD && window.QD.Direct && window.QD.Direct._setPlotBoundary;
    if (setBdy) setBdy(pts, { unbounded: !!unbounded });
  }

    return { recomputeAndRender };
  };
})(typeof window !== 'undefined' ? window : globalThis);
