// =============================================================================
// direct-verify.js -- Verify button for the Direct tab.
//
// Extracted from direct-ui.js by the Phase-3 UI modularization (item E).
// Exposes QD_UI.installDirectVerify(dCtx); direct-ui.js captures runVerify into
// an IIFE-local binding so the output-card handler calls it unchanged.
//
// runVerify dispatches per family: unbounded/bounded weighted families use the
// trusted family identity verifier (or a round-trip through the inverse solver);
// classical QDs use the Fourier boundary-identity diagnostic. sampleAnalyticPhi
// samples the user's input phi on |z|=1 in the mode-appropriate way. Bodies are
// VERBATIM moves. Deps via dCtx: directState + parseComplex. `QD` / `window` /
// `math` are globals.
// =============================================================================
'use strict';

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  global.QD_UI.installDirectVerify = function installDirectVerify(d) {
    const directState  = d.directState;
    const parseComplex = d.parseComplex;

  function runVerify(card) {
    const resBox = card.querySelector('.dir-verify-result');
    const overlayHook = window.QD && window.QD.Direct && window.QD.Direct._setPlotOverlay;
    if (overlayHook) overlayHook(null);                  // clear any stale overlay

    if (!directState.lastH) {
      resBox.style.color = '#b53030';
      resBox.textContent = 'Compute a valid h first.';
      return;
    }

    // UNBOUNDED weighted families (∞∈Ω, Thm 4.3.7): verify the FULL h with the
    // family identity verifier on the built φ (the trusted oracle). For LQD-singular
    // the verifier reads phi.q internally; PQD/LQD non-singular and PQD-singular
    // carry no origin term, so hData (poles + polyPart) suffices.
    if (directState.mode === 'unbounded' && directState.lastWeight && directState.lastWeight !== 'classical' && directState.lastPhi) {
      const fam = QD.Family[directState.lastPhi.family];
      let id;
      try { id = fam.verifyQuadratureIdentity(directState.lastPhi, directState.lastH, {}).maxRelDiff; }
      catch (e) { resBox.style.color = '#b53030'; resBox.textContent = 'Verify error: ' + (e.message || e); return; }
      // Non-singular families hit machine precision (~1e-15); the singular boundary
      // (Blaschke × weight, 0∈Ω) floors near ~1e-6 at the verifier's sample density,
      // so use a slightly looser strong-pass threshold there.
      const okStrong = id < (directState.lastSingular ? 1e-5 : 1e-6);
      const fname = directState.lastPhi.family;
      resBox.style.color = okStrong ? '#2a8f2a' : (id < 1e-2 ? '#b8860b' : '#b53030');
      resBox.innerHTML = 'Quadrature identity (' + fname + ', family verifier): maxRelDiff = <strong>'
        + id.toExponential(2) + '</strong>' + (okStrong ? ' ✓ — h reproduces the unbounded weighted identity.' : '.');
      return;
    }

    // Weighted families: the classical Fourier check (h(φ)−conj(φ) analytic)
    // does NOT apply (the weighted identity is different). Verify by ROUND-TRIP:
    // feed the forward h back to the inverse solver and confirm it reconstructs a
    // univalent Ω whose quadrature identity closes — reusing the trusted inverse
    // identity verifier.
    if (directState.mode === 'bounded' && directState.weight !== 'classical') {
      // SINGULAR (log) forward: verify directly with the family identity verifier
      // on the built φ. (Round-trip is avoided here because the inverse singular-LQD
      // solver doesn't always converge for an arbitrary prescribed origin-residue q
      // — a solver limitation, not a forward-kernel issue; the family verifier IS the
      // trusted identity oracle.)
      if (directState.lastSingular && directState.lastPhi) {
        // Verify the FULL h (finite poles + the origin term) with the family
        // identity verifier — the trusted oracle. For PQD the origin term r₀/w
        // must be added to hData (the verifier reads only hData.poles); for LQD the
        // verifier reads phi.q internally so the finite hData suffices.
        let fam, hVer;
        if (directState.weight === 'power') {
          fam = QD.Family.powerQD_singular;
          const r0 = directState.lastOriginRes;
          hVer = r0 ? { poles: directState.lastH.poles.concat([{ a: { re: 0, im: 0 }, principal: [r0] }]) } : directState.lastH;
        } else {
          fam = QD.Family.boundedLQD_singular; hVer = directState.lastH;
        }
        let id;
        try { id = fam.verifyQuadratureIdentity(directState.lastPhi, hVer, {}).maxRelDiff; }
        catch (e) { resBox.style.color = '#b53030'; resBox.textContent = 'Verify error: ' + (e.message || e); return; }
        const okStrong = id < 1e-6;
        const term = directState.weight === 'power' ? 'r₀/w' : 'q/w';
        resBox.style.color = okStrong ? '#2a8f2a' : (id < 1e-2 ? '#b8860b' : '#b53030');
        resBox.innerHTML = 'Quadrature identity (singular ' + (directState.weight === 'power' ? 'PQD' : 'LQD')
          + ', family verifier): maxRelDiff = <strong>' + id.toExponential(2) + '</strong>'
          + (okStrong ? ' ✓ — h reproduces the weighted identity (incl. the origin term ' + term + ').' : '.');
        return;
      }
      const opts = (directState.weight === 'power')
        ? { alpha: parseFloat(directState.alpha) }
        : { weight: 'log', w0: parseComplex(directState.logW0) };
      let r;
      try { r = QD.solveInverseQD(directState.lastH, opts); }
      catch (e) { resBox.style.color = '#b53030'; resBox.textContent = 'Round-trip solve error: ' + (e.message || e); return; }
      const id = r && r.success && r.primary && r.primary.identity ? r.primary.identity.maxRelDiff : Infinity;
      const univ = r && r.success && r.primary && r.primary.univalent;
      const okStrong = r.success && univ && id < 1e-6;
      resBox.style.color = okStrong ? '#2a8f2a' : (r.success && id < 1e-2 ? '#b8860b' : '#b53030');
      resBox.innerHTML = okStrong
        ? 'Round-trip ✓ — the inverse solver reconstructs a univalent Ω; quadrature identity closes (maxRelDiff = <strong>' + id.toExponential(2) + '</strong>).'
        : (r.success
            ? 'Round-trip weak: solved' + (univ ? '' : ' but NOT univalent') + ', identity maxRelDiff = <strong>' + id.toExponential(2) + '</strong>.'
            : 'Round-trip failed: ' + (r.error || 'no solution') + '.');
      return;
    }

    // Sample φ at N points on |z|=1.
    const N = 500;
    let phiPts;
    try { phiPts = sampleAnalyticPhi(N); }
    catch (e) {
      resBox.style.color = '#b53030';
      resBox.textContent = 'Could not sample φ: ' + (e.message || e);
      return;
    }

    // The correct identity is: h(φ(z)) − conj(φ(z)) is analytic in 𝔻 (after
    // composing with φ), so its Fourier expansion on |z|=1 has only non-
    // negative-frequency terms. We measure the negative-frequency Fourier
    // mass — this should be ≈ 0 for any valid classical QD.
    const v = QD.Direct.verifyBoundaryIdentity(directState.lastH, phiPts);

    // Relative score: negMass normalised by the boundary-data scale.
    const relNeg = v.scale > 0 ? v.negMass / v.scale : v.negMass;
    let color;
    if      (relNeg < 1e-8) color = '#2a8f2a';
    else if (relNeg < 1e-2) color = '#b8860b';
    else                    color = '#b53030';

    resBox.style.color = color;
    resBox.innerHTML =
      'Fourier diagnostic (' + N + ' samples on |z|=1):<br>' +
      '&nbsp;&nbsp;negative-freq mass = <strong>' + v.negMass.toExponential(2) +
      '</strong> (relative ' + relNeg.toExponential(2) + ')<br>' +
      '&nbsp;&nbsp;<span style="color:#888">zero-mode mass = ' + v.zeroMass.toExponential(2) +
      ', positive-freq mass = ' + v.posMass.toExponential(2) + '</span>';
  }

  function sampleAnalyticPhi(N) {
    if (directState.mode === 'bounded') {
      if (directState.coeffsKind === 'rational') {
        const P = directState.coeffsNum.map(parseComplex);
        const Q = directState.coeffsDen.map(parseComplex);
        const pts = new Array(N);
        for (let n = 0; n < N; n++) {
          const t = 2 * Math.PI * n / N;
          const z = { re: Math.cos(t), im: Math.sin(t) };
          const pv = QD.Direct.evalPolyAscending(P, z);
          const qv = QD.Direct.evalPolyAscending(Q, z);
          const d2 = qv.re*qv.re + qv.im*qv.im;
          pts[n] = { re: (pv.re*qv.re + pv.im*qv.im) / d2,
                     im: (pv.im*qv.re - pv.re*qv.im) / d2 };
        }
        return pts;
      }
      const cs = directState.coeffs.map(parseComplex);
      return QD.Direct.sampleBoundaryPolynomial(cs, N);
    } else if (directState.mode === 'unbounded') {
      const c = Number(directState.cValue);
      const F = directState.Fcoeffs.map(parseComplex);
      return QD.Direct.sampleBoundaryLaurent(c, F, N);
    } else {
      // Numerical: re-evaluate the user's expression.
      if (!window.math) throw new Error('math.js not loaded');
      const compiled = window.math.parse(directState.numExpr).compile();
      const pts = new Array(N);
      for (let n = 0; n < N; n++) {
        const theta = 2 * Math.PI * n / N;
        const v = compiled.evaluate({ z: window.math.complex(Math.cos(theta), Math.sin(theta)) });
        if (typeof v === 'number') pts[n] = { re: v, im: 0 };
        else pts[n] = { re: v.re, im: v.im };
      }
      return pts;
    }
  }

    return { runVerify };
  };
})(typeof window !== 'undefined' ? window : globalThis);
