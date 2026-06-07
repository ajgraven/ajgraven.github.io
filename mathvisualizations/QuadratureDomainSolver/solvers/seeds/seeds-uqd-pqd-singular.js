// =============================================================================
// solvers/seeds/seeds-uqd-pqd-singular.js -- Seed strategy for
// Family.unboundedPQD_singular (B3).
//
// Populates QD.Seeds.unboundedPQD_singular = { initialGuess,
// perturbedInitialGuess, diverseInitialGuess }. Self-contained (Complex only).
// solver-uqd-pqd-singular.js aliases these locally.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Complex) {
    throw new Error("seeds-uqd-pqd-singular.js: solver.js / complex.js must be loaded first");
  }

  const Complex = QD.Complex;

  function initialGuess_UPQDS(hData, norm) {
    const c = norm.c, alpha = norm.alpha;
    const phi = {
      family: 'unboundedPQD_singular', unbounded: true, alpha, c,
      w0: undefined, z0: { re: 2, im: 0 }, polyA: [], branches: [],
    };
    let minA = Infinity;
    for (const p of hData.poles) { const m = Complex.abs(p.a); if (m < minA) minA = m; }
    const effC = isFinite(minA) && minA > 0 ? Math.min(c, 0.5 * minA) : c;
    for (const p of hData.poles) {
      let z = Complex.abs2(p.a) < 1e-30 ? { re: 2, im: 0 } : Complex.scale(p.a, 1 / effC);
      const rr = Math.hypot(z.re, z.im);
      if (rr < 1.05) { const s = 1.05 / Math.max(rr, 1e-9); z = { re: z.re * s, im: z.im * s }; }
      const A = [];
      let cPow = 1;
      for (let k = 1; k <= p.principal.length; k++) { cPow *= effC; A.push(Complex.scale(p.principal[k - 1], 1 / cPow)); }
      phi.branches.push({ z, A });
    }
    if (hData.polyPart && hData.polyPart.length > 0) {
      const N = hData.polyPart.length, n = N - 1;
      for (let l = 0; l < N; l++) phi.polyA.push({ re: 0, im: 0 });
      const scale = alpha * Math.pow(c, (n + 1) - alpha);
      phi.polyA[n] = Complex.scale(Complex.conj(hData.polyPart[n]), scale);
    }
    return phi;
  }

  function perturbedInitialGuess_UPQDS(hData, norm, rng, r = 0) {
    const base = initialGuess_UPQDS(hData, norm);
    const sigma = 0.15 + 0.25 * r;
    for (const br of base.branches) {
      br.z = { re: br.z.re + sigma * (rng() - 0.5), im: br.z.im + sigma * (rng() - 0.5) };
      const rr = Math.hypot(br.z.re, br.z.im);
      if (rr < 1.05) { const s = 1.05 / Math.max(rr, 1e-9); br.z.re *= s; br.z.im *= s; }
      for (let k = 0; k < br.A.length; k++) br.A[k] = { re: br.A[k].re * (1 + sigma * (rng() - 0.5)), im: br.A[k].im + sigma * (rng() - 0.5) };
    }
    base.z0 = { re: base.z0.re + sigma * (rng() - 0.5), im: base.z0.im + sigma * (rng() - 0.5) };
    const rz = Math.hypot(base.z0.re, base.z0.im);
    if (rz < 1.05) { const s = 1.05 / Math.max(rz, 1e-9); base.z0.re *= s; base.z0.im *= s; }
    for (let l = 0; l < base.polyA.length; l++) base.polyA[l] = { re: base.polyA[l].re * (1 + sigma * (rng() - 0.5)), im: base.polyA[l].im + sigma * (rng() - 0.5) };
    return base;
  }

  function diverseInitialGuess_UPQDS(hData, norm, rng, r) {
    const base = initialGuess_UPQDS(hData, norm);
    const mag = 1.2 + 3 * ((r % 5) / 4);
    const ang = 2 * Math.PI * rng();
    base.z0 = { re: mag * Math.cos(ang), im: mag * Math.sin(ang) };
    for (const br of base.branches) {
      br.z = { re: br.z.re + 0.5 * (rng() - 0.5), im: br.z.im + 0.5 * (rng() - 0.5) };
      const rr = Math.hypot(br.z.re, br.z.im);
      if (rr < 1.1) { const s = 1.1 / Math.max(rr, 1e-9); br.z.re *= s; br.z.im *= s; }
    }
    return base;
  }

  QD.Seeds = QD.Seeds || {};
  QD.Seeds.unboundedPQD_singular = {
    initialGuess:          initialGuess_UPQDS,
    perturbedInitialGuess: perturbedInitialGuess_UPQDS,
    diverseInitialGuess:   diverseInitialGuess_UPQDS,
  };
})();
