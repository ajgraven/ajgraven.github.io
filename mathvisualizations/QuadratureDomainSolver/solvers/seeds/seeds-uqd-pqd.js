// =============================================================================
// solvers/seeds/seeds-uqd-pqd.js -- Seed strategy for Family.unboundedPQD (B3).
//
// Populates QD.Seeds.unboundedPQD = { initialGuess, perturbedInitialGuess,
// diverseInitialGuess }. Self-contained (Complex only); the G_l (Laurent-at-∞)
// seed uses the monomial closed form. solver-uqd-pqd.js aliases these locally.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Complex) {
    throw new Error("seeds-uqd-pqd.js: solver.js / complex.js must be loaded first");
  }

  const Complex = QD.Complex;

  function initialGuess_UPQD(hData, norm) {
    const c = norm.c, alpha = norm.alpha;
    let minA = Infinity;
    for (const p of hData.poles) {
      const m = Complex.abs(p.a);
      if (m < minA) minA = m;
    }
    const cap = isFinite(minA) && minA > 0 ? 0.5 * minA : Math.min(1, c);
    const effC = Math.min(c, cap);
    const phi = {
      family: 'unboundedPQD', unbounded: true, alpha, c, w0: undefined,
      polyA: [], branches: [],
    };
    for (const p of hData.poles) {
      let z;
      if (Complex.abs2(p.a) < 1e-30) z = { re: 2, im: 0 };
      else z = Complex.scale(p.a, 1 / effC);
      // Keep z_j strictly exterior.
      const rr = Math.hypot(z.re, z.im);
      if (rr < 1.05) { const s = 1.05 / Math.max(rr, 1e-9); z = { re: z.re * s, im: z.im * s }; }
      const A = [];
      let cPow = 1;
      for (let k = 1; k <= p.principal.length; k++) {
        cPow *= effC;
        A.push(Complex.scale(p.principal[k - 1], 1 / cPow));
      }
      phi.branches.push({ z, A });
    }
    // G_l seed. The top coefficient G_{n+1} ≈ α·h_n·c^{(n+1)−α} comes from the
    // monomial closed form (Thm 4.5.3: r# = c^α(1 − γ_k/z^k), G_k = α·h_n·c^{k−α},
    // k = n+1); lower G_l seed to 0. Generalizes the constant-h closed form (n=0).
    if (hData.polyPart && hData.polyPart.length > 0) {
      const N = hData.polyPart.length, n = N - 1;
      for (let l = 0; l < N; l++) phi.polyA.push({ re: 0, im: 0 });
      const scale = alpha * Math.pow(c, (n + 1) - alpha);
      // G_{n+1} ≈ α·conj(h_n)·c^{(n+1)−α} (the reflection form, confirmed for
      // complex coefficients by the identity verifier).
      phi.polyA[n] = Complex.scale(Complex.conj(hData.polyPart[n]), scale);
    }
    return phi;
  }

  function perturbedInitialGuess_UPQD(hData, norm, rng, r = 0) {
    const base = initialGuess_UPQD(hData, norm);
    const sigma = 0.15 + 0.25 * r;
    for (const br of base.branches) {
      br.z = { re: br.z.re + sigma * (rng() - 0.5), im: br.z.im + sigma * (rng() - 0.5) };
      const rr = Math.hypot(br.z.re, br.z.im);
      if (rr < 1.05) { const s = 1.05 / Math.max(rr, 1e-9); br.z.re *= s; br.z.im *= s; }
      for (let k = 0; k < br.A.length; k++) {
        br.A[k] = { re: br.A[k].re * (1 + sigma * (rng() - 0.5)), im: br.A[k].im + sigma * (rng() - 0.5) };
      }
    }
    for (let l = 0; l < base.polyA.length; l++) {
      base.polyA[l] = { re: base.polyA[l].re * (1 + sigma * (rng() - 0.5)), im: base.polyA[l].im + sigma * (rng() - 0.5) };
    }
    return base;
  }

  function diverseInitialGuess_UPQD(hData, norm, rng) {
    const base = initialGuess_UPQD(hData, norm);
    for (const br of base.branches) {
      br.z = { re: br.z.re + 0.4 * (rng() - 0.5), im: br.z.im + 0.4 * (rng() - 0.5) };
      const rr = Math.hypot(br.z.re, br.z.im);
      if (rr < 1.1) { const s = 1.1 / Math.max(rr, 1e-9); br.z.re *= s; br.z.im *= s; }
    }
    return base;
  }

  QD.Seeds = QD.Seeds || {};
  QD.Seeds.unboundedPQD = {
    initialGuess:          initialGuess_UPQD,
    perturbedInitialGuess: perturbedInitialGuess_UPQD,
    diverseInitialGuess:   diverseInitialGuess_UPQD,
  };
})();
