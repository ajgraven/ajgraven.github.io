// =============================================================================
// solvers/seeds/seeds-pqd-singular.js -- Seed strategy for
// Family.powerQD_singular (B3).
//
// Populates QD.Seeds.powerQD_singular = { initialGuess, perturbedInitialGuess,
// diverseInitialGuess }. Reuses the non-singular PQD disk seed
// (QD.diskInitialGuess_PQD, exported by solver-pqd.js, resolved at solve time)
// for z_j / A, then attaches z₀. solver-pqd-singular.js aliases these locally.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Complex) {
    throw new Error("seeds-pqd-singular.js: solver.js / complex.js must be loaded first");
  }

  const Complex = QD.Complex;

  function initialGuess_PQDS(hData, norm) {
    const w0 = norm.w0;
    const alpha = norm.alpha;
    const n = hData.poles.length;
    // Reuse the non-singular PQD disk seed for z_j / A, then attach z₀.
    // With the prefactor Blaschke (b(0)=|z₀|), z₀ is unconstrained in sign;
    // seed it real positive (the canonical example has z₀ = 2/3).
    const base = QD.diskInitialGuess_PQD(hData, w0, alpha);
    return {
      family: 'powerQD_singular',
      alpha,
      unbounded: false,
      w0: Complex.clone(w0),
      z0: { re: 0.5, im: 0 },
      branches: base.branches.map(br => ({
        z: Complex.clone(br.z),
        A: br.A.map(Complex.clone),
      })),
      _n: n,
    };
  }

  function perturbedInitialGuess_PQDS(hData, norm, rng, r) {
    const base = initialGuess_PQDS(hData, norm);
    const sigma = 0.15 + 0.25 * (r || 0);
    for (const br of base.branches) {
      br.z = { re: br.z.re + sigma * (rng() - 0.5), im: br.z.im + sigma * (rng() - 0.5) };
      const rr = Math.hypot(br.z.re, br.z.im);
      if (rr > 0.9) { br.z.re *= 0.85 / rr; br.z.im *= 0.85 / rr; }
      for (let k = 0; k < br.A.length; k++) {
        br.A[k] = { re: br.A[k].re * (1 + sigma * (rng() - 0.5)), im: br.A[k].im + sigma * (rng() - 0.5) };
      }
    }
    // Jitter z₀ inside the disk (kept near the positive real axis).
    const dz = sigma * (rng() - 0.5);
    base.z0 = { re: Math.max(0.05, Math.min(0.95, 0.5 + dz)), im: 0.2 * sigma * (rng() - 0.5) };
    return base;
  }

  function diverseInitialGuess_PQDS(hData, norm, rng, r) {
    const base = initialGuess_PQDS(hData, norm);
    // Sweep |z₀| across the (positive) disk for diversity.
    const mag = 0.2 + 0.7 * ((r % 5) / 4);
    base.z0 = { re: mag, im: 0 };
    for (const br of base.branches) {
      br.z = { re: br.z.re + 0.4 * (rng() - 0.5), im: br.z.im + 0.4 * (rng() - 0.5) };
      const rr = Math.hypot(br.z.re, br.z.im);
      if (rr > 0.9) { br.z.re *= 0.85 / rr; br.z.im *= 0.85 / rr; }
    }
    return base;
  }

  QD.Seeds = QD.Seeds || {};
  QD.Seeds.powerQD_singular = {
    initialGuess:          initialGuess_PQDS,
    perturbedInitialGuess: perturbedInitialGuess_PQDS,
    diverseInitialGuess:   diverseInitialGuess_PQDS,
  };
})();
