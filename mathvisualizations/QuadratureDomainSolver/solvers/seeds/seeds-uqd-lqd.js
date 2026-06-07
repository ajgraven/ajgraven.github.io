// =============================================================================
// solvers/seeds/seeds-uqd-lqd.js -- Seed strategy for Family.unboundedLQD (B3).
//
// Populates QD.Seeds.unboundedLQD = { initialGuess, perturbedInitialGuess,
// diverseInitialGuess }. The β-seed step calls QD.computeTargetF_UQDL (exported
// by solver-uqd-lqd.js, resolved at solve time) and QD.LqdCommon helpers.
// solver-uqd-lqd.js aliases these locally.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Complex) {
    throw new Error("seeds-uqd-lqd.js: solver.js / complex.js must be loaded first");
  }

  const Complex = QD.Complex;

  function initialGuess_UQDL(hData, norm) {
    const c = norm.c;
    // Choose an effective c for the initial guess such that all |z_j| ≥ 1.05;
    // continuation walks it to the target.
    let minA = Infinity;
    for (const p of hData.poles) {
      const m = Complex.abs(p.a);
      if (m > 0 && m < minA) minA = m;
    }
    const cap = isFinite(minA) && minA > 0 ? 0.5 * minA : Math.min(1, c);
    const effC = Math.min(c, cap);

    const branches = hData.poles.map(p => {
      let z;
      if (Complex.abs2(p.a) < 1e-30) {
        // Non-singular guarantees 0 ∉ Ω̄, so a_j ≠ 0; this is a safety net.
        z = { re: 2, im: 0 };
      } else {
        z = Complex.scale(p.a, 1 / effC);
        const r = Complex.abs(z);
        if (r < 1.05) z = Complex.scale(z, 1.05 / Math.max(r, 1e-15));
      }
      const A = [];
      let cPow = 1;
      for (let k = 1; k <= p.principal.length; k++) {
        cPow *= effC;
        const Cjk     = p.principal[k - 1];
        const Cjknext = (k < p.principal.length) ? p.principal[k] : { re: 0, im: 0 };
        const Djk     = Complex.add(Complex.mul(p.a, Cjk), Cjknext);
        A.push(Complex.scale(Djk, 1 / cPow));
      }
      return { z, A };
    });

    // Seed lqdBeta from polyPart by evaluating computeTargetF at this initial
    // φ with β = [0, ..., 0].
    const polyPart = hData.polyPart || [];
    const phiInit = {
      family: 'unboundedLQD',
      unbounded: true,
      c, w0: undefined,
      branches,
      lqdBeta: polyPart.map(() => ({ re: 0, im: 0 })),
    };
    if (polyPart.length > 0) {
      const targetF = QD.computeTargetF_UQDL(phiInit, hData);
      phiInit.lqdBeta = targetF.map(c => ({ re: c.re, im: c.im }));
    }
    return phiInit;
  }

  function perturbedInitialGuess_UQDL(hData, norm, rng, r) {
    const base = initialGuess_UQDL(hData, norm);
    QD.LqdCommon.perturbBranchesInPlace(base.branches, rng, r || 0,
      { side: 'out', zCap: 1.05, zScale: 1.10 });
    return base;
  }

  function diverseInitialGuess_UQDL(hData, norm, rng) {
    const polyPart = hData.polyPart || [];
    const base = {
      family: 'unboundedLQD',
      unbounded: true,
      c: norm.c, w0: undefined,
      branches: QD.LqdCommon.diverseSeedBranches(hData, rng, { zMin: 1.05, zMax: 30 }),
      lqdBeta: polyPart.map(() => ({ re: 0, im: 0 })),
    };
    if (polyPart.length > 0) {
      base.lqdBeta = QD.computeTargetF_UQDL(base, hData).map(c => ({ re: c.re, im: c.im }));
    }
    return base;
  }

  QD.Seeds = QD.Seeds || {};
  QD.Seeds.unboundedLQD = {
    initialGuess:          initialGuess_UQDL,
    perturbedInitialGuess: perturbedInitialGuess_UQDL,
    diverseInitialGuess:   diverseInitialGuess_UQDL,
  };
})();
