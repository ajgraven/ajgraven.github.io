// =============================================================================
// solvers/seeds/seeds-lqd.js -- Seed strategy for Family.boundedLQD (B3).
//
// Populates QD.Seeds.boundedLQD = { initialGuess, perturbedInitialGuess,
// diverseInitialGuess }. Depends on QD.Complex and QD.LqdCommon (both loaded
// before this file). solver-lqd.js aliases these locally.
//
// Math: exp-of-QD bootstrap. z_j ≈ ln(a_j/w₀)/R (capped inside 𝔻), A_{j,k}
// from the residue structure D scaled by R^{-k}.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Complex) {
    throw new Error("seeds-lqd.js: solver.js / complex.js must be loaded first");
  }

  const Complex = QD.Complex;

  function initialGuess_LQD(hData, norm) {
    const w0 = norm.w0;
    const n = hData.poles.length;

    let totalLog = 0;
    for (const p of hData.poles) {
      const ratio = Complex.div(p.a, w0);
      const mag = Math.hypot(Math.log(Complex.abs(ratio)), Math.atan2(ratio.im, ratio.re));
      if (mag > totalLog) totalLog = mag;
    }
    const R = Math.max(totalLog, 0.3);              // floor to keep z_j inside 𝔻
    const cap = 0.85;                               // keep |z_j| safely < 1

    const phi = { family: 'boundedLQD', w0: Complex.clone(w0), branches: [] };

    for (let j = 0; j < n; j++) {
      const p = hData.poles[j];
      const logRatio = {
        re: Math.log(Complex.abs(Complex.div(p.a, w0))),
        im: Math.atan2(p.a.im * w0.re - p.a.re * w0.im,
                       p.a.re * w0.re + p.a.im * w0.im),
      };
      let z = Complex.scale(logRatio, 1 / R);
      const zr = Complex.abs(z);
      if (zr > cap) z = Complex.scale(z, cap / zr);

      const A = [];
      const D = (() => {
        const out = new Array(p.principal.length);
        for (let s = 0; s < p.principal.length; s++) {
          const aC = Complex.mul(p.a, p.principal[s]);
          const next = (s + 1 < p.principal.length) ? p.principal[s + 1] : { re: 0, im: 0 };
          out[s] = Complex.add(aC, next);
        }
        return out;
      })();
      let Rk = 1;
      for (let k = 1; k <= p.principal.length; k++) {
        Rk *= R;
        A.push(Complex.scale(D[k - 1], 1 / Rk));
      }
      phi.branches.push({ z, A });
    }

    return phi;
  }

  function perturbedInitialGuess_LQD(hData, norm, rng, r) {
    const base = initialGuess_LQD(hData, norm);
    QD.LqdCommon.perturbBranchesInPlace(base.branches, rng, r || 0);
    return base;
  }

  function diverseInitialGuess_LQD(hData, norm, rng) {
    return {
      family: 'boundedLQD',
      w0: Complex.clone(norm.w0),
      branches: QD.LqdCommon.diverseSeedBranches(hData, rng),
    };
  }

  QD.Seeds = QD.Seeds || {};
  QD.Seeds.boundedLQD = {
    initialGuess:          initialGuess_LQD,
    perturbedInitialGuess: perturbedInitialGuess_LQD,
    diverseInitialGuess:   diverseInitialGuess_LQD,
  };
})();
