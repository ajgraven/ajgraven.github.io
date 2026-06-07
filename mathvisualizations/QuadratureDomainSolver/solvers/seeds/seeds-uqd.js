// =============================================================================
// solvers/seeds/seeds-uqd.js -- Multistart / continuation seed strategy for
// Family.unboundedQD (B3, completing the A6 seed-split across all families).
//
// Populates QD.Seeds.unboundedQD = { initialGuess, perturbedInitialGuess }.
// solver-uqd.js dispatches through these (and aliases them locally so its
// continuation-in-c loop keeps calling the same names).
//
// Math: exterior-disk seed. z_j placed at p_j / c (clamped so |z_j| > 1),
// A_{j,k} scaled by c^{-k}, Laurent-at-∞ polyA from h's polyPart scaled by c^l.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Complex) {
    throw new Error("seeds-uqd.js: solver.js / complex.js must be loaded first");
  }

  const Complex = QD.Complex;

  function unboundedInitialGuess_UQD(hData, cUser) {
    let minA = Infinity;
    for (const p of hData.poles) {
      const m = Complex.abs(p.a);
      if (m < minA) minA = m;
    }
    const cap = isFinite(minA) && minA > 0 ? 0.5 * minA : Math.min(1, cUser);
    const effC = Math.min(cUser, cap);
    const phi = {
      unbounded: true, c: cUser, w0: undefined, polyA: [], branches: [],
    };
    for (const p of hData.poles) {
      let z;
      if (Complex.abs2(p.a) < 1e-30) {
        z = { re: 2, im: 0 };
      } else {
        z = Complex.scale(p.a, 1 / effC);
      }
      const A = [];
      let cPow = 1;
      for (let k = 1; k <= p.principal.length; k++) {
        cPow *= effC;
        A.push(Complex.scale(p.principal[k - 1], 1 / cPow));
      }
      phi.branches.push({ z, A });
    }
    if (hData.polyPart && hData.polyPart.length > 0) {
      let cPowL = 1;
      for (let l = 0; l < hData.polyPart.length; l++) {
        if (l > 0) cPowL *= cUser;
        phi.polyA.push(Complex.scale(hData.polyPart[l], cPowL));
      }
    }
    return phi;
  }

  function perturbedUnboundedInitialGuess_UQD(hData, c, rng, r = 0) {
    const base = unboundedInitialGuess_UQD(hData, c);
    const sigma = 0.15 + 0.25 * r;
    for (const br of base.branches) {
      br.z = {
        re: br.z.re + sigma * (rng() - 0.5),
        im: br.z.im + sigma * (rng() - 0.5)
      };
      const rr = Math.hypot(br.z.re, br.z.im);
      if (rr < 1.05) {
        const scale = 1.05 / Math.max(rr, 1e-9);
        br.z.re *= scale; br.z.im *= scale;
      }
      for (let k = 0; k < br.A.length; k++) {
        br.A[k] = {
          re: br.A[k].re * (1 + sigma * (rng() - 0.5)),
          im: br.A[k].im + sigma * (rng() - 0.5)
        };
      }
    }
    for (let l = 0; l < base.polyA.length; l++) {
      base.polyA[l] = {
        re: base.polyA[l].re * (1 + sigma * (rng() - 0.5))
              + (Math.abs(base.polyA[l].re) < 1e-9 ? sigma * (rng() - 0.5) : 0),
        im: base.polyA[l].im + sigma * (rng() - 0.5),
      };
    }
    return base;
  }

  QD.Seeds = QD.Seeds || {};
  QD.Seeds.unboundedQD = {
    initialGuess:          unboundedInitialGuess_UQD,
    perturbedInitialGuess: perturbedUnboundedInitialGuess_UQD,
  };
})();
