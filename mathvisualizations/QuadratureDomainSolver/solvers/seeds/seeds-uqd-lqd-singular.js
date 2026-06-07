// =============================================================================
// solvers/seeds/seeds-uqd-lqd-singular.js -- Seed strategy for
// Family.unboundedLQD_singular (B3).
//
// Populates QD.Seeds.unboundedLQD_singular = { initialGuess,
// perturbedInitialGuess, diverseInitialGuess }. The seeds call several helpers
// that solver-uqd-lqd-singular.js exports for this purpose (resolved at solve
// time): QD._finitePolesView_UQDLS, QD._seedLqdGamma_UQDLS,
// QD.computeTargetF_UQDLS, plus QD.LqdCommon.* and QD.solveInverseQD.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Complex) {
    throw new Error("seeds-uqd-lqd-singular.js: solver.js / complex.js must be loaded first");
  }

  const Complex = QD.Complex;

  function initialGuess_UQDLS(hData, norm) {
    const c = norm.c;
    const q = norm.q;

    // Branches correspond 1-to-1 with FINITE poles only.
    const finiteHData = QD._finitePolesView_UQDLS(hData);
    const finitePoles = finiteHData.poles;

    let zj_guess = null, A_guess = null, z0_guess = null;

    // Try companion bootstrap (only if at least one finite pole).
    if (finitePoles.length > 0) {
      try {
        const companion = QD.solveInverseQD(finiteHData, {
          lqd: true, unbounded: true, c,
          identityTol: 1e-3, autoEscalate: false, findAlternates: false,
        });
        if (companion.success && companion.primary && companion.primary.phi) {
          const phiUQDL = companion.primary.phi;
          zj_guess = phiUQDL.branches.map(br => Complex.clone(br.z));
          A_guess  = phiUQDL.branches.map(br => br.A.map(Complex.clone));

          // z₀ = argmin |φ_UQDL(z)| on |z| = 1.01, pushed slightly outward.
          const ring = 1.01;
          let bestZ = null, bestMag = Infinity;
          for (let i = 0; i < 60; i++) {
            const theta = 2 * Math.PI * i / 60;
            const z = { re: ring * Math.cos(theta), im: ring * Math.sin(theta) };
            const mag = Complex.abs(QD.Family.unboundedLQD.evalPhi(z, phiUQDL));
            if (mag < bestMag) { bestMag = mag; bestZ = z; }
          }
          if (bestZ) z0_guess = Complex.scale(bestZ, 1.05);
        }
      } catch (e) { /* fall through to geometric */ }
    }

    // Geometric fallback (also used when no finite poles).
    if (!zj_guess) {
      zj_guess = finitePoles.map(p => {
        let z = Complex.scale(p.a, 1 / c);
        const r = Complex.abs(z);
        if (r < 1.05) z = Complex.scale(z, 1.05 / Math.max(r, 1e-15));
        return z;
      });
      A_guess = finitePoles.map(p => {
        const D = [];
        for (let s = 0; s < p.principal.length; s++) {
          const aC = Complex.mul(p.a, p.principal[s]);
          const next = (s + 1 < p.principal.length) ? p.principal[s + 1] : { re: 0, im: 0 };
          D.push(Complex.add(aC, next));
        }
        let ck = 1;
        const A = [];
        for (let k = 1; k <= p.principal.length; k++) {
          ck *= c;
          A.push(Complex.scale(D[k - 1], 1 / ck));
        }
        return A;
      });
    }
    if (!z0_guess) z0_guess = { re: 2, im: 0 };

    // Seed lqdBeta from polyPart and lqdGamma from a=0 principal.
    const polyPart = hData.polyPart || [];
    const phiInit = {
      family: 'unboundedLQD_singular',
      unbounded: true,
      c, q: Complex.clone(q),
      z0: z0_guess,
      w0: undefined,
      branches: zj_guess.map((z, j) => ({ z, A: A_guess[j].map(Complex.clone) })),
      lqdBeta:  polyPart.map(() => ({ re: 0, im: 0 })),
      lqdGamma: QD._seedLqdGamma_UQDLS(hData),
    };
    if (polyPart.length > 0) {
      const targetF = QD.computeTargetF_UQDLS(phiInit, hData);
      phiInit.lqdBeta = targetF.map(c => ({ re: c.re, im: c.im }));
    }
    return phiInit;
  }

  function perturbedInitialGuess_UQDLS(hData, norm, rng, r) {
    const base = initialGuess_UQDLS(hData, norm);
    QD.LqdCommon.perturbBranchesInPlace(base.branches, rng, r || 0,
      { side: 'out', zCap: 1.05, zScale: 1.10 });
    // Perturb z₀ too, with the same out-side clamp + upper bound.
    const sigma = 0.15 + 0.25 * (r || 0);
    base.z0 = {
      re: base.z0.re + sigma * (rng() - 0.5),
      im: base.z0.im + sigma * (rng() - 0.5),
    };
    const rz0 = Math.hypot(base.z0.re, base.z0.im);
    if (rz0 < 1.05)    { const s = 1.05 / Math.max(rz0, 1e-15); base.z0.re *= s; base.z0.im *= s; }
    else if (rz0 > 50) { const s = 50   / rz0;                  base.z0.re *= s; base.z0.im *= s; }
    // Perturb lqdGamma multiplicatively (real-axis) + additively (im).
    for (let l = 0; l < base.lqdGamma.length; l++) {
      base.lqdGamma[l] = {
        re: base.lqdGamma[l].re * (1 + sigma * (rng() - 0.5)),
        im: base.lqdGamma[l].im + sigma * (rng() - 0.5),
      };
    }
    return base;
  }

  function diverseInitialGuess_UQDLS(hData, norm, rng) {
    const c = norm.c, q = norm.q;
    const mz0 = Math.exp(Math.log(1.05) + rng() * Math.log(30 / 1.05));
    const pz0 = 2 * Math.PI * rng();
    const polyPart = hData.polyPart || [];
    const finiteHData = QD._finitePolesView_UQDLS(hData);
    const base = {
      family: 'unboundedLQD_singular',
      unbounded: true,
      c, q: Complex.clone(q),
      z0: { re: mz0 * Math.cos(pz0), im: mz0 * Math.sin(pz0) },
      w0: undefined,
      branches: QD.LqdCommon.diverseSeedBranches(finiteHData, rng, { zMin: 1.05, zMax: 30 }),
      lqdBeta:  polyPart.map(() => ({ re: 0, im: 0 })),
      lqdGamma: QD._seedLqdGamma_UQDLS(hData),
    };
    if (polyPart.length > 0) {
      base.lqdBeta = QD.computeTargetF_UQDLS(base, hData).map(c => ({ re: c.re, im: c.im }));
    }
    return base;
  }

  QD.Seeds = QD.Seeds || {};
  QD.Seeds.unboundedLQD_singular = {
    initialGuess:          initialGuess_UQDLS,
    perturbedInitialGuess: perturbedInitialGuess_UQDLS,
    diverseInitialGuess:   diverseInitialGuess_UQDLS,
  };
})();
