// =============================================================================
// solvers/seeds/seeds-qd.js -- Multistart / continuation seed strategy for
// Family.boundedQD.
//
// A6 (post-P3) extracts the seed-generation logic out of the math kernel
// (solver-qd.js) so the multistart strategy is tweakable per family
// without touching evalPhi / phiTaylorAt / residual / pack/unpack.
//
// Convention (and the template for the other 5 families):
//   * One seeds-<family>.js per family.
//   * Each populates QD.Seeds.<familyTag> = { initialGuess, perturbedInitialGuess,
//     diverseInitialGuess? }.
//   * solver-<family>.js's Family entries dispatch through QD.Seeds.<tag>.* .
//
// Math used (Theorems 3.2.1 / 3.2.2):
//   diskInitialGuess: seeds with z_j on a circle of radius 1/R times
//     (a_j - w_0); A_{j,k} scaled by R^{-k}. R is chosen large enough that
//     every (a_j - w_0)/R sits inside the unit disk.
//   perturbedInitialGuess: jitter the disk seed by a magnitude that grows
//     with the restart counter r (sigma = 0.15 + 0.25·r).
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Complex) {
    throw new Error("seeds-qd.js: solver.js / complex.js must be loaded first");
  }

  const Complex = QD.Complex;

  // Disk seed: place every z_j inside the unit disk and scale A by R^{-k}.
  // R is chosen as max(sqrt(totalC), 1.5·maxPoleDistance, 1) — big enough
  // that no candidate z_j escapes |z|<1.
  function diskInitialGuess_QD(hData, w0, scale = null) {
    const n = hData.poles.length;

    let totalC = 0;
    for (const p of hData.poles) {
      if (p.principal.length > 0) totalC += Complex.abs(p.principal[0]);
    }
    if (totalC === 0) totalC = 1;
    let R = scale !== null ? scale : Math.sqrt(totalC);
    let maxR = 0;
    for (const p of hData.poles) {
      const d = Complex.abs(Complex.sub(p.a, w0));
      if (d > maxR) maxR = d;
    }
    if (R < 1.5 * maxR) R = 1.5 * maxR;
    if (R === 0) R = 1;

    const phi = { unbounded: false, w0: Complex.clone(w0), c: undefined, branches: [] };
    for (let j = 0; j < n; j++) {
      const z = Complex.scale(Complex.sub(hData.poles[j].a, w0), 1 / R);
      const A = [];
      const mj = hData.poles[j].principal.length;
      let Rk = 1;
      for (let k = 1; k <= mj; k++) {
        Rk *= R;
        A.push(Complex.scale(hData.poles[j].principal[k - 1], 1 / Rk));
      }
      phi.branches.push({ z, A });
    }
    return phi;
  }

  // Perturbed seed for multistart. sigma grows with r so successive
  // restarts cast a wider net. z_j clamped inside |z|<0.9 to avoid the
  // boundary clamp inside newtonSolve.
  function perturbedInitialGuess_QD(hData, w0, rng, r = 0) {
    const base = diskInitialGuess_QD(hData, w0);
    const sigma = 0.15 + 0.25 * r;
    for (const br of base.branches) {
      br.z = {
        re: br.z.re + sigma * (rng() - 0.5),
        im: br.z.im + sigma * (rng() - 0.5),
      };
      const rr = Math.hypot(br.z.re, br.z.im);
      if (rr > 0.9) { br.z.re *= 0.85 / rr; br.z.im *= 0.85 / rr; }
      for (let k = 0; k < br.A.length; k++) {
        br.A[k] = {
          re: br.A[k].re * (1 + sigma * (rng() - 0.5)),
          im: br.A[k].im + sigma * (rng() - 0.5),
        };
      }
    }
    return base;
  }

  QD.Seeds = QD.Seeds || {};
  QD.Seeds.boundedQD = {
    initialGuess:          diskInitialGuess_QD,
    perturbedInitialGuess: perturbedInitialGuess_QD,
  };
})();
