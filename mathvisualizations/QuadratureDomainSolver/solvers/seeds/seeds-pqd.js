// =============================================================================
// solvers/seeds/seeds-pqd.js -- Multistart / continuation seed strategy for
// Family.powerQD (bounded power-weighted quadrature domains, α ≥ 2).
//
// Math reference: Theorem 4.3.2 / Corollary 4.3.1 of Graven 2026. The
// parametrization is
//
//   φ(z) = (R#(z))^(1/α)        (Equation 4.8, bounded case with 0 ∈ Ω, φ_in ≡ 1)
//
// where R# matches the existing boundedQD rational form. The Newton
// variables {z_j, A_{j,k}} are IDENTICAL to boundedQD; the α-power lives
// in the evaluator + residual. So a natural seed strategy is:
//   * Start from boundedQD's disk-seed but with the user-supplied φ(0)
//     raised to α: r0 = w0^α (so R(0) = r0).
//   * Pole-distance-driven R is computed against |a_j^α − r0|
//     rather than |a_j − w0|, so the disk radius scales with the
//     effective "size" of the α-lifted problem.
//   * Initial A's match boundedQD's scaling but with a_j → a_j^α / α
//     (since the LHS of the locator after lifting is a_j^α).
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Complex) {
    throw new Error("seeds-pqd.js: solver.js / complex.js must be loaded first");
  }

  const Complex = QD.Complex;

  // Complex real-power c^α (arbitrary α > 0, principal branch). Delegates to
  // the shared Complex.cpow primitive so seeds work for non-integer α too.
  function cpow(c, alpha) {
    return Complex.cpow(c, alpha);
  }

  // Disk seed: place every z_j inside the unit disk and scale A by R^{-k}.
  // For α-PQDs the natural "lifted" residue is a_j^α at each pole; R is
  // chosen large enough that no candidate z_j escapes |z|<1.
  function diskInitialGuess_PQD(hData, w0, alpha, scale = null) {
    const n = hData.poles.length;
    const r0 = cpow(w0, alpha);

    let totalC = 0;
    for (const p of hData.poles) {
      if (p.principal.length > 0) totalC += Complex.abs(p.principal[0]);
    }
    if (totalC === 0) totalC = 1;
    let R = scale !== null ? scale : Math.sqrt(totalC);
    let maxR = 0;
    for (const p of hData.poles) {
      // Distance from lifted pole-image to r0.
      const aLifted = cpow(p.a, alpha);
      const d = Complex.abs(Complex.sub(aLifted, r0));
      if (d > maxR) maxR = d;
    }
    if (R < 1.5 * maxR) R = 1.5 * maxR;
    if (R === 0) R = 1;

    const phi = {
      family: 'powerQD',
      unbounded: false,
      alpha,
      w0: Complex.clone(w0),
      c: undefined,
      branches: [],
    };
    for (let j = 0; j < n; j++) {
      const aLifted = cpow(hData.poles[j].a, alpha);
      // Direction from r0 to a_j^α, scaled to land inside |z|<1.
      const z = Complex.scale(Complex.sub(aLifted, r0), 1 / R);
      const A = [];
      const mj = hData.poles[j].principal.length;
      let Rk = 1;
      for (let k = 1; k <= mj; k++) {
        Rk *= R;
        // Initial guess for A_{j,k}: the boundedQD-style scaling of
        // C_{j,k} by R^{-k}. The 1/α factor reflects the chain-rule
        // factor when extracting Taylor of (·)^{1/α}.
        A.push(Complex.scale(hData.poles[j].principal[k - 1], 1 / (Rk * alpha)));
      }
      phi.branches.push({ z, A });
    }
    return phi;
  }

  // Perturbed seed for multistart. sigma grows with r so successive
  // restarts cast a wider net. z_j clamped inside |z|<0.85.
  function perturbedInitialGuess_PQD(hData, w0, alpha, rng, r = 0) {
    const base = diskInitialGuess_PQD(hData, w0, alpha);
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
  QD.Seeds.powerQD = {
    initialGuess:          diskInitialGuess_PQD,
    perturbedInitialGuess: perturbedInitialGuess_PQD,
  };
})();
