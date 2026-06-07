// =============================================================================
// solvers/seeds/seeds-lqd-singular.js -- Seed strategy for
// Family.boundedLQD_singular (B3).
//
// Populates QD.Seeds.boundedLQD_singular = { initialGuess,
// perturbedInitialGuess, diverseInitialGuess }. initialGuess bootstraps from a
// companion classical-QD solve (QD.solveInverseQD / QD.evalPhi, resolved at
// call time) then falls back to a geometric disk guess. solver-lqd-singular.js
// aliases these locally.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Complex) {
    throw new Error("seeds-lqd-singular.js: solver.js / complex.js must be loaded first");
  }

  const Complex = QD.Complex;

  function initialGuess_LQDS(hData, norm) {
    const w0 = norm.w0;
    const q = norm.q;
    const n = hData.poles.length;

    // Try companion-QD bootstrap (only if we have at least one finite pole).
    let z0_guess = null;
    let zj_guess = null;
    let A_guess = null;

    if (n > 0) {
      try {
        const companionResult = QD.solveInverseQD(hData, { w0, identityTol: 1e-3, autoEscalate: false });
        if (companionResult.success && companionResult.primary && companionResult.primary.phi) {
          const phiQD = companionResult.primary.phi;
          // Newton-solve φ_QD(z) = 0 starting from z = 0.
          let z = { re: 0, im: 0 };
          for (let it = 0; it < 30; it++) {
            const fz = QD.evalPhi(z, phiQD);
            const fzAbs = Complex.abs(fz);
            if (fzAbs < 1e-10) break;
            const fzPlus = QD.evalPhi({ re: z.re + 1e-7, im: z.im }, phiQD);
            const fzPlusI = QD.evalPhi({ re: z.re, im: z.im + 1e-7 }, phiQD);
            // dF/dz_re and dF/dz_im (treating fz as ℂ ≅ ℝ²)
            const dFdx = { re: (fzPlus.re - fz.re) / 1e-7, im: (fzPlus.im - fz.im) / 1e-7 };
            const dFdy = { re: (fzPlusI.re - fz.re) / 1e-7, im: (fzPlusI.im - fz.im) / 1e-7 };
            // Solve 2x2: [dFdx; dFdy] · [dx, dy] = -fz
            const a = dFdx.re, b = dFdy.re, c = dFdx.im, d = dFdy.im;
            const det = a * d - b * c;
            if (Math.abs(det) < 1e-14) break;
            const rhs1 = -fz.re, rhs2 = -fz.im;
            const dx = ( d * rhs1 - b * rhs2) / det;
            const dy = (-c * rhs1 + a * rhs2) / det;
            z = { re: z.re + dx, im: z.im + dy };
            const r = Math.hypot(z.re, z.im);
            if (r > 0.95) { z.re *= 0.9 / r; z.im *= 0.9 / r; }
          }
          if (Complex.abs(z) > 1e-3 && Complex.abs(z) < 0.99) {
            z0_guess = z;
            zj_guess = phiQD.branches.map(br => Complex.clone(br.z));
            A_guess = phiQD.branches.map(br => br.A.map(Complex.clone));
          }
        }
      } catch (e) {
        // bootstrap failed; fall through to geometric guess
      }
    }

    // Geometric fallback: approximate Ω as a disk centered at w_0 with radius
    // R large enough to contain BOTH the origin and the finite poles. Then
    // φ(z) ≈ w_0 + R·z gives z_0 ≈ -w_0/R (preimage of 0) and z_j ≈ (a_j-w_0)/R.
    if (!z0_guess) {
      // R = 1.2 · max(|w_0|, max_j |a_j - w_0|)  (some headroom)
      let R = Complex.abs(w0);
      for (const p of hData.poles) {
        R = Math.max(R, Complex.abs(Complex.sub(p.a, w0)));
      }
      R = Math.max(R * 1.2, 0.3);

      z0_guess = Complex.scale(w0, -1 / R);
      // Clamp into 𝔻 (with a margin) and off the origin.
      const r0 = Complex.abs(z0_guess);
      if (r0 > 0.9) z0_guess = Complex.scale(z0_guess, 0.85 / r0);
      if (r0 < 1e-2) z0_guess = { re: 0.1, im: 0 };

      zj_guess = hData.poles.map(p => {
        const dir = Complex.sub(p.a, w0);
        const r = Complex.abs(dir);
        if (r < 1e-6) return { re: 0.5, im: 0 };
        return Complex.scale(dir, Math.min(0.85, r / R));
      });

      // A_{j,k} ≈ D_{j,k} / R^k
      A_guess = hData.poles.map(p => {
        const D = [];
        for (let s = 0; s < p.principal.length; s++) {
          const aC = Complex.mul(p.a, p.principal[s]);
          const next = (s + 1 < p.principal.length) ? p.principal[s + 1] : { re: 0, im: 0 };
          D.push(Complex.add(aC, next));
        }
        let Rk = 1;
        const A = [];
        for (let k = 1; k <= p.principal.length; k++) {
          Rk *= R;
          A.push(Complex.scale(D[k - 1], 1 / Rk));
        }
        return A;
      });
    }

    // γ ← w_0 / |z_0|  (from φ_0 equation)
    const absZ0 = Complex.abs(z0_guess);
    const gamma_guess = Complex.scale(w0, 1 / absZ0);

    // If q has a non-trivial imaginary part and our bootstrap produced an
    // essentially-real solution, kick z_0 and γ in the imaginary direction so
    // Newton has a starting point in the complex basin. Without this, real-h
    // + complex-q cases get stuck in the real basin (where the q-equation's
    // imaginary part can't be satisfied).
    const qImAbs = Math.abs(q.im);
    if (qImAbs > 1e-6) {
      const isRealish = Math.abs(z0_guess.im) < 1e-3 && Math.abs(gamma_guess.im) < 1e-3;
      if (isRealish) {
        const kick = Math.min(0.15, qImAbs);
        z0_guess.im += kick * Math.sign(q.im);
        const r0 = Math.hypot(z0_guess.re, z0_guess.im);
        if (r0 > 0.9) { z0_guess.re *= 0.85 / r0; z0_guess.im *= 0.85 / r0; }
        // re-derive γ to keep |γ|·|z_0| ≈ |w_0| (φ_0 constraint magnitude)
        const newAbsZ0 = Math.hypot(z0_guess.re, z0_guess.im);
        const phaseKick = 0.5 * q.im;            // small phase perturbation to γ
        gamma_guess.re = w0.re / newAbsZ0;
        gamma_guess.im = w0.im / newAbsZ0 + phaseKick;
      }
    }

    return {
      family: 'boundedLQD_singular',
      w0: Complex.clone(w0),
      q: Complex.clone(q),
      z0: z0_guess,
      gamma: gamma_guess,
      branches: zj_guess.map((z, j) => ({
        z, A: A_guess[j].map(Complex.clone),
      })),
    };
  }

  function perturbedInitialGuess_LQDS(hData, norm, rng, r) {
    r = r || 0;
    const base = initialGuess_LQDS(hData, norm);
    const sigma = 0.15 + 0.25 * r;

    // perturb z_0
    base.z0 = {
      re: base.z0.re + sigma * 0.3 * (rng() - 0.5),
      im: base.z0.im + sigma * 0.3 * (rng() - 0.5),
    };
    const rz0 = Math.hypot(base.z0.re, base.z0.im);
    if (rz0 > 0.9) { base.z0.re *= 0.85 / rz0; base.z0.im *= 0.85 / rz0; }
    if (rz0 < 0.05) { base.z0.re = 0.1; base.z0.im = 0.05 * (rng() - 0.5); }

    // perturb γ (multiplicatively)
    const dgRe = 1 + sigma * 0.3 * (rng() - 0.5);
    const dgIm = sigma * 0.3 * (rng() - 0.5);
    base.gamma = Complex.mul(base.gamma, { re: dgRe, im: dgIm });

    // perturb z_j and A
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

  function diverseInitialGuess_LQDS(hData, norm, rng) {
    const w0 = norm.w0;
    const q = norm.q;
    // z_0: log-uniform |z_0| ∈ [0.05, 0.9], uniform phase
    const mz0 = Math.exp(Math.log(0.05) + rng() * Math.log(0.9 / 0.05));
    const pz0 = 2 * Math.PI * rng();
    const z0 = { re: mz0 * Math.cos(pz0), im: mz0 * Math.sin(pz0) };
    const gamma = Complex.scale(w0, 1 / Math.max(mz0, 1e-3));

    const branches = [];
    const zMin = 0.05, zMax = 0.95;
    const aMin = 0.1, aMax = 3.0;
    for (const p of hData.poles) {
      const mz = Math.exp(Math.log(zMin) + rng() * Math.log(zMax / zMin));
      const pz = 2 * Math.PI * rng();
      const z = { re: mz * Math.cos(pz), im: mz * Math.sin(pz) };
      const A = [];
      for (let k = 0; k < p.principal.length; k++) {
        const ma = Math.exp(Math.log(aMin) + rng() * Math.log(aMax / aMin));
        const pa = 2 * Math.PI * rng();
        A.push({ re: ma * Math.cos(pa), im: ma * Math.sin(pa) });
      }
      branches.push({ z, A });
    }
    return {
      family: 'boundedLQD_singular',
      w0: Complex.clone(w0), q: Complex.clone(q),
      z0, gamma, branches,
    };
  }

  QD.Seeds = QD.Seeds || {};
  QD.Seeds.boundedLQD_singular = {
    initialGuess:          initialGuess_LQDS,
    perturbedInitialGuess: perturbedInitialGuess_LQDS,
    diverseInitialGuess:   diverseInitialGuess_LQDS,
  };
})();
