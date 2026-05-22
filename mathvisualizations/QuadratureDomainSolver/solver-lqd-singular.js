// =============================================================================
// solver-lqd-singular.js -- Bounded SINGULAR log-weighted quadrature domains
//
// An LQD with 0 ∈ Ω satisfies
//     ∫_Ω f(w) / |w|² dA = ∮_∂Ω f(w) h(w) dw,   ∀ f ∈ L¹_a(Ω; ρ₀)
// with h ∈ Rat(Ω) where h is now allowed an extra simple pole at 0 (inside Ω)
// with residue q ∈ ℂ. Functions f ∈ L¹_a(Ω; ρ₀) are forced to vanish at 0
// (else f/|w|² is not integrable near 0), so the identity has no f(0) term.
//
// Math summary (Chapter V of Graven 2026; Theorems 5.4.1/5.4.2/5.6.2):
//
// • Riemann-map parametrization:
//     φ(z) = γ · b_{z_0}(z) · exp(r#(z))
//   with
//     b_{z_0}(z) = -(conj(z_0)/|z_0|) · (z - z_0) / (1 - conj(z_0) z)   (thesis form)
//     r# ∈ Rat_0(𝔻) with r#(0) = 0 (same shape as non-singular LQD case)
//     z_0 ∈ 𝔻 the unique preimage of 0 (a simple zero by univalence)
//     γ ∈ ℂ* a multiplicative constant
//
//   We defer the degenerate case z_0 = 0 to a later pass.
//
// • Branch-cut-free reformulation:
//     ψ(z) := φ(z) / b_{z_0}(z)   is holomorphic and non-vanishing on 𝔻̄,
//     so log ψ = log γ + r#(z) is single-valued analytic. All inverse-Faber
//     matching is done on log ψ — never on log φ directly.
//
// • Locator (●):  φ(z_j) = a_j        (2n real equations)
//
// • Coefficient (★):  modified-residue Faber match at each finite pole,
//   identical formula to the non-singular LQD case
//   A_{j,k}^target = Σ_{s=k..m_j} (s/k) · D_{j,s} · [t^s] ψ̃_j^k(t)
//   with D_{j,s} = a_j C_{j,s} + C_{j,s+1}. The change vs. non-singular comes
//   from phiTaylorAt at z_j including the b_{z_0} factor in the local data;
//   the formula itself is unchanged.   (2M real equations, M = Σ m_j)
//
// • q-equation (●₀): set G(z) := ln(φ · φ#)(z). The Blaschke identity
//     b_{z_0}(z) · b_{z_0}#(z) ≡ 1  (mod 2πi)
//   collapses the log-singularities at z_0, leaving
//     G(z) = ln|γ|² + r#(z) + r(z),      r(z) = conj(r#(1/conj(z)))
//   which is single-valued analytic near z_0. Then
//     S_0(w) = G(ψ(w)) / w,    Res_{w=0} S_0(w) = G(z_0),
//   so the residue of h at 0 is
//
//        q  =  ln|γ|²  +  r#(z_0)  +  conj(r#(1/conj(z_0)))      … (●₀)
//
//   2 real equations.
//
// • Normalization φ(0) = w_0:  since b_{z_0}(0) = |z_0| and exp(r#(0)) = 1,
//     φ(0) = γ · |z_0|.   2 real equations: γ · |z_0| − w_0 = 0.
//
// • Gauge arg φ'(0) = 0:  φ'(0) = γ · [b_{z_0}'(0) + |z_0| · r#'(0)] where
//     b_{z_0}'(0) = -(conj(z_0)/|z_0|) · (1 − |z_0|²),
//     r#'(0) = Σ_j conj(A_{j,1}).
//   1 real equation: Im φ'(0) = 0 (the sign of Re φ'(0) is resolved in
//   canonicalizePhi).
//
// • Identity test class: monomials f(w) = w^k, k = 1, 2, 3. These are analytic
//   in Ω, vanish at 0, and integrable against ρ₀. RHS is sum of residues at
//   the finite poles a_j only — no contribution from q/w (since w^k · q/w =
//   q·w^{k-1} has no residue at 0 for k ≥ 1).
//
// Unknowns (count): 2n (z_j) + 2M (A_{j,k}) + 2 (z_0) + 2 (γ)  =  2n+2M+4
// Equations (count): 2n (●) + 2M (★) + 2 (●₀) + 2 (φ_0) + 1 (gauge) = 2n+2M+5
// One overdetermined real equation, resolved by Newton's least-squares fallback
// (same handling as the non-singular case).
//
// Registers Family.boundedLQD_singular. Must load AFTER solver-lqd.js so the
// dispatch order is: boundedLQD_singular (most specific) → boundedLQD →
// boundedQD (catch-all).
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family) {
    throw new Error("solver-lqd-singular.js: solver.js must be loaded first");
  }

  // ===========================================================================
  // 0. Blaschke helpers — lifted to QD.LqdCommon (shared with unbounded
  //    singular LQDs since the formula is identical for z_0 ∈ 𝔻 or 𝔻*).
  // ===========================================================================
  const blaschkeEval   = QD.LqdCommon.blaschkeEval;
  const blaschkeTaylor = QD.LqdCommon.blaschkeTaylor;

  // ===========================================================================
  // 1. r# evaluator and (★) Faber loop are shared via QD.LqdCommon.
  // ===========================================================================
  const evalRHash = QD.LqdCommon.evalRHash;
  const rHashTaylorAt = QD.LqdCommon.rHashTaylorAt;

  // ===========================================================================
  // 2. Phi evaluation
  // ===========================================================================

  function evalPhi_LQDS(z, phi) {
    const r = evalRHash(z, phi);
    const ea = Math.exp(r.re);
    const expR = { re: ea * Math.cos(r.im), im: ea * Math.sin(r.im) };
    const b = blaschkeEval(z, phi.z0);
    return Complex.mul(Complex.mul(phi.gamma, b), expR);
  }

  // φ(z_c + t) = γ · b_{z_0}(z_c + t) · exp(r#(z_c + t))
  //            = γ · b_{z_0}(z_c + t) · exp(r#(z_c)) · exp(r#(z_c+t) − r#(z_c))
  function phiTaylorAt_LQDS(zc, phi, L) {
    const bT = blaschkeTaylor(zc, phi.z0, L);                   // Taylor of b
    const rT = rHashTaylorAt(zc, phi, L);                       // Taylor of r#

    // Strip constant from r# to feed Taylor.exp
    const r0 = rT[0];
    const rTilde = rT.slice();
    rTilde[0] = { re: 0, im: 0 };
    const expRTilde = Taylor.exp(rTilde, L);

    // Prefactor γ · exp(r#(z_c))
    const ea = Math.exp(r0.re);
    const expR0 = { re: ea * Math.cos(r0.im), im: ea * Math.sin(r0.im) };
    const scale = Complex.mul(phi.gamma, expR0);

    // φ_T = scale · b · exp(rTilde)
    const bTimesExp = Taylor.mul(bT, expRTilde, L);
    const out = new Array(L + 1);
    for (let l = 0; l <= L; l++) out[l] = Complex.mul(scale, bTimesExp[l]);
    return out;
  }

  // ===========================================================================
  // 3. Targets (★) -- modified residues + standard Faber inversion on
  //                   phiTaylorAt_LQDS at each z_j
  // ===========================================================================

  // (★) Faber loop with modified residues — shared via QD.LqdCommon.
  function computeTargetA_LQDS(phi, hData) {
    return QD.LqdCommon.computeFaberTargetA(phi, hData, phiTaylorAt_LQDS);
  }

  // ===========================================================================
  // 4. Residual
  // ===========================================================================
  //
  // Order:
  //   (●)   2n     locator at finite poles
  //   (★)   2M     coefficient match (Faber)
  //   (●₀)  2      q-equation at origin
  //   (φ_0) 2      γ · |z_0| − w_0 = 0
  //   (G)   1      gauge: Im(φ'(0)) = 0
  function residual_LQDS(phi, hData, options) {
    options = options || {};
    const enforceGauge = options.enforceGauge !== false;
    const enforceQEq   = options.enforceQEq   !== false;
    const out = [];

    // (●) locator
    for (let j = 0; j < hData.poles.length; j++) {
      const zj = phi.branches[j].z;
      const phiZj = evalPhi_LQDS(zj, phi);
      const diff = Complex.sub(phiZj, hData.poles[j].a);
      out.push(diff.re, diff.im);
    }

    // (★) coefficient
    const target = computeTargetA_LQDS(phi, hData);
    for (let j = 0; j < hData.poles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], target[j][k]);
        out.push(diff.re, diff.im);
      }
    }

    // (●₀) q-equation: q − ln|γ|² − r#(z_0) − conj(r#(1/conj(z_0))) = 0
    if (enforceQEq) {
      const z0 = phi.z0;
      const r_at_z0 = evalRHash(z0, phi);
      const absZ02 = Complex.abs2(z0);
      // 1/conj(z_0)
      const oneOverConjZ0 = Complex.scale(z0, 1 / absZ02);
      const r_at_inv = evalRHash(oneOverConjZ0, phi);
      const rConj = Complex.conj(r_at_inv);
      const sum = Complex.add(r_at_z0, rConj);
      const lnGamma2 = Math.log(Complex.abs2(phi.gamma));
      const lhs = { re: sum.re + lnGamma2, im: sum.im };
      const diff = Complex.sub(phi.q, lhs);
      out.push(diff.re, diff.im);
    }

    // (φ_0) normalization: γ · |z_0| = w_0
    const absZ0 = Complex.abs(phi.z0);
    const phi0 = Complex.scale(phi.gamma, absZ0);
    const phi0Diff = Complex.sub(phi0, phi.w0);
    out.push(phi0Diff.re, phi0Diff.im);

    // (G) gauge: Im(φ'(0)) = 0
    //   φ'(0) = γ · [b_{z_0}'(0) + |z_0| · r#'(0)]
    //   b_{z_0}'(0) = -(conj(z_0)/|z_0|) · (1 − |z_0|²)
    //   r#'(0)      = Σ_j conj(A_{j,1})
    if (enforceGauge) {
      const z0 = phi.z0;
      const z0C = Complex.conj(z0);
      const absZ0g = Complex.abs(z0);
      const oneMinusAbs2 = 1 - Complex.abs2(z0);
      // b'_{z_0}(0)
      const bprime0 = Complex.scale(z0C, -oneMinusAbs2 / absZ0g);
      // r#'(0)
      let rprime0 = { re: 0, im: 0 };
      for (const br of phi.branches) {
        if (br.A.length > 0) rprime0 = Complex.add(rprime0, Complex.conj(br.A[0]));
      }
      // bracket
      const bracket = Complex.add(bprime0, Complex.scale(rprime0, absZ0g));
      const phiPrime0 = Complex.mul(phi.gamma, bracket);
      out.push(phiPrime0.im);
    }

    return out;
  }

  // Sign canonicalization: Re(φ'(0)) > 0. If it comes out negative, multiply
  // γ by −1 (which is a Z/2 rotation of the disc, leaving |γ| and z_0 fixed
  // and flipping the sign of every A_{j,k} via the same gauge action that
  // canonicalizePhi uses in the non-singular case). This restores the
  // standard Riemann-map normalization φ'(0) > 0.
  function canonicalizePhi_LQDS(phi) {
    // Compute φ'(0).re sign.
    const z0 = phi.z0;
    const z0C = Complex.conj(z0);
    const absZ0g = Complex.abs(z0);
    const oneMinusAbs2 = 1 - Complex.abs2(z0);
    const bprime0 = Complex.scale(z0C, -oneMinusAbs2 / absZ0g);
    let rprime0 = { re: 0, im: 0 };
    for (const br of phi.branches) {
      if (br.A.length > 0) rprime0 = Complex.add(rprime0, Complex.conj(br.A[0]));
    }
    const bracket = Complex.add(bprime0, Complex.scale(rprime0, absZ0g));
    const phiPrime0 = Complex.mul(phi.gamma, bracket);

    if (phiPrime0.re >= 0) return phi;

    // The Z/2 we have available: z → -z, which acts as
    //   z_0 → -z_0,   z_j → -z_j,   A_{j,k} → (-1)^k A_{j,k},
    //   γ unchanged (φ depends on γ multiplicatively; the rotation acts on the
    //   disc, not the image). Verify φ'(0) sign flips under this:
    //   b_{z_0}(z) → b_{-z_0}(-z) = ... b's response under z → -z gives the
    //   right sign flip because (z - z_0) → -(z + z_0) and (1 - conj(z_0) z)
    //   → (1 + conj(z_0) z), and the leading prefactor likewise flips. The
    //   net effect: φ(z) → -γ · b̃(z) · exp(...) where b̃ is what b becomes
    //   under z → -z; collecting signs gives φ'(0) → -φ'(0).
    return {
      family: 'boundedLQD_singular',
      w0: Complex.clone(phi.w0),
      q:  Complex.clone(phi.q),
      z0: Complex.neg(phi.z0),
      gamma: Complex.clone(phi.gamma),
      branches: phi.branches.map(br => ({
        z: Complex.neg(br.z),
        A: br.A.map((a, k) => (k % 2 === 0 ? Complex.neg(a) : Complex.clone(a))),
      })),
    };
  }

  // ===========================================================================
  // 5. Pack / unpack — schema-driven (R3)
  // ===========================================================================
  // Layout: [{z_j.re, z_j.im}_{j=1..n}, z_0.re, z_0.im, γ.re, γ.im,
  //          {A_{j,k}.re, A_{j,k}.im}]
  // z_j first so Newton's enforceInDisk clamps the right components. z_0
  // gets its own clamp via the schema (both upper |z_0| ≤ 0.9999 and lower
  // |z_0| ≥ 1e-3 — z_0 = 0 is deferred to a later pass).
  const SCHEMA_LQDS = [
    { kind: 'branchesZ', clamp: { side: 'in', cap: 0.9999 } },
    { kind: 'complex', name: 'z0',    clamp: { side: 'in', cap: 0.9999, minR: 1e-3 } },
    { kind: 'complex', name: 'gamma' },
    { kind: 'branchesA' },
  ];

  function packPhi_LQDS(phi)         { return QD.packPhiBySchema(phi, SCHEMA_LQDS); }
  function unpackPhi_LQDS(v, template) {
    return QD.unpackPhiBySchema(v, template, SCHEMA_LQDS);
  }

  // ===========================================================================
  // 6. Initial guess via companion-QD bootstrap
  // ===========================================================================
  //
  // Strategy:
  //   1. Solve the bounded-QD problem with the same finite poles/residues
  //      (ignoring q). Gives φ_QD: 𝔻 → Ω_QD with center w_0_QD.
  //   2. z_0 ← Newton solve of φ_QD(z) = 0, starting from z = 0.
  //   3. z_j and A_{j,k} ← read from φ_QD (treating the QD parametrization
  //      as a starting r# in log space, which is rough but OK for warm-start).
  //   4. γ ← w_0 / |z_0|  (from the φ_0 equation).
  //
  // If bootstrap fails, fall back to a geometric guess.

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

  function diverseInitialGuess_LQDS(hData, norm, rng, r) {
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

  function continuationSolve_LQDS(hData, norm, opts) {
    return { success: false, error: "continuation not implemented for singular LQD", trace: [] };
  }

  // ===========================================================================
  // 7. Identity verification: monomials f(w) = w^k, k = 1, 2, 3.
  // ===========================================================================
  //
  // LHS:  (1/N) Σ_n w_n^k · (ln|w_n|²/w_n) · φ'(z_n) · z_n
  //       = (1/N) Σ_n w_n^{k-1} · ln|w_n|² · φ'(z_n) · z_n
  //   (trapezoidal-rule discretization of ∮_∂Ω f · (ln|w|²/w) dw via the
  //    substitution dw = φ'(z) z · i dθ where z = e^{iθ})
  //
  // Wait — let me re-derive. ∮_∂Ω f(w) g(w) dw with w = φ(e^{iθ}):
  //   dw = φ'(e^{iθ}) · i · e^{iθ} dθ.
  // The Green's-theorem boundary integral is
  //   ∮_∂Ω f(w) · (ln|w|²/w) dw   (no 2πi factor; this is the direct line integral)
  //
  // Numerically: ∫₀^{2π} ... dθ ≈ (2π/N) Σ_n (integrand at θ_n).
  // Substituting g(w) = ln|w|²/w and dw = i·e^{iθ}·φ'(z) dθ:
  //   ∮ f · g · dw = ∫₀^{2π} f(w_n) · ln|w_n|²/w_n · i·z_n · φ'(z_n) dθ
  //                ≈ (2π/N) Σ_n  i z_n φ'(z_n) · f(w_n) · ln|w_n|²/w_n
  //
  // So LHS_test = (2π/N) · i · Σ_n z_n φ'(z_n) f(w_n) ln|w_n|²/w_n.
  //
  // (The non-singular verifier in solver-lqd.js used (1/N) without the 2π i —
  //  that's a different overall scale but the relative residual is the same;
  //  we'll match its convention here for consistency, since identity comparison
  //  is invariant under a global factor as long as LHS and RHS use the same.)
  //
  // RHS = sum of residues of f · h at the finite poles a_j:
  //   For f(w) = w^k and h(w) = Σ_s C_{j,s}/(w-a_j)^s + (other poles) + q/w,
  //   Res_{w=a_j} [w^k · C_{j,s}/(w-a_j)^s] = C_{j,s} · [(d^{s-1}/dw^{s-1}) w^k]/(s-1)! |_{w=a_j}
  //                                          = C_{j,s} · k·(k-1)···(k-s+2) · a_j^{k-s+1} / (s-1)!
  //                                          = C_{j,s} · binom(k, s-1) · a_j^{k-s+1}
  //   (where binom(k, s-1) = k!/((s-1)!(k-s+1)!), valid for k ≥ s-1, else 0).
  //   The 2πi factor on RHS comes from "residue theorem" but in our LHS we
  //   have ∮ directly (not (1/2πi)∮), so RHS = 2πi · Σ_{j,s} Res. The factor
  //   matches our (2π·i/N) scaling on LHS.
  //
  // Final form used below:
  //   LHS = (2π/N) · i · Σ_n z_n · φ'(z_n) · w_n^{k-1} · ln|w_n|²
  //   RHS = 2πi · Σ_{j,s} C_{j,s} · binom(k, s-1) · a_j^{k-s+1}
  //         + 2πi · q · [k=0]  (zero for k ≥ 1, so absent)

  // Refactored to use QD.LqdCommon.verifyIdentityGeneric (R4): the skeleton
  // owns boundary sampling, the LHS trapezoid + i·(2π/N) scaling, the
  // scale-floor heuristic, and the result struct. This family just declares
  // its boundary kernel and its test-function list with RHS values.
  //
  // Test class: monomials f(w) = w^k, k = 1, …, maxDegree. (vanish at 0 as
  // required by f ∈ L¹_a(Ω; ρ₀); analytic in Ω.)
  // Boundary kernel: ln|w|² / w.
  // RHS at degree k: 2πi · Σ_{j, s} C_{j,s} · binom(k, s-1) · a_j^{k-s+1}.
  // (q/w contributes 0 for f = w^k with k ≥ 1.)
  function verifyQuadratureIdentity_LQDS(phi, hData, options) {
    options = options || {};
    const maxOrder = options.maxDegree ?? 3;
    const I = { re: 0, im: 1 };
    const twoPiI = { re: 0, im: 2 * Math.PI };

    return QD.LqdCommon.verifyIdentityGeneric(phi, hData, options, {
      phiTaylorFn: phiTaylorAt_LQDS,
      boundaryKernel(w) {
        const absW2 = Complex.abs2(w);
        if (absW2 < 1e-30) return null;            // skip near origin
        return Complex.scale(Complex.inv(w), Math.log(absW2));
      },
      buildTestFunctions(phi, hData) {
        const tests = [];
        for (let k = 1; k <= maxOrder; k++) {
          // Closed-form RHS = 2πi · Σ residues.
          let rhsSum = { re: 0, im: 0 };
          for (const pole of hData.poles) {
            for (let sIdx = 0; sIdx < pole.principal.length; sIdx++) {
              const s = sIdx + 1;
              const expo = k - s + 1;
              if (expo < 0) continue;
              const coef = QD.binomialCoeff(k, s - 1);
              if (coef === 0) continue;
              const C = pole.principal[sIdx];
              const aPow = Complex.pow(pole.a, expo);
              rhsSum = Complex.add(rhsSum, Complex.scale(Complex.mul(C, aPow), coef));
            }
          }
          const rhs = Complex.mul(twoPiI, rhsSum);
          tests.push({
            label: 'w^' + k,
            f:     (w) => Complex.pow(w, k),
            residueRhs: rhs,
            tag: { k },
          });
        }
        return tests;
      },
      resultFlags: { lqdSingular: true, maxDeg: maxOrder },
    });
  }

  // ===========================================================================
  // 8. Register Family.boundedLQD_singular
  // ===========================================================================
  QD.Family.boundedLQD_singular = {
    name: 'boundedLQD_singular',
    enforceInDisk: true,
    enforceOutDisk: false,
    matches(opts) {
      return !!(opts && opts.lqd && opts.singular && !opts.unbounded);
    },
    normalizeOpts(opts, hData) {
      const w0 = opts.w0;
      if (!w0) throw new Error("Family.boundedLQD_singular: opts.w0 required");
      const q = opts.q || { re: 0, im: 0 };
      // q = 0 reduces to a degenerate singular case (no log pole at origin);
      // we still allow it as the default since the family is dispatched purely
      // by opts.singular.
      return { lqd: true, singular: true, w0: Complex.clone(w0), q: Complex.clone(q) };
    },
    evalPhi: evalPhi_LQDS,
    phiTaylorAt: phiTaylorAt_LQDS,
    computeTargets(phi, hData) {
      return { A: computeTargetA_LQDS(phi, hData), F: null };
    },
    residual: residual_LQDS,
    packPhi: packPhi_LQDS,
    unpackPhi: unpackPhi_LQDS,
    canonicalizePhi: canonicalizePhi_LQDS,
    initialGuess: initialGuess_LQDS,
    perturbedInitialGuess: perturbedInitialGuess_LQDS,
    diverseInitialGuess: diverseInitialGuess_LQDS,
    continuationSolve: continuationSolve_LQDS,
    verifyQuadratureIdentity: verifyQuadratureIdentity_LQDS,
  };

  if (QD.registerFamily) {
    QD.registerFamily('boundedLQD_singular');
  } else {
    throw new Error("solver-lqd-singular.js: QD.registerFamily not found");
  }

})();
