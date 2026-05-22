// =============================================================================
// aqd-bounded.js -- Bounded non-singular Algebraic Quadrature Domain solver.
//
// Family.boundedAQD: Ω ⊂ ℂ bounded, simply connected, with weight
//
//   ρ(w) = |α(w)|²,    α = R′,    R rational
//
// Restrictions for this file (Stage 2):
//   • Ω̄ contains no poles of R          (non-singular)
//   • Ω̄ contains no zeros of R′         (no inner-factor cusps)
//   • h: rational quadrature function with finite poles a_j ∈ Ω.
//   • h has no polynomial part for now (Stage 6 will lift this).
//
// =============================================================================
// STAGE 2a DERIVATION — committed to code as the comment block below.
// -----------------------------------------------------------------------------
// SETUP
//   φ : 𝔻 → Ω is the Riemann map with φ(0) = w₀ (user input).
//   z_j := φ⁻¹(a_j) ∈ 𝔻.
//   Define U(z) := R(w₀) + r#(z) − r#(0)  ─── the ∞-gauge-absorbed form
//   of the inverse-problem RHS. By construction U(0) = R(w₀).
//
//   Parametric form of r# (identical to bounded LQD r#):
//
//     r#(z) = Σ_j Σ_{k=1..m_j} conj(A_{j,k}) · z^k / (1 − conj(z_j) z)^k.
//
//   r# is analytic in 𝔻̄ (poles only at 1/conj(z_j) ∈ 𝔻*).
//
// GOVERNING BOUNDARY IDENTITY (Thm 6.4.1, non-singular reduction)
//   R∘φ(z) = U(z)            for z on |z| = 1.
//
// EQUATION-BUILDING
// (1)  Locator (●) per pole j (complex; 2 real):
//          U(z_j) − R(a_j) = 0.
//      Since r# is finite at z_j, U(z_j) is finite. No R-inversion needed.
//
// (2)  Faber match (★) per pole j (m_j complex; 2 m_j real):
//          A_{j,k} = QD.Faber.inverseFaberAtPole(P_{j,•}, phiTilde_j)
//      where:
//
//      • Residues P_{j,t} are the principal-part coefficients of h/R′ at a_j,
//        computed by Aqd.principalPartsHOverRprime.
//        (LQD analogue used modified residues D_{j,s} = a_j C_{j,s} + C_{j,s+1};
//         AQDs use the Taylor convolution of h with 1/R' at a_j instead.)
//
//      • phiTilde_j = Taylor of (φ(z_j + t) − a_j) in t, with phiTilde[0] = 0,
//        derived from U's Taylor at z_j and R⁻¹'s Taylor at U(z_j):
//           u_l = [t^l] U(z_j + t)         (computed from r# Taylor)
//           w   = R⁻¹(u_0)                 (Newton invert)
//           R̃(h) := R(w + h) − R(w)       (Taylor at w, constant 0)
//           R̃⁻¹                            (Taylor series-inverse)
//           Δ(t) := u_1·t + u_2·t² + …    (deviation of U from u_0)
//           phiTilde_no_const = R̃⁻¹(Δ(t)) (Taylor composition; constant 0)
//           Total Taylor of (φ(z_j + t) − a_j): (w − a_j) + R̃⁻¹(Δ(t)).
//         The Faber primitive treats phiTilde[0] as a "drop the constant"
//         convention — at the SOLUTION w = a_j the constant is 0. Off-solution
//         the constant is nonzero, but Newton drives the locator equation to
//         zero in parallel, which absorbs it.
//
// (3)  Gauge (1 real): rotational disk gauge, same as bounded LQDs:
//          Σ_j Im(A_{j,1}) = 0.
//      Sign canonicalized after Newton converges.
//
// PARAMETER COUNT
//   Per pole j:  2(m_j + 1) real unknowns, 2(m_j + 1) real equations.
//   Plus 1 gauge equation. System is over-determined by 1 and solved in
//   least squares, matching the bounded-LQD pattern.
//
// RECOVERY OF φ (for visualization & identity check, NOT in Newton)
//   At any z ∈ 𝔻̄:
//     • U(z) = R(w₀) + r#(z) − r#(0).
//     • φ(z) = R⁻¹(U(z)) via Aqd.invertR.
//   In evalPhi_AQD we use a small homotopy from z=0 (where φ=w₀) so the seed
//   is robust without per-call caching state.
//
// REUSE
//   • r# evaluator: QD.LqdCommon.evalRHash and rHashTaylorAt — verbatim.
//   • Inverse Faber transform: QD.Faber.inverseFaberAtPole.
//   • Newton driver, multistart, pack/unpack schema: solver.js.
// =============================================================================
'use strict';

(function () {

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family) {
    throw new Error("aqd-bounded.js: solver.js must be loaded first");
  }
  if (!QD.LqdCommon) throw new Error("aqd-bounded.js: solver-lqd-common.js must be loaded first");
  if (!QD.Faber)     throw new Error("aqd-bounded.js: solver-faber.js must be loaded first");
  if (!QD.Aqd)       throw new Error("aqd-bounded.js: aqd-common.js must be loaded first");

  const Aqd          = QD.Aqd;
  const evalRHash    = QD.LqdCommon.evalRHash;
  const rHashTaylorAt= QD.LqdCommon.rHashTaylorAt;

  // ===========================================================================
  // 1. U(z) = R(w₀) + r#(z) − r#(0)   and  its Taylor expansion at z₀.
  // ===========================================================================
  function evalU(z, phi) {
    const rZ = evalRHash(z, phi);
    const rZero = evalRHash({ re: 0, im: 0 }, phi);
    const Rw0 = Aqd.evalR(phi.R, phi.w0);
    // U = Rw0 + rZ - rZero
    return Complex.add(Rw0, Complex.sub(rZ, rZero));
  }

  function uTaylorAt(z0, phi, L) {
    const rT = rHashTaylorAt(z0, phi, L);
    const rZero = evalRHash({ re: 0, im: 0 }, phi);
    const Rw0 = Aqd.evalR(phi.R, phi.w0);
    const out = rT.slice();
    out[0] = Complex.add(out[0], Complex.sub(Rw0, rZero));
    return out;
  }

  // ===========================================================================
  // 2. φ(z) and its Taylor at z₀  via R-inversion of U(z).
  // ===========================================================================
  // For a stateless evalPhi we use an 8-step homotopy along z(t) = t·z₀ to
  // walk the R-inversion seed from φ(0) = w₀ to φ(z₀).
  function evalPhi_AQD(z, phi) {
    return _evalPhiInternal(z, phi, 8);
  }

  function _evalPhiInternal(z, phi, nSteps) {
    let w = Complex.clone(phi.w0);
    for (let k = 1; k <= nSteps; k++) {
      const tk = k / nSteps;
      const zk = { re: tk * z.re, im: tk * z.im };
      const Uzk = evalU(zk, phi);
      const inv = Aqd.invertR(phi.R, Uzk, w);
      if (!inv.converged) {
        // Try with more iterations from origin as a fallback. If still fails,
        // return NaN — solver layer treats this as residual divergence.
        const fallback = Aqd.invertR(phi.R, Uzk, phi.w0, { maxIter: 64 });
        if (!fallback.converged) return { re: NaN, im: NaN };
        w = fallback.w;
      } else {
        w = inv.w;
      }
    }
    return w;
  }

  // Taylor of φ at z₀ up to order L: returns [φ(z₀), φ'(z₀)/1!, ...].
  // This is the visualization / general-purpose evaluator (R-inverts at z₀).
  function phiTaylorAt_AQD(z0, phi, L) {
    const uT = uTaylorAt(z0, phi, L);
    const u0 = uT[0];
    let w;
    const direct = Aqd.invertR(phi.R, u0, phi.w0, { maxIter: 32 });
    if (direct.converged) { w = direct.w; }
    else { w = _evalPhiInternal(z0, phi, 16); }
    if (!isFinite(w.re) || !isFinite(w.im)) {
      const bad = Taylor.zero(L + 1);
      bad[0] = { re: NaN, im: NaN };
      return bad;
    }
    const rT = Aqd.rTaylorAt(phi.R, w, L);
    if (!rT) {
      const bad = Taylor.zero(L + 1);
      bad[0] = { re: NaN, im: NaN };
      return bad;
    }
    const rShifted = Taylor.zero(L + 1);
    for (let i = 1; i <= L; i++) rShifted[i] = rT[i];
    let rShiftedInv;
    try { rShiftedInv = Taylor.invert(rShifted, L); }
    catch (e) {
      const bad = Taylor.zero(L + 1);
      bad[0] = { re: NaN, im: NaN };
      return bad;
    }
    const delta = Taylor.zero(L + 1);
    for (let i = 1; i <= L; i++) delta[i] = uT[i];
    const composed = Taylor.compose(rShiftedInv, delta, L);
    const out = Taylor.zero(L + 1);
    out[0] = Complex.add(w, composed[0] || { re: 0, im: 0 });
    for (let i = 1; i <= L; i++) out[i] = composed[i];
    return out;
  }

  // ===========================================================================
  // 3. (★) Faber target for one h-pole — uses AQD principal-parts of h/R'.
  // ---------------------------------------------------------------------------
  // KEY: to keep the Newton Jacobian numerically clean, we anchor R⁻¹ at the
  // FIXED point a_j (a problem-data constant, not a Newton variable). The
  // phiTilde we feed to the Faber primitive is:
  //
  //   phiTilde[0] = 0   (Faber drops this anyway)
  //   phiTilde[l] for l ≥ 1 = [t^l] R⁻¹(R(a_j) + Δ(t)) where
  //     Δ(t) := u_1·t + u_2·t² + … is U's Taylor at z_j stripped of its
  //     constant, and R⁻¹ is anchored at R(a_j) so its Taylor coefficients
  //     ρ_l = [(R⁻¹)^(l)(R(a_j))/l!] depend ONLY on R and a_j — they're
  //     computed once per solve.
  //
  // At the solution U(z_j) = R(a_j), so this matches the "true" phiTilde
  // exactly. Off-solution, the locator residual U(z_j) − R(a_j) appears
  // separately in (●); Newton drives both to zero. The Faber match remains
  // a smooth analytic function of the Newton variables (z_j, A_{j,k}), with
  // no R-inversion in the inner loop.
  // ===========================================================================
  function computeTargetA_AQD(phi, hData) {
    const target = [];
    for (let j = 0; j < hData.poles.length; j++) {
      const pole = hData.poles[j];
      const m = pole.principal.length;
      if (m === 0) { target.push([]); continue; }

      const P = Aqd.principalPartsHOverRprime(pole, phi.R);
      if (!P) {
        target.push(new Array(m).fill({ re: NaN, im: NaN }));
        continue;
      }

      // R-Taylor at a_j (fixed) → R̃⁻¹ Taylor at 0 (also fixed).
      const rAtA = Aqd.rTaylorAt(phi.R, pole.a, m);
      if (!rAtA) {
        target.push(new Array(m).fill({ re: NaN, im: NaN }));
        continue;
      }
      const rTilde = Taylor.zero(m + 1);
      for (let i = 1; i <= m; i++) rTilde[i] = rAtA[i];
      let rTildeInv;
      try { rTildeInv = Taylor.invert(rTilde, m); }
      catch (e) {
        target.push(new Array(m).fill({ re: NaN, im: NaN }));
        continue;
      }

      // U-Taylor at z_j (depends on Newton vars).
      const zj = phi.branches[j].z;
      const uT = uTaylorAt(zj, phi, m);
      const delta = Taylor.zero(m + 1);
      for (let i = 1; i <= m; i++) delta[i] = uT[i];

      // phiTilde[1..m] = [t^l] (R̃⁻¹ ∘ Δ)(t), constant 0 by convention.
      const composed = Taylor.compose(rTildeInv, delta, m);
      const phiTilde = Taylor.zero(m + 1);
      for (let i = 1; i <= m; i++) phiTilde[i] = composed[i];

      target.push(QD.Faber.inverseFaberAtPole(P, phiTilde));
    }
    return target;
  }

  // ===========================================================================
  // 4. Residual / pack / unpack / canonicalize  (mirrors bounded LQD layout)
  // ===========================================================================
  function residual_AQD(phi, hData, options) {
    options = options || {};
    const enforceGauge = options.enforceGauge !== false;
    const out = [];

    // (●) locator: U(z_j) − R(a_j) = 0
    for (let j = 0; j < hData.poles.length; j++) {
      const zj = phi.branches[j].z;
      const Uzj = evalU(zj, phi);
      const Raj = Aqd.evalR(phi.R, hData.poles[j].a);
      const diff = Complex.sub(Uzj, Raj);
      out.push(diff.re, diff.im);
    }

    // (★) Faber match
    const target = computeTargetA_AQD(phi, hData);
    for (let j = 0; j < hData.poles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], target[j][k]);
        out.push(diff.re, diff.im);
      }
    }

    // Gauge: Σ Im(A_{j,1}) = 0 (same as bounded LQD pattern, up to sign).
    if (enforceGauge) {
      let imSum = 0;
      for (const br of phi.branches) {
        if (br.A.length > 0) imSum += -br.A[0].im;       // Im(conj(a)) = -Im(a)
      }
      out.push(imSum);
    }

    return out;
  }

  function packPhi_AQD(phi) {
    const v = [];
    for (const br of phi.branches) v.push(br.z.re, br.z.im);
    for (const br of phi.branches) for (const a of br.A) v.push(a.re, a.im);
    return v;
  }

  function unpackPhi_AQD(v, template) {
    const phi = {
      family: 'boundedAQD',
      w0: template.w0 ? Complex.clone(template.w0) : { re: 0, im: 0 },
      R: template.R,                                      // share by reference
      branches: [],
    };
    let idx = 0;
    for (let j = 0; j < template.branches.length; j++) {
      phi.branches.push({ z: { re: v[idx], im: v[idx + 1] }, A: [] });
      idx += 2;
    }
    for (let j = 0; j < template.branches.length; j++) {
      const mj = template.branches[j].A.length;
      for (let k = 0; k < mj; k++) {
        phi.branches[j].A.push({ re: v[idx], im: v[idx + 1] });
        idx += 2;
      }
    }
    return phi;
  }

  function canonicalizePhi_AQD(phi) {
    let reSum = 0;
    for (const br of phi.branches) {
      if (br.A.length > 0) reSum += br.A[0].re;
    }
    if (reSum >= 0) return phi;
    return {
      family: 'boundedAQD',
      w0: Complex.clone(phi.w0),
      R: phi.R,
      branches: phi.branches.map(br => ({
        z: Complex.neg(br.z),
        A: br.A.map((a, k) => (k % 2 === 0 ? Complex.neg(a) : Complex.clone(a))),
      })),
    };
  }

  // ===========================================================================
  // 5. Initial guess and seeds
  // ===========================================================================
  // Strategy: assume Ω is near a disk and U(z) ≈ R(w₀) + r#'(0)·z. The
  // governing identity reduces to a LINEAR map in the "R-pulled-back" space:
  //   U(z_j) − R(w₀) = r#(z_j) ≈ r#'(0)·z_j = R(a_j) − R(w₀)    (locator)
  //   r#'(0) = R'(w₀) · φ'(0)
  //
  // So choose a single scalar  s := max_j |(R(a_j) − R(w₀))/R'(w₀)| / 0.7,
  // set z_j = (R(a_j) − R(w₀)) / (R'(w₀) · s), and A_{j,k} = P_{j,k}/s^k
  // (the Faber-match value for the linear ψ̃(t) = t/s).
  //
  // This is the AQD analogue of the existing LQD log-bootstrap, which uses
  // log(a_j/w₀) in the same role for R(w) = log w.
  function initialGuess_AQD(hData, norm) {
    const w0 = norm.w0;
    const R = norm.R;
    const Rp0 = Aqd.evalAlpha(R, w0);                     // R'(w₀), known nonzero
    const Rw0 = Aqd.evalR(R, w0);                         // R(w₀)
    let maxAbs = 0;
    const dRs = [];                                       // dR_j := (R(a_j) - R(w₀))/R'(w₀)
    for (const p of hData.poles) {
      const Raj = Aqd.evalR(R, p.a);
      const dR  = Complex.div(Complex.sub(Raj, Rw0), Rp0);
      dRs.push(dR);
      const m = Complex.abs(dR);
      if (m > maxAbs) maxAbs = m;
    }
    const s = Math.max(maxAbs / 0.7, 0.3);                // keep |z_j| < 0.7
    const phi = {
      family: 'boundedAQD',
      w0: Complex.clone(w0),
      R: R,
      branches: [],
    };
    for (let j = 0; j < hData.poles.length; j++) {
      const p = hData.poles[j];
      let z = Complex.scale(dRs[j], 1 / s);
      const zr = Complex.abs(z);
      if (zr > 0.85) z = Complex.scale(z, 0.85 / zr);
      const P = Aqd.principalPartsHOverRprime(p, R);
      const A = [];
      let sk = 1;
      for (let k = 1; k <= p.principal.length; k++) {
        sk *= s;
        A.push(P && P[k - 1] ? Complex.scale(P[k - 1], 1 / sk) : { re: 0, im: 0 });
      }
      phi.branches.push({ z, A });
    }
    return phi;
  }

  function perturbedInitialGuess_AQD(hData, norm, rng, r) {
    const base = initialGuess_AQD(hData, norm);
    QD.LqdCommon.perturbBranchesInPlace(base.branches, rng, r || 0);
    return base;
  }

  function diverseInitialGuess_AQD(hData, norm, rng, r) {
    return {
      family: 'boundedAQD',
      w0: Complex.clone(norm.w0),
      R: norm.R,
      branches: QD.LqdCommon.diverseSeedBranches(hData, rng),
    };
  }

  // Residue-magnitude continuation: ramp t ∈ (0, 1], scaling each C_{j,s}
  // by t and re-running Newton from the previous t's solution. Helps when
  // the direct initial guess sits outside Newton's basin (e.g. multi-pole
  // or higher-order-pole AQDs). Limitation: for some pole structures the
  // t→0 limit isn't a smooth deformation of the true solution; a future
  // Stage 6 pass will add fixed-point preconditioning for those cases.
  function continuationSolve_AQD(hData, norm, opts) {
    opts = opts || {};
    const newtonOpts = opts.newton || {};
    const numSteps = opts.numSteps ?? 12;
    const trace = [];
    // Build a sequence of scaled-down hDatas with t = 1/numSteps, 2/numSteps, ..., 1.
    let phi = null;
    let lastT = 0;
    for (let i = 1; i <= numSteps; i++) {
      const t = i / numSteps;
      const hStep = {
        poles: hData.poles.map(p => ({
          a: { re: p.a.re, im: p.a.im },
          principal: p.principal.map(C => ({ re: C.re * t, im: C.im * t })),
        })),
        polyPart: hData.polyPart,
      };
      let seed;
      if (!phi) {
        seed = initialGuess_AQD(hStep, norm);
      } else {
        // Rescale phi from t_prev to t: A_{j,k} *= t/t_prev (since A ∝ residues at first order).
        const r = t / lastT;
        seed = {
          family: 'boundedAQD',
          w0: phi.w0,
          R: phi.R,
          branches: phi.branches.map(br => ({
            z: { re: br.z.re, im: br.z.im },
            A: br.A.map((a, k) => ({ re: a.re * r, im: a.im * r })),
          })),
        };
      }
      const res = QD.newtonSolve
        ? QD.newtonSolve(seed, hStep, newtonOpts)
        : null;
      if (!res) return { success: false, error: "continuationSolve: QD.newtonSolve unavailable", trace };
      trace.push({ t, success: res.success, residual: res.residual });
      if (!res.success) {
        return { success: false, error: "continuation Newton failed at t=" + t.toFixed(3), trace, residual: res.residual };
      }
      phi = res.phi;
      lastT = t;
    }
    // Final residual against the true hData
    const finalRes = residual_AQD(phi, hData);
    const norm2 = Math.sqrt(finalRes.reduce((s, x) => s + x * x, 0));
    return { success: norm2 < 1e-6, phi, residual: norm2, trace };
  }

  // ===========================================================================
  // 6. Identity verification — Cauchy test functions f(w) = 1/(w − b)^k.
  // ---------------------------------------------------------------------------
  // CONVENTION (consistent with solver-qd.js and the LQD verifiers): the
  // quadrature identity in this codebase reads
  //
  //   ∫_Ω f · ρ dA  =  (1/(2i)) ∮_∂Ω f · h(w) dw     for f ∈ F(Ω).
  //
  // (The factor 1/(2i) means h = 1/w corresponds to the unit disk under ρ=1,
  //  not h = 2/w.) For AQDs with ρ = |R'|², the LHS area integral becomes a
  //  contour integral via Green's theorem. Since |R'|² = R' · conj(R') and
  //
  //    ∂/∂w̅ (R'(w) · conj(R(w))) = R'(w) · conj(R'(w)) = |R'(w)|²,
  //
  //  we have (for f holomorphic in Ω̄):
  //
  //    ∫_Ω f · |R'|² dA = (1/(2i)) ∮_∂Ω f(w) · R'(w) · conj(R(w)) dw.
  //
  // So the IDENTITY we check is:
  //
  //    (1/(2i)) ∮ f · R' · conj(R) dw  =  (1/(2i)) ∮ f · h dw,
  //
  // equivalently
  //
  //    ∮_∂Ω f · R' · conj(R) dw  =  ∮_∂Ω f · h dw.
  //
  // For numerical efficiency we compute:
  //
  //    LHS = (1/N) Σ_n f(w_n) · R'(w_n) · conj(R(w_n)) · z_n · φ'(z_n)
  //    RHS = Σ_{j,s} C_{j,s} · (-1)^{s-1} · binom(k+s-2, s-1) / (a_j − b)^{k+s-1}
  //
  //  (dropping the common factor 2πi / (2i) = π on each side; this matches the
  //   solver-qd.js normalization so all family verifiers report comparable
  //   relative-error scales).
  //
  // Test functions f = 1/(w − b)^k with b chosen outside Ω̄ AND outside the
  // finite poles of R (so f stays bounded and R · f stays analytic on Ω̄).
  // ===========================================================================
  function verifyQuadratureIdentity_AQD(phi, hData, options) {
    options = options || {};
    const N            = options.numSamples ?? 500;
    const maxOrder     = options.maxDegree ?? 3;
    const numTestPts   = options.numTestPoints ?? 3;

    // Boundary samples (w_n, φ'(z_n) via Taylor of order 1).
    const samples = new Array(N);
    for (let n = 0; n < N; n++) {
      const theta = 2 * Math.PI * n / N;
      const z = { re: Math.cos(theta), im: Math.sin(theta) };
      const T = phiTaylorAt_AQD(z, phi, 1);
      samples[n] = { z, w: T[0], phiPrime: T[1] };
    }

    // Test-point selection: prefer points far from the boundary, also away
    // from R's finite poles/zeros (which would break analyticity of f·|α|²).
    // We compute a bounding box of ∂Ω, then place candidates on a generous
    // ring around it, ranking by min-distance-to-polygon.
    let minRe = Infinity, maxRe = -Infinity, minIm = Infinity, maxIm = -Infinity;
    for (const s of samples) {
      if (s.w.re < minRe) minRe = s.w.re;
      if (s.w.re > maxRe) maxRe = s.w.re;
      if (s.w.im < minIm) minIm = s.w.im;
      if (s.w.im > maxIm) maxIm = s.w.im;
    }
    const spanRe = maxRe - minRe, spanIm = maxIm - minIm;
    const span = Math.max(spanRe, spanIm, 1);
    const cx = 0.5 * (minRe + maxRe), cy = 0.5 * (minIm + maxIm);

    const candidates = [];
    const numCandidates = Math.max(numTestPts * 4, 16);
    for (let i = 0; i < numCandidates; i++) {
      const ang = 2 * Math.PI * i / numCandidates;
      candidates.push({
        re: cx + span * 2.0 * Math.cos(ang),
        im: cy + span * 2.0 * Math.sin(ang),
      });
    }
    // Rank by min-dist-to-boundary, take top numTestPts.
    function minDistToPolygon(p, poly) {
      let dmin = Infinity;
      for (const v of poly) {
        const d = Math.hypot(p.re - v.re, p.im - v.im);
        if (d < dmin) dmin = d;
      }
      return dmin;
    }
    const poly = samples.map(s => s.w);
    candidates.sort((a, b) => minDistToPolygon(b, poly) - minDistToPolygon(a, poly));
    const testPoints = candidates.slice(0, numTestPts);

    // Scale-floor for relative-diff: Σ|C_{j,1}|.
    let scaleRef = 0;
    for (const pole of hData.poles) {
      if (pole.principal.length > 0) scaleRef += Complex.abs(pole.principal[0]);
    }
    if (scaleRef === 0) scaleRef = 1;

    const checks = [];
    let maxRelDiff = 0;
    let maxAbsDiff = 0;

    for (let pIdx = 0; pIdx < testPoints.length; pIdx++) {
      const b = testPoints[pIdx];
      for (let k = 1; k <= maxOrder; k++) {
        // LHS (area integral via Green's theorem; common factor π dropped):
        //   LHS = (1/N) Σ_n f(w_n) · R'(w_n) · conj(R(w_n)) · z_n · φ'(z_n)
        let lhs = { re: 0, im: 0 };
        for (let n = 0; n < N; n++) {
          const s = samples[n];
          if (!isFinite(s.w.re) || !isFinite(s.w.im)) continue;
          const diff = Complex.sub(s.w, b);
          if (Complex.abs2(diff) < 1e-30) continue;
          const fVal = Complex.inv(Complex.pow(diff, k));
          const Rp = Aqd.evalAlpha(phi.R, s.w);          // R'(w_n) = α(w_n)
          const Rv = Aqd.evalR(phi.R, s.w);              // R(w_n)
          if (!isFinite(Rp.re) || !isFinite(Rv.re)) continue;
          let term = Complex.mul(fVal, Rp);
          term = Complex.mul(term, Complex.conj(Rv));
          term = Complex.mul(term, s.z);
          term = Complex.mul(term, s.phiPrime);
          lhs = Complex.add(lhs, term);
        }
        lhs = Complex.scale(lhs, 1 / N);

        // RHS = Σ residues of f·h at a_j inside Ω (same π-dropped convention).
        // For f = 1/(w-b)^k and h = Σ_s C_{j,s}/(w-a_j)^s:
        //   res = Σ_s C_{j,s} · (-1)^{s-1} · binom(k+s-2, s-1) / (a_j-b)^{k+s-1}
        let rhs = { re: 0, im: 0 };
        for (const pole of hData.poles) {
          const aMinusB = Complex.sub(pole.a, b);
          for (let sIdx = 0; sIdx < pole.principal.length; sIdx++) {
            const sOrd = sIdx + 1;
            const Cv = pole.principal[sIdx];
            const sign = (sOrd % 2 === 0) ? -1 : 1;       // (-1)^{s-1}
            const coef = binom(k + sOrd - 2, sOrd - 1);
            const denom = Complex.pow(aMinusB, k + sOrd - 1);
            const term = Complex.div(Cv, denom);
            rhs = Complex.add(rhs, Complex.scale(term, sign * coef));
          }
        }

        const dz = Complex.sub(lhs, rhs);
        const absDiff = Complex.abs(dz);
        const scale = Math.max(Complex.abs(lhs), Complex.abs(rhs), scaleRef);
        const relDiff = absDiff / scale;
        if (relDiff > maxRelDiff) maxRelDiff = relDiff;
        if (absDiff > maxAbsDiff) maxAbsDiff = absDiff;
        checks.push({ bIdx: pIdx, k, lhs, rhs, absDiff, relDiff });
      }
    }

    return {
      checks, maxRelDiff, maxAbsDiff,
      areaScale: scaleRef,
      testPoints,
      maxDeg: maxOrder,
      numSamples: N,
      aqd: true,
    };
  }

  function binom(n, k) {
    if (k < 0 || k > n) return 0;
    let r = 1;
    for (let i = 0; i < k; i++) r = r * (n - i) / (i + 1);
    return r;
  }

  // ===========================================================================
  // 7. Register Family.boundedAQD
  // ===========================================================================
  QD.Family.boundedAQD = {
    name: 'boundedAQD',
    enforceInDisk:  true,
    enforceOutDisk: false,

    // Dispatch: any opts bag with `aqd: true` and no `unbounded`.
    matches(opts) { return !!(opts && opts.aqd && !opts.unbounded); },

    normalizeOpts(opts, hData) {
      if (!opts.R) throw new Error("Family.boundedAQD: opts.R (rational primitive) is required");
      let w0 = opts.w0;
      if (!w0) {
        let sumRe = 0, sumIm = 0;
        for (const p of hData.poles) { sumRe += p.a.re; sumIm += p.a.im; }
        const n = hData.poles.length;
        w0 = n > 0 ? { re: sumRe / n, im: sumIm / n } : { re: 1, im: 0 };
      }
      // Defensive checks: w0 must not be at an R-pole or R'-zero (Stage 2).
      const Rp0 = Aqd.evalAlpha(opts.R, w0);
      if (!isFinite(Rp0.re) || !isFinite(Rp0.im) || Complex.abs(Rp0) < 1e-12) {
        throw new Error("Family.boundedAQD: R'(w₀) is zero or undefined; pick a different φ(0)");
      }
      return { aqd: true, w0, R: opts.R };
    },

    evalPhi: evalPhi_AQD,
    phiTaylorAt: phiTaylorAt_AQD,

    computeTargets(phi, hData) {
      return { A: computeTargetA_AQD(phi, hData), F: null };
    },

    residual: residual_AQD,
    packPhi:  packPhi_AQD,
    unpackPhi: unpackPhi_AQD,
    canonicalizePhi: canonicalizePhi_AQD,

    initialGuess: initialGuess_AQD,
    perturbedInitialGuess: perturbedInitialGuess_AQD,
    diverseInitialGuess: diverseInitialGuess_AQD,
    continuationSolve: continuationSolve_AQD,
    verifyQuadratureIdentity: verifyQuadratureIdentity_AQD,
  };

  if (QD.registerFamily) {
    QD.registerFamily('boundedAQD');
  } else {
    throw new Error("aqd-bounded.js: QD.registerFamily not found");
  }

}());
