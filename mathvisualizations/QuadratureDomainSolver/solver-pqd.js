// =============================================================================
// solver-pqd.js -- Bounded power-weighted quadrature domains (Family.powerQD)
//
// A bounded PQD of order α satisfies
//     ∫_Ω f(w) |w|^{2(α-1)} dA = ∮_∂Ω f(w) h(w) dw,   f ∈ A(Ω) ∩ L¹(Ω; |w|^{2(α-1)})
// with h rational (poles in Ω — the quadrature nodes at finite a_j). This file
// handles ANY real α > 0 with α ≠ 1 (the αth-root power is taken via
// Complex.cpow, so non-integer α is fully supported — QA milestone). α = 1
// recovers classical bounded QDs and routes to Family.boundedQD instead.
//
// Riemann-map characterization (Theorem 4.3.2 / Equation 4.8 of Graven 2026):
//     Ω ∈ QD^α(h) ⇔ φ_out(ξ) extends as the (1/α)th power of a rational.
//
// Parametrization used here (bounded, 0 ∈ Ω so φ_in ≡ 1 per Equation 4.9):
//     φ(ξ) = (R#(ξ))^{1/α}        (αth-root on the branch ANCHORED at φ(0)=w0;
//                                  NOT principal — see QD.PqdCommon.phiAnchored.
//                                  The principal branch fails for off-axis poles
//                                  |arg a| > π/α, where it picks the wrong sheet.)
//     R#(ξ) = r0 + Σ_j Σ_k conj(A_{j,k}) ξ^k / (1 − conj(z_j) ξ)^k
//     r0 = φ(0)^α = w0^α
//
// So the Newton variables {z_j, A_{j,k}} are IDENTICAL to boundedQD; w0 is
// configured (user-supplied), and r0 = w0^α is the constant term of R#.
//
// Locator (●):  φ(z_j) = a_j  ⇔  R#(z_j) = a_j^α  (raised to avoid branch cut)
// Faber-target (★): A_{j,k} = Σ_{s ≥ k} (s/k) · D_{j,s} · [t^s] ψ̃_j^k(t)
//   where ψ̃_j is the local inverse of φ(z_j + ·) − a_j and the modified
//   residues D_{j,s} encode the α-power of the lifted parametrization
//   (see Corollary 4.3.1; for α=1 these reduce to C_{j,s}).
//
// Gauge: Σ_j Im(A_{j,1}) = 0 + canonicalize Re sign (the anchored αth root —
//   anchorArg0 = α·arg(w0) — selects the global sheet; the canonicalizer pins
//   the remaining discrete Z/α rotation).
//
// Companion files:
//   * solvers/seeds/seeds-pqd.js — multistart / continuation seed strategy
//   * solver-qd.js (α=1 catch-all; this file's matches() requires α>0, α≠1)
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family) {
    throw new Error("solver-pqd.js: solver.js must be loaded first");
  }
  if (!QD.Seeds || !QD.Seeds.powerQD) {
    throw new Error("solver-pqd.js: QD.Seeds.powerQD missing — seeds-pqd.js must be loaded first");
  }

  const diskInitialGuess_PQD      = QD.Seeds.powerQD.initialGuess;
  const perturbedInitialGuess_PQD = QD.Seeds.powerQD.perturbedInitialGuess;

  // ===========================================================================
  // Helpers — complex real-power (arbitrary α > 0) and principal αth-root.
  // Both use the principal branch (arg ∈ (-π, π]). cpowA(c, α) = c^α and
  // cprincipalRoot(c, α) = c^{1/α} are inverse operations; both delegate to
  // the shared Complex.cpow primitive (complex base, real exponent).
  // ===========================================================================
  function cpowA(c, alpha) {
    return Complex.cpow(c, alpha);
  }
  // Principal αth root: |c|^{1/α} · exp(i · arg(c) / α), arg ∈ (-π, π].
  function cprincipalRoot(c, alpha) {
    return Complex.cpow(c, 1 / alpha);
  }

  // ===========================================================================
  // 1. R# evaluation (shares structure with boundedQD's evalPhi, but with
  //    r0 = w0^α as constant term).
  // ===========================================================================
  function evalRHash_PQD(z, phi) {
    const r0 = cpowA(phi.w0, phi.alpha);
    let result = Complex.clone(r0);
    for (const br of phi.branches) {
      const zjC = Complex.conj(br.z);
      const denom = Complex.sub(Complex.ONE(), Complex.mul(zjC, z));
      const u = Complex.div(z, denom);
      let uPow = Complex.ONE();
      for (const Ak of br.A) {
        uPow = Complex.mul(uPow, u);
        result = Complex.add(result, Complex.mul(Complex.conj(Ak), uPow));
      }
    }
    return result;
  }

  // φ(z) = (R#(z))^{1/α} on the branch ANCHORED at φ(0)=w0. R# is non-vanishing
  // on 𝔻̄ for a valid PQD, so a single continuous branch exists; we pin it by
  // continuing arg(R#) from argCont(R#(0)) = α·arg(w0) (which makes the root at
  // z=0 equal w0 exactly) along the segment 0→z. This is what makes off-axis
  // poles (|arg a| > π/α) work — the principal root would land on a wrong sheet.
  function evalPhi_PQD(z, phi) {
    return QD.PqdCommon.phiAnchored(
      phi, z, evalRHash_PQD, phi.alpha * Complex.arg(phi.w0), { re: 0, im: 0 });
  }

  // ===========================================================================
  // 2. Taylor of R# at z = z_0 up to order L (identical structure to
  //    boundedQD's phiTaylorAt, plus the r0 constant term).
  // ===========================================================================
  function rHashTaylorAt_PQD(z0, phi, L) {
    const r0 = cpowA(phi.w0, phi.alpha);
    const result = Taylor.zero(L + 1);
    result[0] = Complex.clone(r0);

    for (const br of phi.branches) {
      const zjC = Complex.conj(br.z);
      const alpha_z = Complex.sub(Complex.ONE(), Complex.mul(zjC, z0));
      const alphaInv = Complex.inv(alpha_z);

      const uT = Taylor.zero(L + 1);
      uT[0] = Complex.mul(z0, alphaInv);
      if (L >= 1) {
        let zjcPow = { re: 1, im: 0 };
        let alphaInvPow = Complex.mul(alphaInv, alphaInv);
        for (let l = 1; l <= L; l++) {
          uT[l] = Complex.mul(zjcPow, alphaInvPow);
          zjcPow = Complex.mul(zjcPow, zjC);
          alphaInvPow = Complex.mul(alphaInvPow, alphaInv);
        }
      }

      let uPow = Taylor.truncate(uT, L);
      for (let k = 1; k <= br.A.length; k++) {
        const AkC = Complex.conj(br.A[k - 1]);
        for (let i = 0; i <= L; i++) {
          result[i] = Complex.add(result[i], Complex.mul(AkC, uPow[i]));
        }
        if (k < br.A.length) uPow = Taylor.mul(uPow, uT, L);
      }
    }
    return result;
  }

  // Taylor of φ(z) = (R#(z))^{1/α} at z = z_0 up to order L.
  // Use φ = exp((1/α) · log(R#)). Requires R#(z_0) ≠ 0 (the principal
  // branch is locally analytic away from zeros of R#).
  //
  // When `anchorArg` is supplied it is argCont(R#(z_0)) (the continuous arg on
  // the anchored sheet); we override the constant term of log(R#) with it. This
  // is essential, not cosmetic: Taylor.exp scales the ENTIRE series (slope and
  // all higher coeffs) by e^{(1/α)·logRT[0]}, and inverseFaberAtPole consumes
  // those coefficients — so a wrong-sheet constant rotates the (★) target and
  // Newton cannot converge for off-axis poles. Without `anchorArg` it picks the
  // principal-branch atan2 at z_0 (legacy single-point behaviour).
  function phiTaylorAt_PQD(z0, phi, L, anchorArg) {
    const rT = rHashTaylorAt_PQD(z0, phi, L);
    const logRT = Taylor.log(rT, L);                  // principal-branch log
    if (anchorArg !== undefined) logRT[0] = { re: logRT[0].re, im: anchorArg };
    const scaledLog = Taylor.scaleComplex(logRT, { re: 1 / phi.alpha, im: 0 });
    return Taylor.exp(scaledLog, L);
  }

  // ===========================================================================
  // 2b. Continuous-branch unit-circle sweep — Q1.3.
  // ---------------------------------------------------------------------------
  // R# ∈ A(𝔻̄) is non-vanishing on the closed disk for a valid bounded PQD
  // (since φ = (R#)^{1/α} is the Riemann map and Ω ∌ 0). So arg(R#(e^{iθ}))
  // is a continuous, single-valued function of θ that returns to its starting
  // value after one loop — by the argument principle, winding(R#|∂𝔻 around 0)
  // equals the number of zeros of R# inside 𝔻, which is 0 for a valid PQD.
  // Hence the αth root closes up after one sweep.
  //
  // A naive per-sample principal evaluation (Math.atan2, range cut at the
  // negative real axis) would jump by 2π/α whenever R#(e^{iθ}) crosses that ray
  // (without winding around 0), even though the true arg is continuous. The fix
  // is to unwrap arg(R#) along the θ-sweep AND anchor θ=0 to the interior sheet
  // (anchorSpec below: a walk 0→z=1 from argCont(R#(0))=α·arg(w0)) so the
  // boundary lands on the SAME global branch as evalPhi_PQD / the (★) target —
  // essential for off-axis poles, where the principal sheet is wrong.
  //
  // VALIDITY GUARD: if R# DOES wind around 0 (Newton converged to a spurious
  // fixed point of the (●)+(★)+gauge system that has zeros of R# in 𝔻),
  // sweepUnitCircle_PQD's contArg will accumulate a net non-zero shift over
  // [0, 2π] — meaning φ wouldn't close up. The identity verifier built on
  // top of these samples will then report a huge (~O(1)) maxRelDiff, which
  // causes solveInverseQD's identityOK gate to reject the candidate and
  // multistart continues searching for a valid root. So the spurious
  // fixed points are detected naturally without an extra winding check.
  //
  // sweepUnitCircle_PQD(phi, N) returns N uniformly-spaced samples along ∂𝔻
  // with:
  //   { theta, z, rH, contArg, w, phiPrime }
  // where contArg is the unwrapped arg(R#), w = |rH|^{1/α}·exp(i·contArg/α)
  // is the continuous Riemann map, and phiPrime = R#'(z)/(α·w^{α-1}) is the
  // derivative via the chain rule (with w taken on the same sheet).
  // ===========================================================================
  // Continuous-arg ∂𝔻 sweep via the shared PqdCommon driver. powerQD's φ has
  // NO prefactor (φ = (R#)^{1/α}, so w === root) and φ' = R#'/(α·w^{α-1}),
  // computed on the same continuous sheet the driver selected (wMag^{α-1}·
  // e^{i(α-1)·rootArg}) — arbitrary-real-α safe, no integer-power loop.
  function combine_PQD(ctx) {
    const { rHp, root, rootMag, rootArg, alpha } = ctx;
    const wAm1Mag = Math.pow(rootMag, alpha - 1);
    const wAm1Ang = (alpha - 1) * rootArg;
    const wPowAm1 = { re: wAm1Mag * Math.cos(wAm1Ang), im: wAm1Mag * Math.sin(wAm1Ang) };
    const phiPrime = Complex.div(rHp, Complex.scale(wPowAm1, alpha));
    return { w: root, phiPrime };
  }
  function sweepUnitCircle_PQD(phi, N) {
    // Anchor the sweep's θ=0 sheet to the interior φ(0)=w0 normalization (walk
    // 0→z=1), so the whole boundary is on the single global branch — correct
    // for off-axis poles.
    const anchorSpec = {
      evalRHashFn: evalRHash_PQD,
      anchorPt: { re: 0, im: 0 },
      anchorArg0: phi.alpha * Complex.arg(phi.w0),
    };
    return QD.PqdCommon.sweepUnitCircle(phi, N, rHashTaylorAt_PQD, combine_PQD, anchorSpec);
  }

  // Family.X.sampleBoundary hook (dispatched by sampleBoundaryAdaptive in
  // solver.js for phi.family === 'powerQD'). Continuous-arg sweep + curvature-
  // aware deviation refinement; prefactor = 1 (w = root). See PqdCommon.
  function sampleBoundary_PQD(phi, baseSamples, maxExtra) {
    return QD.PqdCommon.sampleBoundaryViaSweep(
      phi, baseSamples, maxExtra, sweepUnitCircle_PQD,
      (thMid, leftPt, ph) => QD.PqdCommon.boundaryMid(thMid, leftPt, ph, evalRHash_PQD, null),
      0);   // powerQD budget = exactly maxExtra (no N0×fallback)
  }

  // ===========================================================================
  // 3. Target A_{j,k} from (★)  — arbitrary pole order
  //
  // For bounded PQDs (0 ∉ Ω), the explicit closed form for the rational
  // r(z) is (Graven, "Analysis of Log-Weighted Quadrature Domains",
  // §1.1.1 derivation, corrected version):
  //
  //   β_{j,k} = α · (k−1)! · Σ_{n=k}^{n_j} D_{j,n} · B_{n,k}(ψ'(p_j),…)/(n−1)!
  //
  //   D_{j,n} := Σ_{m=0}^{n_j−n} C(1−α, m) · p_j^{1−α−m} · C_{j,n+m}
  //
  // where ψ = φ⁻¹, B_{n,k} is the incomplete Bell polynomial, and
  // C(1−α, m) is the generalized binomial coefficient (equals
  // (1−α)(−α)(−α−1)…(2−α−m)/m! ; for integer α ≥ 1 this is well-defined
  // and non-zero for all m ≥ 0).
  //
  // INSIGHT (and the reason we don't have to implement Bell polynomials
  // by hand): the inner sum
  //
  //   (k−1)! · Σ_{n≥k} D_{j,n} · B_{n,k}(…)/(n−1)!
  //
  // is exactly the order-k residue at z_j of Φ_φ⁻¹[Σ_n D_{j,n}/(w−p_j)^n].
  // The existing primitive `QD.Faber.inverseFaberAtPole(residues, phiTilde)`
  // computes the same quantity via Lagrange-Bürmann inversion of the local
  // Taylor of φ at z_j — equivalent to evaluating Bell polynomials on
  // ψ^{(i)}(p_j) via Faà di Bruno's formula. So the implementation reduces
  // to:
  //
  //   D = modifiedResidues(pole, α);
  //   phiTilde = Taylor of (φ(z_j+t) − p_j) with constant zeroed;
  //   A[k] = α · inverseFaberAtPole(D, phiTilde)[k];
  //
  // SANITY CHECKS:
  //   • α=1: C(0, m) = δ_{m,0} ⇒ D_{j,n} = C_{j,n} · p^0 = C_{j,n}, and
  //     the α=1 prefactor collapses, so the result equals the existing
  //     boundedQD (★) bit-identically. (Family.boundedQD continues to
  //     own α=1 via the dispatcher; this code path is only entered for
  //     α ≥ 2.)
  //   • α=2, simple pole: D_{j,1} = C_{j,1}/p_j and β = α·D/φ'(z_j) =
  //     α·C·p^{1−α}/φ'(z_j) — matches the simple-pole formula and the
  //     existing test (identity verifier 3.86e-15).
  //
  // HIGHER-ORDER POLES: now handled. The realizability constraint for
  // m_j ≥ 2 is more involved than the simple-pole case; the multistart
  // pipeline finds the valid root when one exists.
  // ===========================================================================

  // Generalized binomial coefficient C(x, m) for real x, non-negative integer m.
  //   C(x, 0) = 1
  //   C(x, m) = x(x−1)…(x−m+1)/m!
  function generalizedBinom(x, m) {
    if (m === 0) return 1;
    let r = 1;
    for (let i = 0; i < m; i++) r *= (x - i) / (i + 1);
    return r;
  }

  // Modified residues D_{j,n} for the principal part of h(w)·w^{1−α} at p_j.
  //   D_{j,n} = Σ_{m=0}^{n_j−n} C(1−α, m) · p_j^{1−α−m} · C_{j,n+m}
  function modifiedResidues_PQD(pole, alpha) {
    const C = pole.principal;
    const mj = C.length;
    const p = pole.a;
    const D = new Array(mj);
    // Pre-compute p^{1−α−m} for m = 0..(mj−1); start at the real-power
    // p^{1−α} (valid for arbitrary α via Complex.cpow) and multiply by
    // p^{−1} each step.
    const pInv = Complex.inv(p);
    let pPow = Complex.cpow(p, 1 - alpha);                   // p^{1−α}
    const pPowList = new Array(mj);
    pPowList[0] = pPow;
    for (let i = 1; i < mj; i++) pPowList[i] = Complex.mul(pPowList[i - 1], pInv);

    for (let n = 1; n <= mj; n++) {
      let acc = { re: 0, im: 0 };
      const mMax = mj - n;
      for (let m = 0; m <= mMax; m++) {
        const binom = generalizedBinom(1 - alpha, m);
        if (binom === 0) continue;
        const term = Complex.scale(Complex.mul(pPowList[m], C[n + m - 1]), binom);
        acc = Complex.add(acc, term);
      }
      D[n - 1] = acc;
    }
    return D;
  }

  function computeTargetA_PQD(phi, hData) {
    const alpha = phi.alpha;
    const target = [];
    for (let j = 0; j < hData.poles.length; j++) {
      const pole = hData.poles[j];
      const mj = pole.principal.length;
      const zj = phi.branches[j].z;

      // Taylor of φ at z_j to order m_j. phiTilde[0] = 0 (caller convention),
      // phiTilde[i ≥ 1] = i-th Taylor coefficient of φ at z_j. The Taylor must
      // be built on the ANCHORED sheet (argCont(R#(z_j)) continued from
      // φ(0)=w0), else the (★) target lands on the wrong αth-root sheet for
      // off-axis poles and Newton fails to converge.
      const argZj = QD.PqdCommon.argContAt(
        phi, zj, evalRHash_PQD, alpha * Complex.arg(phi.w0), { re: 0, im: 0 });
      const phiT = phiTaylorAt_PQD(zj, phi, mj, argZj);
      const phiTilde = Taylor.zero(mj + 1);
      for (let i = 1; i <= mj; i++) phiTilde[i] = Complex.clone(phiT[i]);

      // Modified residues D_{j,n} pull in the w^{1−α} factor from Eq. 1.7.
      const D = modifiedResidues_PQD(pole, alpha);

      // β_{j,k}/α = (residue of Φ_φ⁻¹[Σ D_n/(w−p_j)^n] at z_j of order k).
      const A_unscaled = QD.Faber.inverseFaberAtPole(D, phiTilde);
      const A = A_unscaled.map(a => Complex.scale(a, alpha));
      target.push(A);
    }
    return target;
  }

  // ===========================================================================
  // 4. Residual, pack, unpack, canonicalize
  //
  // We enforce the locator in the α-lifted form
  //   R#(z_j) − a_j^α = 0
  // (polynomial, no branch cut). The (★) target equation reconstructs φ via
  // phiTaylorAt_PQD on the ANCHORED αth-root branch (argCont(R#(z_j)) continued
  // from φ(0)=w0 — see computeTargetA_PQD), so the target lands on the correct
  // sheet even for off-axis poles; the canonicalizer pins the global Z/α phase
  // ambiguity.
  // ===========================================================================
  function residual_PQD(phi, hData, options = {}) {
    const enforceGauge = options.enforceGauge !== false;
    const out = [];

    // (●) R#(z_j) = a_j^α   ⇔   (R#(z_j) − a_j^α) = 0
    for (let j = 0; j < hData.poles.length; j++) {
      const rZj = evalRHash_PQD(phi.branches[j].z, phi);
      const aLifted = cpowA(hData.poles[j].a, phi.alpha);
      const diff = Complex.sub(rZj, aLifted);
      out.push(diff.re, diff.im);
    }
    // (★) A_{j,k} = target
    const target = computeTargetA_PQD(phi, hData);
    for (let j = 0; j < hData.poles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], target[j][k]);
        out.push(diff.re, diff.im);
      }
    }
    // Gauge: Σ Im(A_{j,1}) = 0   (same as boundedQD; pins the disk rotation
    // before the αth-root branch is selected by canonicalizePhi).
    if (enforceGauge) {
      let imSum = 0;
      for (const br of phi.branches) if (br.A.length > 0) imSum += br.A[0].im;
      out.push(imSum);
    }
    return out;
  }

  // Schema-driven pack/unpack: same structure as boundedQD but carrying
  // the configured phi.alpha through to unpacked instances.
  const SCHEMA_PQD = [
    { kind: 'branchesZ', clamp: { side: 'in', cap: 0.9999 } },
    { kind: 'branchesA' },
  ];

  function packPhi_PQD(phi) {
    return QD.packPhiBySchema(phi, SCHEMA_PQD);
  }

  function unpackPhi_PQD(v, template) {
    return QD.unpackPhiBySchema(v, template, SCHEMA_PQD, (phi) => {
      phi.family = 'powerQD';
      phi.alpha  = template.alpha;
      phi.unbounded = false;
    });
  }

  // Canonicalize: the αth root is α-fold ambiguous (e^{2πik/α} factors on φ).
  // We pin the global rotation by the same Σ Re(A_{j,1}) > 0 convention as
  // boundedQD (180° rotation z → −z flips the sign of odd-indexed A's). This
  // resolves the Z/2 part; the global αth-root sheet is fixed by the anchored
  // branch (anchorArg0 = α·arg(w0) in QD.PqdCommon.phiAnchored), so φ(0)=w0.
  function canonicalizePhi_PQD(phi) {
    let reSum = 0;
    for (const br of phi.branches) if (br.A.length > 0) reSum += br.A[0].re;
    if (reSum >= 0) return phi;
    return {
      family: 'powerQD',
      alpha: phi.alpha,
      unbounded: false,
      w0: Complex.clone(phi.w0),
      branches: phi.branches.map(br => ({
        z: Complex.neg(br.z),
        A: br.A.map((a, k) => (k % 2 === 0 ? Complex.neg(a) : Complex.clone(a))),
      })),
    };
  }

  // ===========================================================================
  // 5. Continuation — homotopy in RESIDUE STRENGTH s ∈ (0, 1].
  //    Ramp the C_{j,k} residues from a small fraction to full size via
  //    QD.scaleHDataResidues. At small s, R# ≈ w₀^α (a near-disk) with a
  //    near-exact seed and z_j stable, so Newton warm-starts smoothly as s→1.
  //    This is a genuine deformation even for a single pole sitting at w₀
  //    (where the old pole-location homotopy via scaleHDataPoles was vacuous).
  // ===========================================================================
  function continuationSolve_PQD(hData, norm, options = {}) {
    const {
      tStart       = 0.10,
      tStartMin    = 1e-3,
      growFactor   = 1.6,
      shrinkFactor = 0.5,
      minStep      = 5e-4,
      maxSteps     = 80,
      newton       = {},
    } = options;
    const w0 = norm.w0;
    const alpha = norm.alpha;

    const trace = [];
    let s = tStart;
    let phi = diskInitialGuess_PQD(QD.scaleHDataResidues(hData, s), w0, alpha);

    let warmupResult;
    while (true) {
      warmupResult = QD.newtonSolve(phi, QD.scaleHDataResidues(hData, s), newton);
      if (warmupResult.success) { phi = warmupResult.phi; break; }
      s *= shrinkFactor;
      if (s < tStartMin) {
        return {
          success: false,
          error: "continuation: warmup failed even at s=" + s.toExponential(2),
          phi: warmupResult.phi,
          trace,
        };
      }
      phi = diskInitialGuess_PQD(QD.scaleHDataResidues(hData, s), w0, alpha);
    }
    trace.push({ t: s, ok: true, residual: warmupResult.residual });

    let lastSuccessS = s;
    let stepSize = Math.max(s, 0.1);
    for (let step = 0; step < maxSteps; step++) {
      if (lastSuccessS >= 1.0 - 1e-12) break;
      const nextS = Math.min(1.0, lastSuccessS + stepSize);
      const hData_s = QD.scaleHDataResidues(hData, nextS);
      const result = QD.newtonSolve(phi, hData_s, newton);
      if (result.success) {
        phi = result.phi;
        lastSuccessS = nextS;
        trace.push({ t: nextS, ok: true, residual: result.residual });
        stepSize *= growFactor;
      } else {
        stepSize *= shrinkFactor;
        trace.push({ t: nextS, ok: false, residual: result.residual ?? null });
        if (stepSize < minStep) {
          return {
            success: false,
            error: "continuation: step size underflow at s=" + lastSuccessS.toFixed(4),
            phi, trace, lastT: lastSuccessS,
          };
        }
      }
    }
    if (lastSuccessS < 1.0 - 1e-9) {
      return {
        success: false,
        error: "continuation: max steps reached at s=" + lastSuccessS.toFixed(4),
        phi, trace, lastT: lastSuccessS,
      };
    }
    return {
      success: true, phi, iterations: 0,
      residual: trace[trace.length - 1].residual,
      trace, method: "continuation",
    };
  }

  // ===========================================================================
  // 5b. w₀ selection — bootstrap from the classical (α=1) QD.
  //   w₀ = φ(0) is a free gauge (any interior point of Ω works) but MUST land
  //   inside the unknown domain, which the user can't easily guess for a far
  //   pole. The classical bounded-QD solver is robust, so we use it as a
  //   validity oracle: try a few candidate w₀ (pole centroid, then the nodes
  //   themselves — a node is always interior AND realizable, bound = 0), and
  //   return the first for which the classical solve yields a univalent map.
  //   Gated by opts.bootstrapW0 (default on) so per-pixel param-slice sweeps,
  //   which re-solve millions of times, can skip the nested classical solves.
  // ===========================================================================
  function poleCentroid_PQD(hData) {
    let sr = 0, si = 0;
    const n = hData.poles.length;
    for (const p of hData.poles) { sr += p.a.re; si += p.a.im; }
    return n > 0 ? { re: sr / n, im: si / n } : { re: 0, im: 0 };
  }

  function bootstrapW0_PQD(hData, opts) {
    const centroid = poleCentroid_PQD(hData);
    if (!hData.poles || hData.poles.length === 0) return centroid;
    if (opts && opts.bootstrapW0 === false) return centroid;

    // Candidate w₀'s: centroid first, then nodes ordered by descending
    // |residue| (the dominant node is the most robustly-interior point).
    const nodes = hData.poles
      .map(p => ({ a: p.a, w: p.principal.length ? Complex.abs(p.principal[0]) : 0 }))
      .sort((x, y) => y.w - x.w)
      .map(x => x.a);
    const cand = [centroid, ...nodes];

    let fallback = centroid;
    for (const w0c of cand) {
      if (Complex.abs2(w0c) < QD.ZERO_THRESHOLD) continue;   // 0 ∉ Ω required
      if (Complex.abs2(fallback) < QD.ZERO_THRESHOLD) fallback = w0c;
      try {
        const res = QD.solveInverseQD(hData, {
          w0: w0c, findAlternates: false, identityTol: 1e-4,
        });
        if (res && res.success && res.primary && res.primary.univalent &&
            res.primary.phi && res.primary.phi.w0) {
          return Complex.clone(res.primary.phi.w0);
        }
      } catch (e) { /* try next candidate */ }
    }
    return fallback;
  }

  // ===========================================================================
  // 6. Identity verification — power-weighted monomials w^k.
  //
  // ∫_Ω w^k · |w|^{2(α-1)} dA = ∮_∂Ω w^k h(w) dw   for k ≥ 0
  //
  // LHS via Green's theorem with weight |w|^{2(α-1)} = (w·conj(w))^{α-1}:
  //   ∫_Ω w^k · |w|^{2(α-1)} dA
  //     = (1/(2i)) ∮_∂Ω w^k · |w|^{2(α-1)} · conj(w)/α · dw
  //     = (1/α) · (1/(2i)) ∮_∂Ω w^k · (conj(w))^α · w^{α-1} · dw
  // (using d/dw(|w|^{2α}) = α · |w|^{2(α-1)} · conj(w) for the antiderivative).
  //
  // Numerically on z_n = e^{iθ_n}, w_n = φ(z_n):
  //   LHS = (1/(αN)) Σ_n w_n^k · conj(w_n)^α · w_n^{α-1} · φ'(z_n) · z_n
  //
  // RHS via residue calculus at the poles of h:
  //   RHS = Σ_j Σ_s C_{j,s} · binom(k, s-1) · a_j^{k-s+1}
  // (identical to the classical QD case — the weight only enters LHS).
  // ===========================================================================
  function verifyQuadratureIdentity_PQD(phi, hData, options = {}) {
    const N = options.numSamples ?? 500;
    const totalDeg = hData.poles.reduce((s, p) => s + p.principal.length, 0);
    const K = options.maxDegree ?? Math.max(totalDeg, 4);
    const alpha = phi.alpha;

    // Use continuous-arg sweep (Q1.3). Previously each sample called
    // phiTaylorAt_PQD which picks principal-atan2 per-point; the LHS
    // integral computed with such samples picks up spurious O(1)
    // contributions from the artificial branch cut. The swept version
    // tracks unwrapped arg(R#) so w_n and φ'(z_n) are continuous along
    // the contour.
    const samples = sweepUnitCircle_PQD(phi, N);

    // Q1.4 — explicit R# non-vanishing guard. φ = (R#)^{1/α} is a single-valued
    // univalent map only if R# is non-vanishing on 𝔻̄ (winding number 0 about 0;
    // by the argument principle that equals the count of R# zeros inside 𝔻). A
    // spurious Newton root can land on a state whose R# has a zero inside 𝔻
    // (winding ≠ 0) or touching ∂𝔻 (|R#| → 0 on a sample) — there the αth root
    // is multi-valued and the sampled "boundary" is not a valid PQD boundary.
    // Detect it directly and force the identity to fail so the candidate is
    // rejected (previously this was caught only indirectly, via the net phase
    // shift corrupting the LHS integral). `rHashVanishes` is surfaced for
    // diagnostics. Cheap: one extra pass over the existing samples' R# values.
    const rHashGuard = QD.PqdCommon.rHashVanishingGuard(samples);

    let areaScale = 0;
    for (const pole of hData.poles) {
      if (pole.principal.length > 0) areaScale += Complex.abs(pole.principal[0]);
    }
    if (areaScale === 0) areaScale = 1;

    const checks = [];
    let maxRelDiff = 0;
    let maxAbsDiff = 0;

    for (let k = 0; k <= K; k++) {
      // LHS: (1/(αN)) Σ_n w_n^k · |w_n|^{2(α-1)} · conj(w_n) · φ'(z_n) · z_n,
      // accumulated by PqdCommon (single-valued (|w|²)^{α-1} weight — see there).
      // Test monomial f = w^k (analytic in the bounded Ω). No near-zero skip:
      // 0 ∉ Ω̄ for the non-singular family, so |w| > 0 on ∂Ω.
      let lhs = QD.PqdCommon.accumulateWeightedLHS(samples, alpha, (w) => Complex.pow(w, k), 0);
      lhs = Complex.scale(lhs, 1 / (alpha * N));

      // RHS: Σ_j Σ_s C_{j,s} · binom(k, s-1) · a_j^{k-s+1}  (shared bounded form).
      const rhs = QD.PqdCommon.boundedMonomialRHS(hData, k);

      const diff = Complex.sub(lhs, rhs);
      const absDiff = Complex.abs(diff);
      const scale = Math.max(Complex.abs(lhs), Complex.abs(rhs), areaScale);
      const relDiff = absDiff / scale;
      if (relDiff > maxRelDiff) maxRelDiff = relDiff;
      if (absDiff > maxAbsDiff) maxAbsDiff = absDiff;
      checks.push({ k, lhs, rhs, absDiff, relDiff });
    }

    return {
      checks,
      maxRelDiff: rHashGuard.vanishes ? Infinity : maxRelDiff,
      maxAbsDiff,
      areaScale, maxDeg: K, numSamples: N,
      alpha,
      rHashVanishes: rHashGuard.vanishes,
      rHashWinding: rHashGuard.winding,
    };
  }

  // ===========================================================================
  // 7. Register Family.powerQD
  //
  // matches(): selected when opts.alpha is any real > 0, α ≠ 1, AND not
  // unbounded/lqd. α = 1 routes to the existing boundedQD catch-all
  // (bit-identical classical behavior). Non-integer α (incl. the (0,1)
  // LQD-limit regime) is fully supported (QA milestone) — the (★) closed
  // form is α-general and all power ops go through Complex.cpow.
  // ===========================================================================
  QD.Family.powerQD = {
    name: 'powerQD',
    enforceInDisk:  true,
    enforceOutDisk: false,
    matches(opts) {
      const a = opts && opts.alpha;
      return Number.isFinite(a) && a > 0 && a !== 1 && !opts.unbounded && !opts.lqd;
    },

    normalizeOpts(opts, hData) {
      // Accept arbitrary real α > 0 (α ≠ 1). No rounding — non-integer α is
      // a first-class case. α = 1 is classical bounded QD (use that mode);
      // α ≤ 0 is out of scope.
      const alpha = opts.alpha;
      if (!(alpha > 0) || alpha === 1) {
        throw new Error("Family.powerQD: α must be a real number > 0 with α ≠ 1 (α = 1 is classical bounded QD)");
      }
      const userSuppliedW0 = !!opts.w0;
      // w₀ default: bootstrap a guaranteed-interior point from the classical
      // (α=1) QD (see bootstrapW0_PQD). Honor an explicit user w₀ as-is.
      let w0 = userSuppliedW0 ? opts.w0 : bootstrapW0_PQD(hData, opts);
      if (Complex.abs2(w0) < QD.ZERO_THRESHOLD) {
        // 0 ∈ Ω̄ degenerate — the singular case φ = b_{z₀}·(R#)^{1/α} is the
        // separate Family.powerQD_singular (Phase QB), not handled here.
        throw new Error("Family.powerQD: w₀ = φ(0) must be nonzero (0 ∉ Ω; singular 0 ∈ Ω is a separate family)");
      }
      // Realizability guard (single simple pole, user-supplied w₀): a
      // univalent bounded PQD with φ(0)=w₀ exists only if C > (|a^α−w₀^α|/α)².
      // Catch the impossible-domain case early with an actionable message.
      if (userSuppliedW0 && hData.poles.length === 1 && hData.poles[0].principal.length === 1) {
        const p = hData.poles[0];
        const C = Complex.abs(p.principal[0]);
        const diff = Complex.abs(Complex.sub(cpowA(p.a, alpha), cpowA(w0, alpha)));
        const bound = (diff * diff) / (alpha * alpha);
        if (C <= bound) {
          throw new Error("Family.powerQD: no bounded PQD with φ(0)=w₀ for this h "
            + "(realizability: C=" + C.toFixed(4) + " ≤ (|a^α−w₀^α|/α)²=" + bound.toFixed(4)
            + "); choose w₀ nearer the quadrature node a, or omit w₀ to auto-select.");
        }
      }
      return { w0, alpha };
    },

    evalPhi: evalPhi_PQD,
    phiTaylorAt: phiTaylorAt_PQD,
    computeTargets(phi, hData) {
      return { A: computeTargetA_PQD(phi, hData), F: null };
    },
    residual: residual_PQD,
    packPhi: packPhi_PQD,
    unpackPhi: unpackPhi_PQD,
    canonicalizePhi: canonicalizePhi_PQD,
    initialGuess(hData, norm) { return diskInitialGuess_PQD(hData, norm.w0, norm.alpha); },
    perturbedInitialGuess(hData, norm, rng, r) {
      return perturbedInitialGuess_PQD(hData, norm.w0, norm.alpha, rng, r);
    },
    diverseInitialGuess(hData, norm, rng, r) {
      // Reuse the generic diverse strategy, then re-stamp the family + alpha tags.
      const base = QD.diverseInitialGuess(hData, norm, rng, r);
      base.family = 'powerQD';
      base.alpha  = norm.alpha;
      return base;
    },
    continuationSolve(hData, norm, opts) {
      return continuationSolve_PQD(hData, norm, opts);
    },
    verifyQuadratureIdentity: verifyQuadratureIdentity_PQD,
    // Q1.3: optional family hook. sampleBoundaryAdaptive (solver.js)
    // dispatches to this when present to get continuous-arg samples.
    sampleBoundary: sampleBoundary_PQD,
  };
  QD.registerFamily('powerQD');

  // Export helpers used by downstream tabs (Schwarz / Sphere / Direct / tests)
  // and as a load-order sentinel for the dependent families. NOTE: the
  // singular/unbounded PQD families do NOT reuse modifiedResidues_PQD — they
  // build their (★) residues via the branch-cut-free 1/r# convolution
  // (Taylor.reciprocal ∘ Taylor.compose), since their (φ_in∘ψ/w)^α factor is
  // the rational 1/r# rather than the w^{1−α} that modifiedResidues encodes.
  // They depend only on QD.cprincipalRoot (+ the load-order guard on
  // QD.modifiedResidues_PQD confirming solver-pqd.js loaded first).
  QD.evalRHash_PQD        = evalRHash_PQD;
  QD.rHashTaylorAt_PQD    = rHashTaylorAt_PQD;
  QD.diskInitialGuess_PQD = diskInitialGuess_PQD;
  QD.cprincipalRoot       = cprincipalRoot;
  QD.cpowA_PQD            = cpowA;
  QD.modifiedResidues_PQD = modifiedResidues_PQD;
  // NB: QD._rHashVanishingGuard now lives in solver-pqd-common.js (QD.PqdCommon).

})();
