// =============================================================================
// solver-uqd-lqd.js -- Unbounded NON-SINGULAR log-weighted QDs
//                       (Family.unboundedLQD)
//
// Setting: Ω unbounded simply connected with 0 ∉ Ω̄ (so 0 ∈ K, the bounded
// complement), and ∞ ∈ Ω. Identity (with ρ_0 = |w|⁻²):
//     ∫_Ω f(w) / |w|² dA = ∮_∂Ω f(w) h(w) dw,    f ∈ L¹_a(Ω; ρ_0)
// where the test class is analytic-in-Ω + vanishes-at-∞ (the LQD analog of
// A_0(Ω)).
//
// Riemann-map parametrization (Eq. 5.6, non-singular form):
//
//     φ(z) = c · z · exp(r#(z)),      z ∈ 𝔻*,   c > 0
//
// with φ(∞) = ∞ and φ'(∞) = c (real positive). r# has the same shape as in
// bounded LQDs (poles at 1/conj(z_j) ∈ 𝔻 for z_j ∈ 𝔻*, conjugated residues):
//
//     r#(z) = Σ_j Σ_k conj(A_{j,k}) · z^k / (1 - conj(z_j) z)^k.
//
// Gauge (∞-condition): φ'(∞) = c forces exp(r#(∞)) = 1 ⇒ r#(∞) = 0.
// We enforce this by ABSORBING r#(∞) into the parametrization:
//     φ(z) := c · z · exp( r#(z) − r#(∞) ).
// In this gauge r̃# := r# − r#(∞) automatically satisfies r̃#(∞) = 0; we
// retain (z_j, A_{j,k}) as the unknowns and let two of the A_{j,1} degrees
// of freedom be "absorbed" by the r#(∞)-subtraction implicit in evalPhi /
// phiTaylorAt. The system is then SQUARE (no over-determination from the
// ∞-gauge), in contrast to the bounded LQD families where a Z/2 sign
// remains and contributes an extra real equation. Here c > 0 (real
// positive) fully pins the disk-automorphism gauge.
//
// Note: r#(∞) is a closed-form rational evaluation:
//     r#(∞) = Σ_j Σ_k conj(A_{j,k}) · (-1)^k / conj(z_j)^k.
//
// Unknowns / equations:
//   Unknowns: 2n (z_j) + 2M (A_{j,k}), where M = Σ m_j.
//   Equations:
//     (●)   φ(z_j) = a_j                                  2n
//     (★)   modified-residue Faber match at each z_j      2M
//   Total: 2n + 2M; system is square.
//
// Identity test class: f(w) = 1/w^k for k = 1, 2, 3.
//   LHS: ∮_∂Ω (1/w^k)(ln|w|²/w) dw   (same LQD kernel as bounded LQDs)
//   RHS: 2πi · [Σ residues at a_j ∈ Ω of f·h  +  residue at ∞ of f·h]
//     finite-pole residue: C_{j,s} · (-1)^{s-1} · binom(k+s-2, s-1) / a_j^{k+s-1}
//     ∞ residue contribution from polyPart of h: -C_{∞, k-1}  (only if degree k-1 ≤ m_∞)
//
// Polynomial-h support (Andrew Graven derivation, HANDOFF entry #21):
//   The φ parameterization is extended in evalPhi / phiTaylorAt to honor a
//   nonempty phi.lqdBeta via the term B(1/z) = Σ_{l=1..N} β_l/z^l inside the
//   exp argument, with N = polyPart.length = m_∞+1. The (★)_F equations
//   match β against polyPart by applying Φ_φ^{-1} to w·h(w) (not h(w)):
//
//     w·h(w) = Σ_j Σ_s C_{j,s} (1/(w−a_j)^{s−1} + a_j/(w−a_j)^s)
//            + Σ_{i=1..m_∞+1} C_{∞,i-1} w^i
//
//   The finite-pole part of w·h re-indexes to exactly the modified residues
//   D_{j,s} = a_j·C_{j,s} + C_{j,s+1} already used by the (★)_A loop. The
//   (s=1) terms (1/(w−a_j)^0 = constants) accumulate into the augmented
//   polynomial-at-∞:  P̃ = [Σ_j C_{j,1}, C_∞,0, …, C_∞,m_∞].
//
//   Applying QD.Faber.inverseFaberAtInfinity(P̃, f̃, c) gives F̃ of length
//   m_∞+2. The l = 0 output is the gauge-absorbed constant (handled by the
//   r#(z) − r#(∞) normalization) and discarded. For l = 1, …, m_∞+1 we
//   equate F̃_l − (Laurent of r# at ∞ in 1/z^l) = β_l, giving N complex
//   equations matched 1-to-1 against the N free β's. System remains square.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family || !QD.LqdCommon) {
    throw new Error("solver-uqd-lqd.js: solver.js + solver-lqd-common.js must be loaded first");
  }

  // ===========================================================================
  // 1. Phi evaluation: φ(z) = c · z · exp(r#(z) − r#(∞) + B(1/z))
  // ---------------------------------------------------------------------------
  // The B(1/z) addition is the polynomial-h extension (HANDOFF #21–22).
  // Without it (β empty) the parametrization is the "minimal" form
  //   φ = c · z · exp(r#(z) − r#(∞)),
  // matching only finite-pole h. With B(t) = Σ_{l=1..N} β_l t^l (no constant
  // term — that would conflict with the c gauge), the Laurent expansion of φ
  // at ∞ gains N additional free coefficients that are pinned to h's
  // polynomial part via the (★)_F equations (computeTargetF_UQDL below).
  // β lives on phi.lqdBeta as a Complex[] of length N (= polyPart.length).
  // Empty / missing → no polynomial-h extension; this evaluator returns the
  // legacy finite-pole-only φ exactly.
  // ===========================================================================
  function evalPhi_UQDL(z, phi) {
    const r = QD.LqdCommon.evalRHash(z, phi);
    const rInf = rHashAtInfinity(phi);
    const bPart = evalB_OverZ(phi, z);          // B(1/z), zero when lqdBeta empty
    const rEff = Complex.add(Complex.sub(r, rInf), bPart);
    const ea = Math.exp(rEff.re);
    const expR = { re: ea * Math.cos(rEff.im), im: ea * Math.sin(rEff.im) };
    const cz = Complex.scale(z, phi.c);
    return Complex.mul(cz, expR);
  }

  // B(1/z) = Σ_l β_l / z^l. Imported from LqdCommon (HANDOFF #27 dedupe).
  const evalB_OverZ = QD.LqdCommon.evalB_OverZ;

  // ===========================================================================
  // 2. Taylor of φ at z = z_c
  // ---------------------------------------------------------------------------
  // φ(z_c + t) = c (z_c + t) · exp(r#(z_c) − r#(∞) + B(1/z_c))
  //              · exp(r#(z_c+t) − r#(z_c) + B(1/(z_c+t)) − B(1/z_c))
  //            = K · (z_c + t) · exp(rTilde(t) + bTilde(t))
  // with K = c · exp(r#(z_c) − r#(∞) + B(1/z_c)),
  //      rTilde(t) = r#(z_c+t) − r#(z_c),
  //      bTilde(t) = B(1/(z_c+t)) − B(1/z_c).
  // Both rTilde and bTilde have zero constant term, so Taylor.exp converges.
  function phiTaylorAt_UQDL(zc, phi, L) {
    const rT = QD.LqdCommon.rHashTaylorAt(zc, phi, L);   // Taylor of r# at z_c
    const rInf = rHashAtInfinity(phi);
    const bT = bOverZTaylorAt(phi, zc, L);               // Taylor of B(1/z) at z_c
    const b0 = bT[0];                                    // B(1/z_c)
    const r0minusInf = Complex.sub(rT[0], rInf);         // r#(z_c) − r#(∞)
    const expArgConst = Complex.add(r0minusInf, b0);     // full constant in exp

    // Strip constants from both Taylors and sum the variable parts.
    const argTilde = Taylor.zero(L + 1);
    for (let l = 1; l <= L; l++) {
      argTilde[l] = Complex.add(rT[l], bT[l]);
    }
    const expArgTilde = Taylor.exp(argTilde, L);

    const ea = Math.exp(expArgConst.re);
    const expConst = { re: ea * Math.cos(expArgConst.im), im: ea * Math.sin(expArgConst.im) };
    const K = Complex.scale(expConst, phi.c);

    const lin = Taylor.zero(L + 1);
    lin[0] = Complex.clone(zc);
    if (L >= 1) lin[1] = { re: 1, im: 0 };

    const linTimesExp = Taylor.mul(lin, expArgTilde, L);
    const out = new Array(L + 1);
    for (let l = 0; l <= L; l++) out[l] = Complex.mul(K, linTimesExp[l]);
    return out;
  }

  // Taylor expansion of B(1/z) at z = z_c. Imported from LqdCommon
  // (HANDOFF #27 dedupe).
  const bOverZTaylorAt = QD.LqdCommon.bOverZTaylorAt;

  // ===========================================================================
  // 3. r#(∞) — lifted to QD.LqdCommon.rHashAtInfinity (shared with unbounded
  //    singular LQDs, which use the same ∞-gauge-absorption trick).
  // ===========================================================================
  const rHashAtInfinity = QD.LqdCommon.rHashAtInfinity;

  // ===========================================================================
  // 4. Targets — Faber match at finite poles (LQD-common)
  // ===========================================================================
  function computeTargetA_UQDL(phi, hData) {
    return QD.LqdCommon.computeFaberTargetA(phi, hData, phiTaylorAt_UQDL);
  }

  // ===========================================================================
  // 4b. Polynomial-h target  (★)_F  for β
  // ---------------------------------------------------------------------------
  // Derivation in file header. Returns an array of length N = polyPart.length
  // containing the target values for phi.lqdBeta[0], …, phi.lqdBeta[N-1].
  // Returns [] when there is no polynomial part (so the (★)_F block in the
  // residual is silently empty for finite-pole-only h).
  // ===========================================================================
  function computeTargetF_UQDL(phi, hData) {
    const polyPart = hData.polyPart || [];
    const N = polyPart.length;
    if (N === 0) return [];

    // P̃ = [Σ_j C_{j,1},  C_∞,0,  C_∞,1,  …,  C_∞,m_∞]   length N+1.
    let sumRe = 0, sumIm = 0;
    for (const pole of hData.poles) {
      if (pole.principal.length > 0) {
        sumRe += pole.principal[0].re;
        sumIm += pole.principal[0].im;
      }
    }
    const Ptilde = new Array(N + 1);
    Ptilde[0] = { re: sumRe, im: sumIm };
    for (let i = 0; i < N; i++) Ptilde[i + 1] = Complex.clone(polyPart[i]);

    // f̃ = Laurent of φ at ∞ in 1/z^l for l = 0..N-1  (length N).
    const fTilde = QD.LqdCommon.phiLaurentAtInfinity_UQDL(phi, N);

    // F̃ = inverseFaberAtInfinity(P̃, f̃, c)  (length N+1, indices 0..N).
    const Ftilde = QD.Faber.inverseFaberAtInfinity(Ptilde, fTilde, phi.c);

    // (★)_F target: β_l = F̃_l for l = 1, …, N. (The F̃_0 output is the
    // gauge-absorbed constant of r at ∞, discarded.) Empirically — via a
    // direct β-sweep at fixed (z_j, A) — this no-subtraction form is what
    // makes the LQD boundary identity hold. The classical-QD−style
    // subtraction of the r# Laurent (which would be the natural-looking
    // formula) is not the correct one for the LQD kernel.
    const targets = new Array(N);
    for (let l = 1; l <= N; l++) {
      targets[l - 1] = Complex.clone(Ftilde[l]);
    }
    return targets;
  }

  // ===========================================================================
  // 5. Residual: 2n + 2M + 2N (square).
  // ---------------------------------------------------------------------------
  //   (●)    2n     φ(z_j) = a_j                  locator
  //   (★)_A  2M     A_{j,k} = target_A             modified-residue Faber
  //   (★)_F  2N     β_l    = target_F              poly-h-at-∞ Faber
  // The ∞-gauge r#(∞) = 0 is built into evalPhi / phiTaylorAt via the
  // r#(z) − r#(∞) subtraction; no explicit equation here. The (★)_F block
  // is empty (N = 0) when hData.polyPart is empty / absent, so finite-pole-
  // only callers see the exact same system as before.
  function residual_UQDL(phi, hData, options) {
    options = options || {};
    const out = [];
    for (let j = 0; j < hData.poles.length; j++) {
      const phiZj = evalPhi_UQDL(phi.branches[j].z, phi);
      const diff = Complex.sub(phiZj, hData.poles[j].a);
      out.push(diff.re, diff.im);
    }
    const target = computeTargetA_UQDL(phi, hData);
    for (let j = 0; j < hData.poles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], target[j][k]);
        out.push(diff.re, diff.im);
      }
    }
    // (★)_F polynomial-h block. Inert when phi.lqdBeta is empty.
    if (phi.lqdBeta && phi.lqdBeta.length > 0) {
      const targetF = computeTargetF_UQDL(phi, hData);
      const N = phi.lqdBeta.length;
      for (let l = 0; l < N; l++) {
        const tl = (l < targetF.length) ? targetF[l] : { re: 0, im: 0 };
        const diff = Complex.sub(phi.lqdBeta[l], tl);
        out.push(diff.re, diff.im);
      }
    }
    return out;
  }

  // No Z/2 ambiguity in the unbounded case — c > 0 already pins everything.
  function canonicalizePhi_UQDL(phi) { return phi; }

  // ===========================================================================
  // 6. Pack / unpack — schema-driven
  // ===========================================================================
  // Layout: {z_j} (n) + {A_{j,k}} (M) + {β_l} (N). z_j ∈ 𝔻*, clamped to
  // |z| ≥ 1.0001. phi.lqdBeta length = template.lqdBeta.length = polyPart
  // .length = N (= m_∞ + 1); zero-length when polyPart is empty/absent so
  // the slot collapses cleanly for finite-pole-only h.
  const SCHEMA_UQDL = [
    { kind: 'branchesZ', clamp: { side: 'out', cap: QD.DISK_CLAMP_OUT } },
    { kind: 'branchesA' },
    { kind: 'complexList', name: 'lqdBeta' },
  ];
  function packPhi_UQDL(phi)         { return QD.packPhiBySchema(phi, SCHEMA_UQDL); }
  function unpackPhi_UQDL(v, template) {
    const phi = QD.unpackPhiBySchema(v, template, SCHEMA_UQDL);
    phi.family = 'unboundedLQD';
    phi.unbounded = true;
    phi.c = template.c;
    return phi;
  }

  // ===========================================================================
  // 7. Initial guess — trivial-LQD bootstrap (r# = 0 ⇒ z_j = a_j/c, A = D/c^k)
  // ===========================================================================
  // For non-singular unbounded LQD the trivial case r# ≡ 0 gives φ(z) = c·z
  // (the exterior of the disk of radius c). For general h:
  //   z_j ≈ a_j / c   (preimage of a_j under φ ≈ cz)
  //   A_{j,k} ≈ D_{j,k} / c^k   (Faber match for ψ̃(t) = t/c)
  // where D_{j,s} = a_j · C_{j,s} + C_{j,s+1} are the modified residues
  // (Theorem 5.4.2). At small c the preimages have |z_j| ≫ 1 (well outside
  // 𝔻̄), making this an excellent warm-start; for larger c, Newton walks
  // the solution via continuation in c.
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
    // φ with β = [0, ..., 0]. In the trivial r# ≈ 0, β ≈ 0 limit this reduces
    // to β_l ≈ c^l · conj(polyPart[l-1]); using the full target instead is
    // marginally more accurate (the existing z_j ≠ 0 contributes via
    // rHashLaurentAtInfinity) and reuses the same code path.
    const polyPart = hData.polyPart || [];
    const phiInit = {
      family: 'unboundedLQD',
      unbounded: true,
      c, w0: undefined,
      branches,
      lqdBeta: polyPart.map(() => ({ re: 0, im: 0 })),
    };
    if (polyPart.length > 0) {
      const targetF = computeTargetF_UQDL(phiInit, hData);
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

  function diverseInitialGuess_UQDL(hData, norm, rng, r) {
    const polyPart = hData.polyPart || [];
    const base = {
      family: 'unboundedLQD',
      unbounded: true,
      c: norm.c, w0: undefined,
      branches: QD.LqdCommon.diverseSeedBranches(hData, rng, { zMin: 1.05, zMax: 30 }),
      lqdBeta: polyPart.map(() => ({ re: 0, im: 0 })),
    };
    if (polyPart.length > 0) {
      base.lqdBeta = computeTargetF_UQDL(base, hData).map(c => ({ re: c.re, im: c.im }));
    }
    return base;
  }

  // ===========================================================================
  // 8. Continuation in c
  // ===========================================================================
  // Analogous to continuationInC_UQD (unbounded QD): walk c from a small
  // starting value up to the user's target, warm-starting Newton each step.
  function continuationSolve_UQDL(hData, norm, options) {
    options = options || {};
    const {
      cStart       = null,
      growFactor   = 1.6,
      shrinkFactor = 0.5,
      minStep      = 1e-4,
      maxSteps     = 80,
      newton       = {},
    } = options;
    const cTarget = norm.c;

    let minA = Infinity;
    for (const p of hData.poles) {
      const m = Complex.abs(p.a);
      if (m > 0 && m < minA) minA = m;
    }
    const startGuess = cStart ?? Math.min(cTarget, isFinite(minA) ? 0.25 * minA : 0.25);
    if (startGuess <= 0) {
      return { success: false, error: "continuationInC (LQD): invalid starting c", trace: [] };
    }

    const trace = [];
    let c = startGuess;
    let phi = initialGuess_UQDL(hData, { c });

    let warmup;
    while (true) {
      warmup = QD.newtonSolve(phi, hData, newton);
      if (warmup.success) { phi = warmup.phi; break; }
      c *= shrinkFactor;
      if (c < minStep) {
        return {
          success: false,
          error: "continuationInC (LQD): warmup failed even at c=" + c.toExponential(2),
          phi: warmup.phi, trace,
        };
      }
      phi = initialGuess_UQDL(hData, { c });
    }
    trace.push({ c, ok: true, residual: warmup.residual });

    if (c >= cTarget - 1e-12) {
      return { success: true, phi, iterations: 0, residual: warmup.residual,
               trace, method: "continuation-in-c-lqd" };
    }

    let lastSuccessC = c;
    let stepSize = Math.max((cTarget - c) * 0.4, minStep);
    for (let step = 0; step < maxSteps; step++) {
      if (lastSuccessC >= cTarget - 1e-12) break;
      const nextC = Math.min(cTarget, lastSuccessC + stepSize);
      const phiNext = QD.clonePhi(phi);
      phiNext.c = nextC;
      const result = QD.newtonSolve(phiNext, hData, newton);
      if (result.success) {
        phi = result.phi;
        lastSuccessC = nextC;
        trace.push({ c: nextC, ok: true, residual: result.residual });
        stepSize *= growFactor;
      } else {
        stepSize *= shrinkFactor;
        trace.push({ c: nextC, ok: false, residual: result.residual ?? null });
        if (stepSize < minStep) {
          return {
            success: false,
            error: "continuationInC (LQD): step underflow at c=" + lastSuccessC.toFixed(4),
            phi, trace, lastC: lastSuccessC,
          };
        }
      }
    }
    if (lastSuccessC < cTarget - 1e-9) {
      return {
        success: false,
        error: "continuationInC (LQD): max steps reached at c=" + lastSuccessC.toFixed(4),
        phi, trace, lastC: lastSuccessC,
      };
    }
    return { success: true, phi, iterations: 0,
             residual: trace[trace.length - 1].residual,
             trace, method: "continuation-in-c-lqd" };
  }

  // ===========================================================================
  // 9. Identity verifier — test class f(w) = 1/w^k via LqdCommon skeleton
  // ===========================================================================
  function verifyQuadratureIdentity_UQDL(phi, hData, options) {
    options = options || {};
    const maxOrder = options.maxDegree ?? 3;
    // For unbounded Ω, the parametrization z = e^{iθ} with θ increasing
    // gives ∂Ω traversed CW around Ω (≡ CCW around K). The residue theorem
    // takes ∂Ω CCW around Ω, so the trapezoidal-rule integral picks up a
    // sign of −1 relative to the residue sum. We pre-multiply RHS by −2πi
    // (instead of +2πi) to absorb this orientation difference.
    const minusTwoPiI = { re: 0, im: -2 * Math.PI };

    return QD.LqdCommon.verifyIdentityGeneric(phi, hData, options, {
      phiTaylorFn: phiTaylorAt_UQDL,
      // Same boundary kernel as bounded LQDs: ln|w|² / w.
      boundaryKernel(w) {
        const absW2 = Complex.abs2(w);
        if (absW2 < 1e-30) return null;
        return Complex.scale(Complex.inv(w), Math.log(absW2));
      },
      buildTestFunctions(phi, hData) {
        const tests = [];
        const polyPart = hData.polyPart || [];
        const m_inf = polyPart.length - 1;
        for (let k = 1; k <= maxOrder; k++) {
          // RHS = 2πi · [Σ Res at finite poles + Res at ∞]
          // Finite-pole residue (b = 0 in the standard 1/(w-b)^k formula):
          //   C_{j,s} · (-1)^{s-1} · binom(k+s-2, s-1) / a_j^{k+s-1}
          let rhsSum = { re: 0, im: 0 };
          for (const pole of hData.poles) {
            for (let sIdx = 0; sIdx < pole.principal.length; sIdx++) {
              const s = sIdx + 1;
              const C = pole.principal[sIdx];
              const sign = (s % 2 === 0) ? -1 : 1;
              const coef = QD.binomialCoeff(k + s - 2, s - 1);
              if (coef === 0) continue;
              const denom = Complex.pow(pole.a, k + s - 1);
              const term = Complex.div(C, denom);
              rhsSum = Complex.add(rhsSum, Complex.scale(term, sign * coef));
            }
          }
          // Residue at ∞ contribution from polyPart at b = 0 (matches the
          // formula in verifyQuadratureIdentity_UQD at b = 0): only l = k-1
          // contributes, giving −C_{∞, k-1}.
          if (m_inf >= 0) {
            const l = k - 1;
            if (l >= 0 && l <= m_inf) {
              rhsSum = Complex.add(rhsSum, Complex.scale(polyPart[l], -1));
            }
          }
          const rhs = Complex.mul(minusTwoPiI, rhsSum);
          tests.push({
            label: '1/w^' + k,
            f: (w) => Complex.inv(Complex.pow(w, k)),
            residueRhs: rhs,
            tag: { k },
          });
        }
        return tests;
      },
      // Flag for the UI's describeTestClass.
      resultFlags: { lqdUnbounded: true, unbounded: true, maxDeg: maxOrder },
    });
  }

  // ===========================================================================
  // 10. Register Family.unboundedLQD
  // ===========================================================================
  QD.Family.unboundedLQD = {
    name: 'unboundedLQD',
    enforceInDisk:  false,
    enforceOutDisk: true,
    matches(opts) { return !!(opts && opts.lqd && opts.unbounded && !opts.singular); },
    normalizeOpts(opts, hData) {
      const c = opts.c;
      if (typeof c !== 'number' || !(c > 0)) {
        throw new Error("Family.unboundedLQD: opts.c must be a positive number");
      }
      return { lqd: true, unbounded: true, c };
    },
    evalPhi: evalPhi_UQDL,
    phiTaylorAt: phiTaylorAt_UQDL,
    computeTargets(phi, hData) {
      return { A: computeTargetA_UQDL(phi, hData), F: computeTargetF_UQDL(phi, hData) };
    },
    residual: residual_UQDL,
    packPhi: packPhi_UQDL,
    unpackPhi: unpackPhi_UQDL,
    canonicalizePhi: canonicalizePhi_UQDL,
    initialGuess: initialGuess_UQDL,
    perturbedInitialGuess: perturbedInitialGuess_UQDL,
    diverseInitialGuess: diverseInitialGuess_UQDL,
    continuationSolve: continuationSolve_UQDL,
    verifyQuadratureIdentity: verifyQuadratureIdentity_UQDL,
  };
  QD.registerFamily('unboundedLQD');

})();
