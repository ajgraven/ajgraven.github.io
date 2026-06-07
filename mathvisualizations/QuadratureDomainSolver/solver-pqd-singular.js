// =============================================================================
// solver-pqd-singular.js -- Bounded SINGULAR power-weighted quadrature domains
// (Family.powerQD_singular). The singular case is 0 ∈ Ω: the origin is an
// interior point of the domain, the image of a preimage z₀ ∈ 𝔻 (φ(z₀)=0).
//
// A bounded singular PQD of order α > 0 satisfies
//     ∫_Ω f(w) |w|^{2(α-1)} dA = ∮_∂Ω f(w) h(w) dw
// for analytic test functions f. Unlike singular LQDs there is NO point
// charge q at the origin (user-confirmed: the α-weight |w|^{2(α-1)} makes
// the quadrature data unique).
//
// Riemann-map characterization (thesis §4.3, Equation 1.2 of the
// "Analysis of LQDs" derivation; Corollary 1.1 singular case):
//     φ(z) = b_{z₀}(z) · (R#(z))^{1/α}
//   with the αth-root on the branch ANCHORED at the φ(0) normalization (QD.PqdCommon.phiAnchored;
  //   not principal — that fails for off-axis poles), where
//     b_{z₀}(z) = −(conj(z₀)/|z₀|)·(z − z₀)/(1 − conj(z₀) z)    (b(0) = |z₀|)
//     R#(z) = w₀^α/|z₀|^α + Σ_j Σ_k conj(A_{j,k}) z^k/(1 − conj(z_j) z)^k
//   (unimodular-prefactor Blaschke, shared with the singular-LQD family.)
//   The constant term w₀^α/|z₀|^α is HARDWIRED (the Faber-transform part of
//   r vanishes at ∞, so R#(0) = conj(r(∞)) = w₀^α/|z₀|^α). With b(0)=|z₀|,
//     φ(0) = |z₀| · (w₀^α/|z₀|^α)^{1/α} = |z₀| · w₀/|z₀| = w₀,
//   which holds for ANY z₀ — so φ(0)=w₀ does NOT constrain z₀ (see the
//   |z₀|-closing note below). arg(z₀) and |z₀| are both pinned by the
//   quadrature data: arg(z₀) via the locator/(★)/gauge, |z₀| via the mass
//   constraint (M).
//
// Newton unknowns: {z_j (2n), A_{j,k} (2M), z₀ (2)}. No γ (constant
// hardwired), no q (no point charge).
//
// Residual:
//   (●)   2n  locator: R#(z_j) = (a_j / b_{z₀}(z_j))^α       (α-lifted form)
//   (★)   2M  modified-residue Faber match (same D_{j,n} as non-singular
//             powerQD; phiTaylorAt now carries the Blaschke factor)
//   (φ₀)  2   φ(0) = w₀  (pins the discrete αth-root branch + off-axis w₀)
//   (G)   1   rotation gauge: Im φ'(0) = 0  (φ'(0) > 0 fixed in canonicalize)
//   (M)   2   MASS / AREA: the f=1 quadrature identity,
//                 ∫_Ω |w|^{2(α-1)} dA = ∮_∂Ω h dw   (= Σ residues).
//
// WHY (M) IS REQUIRED — the |z₀|-closing constraint:
//   The hardwired constant R#(0) = w₀^α/|z₀|^α was DERIVED assuming φ(0)=w₀,
//   so substituting it back makes φ(0)=w₀ identically true for EVERY z₀
//   (the |z₀| cancels: φ(0) = |z₀|·(w₀^α/|z₀|^α)^{1/α} = w₀, independent of
//   |z₀| and arg z₀). Hence (φ₀) contributes NO live constraint on z₀, and
//   (●)+(★)+(G) pin everything except |z₀| — one real DOF. A rank diagnostic
//   confirms this: the (●)+(★)+(φ₀)+(G) Jacobian is rank-deficient by exactly
//   1 (smallest σ ~5e-9, gap ~1e8 to the next), with the null direction in
//   z₀+A. The Riemann map is unique only up to a disk automorphism; with two
//   marked points (z=0↦w₀ and z=z₀↦0, the latter automatic since b_{z₀}(z₀)=0)
//   plus a rotation we over-specify the 3-DOF automorphism group, and one
//   combination (|z₀|) goes free.
//   The f=1 (mass) case of the quadrature identity closes it: the weighted
//   area ∫_Ω|w|^{2(α-1)}dA depends sensitively & monotonically on |z₀|, while
//   the RHS Σ residues is fixed data. Adding (M) makes the Jacobian full rank
//   (cond ~23) and the full weighted identity then passes to machine
//   precision for all test monomials. (The (★) pole-matching enforces only the
//   principal-part / local balance; (M) supplies the global mass balance.)
//
// Template: Family.boundedLQD_singular (solver-lqd-singular.js) for the
// Blaschke / z₀ schema; Family.powerQD (solver-pqd.js) for the R#/αth-root
// machinery (the anchored-branch root via QD.PqdCommon.phiAnchored). The (★)
// here uses the branch-cut-free 1/R# convolution (NOT powerQD's modifiedResidues
// — the singular factor (b_{z₀}∘ψ/w)^α = 1/R# is rational), so it shares the
// anchored αth-root + inverseFaberAtPole, not the generalized-binomial residue
// helper.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family) {
    throw new Error("solver-pqd-singular.js: solver.js must be loaded first");
  }
  if (!QD.modifiedResidues_PQD || !QD.cprincipalRoot) {
    throw new Error("solver-pqd-singular.js: solver-pqd.js must be loaded first");
  }

  // (φ's αth-root now uses the anchored continuous branch via
  // QD.PqdCommon.phiAnchored, not QD.cprincipalRoot — see evalPhi_PQDS.)
  // Trapezoidal-sweep resolution for the (M) mass/area residual. The boundary
  // integrand is analytic & periodic, so the sweep converges spectrally —
  // 512 samples is machine-precision on a smooth ∂Ω while keeping the
  // per-Newton-step cost modest (this is the only quadrature in the residual).
  const MASS_RESIDUAL_SAMPLES = 512;
  // NOTE: the singular (★) builds its modified residues inline (w·h base
  // convolved with the Taylor of 1/R#(ψ(w))); QD.modifiedResidues_PQD is the
  // non-singular helper and is intentionally not reused here.

  // ===========================================================================
  // 0. Blaschke factor — the unimodular-prefactor (thesis) convention,
  //    shared with the singular-LQD family:
  //      b_{z₀}(z) = −(conj(z₀)/|z₀|)·(z − z₀)/(1 − conj(z₀) z),   b_{z₀}(0) = |z₀|.
  //    b(0) = |z₀| (real positive) is CONSISTENT with the hardwired constant
  //    R#(0) = w₀^α/|z₀|^α: then φ(0) = |z₀|·(w₀^α/|z₀|^α)^{1/α} = w₀ for any
  //    z₀ (so z₀ is real positive in the canonical example, e.g. z₀ = 2/3).
  //    The rotation gauge is fixed by (G) Im φ'(0) = 0; |z₀| by the mass
  //    constraint (M). See the header for the DOF accounting.
  // ===========================================================================
  const bEval   = QD.LqdCommon.blaschkeEval;     // b_{z₀}(z)
  const bTaylor = QD.LqdCommon.blaschkeTaylor;    // Taylor of b at zc

  // ===========================================================================
  // 1. R# (constant w₀^α/|z₀|^α) + its Taylor. Branch structure identical to
  //    the non-singular powerQD; only the constant term differs.
  // ===========================================================================
  function r0Const(phi) {
    // w₀^α / |z₀|^α  (|z₀|^α is the real-power of |z₀|²).
    const w0a = Complex.cpow(phi.w0, phi.alpha);
    const z0abs_pow = Math.pow(Complex.abs2(phi.z0), 0.5 * phi.alpha);  // |z₀|^α
    return Complex.scale(w0a, 1 / z0abs_pow);
  }

  function evalRHash_PQDS(z, phi) {
    let result = r0Const(phi);
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

  function rHashTaylorAt_PQDS(z0pt, phi, L) {
    const result = Taylor.zero(L + 1);
    result[0] = r0Const(phi);
    for (const br of phi.branches) {
      const zjC = Complex.conj(br.z);
      const alpha_z = Complex.sub(Complex.ONE(), Complex.mul(zjC, z0pt));
      const alphaInv = Complex.inv(alpha_z);
      const uT = Taylor.zero(L + 1);
      uT[0] = Complex.mul(z0pt, alphaInv);
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

  // ===========================================================================
  // 2. φ evaluation:  φ(z) = b_{z₀}(z) · (R#(z))^{1/α}.
  // ===========================================================================
  // φ = b_{z₀}·(R#)^{1/α} with the αth-root on the branch ANCHORED at φ(0)=w0.
  // At z=0, R#(0)=w0^α/|z₀|^α has arg α·arg(w0) (|z₀|^α∈ℝ⁺) and b(0)=|z₀|∈ℝ⁺,
  // so anchoring argCont(R#(0))=α·arg(w0) yields φ(0)=w0 exactly and makes
  // off-axis poles work (the principal root would land on a wrong sheet).
  function evalPhi_PQDS(z, phi) {
    const root = QD.PqdCommon.phiAnchored(
      phi, z, evalRHash_PQDS, phi.alpha * Complex.arg(phi.w0), { re: 0, im: 0 });
    return Complex.mul(bEval(z, phi.z0), root);
  }

  // Taylor of φ at zc to order L. φ = b · exp((1/α)·log(R#)). When `anchorArg`
  // is supplied it overrides the constant term of log(R#) with the anchored
  // continuous arg — essential so the (★) target (which Taylor.exp scales by
  // e^{(1/α)·logRT[0]}) lands on the correct sheet for off-axis poles.
  function phiTaylorAt_PQDS(zc, phi, L, anchorArg) {
    const rT = rHashTaylorAt_PQDS(zc, phi, L);
    const logRT = Taylor.log(rT, L);
    if (anchorArg !== undefined) logRT[0] = { re: logRT[0].re, im: anchorArg };
    const scaledLog = Taylor.scaleComplex(logRT, { re: 1 / phi.alpha, im: 0 });
    const rootT = Taylor.exp(scaledLog, L);
    const bT = bTaylor(zc, phi.z0, L);
    return Taylor.mul(bT, rootT, L);
  }

  // ===========================================================================
  // 3. (★) target — branch-cut-free form (per the author's correction).
  //
  // Theorem 1.1: r = α·Φ_φ⁻¹(C_Ω*[ w·h(w)·(φ_in∘ψ(w)/w)^α ]) + C. Crucially,
  // the combined factor is branch-cut-free because
  //   (φ_in∘ψ(w)/w)^α = (b_{z₀}(z)/φ(z))^α = ((R#(z))^{−1/α})^α = 1/R#(z),
  // which is RATIONAL (z = ψ(w)). Splitting it as b_{z₀}(ψ)^α · w^{−α} would
  // reintroduce a spurious branch mismatch for non-integer α; we avoid that
  // by working with 1/R# directly.
  //
  // So the (★) integrand factors as:
  //   • base  w·h(w):  principal part at a_j has residues
  //       D̃_{j,k} = a_j·C_{j,k} + C_{j,k+1}        (C_{j,m_j+1} ≡ 0)
  //   • analytic factor 1/R#(ψ(w)): Taylor F_{j,i} at a_j (computed as the
  //       Taylor of 1/R#(z) at z_j, composed with the local inverse ψ̃).
  // Convolve:
  //       D'_{j,k} = Σ_{i=0}^{m_j−k} D̃_{j,k+i} · F_{j,i}
  // then A_{j,k} = α · inverseFaberAtPole(D'_{j,·}, phiTilde)[k], with
  // phiTilde the Taylor of (φ − a_j) at z_j (Blaschke included).
  //
  // Reduces to the non-singular powerQD (★): for a simple pole F_0 = 1/R#(z_j)
  // = 1/a_j^α and D̃_1 = a_j·C_1, giving D'_1 = C_1·a_j^{1−α} as expected.
  // ===========================================================================
  function computeTargetA_PQDS(phi, hData) {
    const alpha = phi.alpha;
    const target = [];
    for (let j = 0; j < hData.poles.length; j++) {
      const pole = hData.poles[j];
      const C = pole.principal;
      const mj = C.length;
      const aj = pole.a;
      const zj = phi.branches[j].z;

      // φ Taylor at z_j (Blaschke-included), constant zeroed → phiTilde. Built
      // on the anchored sheet (argCont(R#(z_j)) from φ(0)=w0) so the (★) target
      // is on the right αth-root sheet for off-axis poles.
      const argZj = QD.PqdCommon.argContAt(
        phi, zj, evalRHash_PQDS, alpha * Complex.arg(phi.w0), { re: 0, im: 0 });
      const phiT = phiTaylorAt_PQDS(zj, phi, mj, argZj);
      const phiTilde = Taylor.zero(mj + 1);
      for (let i = 1; i <= mj; i++) phiTilde[i] = Complex.clone(phiT[i]);

      // Base residues from w·h(w):  D̃_k = a_j·C_k + C_{k+1}.
      const Dtilde = new Array(mj);
      for (let k = 1; k <= mj; k++) {
        const Ck = C[k - 1];
        const Ck1 = (k < mj) ? C[k] : { re: 0, im: 0 };
        Dtilde[k - 1] = Complex.add(Complex.mul(aj, Ck), Ck1);
      }

      // Analytic factor 1/R#(ψ(w)): Taylor of 1/R#(z) at z_j, composed with ψ̃.
      const rT = rHashTaylorAt_PQDS(zj, phi, mj);
      const recipR = Taylor.reciprocal(rT, mj);                // 1/R# in (z−z_j)
      const psiTilde = Taylor.invert(phiTilde, mj);            // ψ̃ in s=(w−a_j)
      const F = Taylor.compose(recipR, psiTilde, mj);          // 1/R#(ψ(w)) in s

      // Convolve: D'_{j,k} = Σ_{i=0}^{m_j−k} D̃_{k+i} · F_i.
      const Dp = new Array(mj);
      for (let k = 1; k <= mj; k++) {
        let acc = { re: 0, im: 0 };
        for (let i = 0; i <= mj - k; i++) {
          acc = Complex.add(acc, Complex.mul(Dtilde[k + i - 1], F[i]));
        }
        Dp[k - 1] = acc;
      }

      const A_unscaled = QD.Faber.inverseFaberAtPole(Dp, phiTilde);
      target.push(A_unscaled.map(a => Complex.scale(a, alpha)));
    }
    return target;
  }

  // ===========================================================================
  // 4. Residual:  (●) locator, (★) coefficient, (φ₀) φ(0)=w₀.
  // ===========================================================================
  function residual_PQDS(phi, hData /*, options */) {
    const out = [];
    const alpha = phi.alpha;

    // (●) locator in the α-lifted, branch-cut-free form:
    //   φ(z_j) = a_j  ⇔  R#(z_j) = (a_j / b_{z₀}(z_j))^α.
    for (let j = 0; j < hData.poles.length; j++) {
      const zj = phi.branches[j].z;
      const rZj = evalRHash_PQDS(zj, phi);
      const bj  = bEval(zj, phi.z0);
      const lifted = Complex.cpow(Complex.div(hData.poles[j].a, bj), alpha);
      const diff = Complex.sub(rZj, lifted);
      out.push(diff.re, diff.im);
    }

    // (★) coefficient match.
    const target = computeTargetA_PQDS(phi, hData);
    for (let j = 0; j < hData.poles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], target[j][k]);
        out.push(diff.re, diff.im);
      }
    }

    // (φ₀) centering: φ(0) = w₀. With b(0)=|z₀|, R#(0)=w₀^α/|z₀|^α and the
    // αth-root anchored at argCont(R#(0))=α·arg(w₀), evalPhi_PQDS reproduces
    // φ(0)=w₀ by construction (for any w₀, off-axis included), so this residual
    // is ~0; we keep it as an explicit, harmless redundancy that also pins the
    // αth-root discrete branch.
    const phi0 = evalPhi_PQDS({ re: 0, im: 0 }, phi);
    const d0 = Complex.sub(phi0, phi.w0);
    out.push(d0.re, d0.im);

    // (G) rotation gauge: Im(φ'(0)) = 0 (standard Riemann normalization
    // φ'(0) > 0; sign resolved in canonicalizePhi). In the SINGULAR case
    // φ'(0) carries the Blaschke term, so Σ Im(A_{j,1})=0 (the non-singular
    // gauge) would fix the wrong direction. We read φ'(0) straight off the
    // order-1 Taylor of φ at 0 (which includes b'_{z₀}(0)).
    const phip0 = phiTaylorAt_PQDS({ re: 0, im: 0 }, phi, 1)[1];
    out.push(phip0.im);

    // (M) mass / area constraint — the |z₀|-closing equation (see header).
    const mass = massResidual_PQDS(phi, hData, MASS_RESIDUAL_SAMPLES);
    out.push(mass[0], mass[1]);

    return out;
  }

  // ===========================================================================
  // 5. Schema-driven pack/unpack. {z_j, z₀, A_{j,k}} — no γ.
  // ===========================================================================
  const SCHEMA_PQDS = [
    { kind: 'branchesZ', clamp: { side: 'in', cap: 0.9999 } },
    { kind: 'complex', name: 'z0', clamp: { side: 'in', cap: 0.9999, minR: 1e-3 } },
    { kind: 'branchesA' },
  ];
  function packPhi_PQDS(phi) { return QD.packPhiBySchema(phi, SCHEMA_PQDS); }
  function unpackPhi_PQDS(v, template) {
    return QD.unpackPhiBySchema(v, template, SCHEMA_PQDS, (phi) => {
      phi.family = 'powerQD_singular';
      phi.alpha  = template.alpha;
      phi.unbounded = false;
    });
  }

  // Canonicalize: the αth root is multi-valued; (●)/(★)/(G)/(M) pin z_j, A,
  // arg(z₀) and |z₀|, so the only residual ambiguity is the Z/2 sign resolved
  // exactly as boundedQD/powerQD (z → −z flips odd A's, negates z_j and z₀).
  // Enforce the standard Riemann normalization φ'(0) > 0.
  function canonicalizePhi_PQDS(phi) {
    // Standard Riemann normalization: φ'(0) > 0. If Re(φ'(0)) < 0, apply the
    // Z/2 rotation z → −z (flips z_j, z₀, and odd-indexed A's).
    const phip0 = phiTaylorAt_PQDS({ re: 0, im: 0 }, phi, 1)[1];
    if (phip0.re >= 0) return phi;
    return {
      family: 'powerQD_singular',
      alpha: phi.alpha,
      unbounded: false,
      w0: Complex.clone(phi.w0),
      z0: Complex.neg(phi.z0),
      branches: phi.branches.map(br => ({
        z: Complex.neg(br.z),
        A: br.A.map((a, k) => (k % 2 === 0 ? Complex.neg(a) : Complex.clone(a))),
      })),
    };
  }

  // ===========================================================================
  // 6. Seeds. With the prefactor Blaschke (b(0)=|z₀|), φ(0)=w₀ holds for any
  //    z₀ and does NOT constrain its sign; |z₀| is pinned by the mass
  //    constraint (M), arg(z₀) by the locator/(★)/gauge. Seed z₀ real positive
  //    (the canonical example lands at z₀ = 2/3, modulo the canonicalize sign
  //    flip → −2/3). z_j and A_{j,k} mirror the non-singular powerQD disk seed.
  // ===========================================================================
  // Seed strategy extracted to solvers/seeds/seeds-pqd-singular.js (B3).
  // Aliased locally so the Family entry + any internal callers keep their names.
  if (!QD.Seeds || !QD.Seeds.powerQD_singular) {
    throw new Error("solver-pqd-singular.js: QD.Seeds.powerQD_singular missing — seeds-pqd-singular.js must be loaded first");
  }
  const initialGuess_PQDS          = QD.Seeds.powerQD_singular.initialGuess;
  const perturbedInitialGuess_PQDS = QD.Seeds.powerQD_singular.perturbedInitialGuess;
  const diverseInitialGuess_PQDS   = QD.Seeds.powerQD_singular.diverseInitialGuess;

  // ===========================================================================
  // 6b. Continuous-branch unit-circle sweep (Q1.3 analog for the singular
  //     family). φ = b_{z₀}·(R#)^{1/α}: the Blaschke is single-valued, but
  //     the αth root needs the arg(R#) unwrap to avoid the atan2 cut. Returns
  //     { theta, z, w, phiPrime } with w/φ' on the continuous sheet.
  // ===========================================================================
  // φ = b_{z₀}(z)·(R#)^{1/α}: prefactor is the Blaschke factor b; chain rule
  // φ' = b'·root + b·root·R#'/(α·R#). Sweep driver in PqdCommon.
  function combine_PQDS(ctx) {
    const { z, rH, rHp, root, invAlpha, phi } = ctx;
    const bT = bTaylor(z, phi.z0, 1);
    const b = bT[0], bp = bT[1];
    const w = Complex.mul(b, root);
    const rootOverR = Complex.div(root, rH);             // (R#)^{1/α-1}
    const term2 = Complex.mul(Complex.mul(b, rootOverR), Complex.scale(rHp, invAlpha));
    const phiPrime = Complex.add(Complex.mul(bp, root), term2);
    return { w, phiPrime };
  }
  function sweepUnitCircle_PQDS(phi, N) {
    const anchorSpec = {
      evalRHashFn: evalRHash_PQDS,
      anchorPt: { re: 0, im: 0 },
      anchorArg0: phi.alpha * Complex.arg(phi.w0),
    };
    return QD.PqdCommon.sweepUnitCircle(phi, N, rHashTaylorAt_PQDS, combine_PQDS, anchorSpec);
  }

  // (M) Mass / area residual: the f=1 case of the weighted quadrature
  // identity,  ∫_Ω |w|^{2(α-1)} dA − ∮_∂Ω h dw = 0. The LHS is the k=0 term of
  // the verifier (boundary-integral / Green's form on the continuous sheet);
  // the RHS is Σ residues (= Σ first-principal-coefficients). This is the one
  // global constraint that pins the free |z₀| DOF (see header). N is large so
  // the trapezoidal sweep is at machine precision on the smooth ∂Ω.
  function massResidual_PQDS(phi, hData, N) {
    const alpha = phi.alpha;
    const samples = sweepUnitCircle_PQDS(phi, N);
    let lhs = { re: 0, im: 0 };
    for (let i = 0; i < N; i++) {
      const s = samples[i];
      const w2 = s.w.re * s.w.re + s.w.im * s.w.im;
      if (w2 < 1e-30) continue;
      const weight = Math.pow(w2, alpha - 1);              // |w|^{2(α-1)}
      let term = Complex.scale(Complex.conj(s.w), weight); // conj(w)·|w|^{2(α-1)}  (w^0)
      term = Complex.mul(term, s.phiPrime);
      term = Complex.mul(term, s.z);
      lhs = Complex.add(lhs, term);
    }
    lhs = Complex.scale(lhs, 1 / (alpha * N));
    let rhs = { re: 0, im: 0 };
    for (const pole of hData.poles) {
      if (pole.principal.length > 0) rhs = Complex.add(rhs, pole.principal[0]);
    }
    const d = Complex.sub(lhs, rhs);
    return [d.re, d.im];
  }

  // Family hook: continuous-arg boundary sampler (dispatched by
  // sampleBoundaryAdaptive in solver.js). φ = b_{z₀}(z)·(R#)^{1/α}; the
  // Blaschke prefactor is single-valued, the αth root unwraps from the left
  // neighbour's contArg. Budget densifies 3× near the higher-curvature origin.
  function sampleBoundary_PQDS(phi, baseSamples, maxExtra) {
    return QD.PqdCommon.sampleBoundaryViaSweep(
      phi, baseSamples, maxExtra, sweepUnitCircle_PQDS,
      (thMid, leftPt, ph) => QD.PqdCommon.boundaryMid(
        thMid, leftPt, ph, evalRHash_PQDS, (z, p) => bEval(z, p.z0)),
      3);
  }

  // ===========================================================================
  // 7. Identity verification — power-weighted monomials, single-valued form.
  //    Same as non-singular powerQD: LHS = (1/(αN)) Σ w^k·|w|^{2(α-1)}·conj(w)·φ'·z.
  //    Test monomials w^k (analytic in Ω). Origin ∈ Ω, but w^k is regular
  //    there, so no special treatment.
  // ===========================================================================
  function verifyQuadratureIdentity_PQDS(phi, hData, options = {}) {
    const N = options.numSamples ?? 500;
    const totalDeg = hData.poles.reduce((s, p) => s + p.principal.length, 0);
    const K = options.maxDegree ?? Math.max(totalDeg, 4);
    const alpha = phi.alpha;

    // Continuous-arg sweep (Q1.3): per-sample principal Taylor.log mis-
    // measures when arg(R#) crosses the atan2 cut; the unwrapped sweep gives
    // w and φ' on a consistent sheet.
    const samples = sweepUnitCircle_PQDS(phi, N);

    // Q1.4 — R# non-vanishing guard (shared with the bounded PQD family). For
    // φ = b_{z₀}·(R#)^{1/α} the αth root is single-valued only if R# has no
    // zero in 𝔻̄; a spurious root with winding(R#) ≠ 0 (or |R#| → 0 on ∂𝔻) is
    // rejected by forcing the identity to fail. The Blaschke factor is separate
    // and does not affect this condition.
    const rHashGuard = QD.PqdCommon.rHashVanishingGuard(samples);

    let areaScale = 0;
    for (const pole of hData.poles) {
      if (pole.principal.length > 0) areaScale += Complex.abs(pole.principal[0]);
    }
    if (areaScale === 0) areaScale = 1;

    const checks = [];
    let maxRelDiff = 0, maxAbsDiff = 0;
    for (let k = 0; k <= K; k++) {
      // LHS via PqdCommon (same single-valued weighted sum as powerQD), test
      // monomial f = w^k. skipNearZeroW2 = 1e-30: 0 ∈ Ω here, so drop samples
      // where the boundary nearly touches the origin. Scale +1/(αN).
      let lhs = QD.PqdCommon.accumulateWeightedLHS(samples, alpha, (w) => Complex.pow(w, k), 1e-30);
      lhs = Complex.scale(lhs, 1 / (alpha * N));

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
      maxAbsDiff, areaScale, maxDeg: K, numSamples: N, alpha, singular: true,
      rHashVanishes: rHashGuard.vanishes, rHashWinding: rHashGuard.winding,
    };
  }

  // ===========================================================================
  // 8. Register Family.powerQD_singular (more specific than powerQD).
  // ===========================================================================
  QD.Family.powerQD_singular = {
    name: 'powerQD_singular',
    enforceInDisk:  true,
    enforceOutDisk: false,
    matches(opts) {
      const a = opts && opts.alpha;
      return Number.isFinite(a) && a > 0 && a !== 1
          && !!opts.singular && !opts.unbounded && !opts.lqd;
    },
    normalizeOpts(opts, hData) {
      let w0 = opts.w0;
      if (!w0) {
        let sumRe = 0, sumIm = 0;
        for (const p of hData.poles) { sumRe += p.a.re; sumIm += p.a.im; }
        const n = hData.poles.length;
        w0 = n > 0 ? { re: sumRe / n, im: sumIm / n } : { re: 1, im: 0 };
      }
      const alpha = opts.alpha;
      if (!(alpha > 0) || alpha === 1) {
        throw new Error("Family.powerQD_singular: α must be real > 0, α ≠ 1");
      }
      if (Complex.abs2(w0) < QD.ZERO_THRESHOLD) {
        throw new Error("Family.powerQD_singular: w₀ = φ(0) must be nonzero");
      }
      return { w0, alpha, singular: true };
    },
    evalPhi: evalPhi_PQDS,
    phiTaylorAt: phiTaylorAt_PQDS,
    computeTargets(phi, hData) { return { A: computeTargetA_PQDS(phi, hData), F: null }; },
    residual: residual_PQDS,
    packPhi: packPhi_PQDS,
    unpackPhi: unpackPhi_PQDS,
    canonicalizePhi: canonicalizePhi_PQDS,
    initialGuess: initialGuess_PQDS,
    perturbedInitialGuess: perturbedInitialGuess_PQDS,
    diverseInitialGuess: diverseInitialGuess_PQDS,
    // Continuation in α from the classical limit (residue-/c-homotopies break
    // 0∈Ω for singular families). See QD.PqdCommon.continuationInAlpha.
    continuationSolve(hData, norm, options = {}) {
      return QD.PqdCommon.continuationInAlpha(hData, norm, options);
    },
    verifyQuadratureIdentity: verifyQuadratureIdentity_PQDS,
    sampleBoundary: sampleBoundary_PQDS,
  };
  QD.registerFamily('powerQD_singular');

  QD.evalRHash_PQDS     = evalRHash_PQDS;
  QD.rHashTaylorAt_PQDS = rHashTaylorAt_PQDS;
  QD.blaschkeNoPrefix_PQDS = bEval;

})();
