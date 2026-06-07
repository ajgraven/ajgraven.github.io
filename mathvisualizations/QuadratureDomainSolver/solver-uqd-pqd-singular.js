// =============================================================================
// solver-uqd-pqd-singular.js -- Unbounded SINGULAR power-weighted QDs
// (Family.unboundedPQD_singular). 0 ∈ Ω (origin interior), ∞ ∈ Ω; α > 0, α ≠ 1.
//
//     ∫_Ω f(w) |w|^{2(α-1)} dA = ∮_∂Ω f(w) h(w) dw,   f ∈ A₀(Ω)  (vanishing at ∞)
//
// Riemann-map characterization (thesis Eq 4.8/4.10, singular unbounded case;
// Theorem 4.4.3 / 4.5.2), φ : 𝔻* → Ω, φ(∞)=∞, φ'(∞)=c>0:
//     φ(z) = z · b_{z₀}(z) · (r#(z))^{1/α}            (principal αth root)
//     r#(z) = |c·z₀|^α + Σ_l G_l/z^l + Σ_j Σ_k conj(A_{j,k}) z^k/(1−conj(z_j)z)^k
// where z₀ ∈ 𝔻* is the unique preimage of the origin (φ(z₀)=0, via b_{z₀}(z₀)=0).
// The constant r#(∞)=|cz₀|^α is HARDWIRED: with the prefactor Blaschke
// (b_{z₀}(∞)=1/|z₀|), φ(z) ~ z·(1/|z₀|)·(|cz₀|^α)^{1/α} = c·z, so φ'(∞)=c.
// NO point charge q (the α-weight makes the quadrature data unique).
//
// Template: solver-uqd-pqd.js (non-singular UPQD — r#/φ/Laurent-at-∞ machinery)
// × solver-uqd-lqd-singular.js (Blaschke b_{z₀}, z₀∈𝔻* schema).
//
// Unknowns: {z_j (2n), A_{j,k} (2M), G_l (2N_G), z₀ (2)}.  No q.
// Residual: (●) 2n + (★)_A 2M + (★)_F 2N_G + (●_{z₀}) 2  →  square.
//
// z₀-CLOSURE (●_{z₀}) — r(z₀) = 0 (thesis Proposition 4.6.3): when h has no
// pole at 0, the rational r (= reflection of r#) must have a ROOT at the
// origin-preimage z₀. This supplies the 2 real equations that pin z₀; without
// it the system is rank-deficient by 2 (a Jacobian SV diagnostic on the
// singular monomial gives two zero singular values — hardwiring r#(∞)=|cz₀|^α
// makes the leading match z₀-independent, the same vacuousness as QB's
// φ(0)=w₀). With r(z₀)=0 the system is full rank and recovers the thesis
// monomial ground truth z₀ = γ^{1/(2α−1)} (Thm 4.5.2) to machine precision.
// (General case: h with a pole at 0 of order m₀ → r has an m₀-related root/
// pole at z₀ per Prop 4.6.3 — TODO, not yet handled; standard PQDs have h
// analytic at 0.)
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family || !QD.LqdCommon) {
    throw new Error("solver-uqd-pqd-singular.js: solver.js + solver-lqd-common.js must be loaded first");
  }
  if (!QD.cprincipalRoot) {
    throw new Error("solver-uqd-pqd-singular.js: solver-pqd.js must be loaded first");
  }

  const cprincipalRoot = QD.cprincipalRoot;
  const cpowA = (c, a) => Complex.cpow(c, a);
  const bEval = QD.LqdCommon.blaschkeEval;       // b_{z₀}(z), prefactor (b(z₀)=0)
  const bTaylor = QD.LqdCommon.blaschkeTaylor;   // Taylor of b at zc

  // r#(∞) = |c·z₀|^α  (real positive; c > 0).
  function r0Const(phi) {
    const cz0 = phi.c * Math.hypot(phi.z0.re, phi.z0.im);
    return { re: Math.pow(cz0, phi.alpha), im: 0 };
  }

  // ===========================================================================
  // 1. r# evaluation + Taylor  (identical structure to UPQD; only r0 differs).
  // ===========================================================================
  function evalRHash_UPQDS(z, phi) {
    let result = r0Const(phi);
    if (phi.polyA && phi.polyA.length > 0) {
      let zPow = Complex.clone(z);
      for (let l = 1; l <= phi.polyA.length; l++) {
        result = Complex.add(result, Complex.div(phi.polyA[l - 1], zPow));
        if (l < phi.polyA.length) zPow = Complex.mul(zPow, z);
      }
    }
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

  function rHashTaylorAt_UPQDS(z0pt, phi, L) {
    const result = Taylor.zero(L + 1);
    result[0] = r0Const(phi);
    if (phi.polyA && phi.polyA.length > 0) {
      const z0inv = Complex.inv(z0pt);
      let z0invPowL = Complex.clone(z0inv);
      for (let l = 1; l <= phi.polyA.length; l++) {
        let z0invPowLM = Complex.clone(z0invPowL);
        let binom = 1;
        for (let m = 0; m <= L; m++) {
          const sign = (m % 2 === 0) ? 1 : -1;
          const coef = Complex.scale(Complex.mul(phi.polyA[l - 1], z0invPowLM), sign * binom);
          result[m] = Complex.add(result[m], coef);
          if (m < L) { z0invPowLM = Complex.mul(z0invPowLM, z0inv); binom = binom * (l + m) / (m + 1); }
        }
        if (l < phi.polyA.length) z0invPowL = Complex.mul(z0invPowL, z0inv);
      }
    }
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
        for (let i = 0; i <= L; i++) result[i] = Complex.add(result[i], Complex.mul(AkC, uPow[i]));
        if (k < br.A.length) uPow = Taylor.mul(uPow, uT, L);
      }
    }
    return result;
  }

  // r(z) — the Schwarz companion / reflection of r# (r = conj(r#(1/conj(z)))):
  //   r(z) = |cz₀|^α + Σ_l conj(G_l) z^l + Σ_j Σ_k A_{j,k}/(z − z_j)^k.
  // Used by the z₀-closure (●_{z₀}): r(z₀) = 0 (see residual_UPQDS).
  function evalRCompanion_UPQDS(z, phi) {
    let result = r0Const(phi);                              // |cz₀|^α (real)
    if (phi.polyA && phi.polyA.length > 0) {
      let zPow = Complex.clone(z);                          // z^1
      for (let l = 1; l <= phi.polyA.length; l++) {
        result = Complex.add(result, Complex.mul(Complex.conj(phi.polyA[l - 1]), zPow));
        if (l < phi.polyA.length) zPow = Complex.mul(zPow, z);
      }
    }
    for (const br of phi.branches) {
      const diff = Complex.sub(z, br.z);
      let dPow = Complex.clone(diff);                        // (z − z_j)^1
      for (let k = 1; k <= br.A.length; k++) {
        result = Complex.add(result, Complex.div(br.A[k - 1], dPow));
        if (k < br.A.length) dPow = Complex.mul(dPow, diff);
      }
    }
    return result;
  }

  // ===========================================================================
  // 2. φ = z · b_{z₀}(z) · (r#(z))^{1/α}.
  // ===========================================================================
  function evalPhi_UPQDS(z, phi) {
    const root = cprincipalRoot(evalRHash_UPQDS(z, phi), phi.alpha);
    return Complex.mul(Complex.mul(z, bEval(z, phi.z0)), root);
  }

  // Taylor of φ at z0pt:  (z0pt+t) · b_{z₀}(z0pt+t) · exp((1/α) log r#).
  function phiTaylorAt_UPQDS(z0pt, phi, L) {
    const rT = rHashTaylorAt_UPQDS(z0pt, phi, L);
    const rootT = Taylor.exp(Taylor.scaleComplex(Taylor.log(rT, L), { re: 1 / phi.alpha, im: 0 }), L);
    const linT = Taylor.zero(L + 1);
    linT[0] = Complex.clone(z0pt);
    if (L >= 1) linT[1] = { re: 1, im: 0 };
    const bT = bTaylor(z0pt, phi.z0, L);
    return Taylor.mul(Taylor.mul(linT, bT, L), rootT, L);
  }

  // ===========================================================================
  // 3. (★)_A — finite poles. Factor (φ_in∘ψ/w)^α = (z·b/w)^α = 1/r#. Same
  //    branch-cut-free convolution as powerQD_singular / UPQD; phiTilde now
  //    carries the Blaschke (via phiTaylorAt_UPQDS).
  // ===========================================================================
  function computeTargetA_UPQDS(phi, hData) {
    const alpha = phi.alpha;
    const target = [];
    for (let j = 0; j < hData.poles.length; j++) {
      const pole = hData.poles[j];
      const C = pole.principal;
      const mj = C.length;
      const aj = pole.a;
      const zj = phi.branches[j].z;
      const phiT = phiTaylorAt_UPQDS(zj, phi, mj);
      const phiTilde = Taylor.zero(mj + 1);
      for (let i = 1; i <= mj; i++) phiTilde[i] = Complex.clone(phiT[i]);
      const Dtilde = new Array(mj);
      for (let k = 1; k <= mj; k++) {
        const Ck = C[k - 1];
        const Ck1 = (k < mj) ? C[k] : { re: 0, im: 0 };
        Dtilde[k - 1] = Complex.add(Complex.mul(aj, Ck), Ck1);
      }
      const rT = rHashTaylorAt_UPQDS(zj, phi, mj);
      const recipR = Taylor.reciprocal(rT, mj);
      const psiTilde = Taylor.invert(phiTilde, mj);
      const F = Taylor.compose(recipR, psiTilde, mj);
      const Dp = new Array(mj);
      for (let k = 1; k <= mj; k++) {
        let acc = { re: 0, im: 0 };
        for (let i = 0; i <= mj - k; i++) acc = Complex.add(acc, Complex.mul(Dtilde[k + i - 1], F[i]));
        Dp[k - 1] = acc;
      }
      const A_unscaled = QD.Faber.inverseFaberAtPole(Dp, phiTilde);
      target.push(A_unscaled.map(a => Complex.scale(a, alpha)));
    }
    return target;
  }

  // ===========================================================================
  // 3b. (★)_F — ∞-pole Laurent-matching residual (mirrors UPQD), but φ now
  //     carries the Blaschke b_{z₀}. M·sⁿ coeffs [0..n] = 0, n = polyPart−1.
  //     φ(s) = (1/s)·b_{z₀}(1/s)·(r#(s))^{1/α}.
  // ===========================================================================
  function shiftUp(T, j, L) {
    const out = Taylor.zero(L + 1);
    for (let i = 0; i + j <= L && i < T.length; i++) out[i + j] = Complex.clone(T[i]);
    return out;
  }
  // r#(s) Taylor in s = 1/z to order L (constant |cz₀|^α + G_l s^l + branch).
  function rHashInS_UPQDS(phi, L) {
    const out = Taylor.zero(L + 1);
    out[0] = r0Const(phi);
    const G = phi.polyA || [];
    for (let l = 1; l <= G.length && l <= L; l++) out[l] = Complex.add(out[l], G[l - 1]);
    for (const br of phi.branches) {
      const zjCinv = Complex.inv(Complex.conj(br.z));
      const uT = Taylor.zero(L + 1);
      let pw = Complex.clone(zjCinv);
      for (let m = 0; m <= L; m++) { uT[m] = Complex.neg(pw); pw = Complex.mul(pw, zjCinv); }
      let uPow = Taylor.truncate(uT, L);
      for (let k = 1; k <= br.A.length; k++) {
        const AkC = Complex.conj(br.A[k - 1]);
        for (let i = 0; i <= L; i++) out[i] = Complex.add(out[i], Complex.mul(AkC, uPow[i]));
        if (k < br.A.length) uPow = Taylor.mul(uPow, uT, L);
      }
    }
    return out;
  }
  // b_{z₀}(1/s) as a Taylor in s to order L:
  //   b(1/s) = −(conj(z₀)/|z₀|)·(1 − z₀ s)/(s − conj(z₀)).
  function blaschkeInS_UPQDS(phi, L) {
    const z0 = phi.z0;
    const z0C = Complex.conj(z0);
    const absZ0 = Math.hypot(z0.re, z0.im);
    const pref = Complex.scale(z0C, -1 / absZ0);             // −conj(z₀)/|z₀|
    // numerator (1 − z₀ s): Taylor [1, −z₀, 0, ...].
    const num = Taylor.zero(L + 1);
    num[0] = { re: 1, im: 0 };
    if (L >= 1) num[1] = Complex.neg(z0);
    // 1/(s − conj(z₀)) = −(1/conj(z₀))·1/(1 − s/conj(z₀)) = −Σ s^m / conj(z₀)^{m+1}.
    const denomInv = Taylor.zero(L + 1);
    const z0Cinv = Complex.inv(z0C);
    let pw = Complex.clone(z0Cinv);
    for (let m = 0; m <= L; m++) { denomInv[m] = Complex.neg(pw); pw = Complex.mul(pw, z0Cinv); }
    return Taylor.scaleComplex(Taylor.mul(num, denomInv, L), pref);
  }
  function laurentMatchAtInfinity_UPQDS(phi, hData) {
    const polyPart = hData.polyPart || [];
    const N = polyPart.length;
    if (N === 0) return [];
    const n = N - 1, alpha = phi.alpha;
    const rH = rHashInS_UPQDS(phi, n);
    const rho = Taylor.exp(Taylor.scaleComplex(Taylor.log(rH, n), { re: 1 / alpha, im: 0 }), n);
    const bS = blaschkeInS_UPQDS(phi, n);                    // b_{z₀}(1/s), order n
    // φ·s = b(1/s)·ρ(s)  (since φ = (1/s)·b·ρ ⇒ φ·s = b·ρ). So 1/φ = s/(b·ρ).
    const phiS = Taylor.mul(bS, rho, n);                     // = φ·s
    // r·s^{n+1} mod s^{n+1} = Σ_{l=1}^{n+1} conj(G_l) s^{n+1−l}.
    const G = phi.polyA || [];
    const rTimesS = Taylor.zero(n + 1);
    for (let l = 1; l <= N; l++) rTimesS[n + 1 - l] = Complex.conj(G[l - 1]);
    // firstTerm = (r·s^{n+1})·r# / (α · φ·s) = (r·r#·s^{n+1})/(α φ s)
    //           = sⁿ·(r·r#)/(αφ).   [since (r·s^{n+1})/(φ·s) = r·sⁿ/φ]
    const denom = Taylor.scaleComplex(phiS, { re: alpha, im: 0 });
    const recipDenom = Taylor.reciprocal(denom, n);
    let firstTerm = Taylor.mul(rTimesS, rH, n);
    firstTerm = Taylor.mul(firstTerm, recipDenom, n);
    // hTerm = Σ_{m=0}^n h_m · φ^m · sⁿ ;  φ^m·sⁿ = (φ·s)^m · s^{n−m} = phiS^m·s^{n−m}.
    let hTerm = Taylor.zero(n + 1);
    let phiSPow = Taylor.constant({ re: 1, im: 0 }, n + 1);  // (φ·s)^0
    for (let m = 0; m <= n; m++) {
      if (m > 0) phiSPow = Taylor.mul(phiSPow, phiS, n);
      hTerm = Taylor.add(hTerm, shiftUp(Taylor.scaleComplex(phiSPow, polyPart[m]), n - m, n));
    }
    const M = Taylor.sub(firstTerm, hTerm);
    const res = [];
    for (let i = 0; i <= n; i++) res.push(M[i] || { re: 0, im: 0 });
    return res;
  }

  // ===========================================================================
  // 4. Residual:  (●) locator, (★)_A, (★)_F.  z₀-closure appended below once
  //    the rank diagnostic identifies what (if anything) it needs.
  // ===========================================================================
  function residual_UPQDS(phi, hData, options) {
    options = options || {};
    const out = [];
    const alpha = phi.alpha;

    // (●) locator: r#(z_j) = (a_j / (z_j·b_{z₀}(z_j)))^α.
    for (let j = 0; j < hData.poles.length; j++) {
      const zj = phi.branches[j].z;
      const rZj = evalRHash_UPQDS(zj, phi);
      const base = Complex.mul(zj, bEval(zj, phi.z0));
      const lifted = cpowA(Complex.div(hData.poles[j].a, base), alpha);
      const diff = Complex.sub(rZj, lifted);
      out.push(diff.re, diff.im);
    }

    // (★)_A.
    const targetA = computeTargetA_UPQDS(phi, hData);
    for (let j = 0; j < hData.poles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], targetA[j][k]);
        out.push(diff.re, diff.im);
      }
    }

    // (★)_F ∞-pole Laurent matching.
    if (phi.polyA && phi.polyA.length > 0) {
      const mInf = laurentMatchAtInfinity_UPQDS(phi, hData);
      for (const coeff of mInf) out.push(coeff.re, coeff.im);
    }

    // (●_{z₀}) z₀-CLOSURE: r(z₀) = 0 (Proposition 4.6.3). When h has no pole at
    // 0, the rational r must have a root at the origin-preimage z₀. This is the
    // 2 real equations that pin z₀ (the (●)+(★)_A+(★)_F system leaves z₀ free —
    // hardwiring r#(∞)=|cz₀|^α makes the leading match z₀-independent). For the
    // singular monomial it forces z₀ real (= γ^{1/(2α−1)}), matching Thm 4.5.2.
    // TODO: h with a pole at 0 of order m₀ → r has an order-(m₀-related) root/
    // pole at z₀ (Prop 4.6.3, general case); not yet handled.
    const rZ0 = evalRCompanion_UPQDS(phi.z0, phi);
    out.push(rZ0.re, rZ0.im);

    return out;
  }

  // ===========================================================================
  // 5. Pack / unpack — {z_j, z₀, A_{j,k}, G_l}. z_j and z₀ both in 𝔻*.
  // ===========================================================================
  const SCHEMA_UPQDS = [
    { kind: 'branchesZ', clamp: { side: 'out', cap: QD.DISK_CLAMP_OUT } },
    { kind: 'complex', name: 'z0', clamp: { side: 'out', cap: QD.DISK_CLAMP_OUT, maxR: QD.Z0_MAX_RADIUS } },
    { kind: 'branchesA' },
    { kind: 'complexList', name: 'polyA' },
  ];
  function packPhi_UPQDS(phi) { return QD.packPhiBySchema(phi, SCHEMA_UPQDS); }
  function unpackPhi_UPQDS(v, template) {
    const phi = QD.unpackPhiBySchema(v, template, SCHEMA_UPQDS);
    phi.family = 'unboundedPQD_singular';
    phi.unbounded = true;
    phi.alpha = template.alpha;
    phi.c = template.c;
    return phi;
  }
  function canonicalizePhi_UPQDS(phi) { return phi; }

  // ===========================================================================
  // 6. Seeds. z_j ≈ a_j/c (exterior); z₀ from the monomial closed form when
  //    available (Thm 4.5.2: z₀ ≈ γ^{1/(2α−1)}, real); G_l seed as UPQD.
  // ===========================================================================
  // Seed strategy extracted to solvers/seeds/seeds-uqd-pqd-singular.js (B3).
  // Aliased locally so the continuation loop + Family entry keep their names.
  if (!QD.Seeds || !QD.Seeds.unboundedPQD_singular) {
    throw new Error("solver-uqd-pqd-singular.js: QD.Seeds.unboundedPQD_singular missing — seeds-uqd-pqd-singular.js must be loaded first");
  }
  const initialGuess_UPQDS          = QD.Seeds.unboundedPQD_singular.initialGuess;
  const perturbedInitialGuess_UPQDS = QD.Seeds.unboundedPQD_singular.perturbedInitialGuess;
  const diverseInitialGuess_UPQDS   = QD.Seeds.unboundedPQD_singular.diverseInitialGuess;

  // ===========================================================================
  // 7. Continuous-branch sweep + identity verifier (test class A₀, f=1/(w−b)^k).
  // ===========================================================================
  // φ = z·b_{z₀}(z)·(r#)^{1/α}: prefactor is z·b; chain rule
  // φ' = root·(b + z·b') + z·b·root·r#'/(α·r#). Sweep driver in PqdCommon.
  function combine_UPQDS(ctx) {
    const { z, rH, rHp, root, alpha, phi } = ctx;
    const bT = bTaylor(z, phi.z0, 1);
    const b = bT[0], bp = bT[1];
    const zb = Complex.mul(z, b);
    const w = Complex.mul(zb, root);
    const dzb = Complex.add(b, Complex.mul(z, bp));        // (z·b)'
    const term1 = Complex.mul(root, dzb);
    const term2 = Complex.mul(zb, Complex.mul(root, Complex.div(rHp, Complex.scale(rH, alpha))));
    const phiPrime = Complex.add(term1, term2);
    return { w, phiPrime };
  }
  function sweepUnitCircle_UPQDS(phi, N) {
    return QD.PqdCommon.sweepUnitCircle(phi, N, rHashTaylorAt_UPQDS, combine_UPQDS);
  }

  // Family hook: continuous-arg boundary sampler. φ = z·b_{z₀}(z)·(r#)^{1/α};
  // the z·Blaschke prefactor is single-valued, the αth root unwraps from the
  // left neighbour's contArg. Budget densifies 3× near the origin.
  function sampleBoundary_UPQDS(phi, baseSamples, maxExtra) {
    return QD.PqdCommon.sampleBoundaryViaSweep(
      phi, baseSamples, maxExtra, sweepUnitCircle_UPQDS,
      (thMid, leftPt, ph) => QD.PqdCommon.boundaryMid(
        thMid, leftPt, ph, evalRHash_UPQDS, (z, p) => Complex.mul(z, bEval(z, p.z0))),
      3);
  }
  function verifyQuadratureIdentity_UPQDS(phi, hData, options = {}) {
    // The singular boundary (Blaschke × αth-root, with 0 ∈ Ω) needs a denser
    // trapezoidal sweep than the non-singular families to resolve the
    // higher-curvature region near the origin (spectral convergence: 600→4e-4,
    // 1200→1e-10, 2400→4e-14 on the one-pole α=2 case). Enforce a MINIMUM of
    // 2000 even when the solve pipeline passes a smaller univalenceSamples,
    // so the reported identity reflects the true (machine-precision) accuracy.
    const N = Math.max(options.numSamples || 0, 2000);
    const maxOrder = options.maxDegree ?? 3;
    const numTestPoints = options.numTestPoints ?? 3;
    const alpha = phi.alpha;
    const samples = sweepUnitCircle_UPQDS(phi, N);
    let cx = 0, cy = 0;
    for (const s of samples) { cx += s.w.re; cy += s.w.im; }
    cx /= N; cy /= N;
    let maxDev = 0;
    for (const s of samples) { const d = Math.hypot(s.w.re - cx, s.w.im - cy); if (d > maxDev) maxDev = d; }
    // Test points in K (the bounded complement), avoiding the origin (∈ Ω).
    const polygonPts = samples.map(s => s.w);
    const inside = (x, y) => {
      let cr = 0;
      for (let i = 0; i < polygonPts.length; i++) {
        const j = (i + 1) % polygonPts.length;
        const yi = polygonPts[i].im, yj = polygonPts[j].im;
        if ((yi > y) !== (yj > y)) { const t = (y - yi) / (yj - yi); if (polygonPts[i].re + t * (polygonPts[j].re - polygonPts[i].re) > x) cr++; }
      }
      return (cr % 2) === 1;
    };
    const cand = [{ re: cx, im: cy }];
    for (const frac of [0.15, 0.3, 0.45]) for (let i = 0; i < 8; i++) {
      const a = 2 * Math.PI * i / 8, r = frac * maxDev;
      cand.push({ re: cx + r * Math.cos(a), im: cy + r * Math.sin(a) });
    }
    const testPoints = [];
    for (const b of cand) {
      if (Math.hypot(b.re, b.im) < 1e-2) continue;
      if (inside(b.re, b.im)) testPoints.push(b);
      if (testPoints.length >= numTestPoints) break;
    }
    let areaScale = 0;
    for (const pole of hData.poles) if (pole.principal.length > 0) areaScale += Complex.abs(pole.principal[0]);
    for (const cc of (hData.polyPart || [])) areaScale += Complex.abs(cc);
    if (areaScale === 0) areaScale = 1;
    const checks = [];
    let maxRelDiff = 0, maxAbsDiff = 0;
    for (let pIdx = 0; pIdx < testPoints.length; pIdx++) {
      const b = testPoints[pIdx];
      for (let k = 1; k <= maxOrder; k++) {
        // LHS via PqdCommon, f = 1/(w−b)^k. skipNearZeroW2 = 1e-30 (0 ∈ Ω here).
        // Scale −1/(αN). RHS = shared unbounded finite-pole + polyPart residues.
        let lhs = QD.PqdCommon.accumulateWeightedLHS(
          samples, alpha, (w) => Complex.inv(Complex.pow(Complex.sub(w, b), k)), 1e-30);
        lhs = Complex.scale(lhs, -1 / (alpha * N));
        const rhs = QD.PqdCommon.unboundedTestPointRHS(hData, k, b);
        const dd = Complex.sub(lhs, rhs);
        const absDiff = Complex.abs(dd);
        const scale = Math.max(Complex.abs(lhs), Complex.abs(rhs), areaScale);
        const relDiff = absDiff / scale;
        if (relDiff > maxRelDiff) maxRelDiff = relDiff;
        if (absDiff > maxAbsDiff) maxAbsDiff = absDiff;
        checks.push({ bIdx: pIdx, k, lhs, rhs, absDiff, relDiff });
      }
    }
    return { checks, maxRelDiff, maxAbsDiff, areaScale, testPoints, maxDeg: maxOrder, numSamples: N, alpha, unbounded: true, singular: true };
  }

  // ===========================================================================
  // 8. Register Family.unboundedPQD_singular.
  // ===========================================================================
  QD.Family.unboundedPQD_singular = {
    name: 'unboundedPQD_singular',
    enforceInDisk: false,
    enforceOutDisk: true,
    matches(opts) {
      const a = opts && opts.alpha;
      return Number.isFinite(a) && a > 0 && a !== 1
          && !!opts.unbounded && !!opts.singular && !opts.lqd;
    },
    normalizeOpts(opts, hData) {
      const c = opts.c;
      if (typeof c !== 'number' || !(c > 0)) throw new Error("Family.unboundedPQD_singular: opts.c must be positive");
      const alpha = opts.alpha;
      if (!(alpha > 0) || alpha === 1) throw new Error("Family.unboundedPQD_singular: α must be real > 0, α ≠ 1");
      const nPoles = (hData && hData.poles && hData.poles.length) || 0;
      const nPoly  = (hData && hData.polyPart && hData.polyPart.length) || 0;
      if (nPoles === 0 && nPoly === 0) {
        throw new Error("Family.unboundedPQD_singular: no quadrature data — h needs a finite pole or a polynomial part");
      }
      return { unbounded: true, singular: true, alpha, c };
    },
    evalPhi: evalPhi_UPQDS,
    phiTaylorAt: phiTaylorAt_UPQDS,
    computeTargets(phi, hData) { return { A: computeTargetA_UPQDS(phi, hData), F: laurentMatchAtInfinity_UPQDS(phi, hData) }; },
    residual: residual_UPQDS,
    packPhi: packPhi_UPQDS,
    unpackPhi: unpackPhi_UPQDS,
    canonicalizePhi: canonicalizePhi_UPQDS,
    initialGuess: initialGuess_UPQDS,
    perturbedInitialGuess: perturbedInitialGuess_UPQDS,
    diverseInitialGuess: diverseInitialGuess_UPQDS,
    // Continuation in α from the classical limit (residue-/c-homotopies break
    // the unbounded singular ansatz). See QD.PqdCommon.continuationInAlpha.
    continuationSolve(hData, norm, options = {}) {
      return QD.PqdCommon.continuationInAlpha(hData, norm, options);
    },
    verifyQuadratureIdentity: verifyQuadratureIdentity_UPQDS,
    sampleBoundary: sampleBoundary_UPQDS,
  };
  QD.registerFamily('unboundedPQD_singular');

  QD.evalRHash_UPQDS = evalRHash_UPQDS;
  QD.laurentMatchAtInfinity_UPQDS = laurentMatchAtInfinity_UPQDS;

})();
