// =============================================================================
// solver-uqd-pqd.js -- Unbounded power-weighted quadrature domains
// (Family.unboundedPQD). 0 ∉ Ω (non-singular); α > 0, α ≠ 1.
//
// An unbounded simply-connected Ω with bounded complement K, satisfying
//     ∫_Ω f(w) |w|^{2(α-1)} dA = ∮_∂Ω f(w) h(w) dw,   f ∈ A₀(Ω)
// where A₀(Ω) = analytic on Ω AND vanishing at ∞. The decay of A₀ test
// functions makes the weighted integral converge for ALL α > 0 even though
// Ω is unbounded. α = 1 recovers classical unbounded QDs (Family.unboundedQD);
// this file handles α ≠ 1.
//
// Riemann-map characterization (thesis Theorems 4.3.1/4.3.3/4.3.6, Eq 4.8/4.10,
// Corollary 4.3.3), φ : 𝔻* → Ω, φ(∞) = ∞, φ'(∞) = c > 0:
//     φ(z) = z · (r#(z))^{1/α}                         (principal αth root)
//     r#(z) = c^α + Σ_l G_l/z^l + Σ_j Σ_k conj(A_{j,k}) z^k/(1 − conj(z_j) z)^k
// r# is analytic in 𝔻* (incl. ∞) with r#(∞) = c^α — the conformal-radius
// normalization (φ ~ c·z at ∞), which we HARDWIRE (user-confirmed). The
// Laurent block Σ_l G_l/z^l (l = 1..|polyPart|) is the reflection of the
// degree-n polynomial p in Cor 4.3.3 and is present only when h has a
// polynomial part (pole at ∞ / nonzero h(∞)); G_l absent when h(∞) = 0.
//
// Newton unknowns: {z_j (2n, exterior |z_j| > 1), A_{j,k} (2M), G_l (2N_G)}.
// c is a user input (as Family.unboundedQD). No gauge (c > 0 fixes rotation).
//
// Residual:
//   (●)   2n   locator: r#(z_j) = (a_j / z_j)^α       (α-lifted, branch-cut-free)
//   (★)_A 2M   modified-residue Faber match — IDENTICAL skeleton to
//              computeTargetA_PQDS: base w·h convolved with Taylor of
//              1/r#(ψ(w)) (since (φ_in∘ψ/w)^α = (ψ/w)^α = 1/r#, rational),
//              then α·inverseFaberAtPole.
//   (★)_F 2N_G ∞-pole block (polynomial-h). Hardwired r#(∞)=c^α means the
//              constant is NOT a (★) unknown.
//
// Template: solver-uqd.js (exterior schema, Laurent-at-∞, verifier RHS,
// continuation-in-c) × solver-pqd.js / solver-pqd-singular.js (cpow, αth-root
// Taylor, branch-cut-free (★), continuous-branch sweep).
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family) {
    throw new Error("solver-uqd-pqd.js: solver.js must be loaded first");
  }
  if (!QD.cprincipalRoot) {
    throw new Error("solver-uqd-pqd.js: solver-pqd.js must be loaded first");
  }

  const cprincipalRoot = QD.cprincipalRoot;       // c^{1/α}, principal branch
  const cpowA = (c, a) => Complex.cpow(c, a);     // c^α, principal branch

  // ===========================================================================
  // 1. r# evaluation + Taylor.  r#(z) = c^α + Σ_l G_l/z^l + Σ_j Σ_k conj(A)·u^k,
  //    u = z/(1 − conj(z_j) z).  Constant c^α hardwired from phi.c.
  // ===========================================================================
  function r0Const(phi) { return cpowA({ re: phi.c, im: 0 }, phi.alpha); }   // (c^α; c real > 0)

  function evalRHash_UPQD(z, phi) {
    let result = r0Const(phi);
    // Laurent-at-∞ block: Σ_l G_l / z^l  (G_l = phi.polyA[l-1], l = 1..N_G).
    if (phi.polyA && phi.polyA.length > 0) {
      let zPow = Complex.clone(z);                            // z^1
      for (let l = 1; l <= phi.polyA.length; l++) {
        result = Complex.add(result, Complex.div(phi.polyA[l - 1], zPow));
        if (l < phi.polyA.length) zPow = Complex.mul(zPow, z);
      }
    }
    // Finite-pole branch terms (exterior z_j).
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

  // Taylor of r# at z = z0 up to order L.
  function rHashTaylorAt_UPQD(z0, phi, L) {
    const result = Taylor.zero(L + 1);
    result[0] = r0Const(phi);

    // Laurent block G_l/z^l Taylor at z0 (same expansion as solver-uqd.js's
    // polyA, but starting at l = 1 — no z^0 term, that is the hardwired c^α).
    if (phi.polyA && phi.polyA.length > 0) {
      const z0inv = Complex.inv(z0);
      let z0invPowL = Complex.clone(z0inv);                   // z0^{-l}, l starts at 1
      for (let l = 1; l <= phi.polyA.length; l++) {
        let z0invPowLM = Complex.clone(z0invPowL);
        let binom = 1;
        for (let m = 0; m <= L; m++) {
          const sign = (m % 2 === 0) ? 1 : -1;
          const coef = Complex.scale(Complex.mul(phi.polyA[l - 1], z0invPowLM), sign * binom);
          result[m] = Complex.add(result[m], coef);
          if (m < L) {
            z0invPowLM = Complex.mul(z0invPowLM, z0inv);
            binom = binom * (l + m) / (m + 1);
          }
        }
        if (l < phi.polyA.length) z0invPowL = Complex.mul(z0invPowL, z0inv);
      }
    }

    // Finite-pole contributions.
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

  // ===========================================================================
  // 2. φ evaluation:  φ(z) = z · (r#(z))^{1/α}.
  // ===========================================================================
  function evalPhi_UPQD(z, phi) {
    const root = cprincipalRoot(evalRHash_UPQD(z, phi), phi.alpha);
    return Complex.mul(z, root);
  }

  // Taylor of φ at z0:  φ = (z0 + t) · exp((1/α) log r#).
  function phiTaylorAt_UPQD(z0, phi, L) {
    const rT = rHashTaylorAt_UPQD(z0, phi, L);
    const logRT = Taylor.log(rT, L);
    const scaledLog = Taylor.scaleComplex(logRT, { re: 1 / phi.alpha, im: 0 });
    const rootT = Taylor.exp(scaledLog, L);
    const linT = Taylor.zero(L + 1);
    linT[0] = Complex.clone(z0);
    if (L >= 1) linT[1] = { re: 1, im: 0 };
    return Taylor.mul(linT, rootT, L);
  }

  // ===========================================================================
  // 3. (★)_A — finite poles. Branch-cut-free form: factor (φ_in∘ψ/w)^α =
  //    (ψ(w)/w)^α = 1/r#(z). Verbatim port of computeTargetA_PQDS.
  // ===========================================================================
  function computeTargetA_UPQD(phi, hData) {
    const alpha = phi.alpha;
    const target = [];
    for (let j = 0; j < hData.poles.length; j++) {
      const pole = hData.poles[j];
      const C = pole.principal;
      const mj = C.length;
      const aj = pole.a;
      const zj = phi.branches[j].z;

      const phiT = phiTaylorAt_UPQD(zj, phi, mj);
      const phiTilde = Taylor.zero(mj + 1);
      for (let i = 1; i <= mj; i++) phiTilde[i] = Complex.clone(phiT[i]);

      // Base residues from w·h(w):  D̃_k = a_j·C_k + C_{k+1}.
      const Dtilde = new Array(mj);
      for (let k = 1; k <= mj; k++) {
        const Ck = C[k - 1];
        const Ck1 = (k < mj) ? C[k] : { re: 0, im: 0 };
        Dtilde[k - 1] = Complex.add(Complex.mul(aj, Ck), Ck1);
      }

      // Analytic factor 1/r#(ψ(w)): Taylor of 1/r#(z) at z_j composed with ψ̃.
      const rT = rHashTaylorAt_UPQD(zj, phi, mj);
      const recipR = Taylor.reciprocal(rT, mj);
      const psiTilde = Taylor.invert(phiTilde, mj);
      const F = Taylor.compose(recipR, psiTilde, mj);

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
  // 3b. (★)_F — ∞-pole block (polynomial-h), the Laurent-at-∞ matching residual.
  //
  // The generalized Schwarz function S_α(w) = (r·r#)(ψ(w))/(α·w) splits as
  // S_α = h + G with G ∈ A₀(Ω) (analytic in Ω, → 0 at ∞). Since the pole at ∞
  // is in Ω, h's polynomial part = the polynomial-at-∞ part of S_α. With
  // s = 1/z this is the condition M(z) := (r·r#)/(α·φ) − h_poly(φ) = O(1/z),
  // i.e. the s⁰…sⁿ coeffs of M·sⁿ vanish — n+1 complex equations for the n+1
  // unknown G_l (l = 1..n+1, n = degree of h's polynomial part).
  //
  // M·sⁿ is a genuine power series in s (built below with the Taylor API). The
  // r·s^{n+1} factor reduces to Σ_{l=1}^{n+1} conj(G_l) s^{n+1−l} (the c^α +
  // branch parts of r = c^α + Σ conj(G_l)z^l + R## land at orders ≥ n+1 → out
  // of range). The conjugation is FORCED by r being the reflection of r#; it is
  // confirmed empirically by the identity verifier (resolves the complex-h sign).
  //
  // Reduces to the n=0 closed form G_1 = α·h_0·c^{1−α} (kept as the seed in
  // initialGuess_UPQD). This is a residual (nonlinear in G_l), not a target.
  // ===========================================================================

  // Multiply a Taylor series T (coeffs of s^0..s^L) by s^j: prepend j zeros,
  // truncate back to order L.
  function shiftUp(T, j, L) {
    const out = Taylor.zero(L + 1);
    for (let i = 0; i + j <= L && i < T.length; i++) out[i + j] = Complex.clone(T[i]);
    return out;
  }

  // r#(s) as a Taylor in s = 1/z to order L: c^α + Σ_l G_l s^l + Σ branch u(s)^k,
  // u(s) = 1/(s − conj(z_j)),  u(s)[m] = −(1/conj(z_j))^{m+1}.
  function rHashInS_UPQD(phi, L) {
    const out = Taylor.zero(L + 1);
    out[0] = r0Const(phi);
    const G = phi.polyA || [];
    for (let l = 1; l <= G.length && l <= L; l++) out[l] = Complex.add(out[l], G[l - 1]);
    for (const br of phi.branches) {
      const zjCinv = Complex.inv(Complex.conj(br.z));  // 1/conj(z_j)
      const uT = Taylor.zero(L + 1);
      let pw = Complex.clone(zjCinv);                  // (1/conj(z_j))^{m+1}, m=0 → ^1
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

  // (★)_F residual coefficients: the s⁰..sⁿ coeffs of M·sⁿ (n = polyPart.length−1).
  function laurentMatchAtInfinity_UPQD(phi, hData) {
    const polyPart = hData.polyPart || [];
    const N = polyPart.length;
    if (N === 0) return [];
    const n = N - 1;                                   // degree of h's poly part
    const alpha = phi.alpha;

    const rH = rHashInS_UPQD(phi, n);                  // r#(s), order n
    const logRH = Taylor.log(rH, n);
    const rho = Taylor.exp(Taylor.scaleComplex(logRH, { re: 1 / alpha, im: 0 }), n); // ρ=(r#)^{1/α}

    // r·s^{n+1} mod s^{n+1} = Σ_{l=1}^{n+1} conj(G_l) s^{n+1−l}  (coeff[n+1−l]).
    const G = phi.polyA || [];
    const rTimesS = Taylor.zero(n + 1);
    for (let l = 1; l <= N; l++) rTimesS[n + 1 - l] = Complex.conj(G[l - 1]);

    // firstTerm = r·s^{n+1} · r# / (α·ρ).
    const denom = Taylor.scaleComplex(rho, { re: alpha, im: 0 });
    const recipDenom = Taylor.reciprocal(denom, n);
    let firstTerm = Taylor.mul(rTimesS, rH, n);
    firstTerm = Taylor.mul(firstTerm, recipDenom, n);

    // hTerm = Σ_{m=0}^n h_m · ρ^m · s^{n−m}.
    let hTerm = Taylor.zero(n + 1);
    let rhoPow = Taylor.constant({ re: 1, im: 0 }, n + 1);   // ρ^0
    for (let m = 0; m <= n; m++) {
      if (m > 0) rhoPow = Taylor.mul(rhoPow, rho, n);
      const term = shiftUp(Taylor.scaleComplex(rhoPow, polyPart[m]), n - m, n);
      hTerm = Taylor.add(hTerm, term);
    }

    const M = Taylor.sub(firstTerm, hTerm);            // M·sⁿ, coeffs [0..n] = 0.
    const res = [];
    for (let i = 0; i <= n; i++) res.push(M[i] || { re: 0, im: 0 });
    return res;
  }

  // ===========================================================================
  // 4. Residual:  (●) locator, (★)_A, (★)_F. No gauge (c > 0 fixes rotation).
  // ===========================================================================
  function residual_UPQD(phi, hData /*, options */) {
    const out = [];
    const alpha = phi.alpha;

    // (●) r#(z_j) = (a_j / z_j)^α.
    for (let j = 0; j < hData.poles.length; j++) {
      const zj = phi.branches[j].z;
      const rZj = evalRHash_UPQD(zj, phi);
      const lifted = cpowA(Complex.div(hData.poles[j].a, zj), alpha);
      const diff = Complex.sub(rZj, lifted);
      out.push(diff.re, diff.im);
    }

    // (★)_A coefficient match.
    const targetA = computeTargetA_UPQD(phi, hData);
    for (let j = 0; j < hData.poles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], targetA[j][k]);
        out.push(diff.re, diff.im);
      }
    }

    // (★)_F ∞-pole / Laurent-at-∞ matching residual.
    if (phi.polyA && phi.polyA.length > 0) {
      const mInf = laurentMatchAtInfinity_UPQD(phi, hData);
      for (const coeff of mInf) out.push(coeff.re, coeff.im);
    }
    return out;
  }

  // ===========================================================================
  // 5. Pack / unpack (manual, mirroring solver-uqd.js) + canonicalize.
  // ===========================================================================
  function packPhi_UPQD(phi) {
    const v = [];
    for (const br of phi.branches) v.push(br.z.re, br.z.im);
    for (const br of phi.branches) for (const a of br.A) v.push(a.re, a.im);
    if (phi.polyA) for (const G of phi.polyA) v.push(G.re, G.im);
    return v;
  }

  function unpackPhi_UPQD(v, template) {
    const phi = {
      family: 'unboundedPQD',
      unbounded: true,
      alpha: template.alpha,
      c: template.c,
      w0: undefined,
      polyA: [],
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
    if (template.polyA) {
      for (let l = 0; l < template.polyA.length; l++) {
        phi.polyA.push({ re: v[idx], im: v[idx + 1] });
        idx += 2;
      }
    }
    return phi;
  }

  // No gauge ambiguity (c > 0 pins the rotation); the αth-root Z/α sheet is
  // fixed by the principal branch in cprincipalRoot.
  function canonicalizePhi_UPQD(phi) { return phi; }

  // ===========================================================================
  // 6. Seeds — reuse the unboundedQD disk seed (z_j ≈ a_j/c, exterior),
  //    with α-power applied to the residues, then attach the G_l Laurent seed.
  // ===========================================================================
  // Seed strategy extracted to solvers/seeds/seeds-uqd-pqd.js (B3). Aliased
  // locally so the continuation loop + Family entry keep their names.
  if (!QD.Seeds || !QD.Seeds.unboundedPQD) {
    throw new Error("solver-uqd-pqd.js: QD.Seeds.unboundedPQD missing — seeds-uqd-pqd.js must be loaded first");
  }
  const initialGuess_UPQD          = QD.Seeds.unboundedPQD.initialGuess;
  const perturbedInitialGuess_UPQD = QD.Seeds.unboundedPQD.perturbedInitialGuess;
  const diverseInitialGuess_UPQD   = QD.Seeds.unboundedPQD.diverseInitialGuess;

  // ===========================================================================
  // 7. Continuous-branch unit-circle sweep + identity verifier.
  //    φ = z·(r#)^{1/α}: track unwrapped arg(r#) to avoid the atan2 cut.
  // ===========================================================================
  // φ = z·(r#)^{1/α}: prefactor is the leading z; chain rule
  // φ' = root·(1 + z·r#'/(α·r#)). Sweep driver in PqdCommon.
  function combine_UPQD(ctx) {
    const { z, rH, rHp, root, alpha } = ctx;
    const w = Complex.mul(z, root);
    const inner = Complex.add(Complex.ONE(), Complex.div(Complex.mul(z, rHp), Complex.scale(rH, alpha)));
    const phiPrime = Complex.mul(root, inner);
    return { w, phiPrime };
  }
  function sweepUnitCircle_UPQD(phi, N) {
    return QD.PqdCommon.sweepUnitCircle(phi, N, rHashTaylorAt_UPQD, combine_UPQD);
  }

  // Family hook: continuous-arg boundary sampler. φ = z·(r#)^{1/α}; the leading
  // z prefactor is single-valued, the αth root unwraps from the left
  // neighbour's contArg. Budget densifies 3× as a fallback.
  function sampleBoundary_UPQD(phi, baseSamples, maxExtra) {
    return QD.PqdCommon.sampleBoundaryViaSweep(
      phi, baseSamples, maxExtra, sweepUnitCircle_UPQD,
      (thMid, leftPt, ph) => QD.PqdCommon.boundaryMid(
        thMid, leftPt, ph, evalRHash_UPQD, (z) => z),
      3);
  }

  // Identity verifier — test functions f(w) = 1/(w − b)^k for b ∈ K (the
  // bounded complement). Mirrors verifyQuadratureIdentity_UQD with the
  // |w|^{2(α-1)} weight and 1/α scaling on the LHS (the generalized Schwarz
  // function S_α = (1/α) w |w|^{2(α-1)} on ∂Ω).
  function verifyQuadratureIdentity_UPQD(phi, hData, options = {}) {
    // ≥1500-node contour integral (see verifyQuadratureIdentity_UQD): the peaked
    // 1/(w−b)^k integrand is grossly under-resolved at the old 600 nodes as c grows.
    const N = Math.max(options.numSamples ?? 0, 1500);
    const maxOrder = options.maxDegree ?? 3;
    const numTestPoints = options.numTestPoints ?? 3;
    const alpha = phi.alpha;

    const samples = sweepUnitCircle_UPQD(phi, N);

    // Test points b ∈ K, ranked by clearance from ∂Ω and the poles of h (shared
    // QD.chooseHoleTestPoints) — replaces the geometry-blind centroid+0.18·maxDev
    // placement that drifted onto a pole/the boundary at large c (the c* bug).
    const testPoints = QD.chooseHoleTestPoints(samples.map(s => s.w), hData.poles, { numTestPoints });

    let areaScale = 0;
    for (const pole of hData.poles) {
      if (pole.principal.length > 0) areaScale += Complex.abs(pole.principal[0]);
    }
    for (const cc of (hData.polyPart || [])) areaScale += Complex.abs(cc);
    if (areaScale === 0) areaScale = 1;

    if (testPoints.length === 0) {
      return { checks: [], maxRelDiff: Infinity, maxAbsDiff: Infinity, areaScale,
               testPoints: [], maxDeg: maxOrder, numSamples: N, alpha, unbounded: true,
               warning: 'no test points clear of ∂Ω/poles' };
    }

    const checks = [];
    let maxRelDiff = 0, maxAbsDiff = 0;
    for (let pIdx = 0; pIdx < testPoints.length; pIdx++) {
      const b = testPoints[pIdx];
      for (let k = 1; k <= maxOrder; k++) {
        // LHS = (−1/(αN)) Σ f(w)·(|w|²)^{α-1}·conj(w)·φ'·z, accumulated by
        // PqdCommon. Test function f = 1/(w−b)^k (analytic in the unbounded Ω).
        // No near-zero skip (0 ∉ Ω̄ for the non-singular family).
        let lhs = QD.PqdCommon.accumulateWeightedLHS(
          samples, alpha, (w) => Complex.inv(Complex.pow(Complex.sub(w, b), k)), 0);
        lhs = Complex.scale(lhs, -1 / (alpha * N));

        // RHS = residues of f(w)·h(w): finite poles + polyPart at ∞ (shared form).
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
    return { checks, maxRelDiff, maxAbsDiff, areaScale, testPoints, maxDeg: maxOrder, numSamples: N, alpha, unbounded: true };
  }

  // ===========================================================================
  // 8. Register Family.unboundedPQD (more specific than unboundedQD; α ≠ 1).
  // ===========================================================================
  QD.Family.unboundedPQD = {
    name: 'unboundedPQD',
    enforceInDisk: false,
    enforceOutDisk: true,
    matches(opts) {
      const a = opts && opts.alpha;
      return Number.isFinite(a) && a > 0 && a !== 1
          && !!opts.unbounded && !opts.lqd && !opts.singular;
    },
    normalizeOpts(opts, hData) {
      const c = opts.c;
      if (typeof c !== 'number' || !(c > 0)) {
        throw new Error("Family.unboundedPQD: opts.c must be a positive number");
      }
      const alpha = opts.alpha;
      if (!(alpha > 0) || alpha === 1) {
        throw new Error("Family.unboundedPQD: α must be real > 0, α ≠ 1 (α = 1 is classical unbounded QD)");
      }
      const nPoles = (hData && hData.poles && hData.poles.length) || 0;
      const nPoly  = (hData && hData.polyPart && hData.polyPart.length) || 0;
      if (nPoles === 0 && nPoly === 0) {
        throw new Error("Family.unboundedPQD: no quadrature data — h needs a finite pole or a polynomial part");
      }
      return { unbounded: true, alpha, c };
    },
    evalPhi: evalPhi_UPQD,
    phiTaylorAt: phiTaylorAt_UPQD,
    computeTargets(phi, hData) {
      return { A: computeTargetA_UPQD(phi, hData), F: laurentMatchAtInfinity_UPQD(phi, hData) };
    },
    residual: residual_UPQD,
    packPhi: packPhi_UPQD,
    unpackPhi: unpackPhi_UPQD,
    canonicalizePhi: canonicalizePhi_UPQD,
    initialGuess: initialGuess_UPQD,
    perturbedInitialGuess: perturbedInitialGuess_UPQD,
    diverseInitialGuess: diverseInitialGuess_UPQD,
    // Continuation in α from the classical limit (a small conformal radius
    // blows the unbounded seed up). See QD.PqdCommon.continuationInAlpha.
    continuationSolve(hData, norm, options = {}) {
      return QD.PqdCommon.continuationInAlpha(hData, norm, options);
    },
    verifyQuadratureIdentity: verifyQuadratureIdentity_UPQD,
    sampleBoundary: sampleBoundary_UPQD,
  };
  QD.registerFamily('unboundedPQD');

  QD.evalRHash_UPQD     = evalRHash_UPQD;
  QD.rHashTaylorAt_UPQD = rHashTaylorAt_UPQD;

})();
