// =============================================================================
// solver-uqd-lqd-singular.js -- Unbounded SINGULAR log-weighted QDs
//                                (Family.unboundedLQD_singular)
//
// Setting: Ω unbounded simply connected with 0 ∈ Ω AND ∞ ∈ Ω. Identity:
//     ∫_Ω f(w)/|w|² dA = ∮_∂Ω f(w) h(w) dw,   f ∈ L¹_a(Ω; ρ₀)
// where the test class is analytic in Ω with f(0) = 0 AND f(∞) = 0 (the
// LQD-singular analog of A₀_*(Ω)).
//
// Riemann-map parametrization (Eq. 5.6, singular form):
//
//     φ(z) = c · |z₀| · z · b_{z₀}(z) · exp( r#(z) − r#(∞) ),   z ∈ 𝔻*,
//
// with:
//   • c > 0 user input, the conformal radius:  φ'(∞) = c exactly via the
//     (− r#(∞)) absorption (mirroring unboundedLQD)
//   • z₀ ∈ 𝔻* unknown, the preimage of 0
//   • b_{z₀} thesis form (shared with bounded singular via LqdCommon)
//   • r# rational, same form as other LQDs
//
// q-equation (●₀):  q (residue of h at 0) is user input. Derivation parallels
// the bounded-singular case:
//     G(z) := ln(φ · φ#)(z) = ln(c²|z₀|²) + r̃#(z) + r(z),
//     r̃# := r# − r#(∞),  r(z) := conj(r̃#(1/conj(z))),
// using the Blaschke identity b·b# ≡ 1 (holds for z₀ ∈ 𝔻 OR 𝔻*). Then
//     q = Res_{w=0} S₀(w) = G(z₀) = ln(c²|z₀|²) + r̃#(z₀) + conj(r̃#(1/conj(z₀))).
//
// Unknowns / equations — square system:
//   Unknowns:  {z_j} (2n) + {A_{j,k}} (2M) + z₀ (2)
//   Equations: (●) 2n + (★) 2M + (●₀) 2
//
// Identity test class:  f(w) = w/(w−b)^k for k = 2, 3, 4, with b ∈ K (the
// bounded complement). f vanishes at 0 (since b ≠ 0) and at ∞ (k ≥ 2). K may
// be non-convex, so candidate b's are ray-cast against the boundary polygon
// and rejected if not actually inside K.
//
// Polynomial-h support (HANDOFF #22):
//   The (●₀) q-equation is the FULL Andrew Graven q-formula
//     q = ln(c²|z₀|²) + R(z₀) + R#(z₀)
//   where R(z) := r̃#(z) + B(1/z) is the full exponent in φ (including the
//   polynomial-h β-correction via B(1/z) = Σ_{l=1..N} β_l/z^l) and
//   R#(z) := conj(R(1/conj(z))) is its Schwarz reflection. Expanding:
//     q = ln(c²|z₀|²)
//       + r̃#(z₀) + conj(r̃#(1/conj(z₀)))             (existing β = 0 case)
//       + B(1/z₀) + conj(B(conj(z₀)))                 (NEW β-correction)
//   This is derived from S₀(w) = ln(φ·φ#)(ψ(w))/w (the logarithmic
//   generalized Schwarz function), with the Blaschke identity b·b# ≡ 1
//   making the numerator finite at w = 0. The (★)_F equations match β to
//   h's polyPart via β_l = F̃_l with F̃ = inverseFaberAtInfinity(P̃, f̃, c),
//   same form as the non-singular family but using the Blaschke-aware
//   φ-Laurent helper phiLaurentAtInfinity_UQDLS in LqdCommon.
//
// Case (a) — no finite poles + polyPart (HANDOFF #23, shipped): with the
// HANDOFF #22 polynomial-h support, the `h = q/w` rejection now applies
// only when polyPart is ALSO empty. polyPart provides the structure that
// pins φ when there are no finite-pole landmarks.
//
// Case (b) — higher-order pole at 0 in h (HANDOFF #24, shipped):
//   An order-(m₀+1) pole at w = 0 in h is parametrized by a SYNTHETIC
//   BRANCH added to r̃# anchored at z = z₀ (the Blaschke zero, NOT z = 0):
//      r̃#_syn(z) = Σ_{l=1..m₀} conj(c_l) · z^l / (1 − conj(z₀) · z)^l
//   Same Möbius form as existing finite-pole branches at z_j ∈ 𝔻*, just
//   anchored at z₀. Has a pole at z = 1/conj(z₀) of order m₀. Under the
//   Schwarz reflection R#(z) = conj(R(1/conj(z))) the pole reflects to
//   z = z₀, and since ψ(0) = z₀ this propagates to a genuine order-m₀
//   pole of R#(ψ(w)) at w = 0 — giving S₀(w) an order-(m₀+1) pole at
//   w = 0 with the literal `q_k = coefficient of 1/w^k in S₀(w)` formula.
//
//   Data convention: hData.poles holds a single entry with
//   `a = {re:0,im:0}` and `principal = [q_2, …, q_{m₀+1}]` (length m₀).
//   The simple-pole residue q_1 stays in `opts.q` / `phi.q` as before.
//
//   Implementation: phi.lqdGamma = [c_1, …, c_{m₀}] is the new Newton-
//   vector slot. The `_phiWithSyntheticBranch` helper merges γ into
//   phi.branches as a `{ z: phi.z0, A: phi.lqdGamma }` entry on the fly
//   for r̃# evaluations, so the existing Möbius-form r̃# machinery in
//   LqdCommon picks it up transparently. The new (★)_A block in
//   residual_UQDLS matches phi.lqdGamma to principal directly via
//   `inverseFaberAtPole(principal, phiTilde_at_z0)`. The (●₀) q-equation
//   gains two closed-form γ-correction sums (synth eval at z₀, plus the
//   regularized analytic part of synth's Schwarz reflection at 1/conj(z₀)).
//
//   See HANDOFF.md §10 / entry #24 for the full derivation and history.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family || !QD.LqdCommon) {
    throw new Error("solver-uqd-lqd-singular.js: solver.js + solver-lqd-common.js must be loaded first");
  }

  const blaschkeEval        = QD.LqdCommon.blaschkeEval;
  const blaschkeTaylor      = QD.LqdCommon.blaschkeTaylor;
  const evalRHash           = QD.LqdCommon.evalRHash;
  const rHashTaylorAt       = QD.LqdCommon.rHashTaylorAt;
  const rHashAtInfinity     = QD.LqdCommon.rHashAtInfinity;
  const computeFaberTargetA = QD.LqdCommon.computeFaberTargetA;
  const perturbBranchesInPlace = QD.LqdCommon.perturbBranchesInPlace;
  const diverseSeedBranches    = QD.LqdCommon.diverseSeedBranches;
  const sampleBoundaryWithDerivative = QD.LqdCommon.sampleBoundaryWithDerivative;
  const verifyIdentityGeneric  = QD.LqdCommon.verifyIdentityGeneric;
  const residueScaleFloor      = QD.LqdCommon.residueScaleFloor;

  // ===========================================================================
  // 1. φ evaluation: c · |z₀| · z · b_{z₀}(z) · exp(r̃#(z) + B(1/z))
  // ---------------------------------------------------------------------------
  // B(1/z) = Σ_{l=1..N} β_l / z^l is the polynomial-h extension (HANDOFF #22).
  // The synthetic γ-branch (HANDOFF #24, case (b) higher-order pole at 0)
  // is folded INTO r̃# via _phiWithSyntheticBranch: when phi.lqdGamma is
  // nonempty, evalRHash/rHashAtInfinity see one extra Möbius branch at
  // z = z₀ contributing r̃#_syn(z) = Σ conj(c_l)·z^l/(1−conj(z₀)·z)^l.
  // ===========================================================================
  function evalPhi_UQDLS(z, phi) {
    const phiX = _phiWithSyntheticBranch(phi);
    const r = evalRHash(z, phiX);
    const rInf = rHashAtInfinity(phiX);
    const bPart = evalB_OverZ(phi, z);          // B(1/z), zero when lqdBeta empty
    const rEff = Complex.add(Complex.sub(r, rInf), bPart);
    const ea = Math.exp(rEff.re);
    const expR = { re: ea * Math.cos(rEff.im), im: ea * Math.sin(rEff.im) };
    const b = blaschkeEval(z, phi.z0);
    const absZ0 = Complex.abs(phi.z0);
    // c · |z_0| · z · b · expR
    const scale = Complex.scale(z, phi.c * absZ0);
    return Complex.mul(Complex.mul(scale, b), expR);
  }

  // Merge phi.lqdGamma as a synthetic branch anchored at z₀, so the
  // existing Möbius-form r̃# machinery picks up the new contribution
  // transparently. No-op when lqdGamma is empty.
  //
  // Memoized via a `_cachedSyntheticMerge` property stashed on the input
  // phi (HANDOFF #27). One residual evaluation makes O(n) calls into
  // evalPhi_UQDLS / phiTaylorAt_UQDLS / computeTarget*, all of which
  // would otherwise rebuild the merged branches list each call. The cache
  // is naturally invalidated when Newton rebuilds phi from a new packed
  // vector each iteration (a fresh phi object → no cache).
  function _phiWithSyntheticBranch(phi) {
    const gamma = phi.lqdGamma || [];
    if (gamma.length === 0) return phi;
    if (phi._cachedSyntheticMerge) return phi._cachedSyntheticMerge;
    const merged = Object.assign({}, phi, {
      branches: phi.branches.concat([{ z: phi.z0, A: gamma }]),
    });
    phi._cachedSyntheticMerge = merged;
    return merged;
  }

  // B(1/z) = Σ_l β_l / z^l. Imported from LqdCommon (HANDOFF #27 dedupe).
  const evalB_OverZ = QD.LqdCommon.evalB_OverZ;

  // ===========================================================================
  // 2. φ-Taylor at z = z_c
  // ===========================================================================
  // φ(z_c + t) = c·|z₀|·(z_c+t)·b_{z₀}(z_c+t)·exp(r̃#(z_c+t) + B(1/(z_c+t)))
  //            = K · lin(t) · bT(t) · exp(argTilde(t))
  // with
  //   K          = c·|z₀|·exp(r̃#(z_c) + B(1/z_c))
  //   lin(t)     = z_c + t                                (Taylor [z_c, 1, 0, ...])
  //   bT(t)      = b_{z₀}(z_c + t)                        (blaschkeTaylor)
  //   argTilde   = rTilde(t) + bTilde(t)                  (both with zero constant)
  //
  // The synthetic γ-branch (HANDOFF #24) is merged into phi.branches via
  // _phiWithSyntheticBranch so rHashTaylorAt/rHashAtInfinity see it
  // automatically. Caveat: when zc = z₀ AND lqdGamma ≠ [], the merged
  // branch contributes a Taylor expansion of its own form around z₀ — its
  // pole at 1/conj(z₀) is well separated from z₀ for |z₀| > 1, so the
  // Taylor is finite. The (★)_A synthetic-branch matching equations then
  // pin lqdGamma against the principal at z₀ from h.
  function phiTaylorAt_UQDLS(zc, phi, L) {
    const phiX = _phiWithSyntheticBranch(phi);
    const rT = rHashTaylorAt(zc, phiX, L);
    const rInf = rHashAtInfinity(phiX);
    const bT_OverZ = bOverZTaylorAt(phi, zc, L);          // Taylor of B(1/z) at z_c
    const b0 = bT_OverZ[0];
    const r0minusInf = Complex.sub(rT[0], rInf);
    const expArgConst = Complex.add(r0minusInf, b0);

    const argTilde = Taylor.zero(L + 1);
    for (let l = 1; l <= L; l++) {
      argTilde[l] = Complex.add(rT[l], bT_OverZ[l]);
    }
    const expArgTilde = Taylor.exp(argTilde, L);

    const ea = Math.exp(expArgConst.re);
    const expConst = { re: ea * Math.cos(expArgConst.im), im: ea * Math.sin(expArgConst.im) };
    const absZ0 = Complex.abs(phi.z0);
    const K = Complex.scale(expConst, phi.c * absZ0);

    const lin = Taylor.zero(L + 1);
    lin[0] = Complex.clone(zc);
    if (L >= 1) lin[1] = { re: 1, im: 0 };

    const bT = blaschkeTaylor(zc, phi.z0, L);

    const t1 = Taylor.mul(lin, bT, L);
    const t2 = Taylor.mul(t1, expArgTilde, L);
    const out = new Array(L + 1);
    for (let l = 0; l <= L; l++) out[l] = Complex.mul(K, t2[l]);
    return out;
  }

  // Taylor expansion of B(1/z) at z = z_c. Imported from LqdCommon
  // (HANDOFF #27 dedupe).
  const bOverZTaylorAt = QD.LqdCommon.bOverZTaylorAt;

  // Filter out the a = 0 entry (higher-order pole at origin) from
  // hData.poles for the standard finite-pole loops. The a = 0 entry is
  // handled separately via the synthetic-branch (★)_A block.
  function _finitePolesView(hData) {
    const polesAll = hData.poles || [];
    let hasA0 = false;
    const finite = [];
    for (const p of polesAll) {
      if (Complex.abs2(p.a) < QD.ZERO_THRESHOLD) hasA0 = true;
      else finite.push(p);
    }
    return hasA0
      ? Object.assign({}, hData, { poles: finite })
      : hData;
  }

  // Locate the a = 0 entry if present, else null.
  function _findA0Pole(hData) {
    for (const p of hData.poles || []) {
      if (Complex.abs2(p.a) < QD.ZERO_THRESHOLD) return p;
    }
    return null;
  }

  function computeTargetA_UQDLS(phi, hData) {
    return computeFaberTargetA(phi, _finitePolesView(hData), phiTaylorAt_UQDLS);
  }

  // (Earlier attempts at the higher-order q-equations were removed in
  // HANDOFF #23 attempt-1 cleanup; see HANDOFF #24 for the correct
  // synthetic-branch approach now in use.)

  // ===========================================================================
  // (★)_F  Polynomial-h target for β (UQDLS).
  // ---------------------------------------------------------------------------
  // Same algebraic structure as computeTargetF_UQDL (solver-uqd-lqd.js
  // §4b) but using the Blaschke-aware φ-Laurent helper. The augmented
  // polynomial P̃ is identical to the non-singular case — the s=1 Σⱼ
  // C_{j,1} constants from each finite pole accumulate alongside h's
  // shifted polyPart. (Note: the q/w pole at the origin contributes the
  // constant q to P̃_0 as well, but P̃_0 only feeds the discarded F̃_0
  // output of inverseFaberAtInfinity, so we omit it here too — same
  // reason it's omitted from the UQDL implementation.)
  // ===========================================================================
  function computeTargetF_UQDLS(phi, hData) {
    const polyPart = hData.polyPart || [];
    const N = polyPart.length;
    if (N === 0) return [];

    // Sum simple-pole C_{j,1} contributions from FINITE poles only. The
    // a = 0 entry's principal[0] is q_2 (not a simple-pole residue) and
    // doesn't belong in this accumulator. P̃_0 only feeds the discarded
    // F̃_0 output anyway, so this is mostly a hygiene fix.
    const finiteHData = _finitePolesView(hData);
    let sumRe = 0, sumIm = 0;
    for (const pole of finiteHData.poles) {
      if (pole.principal.length > 0) {
        sumRe += pole.principal[0].re;
        sumIm += pole.principal[0].im;
      }
    }
    const Ptilde = new Array(N + 1);
    Ptilde[0] = { re: sumRe, im: sumIm };
    for (let i = 0; i < N; i++) Ptilde[i + 1] = Complex.clone(polyPart[i]);

    // Merge γ into phi.branches so rHashLaurentAtInfinity picks up the
    // synthetic-branch contribution at ∞.
    const phiX = _phiWithSyntheticBranch(phi);
    const fTilde = QD.LqdCommon.phiLaurentAtInfinity_UQDLS(phiX, N);
    const Ftilde = QD.Faber.inverseFaberAtInfinity(Ptilde, fTilde, phi.c);

    const targets = new Array(N);
    for (let l = 1; l <= N; l++) {
      targets[l - 1] = Complex.clone(Ftilde[l]);
    }
    return targets;
  }

  // (★)_A target for the synthetic γ-branch at z = z₀: simply
  // inverseFaberAtPole(principal, phiTilde_at_z0), Option A confirmed by
  // Andrew (no modified-residue shift at a = 0).
  function computeTargetGamma_UQDLS(phi, hData) {
    const a0Pole = _findA0Pole(hData);
    if (!a0Pole || a0Pole.principal.length === 0) return [];
    const m0 = a0Pole.principal.length;
    const phiTilde = phiTaylorAt_UQDLS(phi.z0, phi, m0 + 1);
    phiTilde[0] = { re: 0, im: 0 };          // φ(z₀) = 0 by Blaschke; absorb FP noise
    return QD.Faber.inverseFaberAtPole(a0Pole.principal, phiTilde);
  }

  // B(1/z₀) = Σ_l β_l / z₀^l — the β contribution to the (●₀) q-equation
  // (HANDOFF #22 polynomial-h β-correction).
  function _evalB_OverZ0(phi) {
    const beta = phi.lqdBeta || [];
    if (beta.length === 0) return { re: 0, im: 0 };
    const zInv = Complex.inv(phi.z0);
    let pow = Complex.clone(zInv);
    let acc = { re: 0, im: 0 };
    for (let l = 0; l < beta.length; l++) {
      acc = Complex.add(acc, Complex.mul(beta[l], pow));
      if (l + 1 < beta.length) pow = Complex.mul(pow, zInv);
    }
    return acc;
  }

  // conj(B(conj(z₀))) = Σ conj(β_l) · z₀^l — the Schwarz reflection of B
  // contributing to the (●₀) q-equation.
  function _evalConjBConjZ0(phi) {
    const beta = phi.lqdBeta || [];
    if (beta.length === 0) return { re: 0, im: 0 };
    let pow = Complex.clone(phi.z0);
    let acc = { re: 0, im: 0 };
    for (let l = 0; l < beta.length; l++) {
      acc = Complex.add(acc, Complex.mul(Complex.conj(beta[l]), pow));
      if (l + 1 < beta.length) pow = Complex.mul(pow, phi.z0);
    }
    return acc;
  }

  // γ-correction sum #1 (HANDOFF #24): synthetic-branch r̃#_syn evaluated at
  // z = z₀ contributes  Σ_l conj(c_l) · z₀^l / (1 − |z₀|²)^l  to r̃#(z₀).
  // The base-r̃# (finite-pole branches only) is evaluated through evalRHash
  // on the un-merged phi, so this sum captures the synth's contribution that
  // would otherwise be picked up by evalRHash(z₀, merged-phi) — kept separate
  // here because (●₀) uses evalRHash(z₀, phi) on the BASE phi (no merge)
  // so β-only and γ-only contributions stay distinct.
  function _evalSyntheticAtZ0(phi) {
    const gamma = phi.lqdGamma || [];
    if (gamma.length === 0) return { re: 0, im: 0 };
    const z0 = phi.z0;
    const denom = 1 - Complex.abs2(z0);            // real; nonzero for z0 in 𝔻*
    let acc = { re: 0, im: 0 };
    let z0Pow = { re: 1, im: 0 };
    let denomPow = 1;
    for (let l = 1; l <= gamma.length; l++) {
      z0Pow = Complex.mul(z0Pow, z0);              // z0^l
      denomPow *= denom;                            // (1 - |z0|²)^l (signed; negative for |z0| > 1)
      const term = Complex.scale(
        Complex.mul(Complex.conj(gamma[l - 1]), z0Pow),
        1 / denomPow
      );
      acc = Complex.add(acc, term);
    }
    return acc;
  }

  // γ-correction sum #2 (HANDOFF #24): the analytic part of
  // conj(r̃#_syn(1/conj(z₀))), needed because r̃#_syn has a literal pole at
  // 1/conj(z₀) and the naive `conj(r̃#(1/conj(z₀)))` term in the q-equation
  // diverges. Subtracting the principal part leaves only the analytic
  // constant
  //     Σ_l c_l · (−1)^l / z₀^l
  // (derivation: at z = 1/conj(z₀) + u for small u, the l-th term of r̃#_syn
  // expands to give (−1)^l / conj(z₀)^l · conj(c_l) as its u^0 coefficient;
  // conjugating gives c_l · (−1)^l / z₀^l. The negative-power terms in u are
  // the principal part, dropped in the regularized formula.)
  function _evalRegularizedSyntheticAt1OverConjZ0(phi) {
    const gamma = phi.lqdGamma || [];
    if (gamma.length === 0) return { re: 0, im: 0 };
    const z0 = phi.z0;
    let acc = { re: 0, im: 0 };
    let z0Pow = { re: 1, im: 0 };
    let sign = 1;
    for (let l = 1; l <= gamma.length; l++) {
      z0Pow = Complex.mul(z0Pow, z0);              // z0^l
      sign = -sign;                                  // (-1)^l
      const inv = Complex.inv(z0Pow);
      const term = Complex.scale(Complex.mul(gamma[l - 1], inv), sign);
      acc = Complex.add(acc, term);
    }
    return acc;
  }

  // ===========================================================================
  // 3. Residual: (●) 2n + (★)_A 2M + (●₀) 2 + (★)_F 2N + (★)_Γ 2m₀
  // ---------------------------------------------------------------------------
  // Standard finite-pole loops (●, ★)_A run over the FINITE-poles view of
  // hData (a = 0 entry filtered out — see _finitePolesView). The a = 0
  // entry's m₀ principal residues [q_2, …, q_{m₀+1}] are matched via a
  // separate (★)_Γ block at the bottom using inverseFaberAtPole(principal,
  // phiTilde_at_z0) — same Faber machinery, anchored at the synthetic
  // branch's host point z = z₀.
  //
  // (●₀) q-equation (HANDOFF #22 β-correction + HANDOFF #24 γ-correction):
  //   q_1 = ln(c²|z₀|²)
  //       + r̃#_base(z₀) + conj(r̃#_base(1/conj(z₀)))   existing (β = γ = 0)
  //       + B(1/z₀) + conj(B(conj(z₀)))                 β-correction (HANDOFF #22)
  //       + Σ_l conj(c_l) z₀^l / (1 − |z₀|²)^l          γ #1 (synth at z₀)
  //       + Σ_l c_l (−1)^l / z₀^l                       γ #2 (regularized at 1/conj(z₀))
  // Note: r̃#_base uses evalRHash on the BASE phi (no synthetic-branch
  // merge), since synth's contributions enter through the two closed-form
  // γ-correction sums. This keeps the β-only and γ-only pieces algebraically
  // separated and avoids evaluating evalRHash at 1/conj(z₀) where the
  // synthetic branch has its literal pole.
  //
  // The (★)_F block matches phi.lqdBeta to h's polyPart via β_l = F̃_l using
  // phiLaurentAtInfinity_UQDLS on the merged phi (so the Laurent picks up
  // the synthetic-branch contribution at ∞).
  // ===========================================================================
  function residual_UQDLS(phi, hData, options) {
    options = options || {};
    const out = [];

    const finiteHData = _finitePolesView(hData);
    const finitePoles = finiteHData.poles;

    // (●) locator: φ(z_j) = a_j on FINITE poles (1-to-1 with phi.branches).
    for (let j = 0; j < finitePoles.length; j++) {
      const phiZj = evalPhi_UQDLS(phi.branches[j].z, phi);
      const diff = Complex.sub(phiZj, finitePoles[j].a);
      out.push(diff.re, diff.im);
    }

    // (★)_A modified-residue Faber match on FINITE poles.
    const target = computeTargetA_UQDLS(phi, hData);
    for (let j = 0; j < finitePoles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], target[j][k]);
        out.push(diff.re, diff.im);
      }
    }

    // (●₀) q-equation. Uses BASE phi for evalRHash (γ contributes via
    // the closed-form correction sums, not through evalRHash).
    const rInf = rHashAtInfinity(phi);
    const rZ0 = evalRHash(phi.z0, phi);
    const absZ02 = Complex.abs2(phi.z0);
    const oneOverConjZ0 = Complex.scale(phi.z0, 1 / absZ02);    // 1/conj(z₀)
    const rInvZ0 = evalRHash(oneOverConjZ0, phi);
    const rTildeZ0  = Complex.sub(rZ0, rInf);
    const rTildeInv = Complex.sub(rInvZ0, rInf);
    let sum = Complex.add(rTildeZ0, Complex.conj(rTildeInv));
    sum = Complex.add(sum, _evalB_OverZ0(phi));
    sum = Complex.add(sum, _evalConjBConjZ0(phi));
    sum = Complex.add(sum, _evalSyntheticAtZ0(phi));
    sum = Complex.add(sum, _evalRegularizedSyntheticAt1OverConjZ0(phi));
    const lnFactor = Math.log(phi.c * phi.c * absZ02);
    const lhs = { re: sum.re + lnFactor, im: sum.im };
    const dq1 = Complex.sub(phi.q, lhs);
    out.push(dq1.re, dq1.im);

    // (★)_F polynomial-h block. Inert when phi.lqdBeta is empty.
    if (phi.lqdBeta && phi.lqdBeta.length > 0) {
      const targetF = computeTargetF_UQDLS(phi, hData);
      const Nbeta = phi.lqdBeta.length;
      for (let l = 0; l < Nbeta; l++) {
        const tl = (l < targetF.length) ? targetF[l] : { re: 0, im: 0 };
        const dF = Complex.sub(phi.lqdBeta[l], tl);
        out.push(dF.re, dF.im);
      }
    }

    // (★)_Γ synthetic-branch (★)_A block — matches phi.lqdGamma to the
    // a = 0 entry's principal via inverseFaberAtPole at z = z₀. Inert when
    // phi.lqdGamma is empty (no higher-order pole at origin).
    if (phi.lqdGamma && phi.lqdGamma.length > 0) {
      const targetG = computeTargetGamma_UQDLS(phi, hData);
      const mG = phi.lqdGamma.length;
      for (let l = 0; l < mG; l++) {
        const tl = (l < targetG.length) ? targetG[l] : { re: 0, im: 0 };
        const dG = Complex.sub(phi.lqdGamma[l], tl);
        out.push(dG.re, dG.im);
      }
    }

    return out;
  }

  // c > 0 pins the rotation gauge; no Z/2 ambiguity.
  function canonicalizePhi_UQDLS(phi) { return phi; }

  // ===========================================================================
  // 4. Pack / unpack — schema-driven
  // ===========================================================================
  // Layout: [{z_j}_{j=1..n}, z₀, {A_{j,k}}, {β_l}, {c_l}]. z_j and z₀
  // both clamped to 𝔻*; z₀ additionally bounded above (maxR=1000) to keep
  // Newton out of the deferred z₀ → ∞ degeneracy. β slot (polynomial-h at
  // ∞) collapses to length 0 when template.lqdBeta = [] (finite-pole-only
  // h). γ slot (higher-order pole at origin, HANDOFF #24) collapses to
  // length 0 when template.lqdGamma = [] (no a = 0 entry in hData).
  const SCHEMA_UQDLS = [
    { kind: 'branchesZ', clamp: { side: 'out', cap: QD.DISK_CLAMP_OUT } },
    { kind: 'complex',   name: 'z0', clamp: { side: 'out', cap: QD.DISK_CLAMP_OUT, maxR: QD.Z0_MAX_RADIUS } },
    { kind: 'branchesA' },
    { kind: 'complexList', name: 'lqdBeta'  },
    { kind: 'complexList', name: 'lqdGamma' },
  ];
  function packPhi_UQDLS(phi) { return QD.packPhiBySchema(phi, SCHEMA_UQDLS); }
  function unpackPhi_UQDLS(v, template) {
    const phi = QD.unpackPhiBySchema(v, template, SCHEMA_UQDLS);
    phi.family = 'unboundedLQD_singular';
    phi.unbounded = true;
    phi.c = template.c;
    return phi;
  }

  // ===========================================================================
  // 5. Initial guess — companion unbounded non-singular LQD + argmin |φ|
  // ===========================================================================
  // Seed lqdGamma from the a=0 entry's principal (simple direct copy —
  // motivated by inverseFaberAtPole's leading-order behavior; Newton
  // refines from there).
  function _seedLqdGamma(hData) {
    const a0Pole = _findA0Pole(hData);
    return a0Pole
      ? a0Pole.principal.map(p => ({ re: p.re, im: p.im }))
      : [];
  }

  function initialGuess_UQDLS(hData, norm) {
    const c = norm.c;
    const q = norm.q;

    // Branches correspond 1-to-1 with FINITE poles only.
    const finiteHData = _finitePolesView(hData);
    const finitePoles = finiteHData.poles;

    let zj_guess = null, A_guess = null, z0_guess = null;

    // Try companion bootstrap (only if at least one finite pole).
    if (finitePoles.length > 0) {
      try {
        const companion = QD.solveInverseQD(finiteHData, {
          lqd: true, unbounded: true, c,
          identityTol: 1e-3, autoEscalate: false, findAlternates: false,
        });
        if (companion.success && companion.primary && companion.primary.phi) {
          const phiUQDL = companion.primary.phi;
          zj_guess = phiUQDL.branches.map(br => Complex.clone(br.z));
          A_guess  = phiUQDL.branches.map(br => br.A.map(Complex.clone));

          // z₀ = argmin |φ_UQDL(z)| on |z| = 1.01, pushed slightly outward.
          const ring = 1.01;
          let bestZ = null, bestMag = Infinity;
          for (let i = 0; i < 60; i++) {
            const theta = 2 * Math.PI * i / 60;
            const z = { re: ring * Math.cos(theta), im: ring * Math.sin(theta) };
            const mag = Complex.abs(QD.Family.unboundedLQD.evalPhi(z, phiUQDL));
            if (mag < bestMag) { bestMag = mag; bestZ = z; }
          }
          if (bestZ) z0_guess = Complex.scale(bestZ, 1.05);
        }
      } catch (e) { /* fall through to geometric */ }
    }

    // Geometric fallback (also used when no finite poles).
    if (!zj_guess) {
      zj_guess = finitePoles.map(p => {
        let z = Complex.scale(p.a, 1 / c);
        const r = Complex.abs(z);
        if (r < 1.05) z = Complex.scale(z, 1.05 / Math.max(r, 1e-15));
        return z;
      });
      A_guess = finitePoles.map(p => {
        const D = [];
        for (let s = 0; s < p.principal.length; s++) {
          const aC = Complex.mul(p.a, p.principal[s]);
          const next = (s + 1 < p.principal.length) ? p.principal[s + 1] : { re: 0, im: 0 };
          D.push(Complex.add(aC, next));
        }
        let ck = 1;
        const A = [];
        for (let k = 1; k <= p.principal.length; k++) {
          ck *= c;
          A.push(Complex.scale(D[k - 1], 1 / ck));
        }
        return A;
      });
    }
    if (!z0_guess) z0_guess = { re: 2, im: 0 };

    // Seed lqdBeta from polyPart and lqdGamma from a=0 principal.
    const polyPart = hData.polyPart || [];
    const phiInit = {
      family: 'unboundedLQD_singular',
      unbounded: true,
      c, q: Complex.clone(q),
      z0: z0_guess,
      w0: undefined,
      branches: zj_guess.map((z, j) => ({ z, A: A_guess[j].map(Complex.clone) })),
      lqdBeta:  polyPart.map(() => ({ re: 0, im: 0 })),
      lqdGamma: _seedLqdGamma(hData),
    };
    if (polyPart.length > 0) {
      const targetF = computeTargetF_UQDLS(phiInit, hData);
      phiInit.lqdBeta = targetF.map(c => ({ re: c.re, im: c.im }));
    }
    return phiInit;
  }

  function perturbedInitialGuess_UQDLS(hData, norm, rng, r) {
    const base = initialGuess_UQDLS(hData, norm);
    perturbBranchesInPlace(base.branches, rng, r || 0,
      { side: 'out', zCap: 1.05, zScale: 1.10 });
    // Perturb z₀ too, with the same out-side clamp + upper bound.
    const sigma = 0.15 + 0.25 * (r || 0);
    base.z0 = {
      re: base.z0.re + sigma * (rng() - 0.5),
      im: base.z0.im + sigma * (rng() - 0.5),
    };
    const rz0 = Math.hypot(base.z0.re, base.z0.im);
    if (rz0 < 1.05)    { const s = 1.05 / Math.max(rz0, 1e-15); base.z0.re *= s; base.z0.im *= s; }
    else if (rz0 > 50) { const s = 50   / rz0;                  base.z0.re *= s; base.z0.im *= s; }
    // Perturb lqdGamma multiplicatively (real-axis) + additively (im).
    for (let l = 0; l < base.lqdGamma.length; l++) {
      base.lqdGamma[l] = {
        re: base.lqdGamma[l].re * (1 + sigma * (rng() - 0.5)),
        im: base.lqdGamma[l].im + sigma * (rng() - 0.5),
      };
    }
    return base;
  }

  function diverseInitialGuess_UQDLS(hData, norm, rng, r) {
    const c = norm.c, q = norm.q;
    const mz0 = Math.exp(Math.log(1.05) + rng() * Math.log(30 / 1.05));
    const pz0 = 2 * Math.PI * rng();
    const polyPart = hData.polyPart || [];
    const finiteHData = _finitePolesView(hData);
    const base = {
      family: 'unboundedLQD_singular',
      unbounded: true,
      c, q: Complex.clone(q),
      z0: { re: mz0 * Math.cos(pz0), im: mz0 * Math.sin(pz0) },
      w0: undefined,
      branches: diverseSeedBranches(finiteHData, rng, { zMin: 1.05, zMax: 30 }),
      lqdBeta:  polyPart.map(() => ({ re: 0, im: 0 })),
      lqdGamma: _seedLqdGamma(hData),
    };
    if (polyPart.length > 0) {
      base.lqdBeta = computeTargetF_UQDLS(base, hData).map(c => ({ re: c.re, im: c.im }));
    }
    return base;
  }

  // ===========================================================================
  // 6. Continuation in c
  // ===========================================================================
  function continuationSolve_UQDLS(hData, norm, options) {
    options = options || {};
    const { cStart = null, growFactor = 1.6, shrinkFactor = 0.5,
            minStep = 1e-4, maxSteps = 80, newton = {} } = options;
    const cTarget = norm.c;

    let minA = Infinity;
    for (const p of hData.poles) {
      const m = Complex.abs(p.a);
      if (m > 0 && m < minA) minA = m;
    }
    const startGuess = cStart ?? Math.min(cTarget, isFinite(minA) ? 0.25 * minA : 0.25);
    if (startGuess <= 0) {
      return { success: false, error: "continuationInC (UQDLS): invalid starting c", trace: [] };
    }

    const trace = [];
    let c = startGuess;
    let phi = initialGuess_UQDLS(hData, { c, q: norm.q });

    let warmup;
    while (true) {
      warmup = QD.newtonSolve(phi, hData, newton);
      if (warmup.success) { phi = warmup.phi; break; }
      c *= shrinkFactor;
      if (c < minStep) {
        return {
          success: false,
          error: "continuationInC (UQDLS): warmup failed even at c=" + c.toExponential(2),
          phi: warmup.phi, trace,
        };
      }
      phi = initialGuess_UQDLS(hData, { c, q: norm.q });
    }
    trace.push({ c, ok: true, residual: warmup.residual });

    if (c >= cTarget - 1e-12) {
      return { success: true, phi, iterations: 0, residual: warmup.residual,
               trace, method: "continuation-in-c-uqdls" };
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
            error: "continuationInC (UQDLS): step underflow at c=" + lastSuccessC.toFixed(4),
            phi, trace, lastC: lastSuccessC,
          };
        }
      }
    }
    if (lastSuccessC < cTarget - 1e-9) {
      return {
        success: false,
        error: "continuationInC (UQDLS): max steps reached at c=" + lastSuccessC.toFixed(4),
        phi, trace, lastC: lastSuccessC,
      };
    }
    return { success: true, phi, iterations: 0,
             residual: trace[trace.length - 1].residual,
             trace, method: "continuation-in-c-uqdls" };
  }

  // ===========================================================================
  // 7. Identity verifier — test class w/(w−b)^k via LqdCommon skeleton
  // ===========================================================================
  // Test points b ∈ K verified by ray-casting against the boundary polygon
  // (K may be non-convex, so we can't trust the centroid heuristic alone).
  // q/w pole at 0 contributes 0 to RHS (f·h = q/(w−b)^k near 0 has no
  // residue at 0); q is checked algebraically via the q-equation in
  // residual_UQDLS instead.
  function verifyQuadratureIdentity_UQDLS(phi, hData, options) {
    options = options || {};
    const maxOrder = options.maxDegree ?? 4;
    const minOrder = 2;                                  // f vanishes at ∞ ⇒ k ≥ 2
    const N = options.numSamples ?? 500;
    const minusTwoPiI = { re: 0, im: -2 * Math.PI };

    const samples = sampleBoundaryWithDerivative(phi, N, phiTaylorAt_UQDLS);
    const polygonPts = samples.map(s => s.w);

    // K-bounding-box for offset sizing + a centroid for candidate b's.
    let cx = 0, cy = 0;
    for (const s of samples) { cx += s.w.re; cy += s.w.im; }
    cx /= N; cy /= N;
    let maxDev = 0;
    for (const s of samples) {
      const d = Math.hypot(s.w.re - cx, s.w.im - cy);
      if (d > maxDev) maxDev = d;
    }

    // Ray-cast: inside polygon (= K, since boundary traverses CCW around K
    // for our z = e^{iθ} parametrization of an unbounded Ω).
    const insidePolygon = (x, y) => {
      let cross = 0;
      for (let i = 0; i < polygonPts.length; i++) {
        const j = (i + 1) % polygonPts.length;
        const yi = polygonPts[i].im, yj = polygonPts[j].im;
        if ((yi > y) !== (yj > y)) {
          const t = (y - yi) / (yj - yi);
          if (polygonPts[i].re + t * (polygonPts[j].re - polygonPts[i].re) > x) cross++;
        }
      }
      return (cross % 2) === 1;
    };

    // Minimum distance from a point to the boundary polygon — to ensure
    // candidate b's are well INSIDE K, not just barely inside. If b is too
    // close to ∂Ω, the integrand 1/(w−b)^k near boundary samples blows up
    // and the trapezoidal rule loses accuracy.
    const minDistToBoundary = (x, y) => {
      let m = Infinity;
      for (const p of polygonPts) {
        const d = Math.hypot(p.re - x, p.im - y);
        if (d < m) m = d;
      }
      return m;
    };

    // Generate candidate b's: centroid + offsets at multiple radii. Filter
    // by (i) inside polygon (= inside K), (ii) not at the origin (∈ Ω),
    // (iii) far enough from ∂Ω. Take the FAREST-from-boundary candidates.
    const candidates = [{ re: cx, im: cy }];
    for (const frac of [0.1, 0.2, 0.3, 0.45, 0.6]) {
      for (let i = 0; i < 12; i++) {
        const ang = 2 * Math.PI * i / 12;
        const r = frac * maxDev;
        candidates.push({ re: cx + r * Math.cos(ang), im: cy + r * Math.sin(ang) });
      }
    }
    const ranked = [];
    for (const b of candidates) {
      if (Math.hypot(b.re, b.im) < 1e-3) continue;          // avoid w = 0 ∈ Ω
      if (!insidePolygon(b.re, b.im)) continue;
      const d = minDistToBoundary(b.re, b.im);
      ranked.push({ b, d });
    }
    // Take the 3 candidates that are farthest from ∂Ω.
    ranked.sort((p, q) => q.d - p.d);
    const testPoints = ranked.slice(0, 3).map(r => r.b);
    if (testPoints.length === 0) {
      // K too thin / origin too central / boundary near-degenerate
      return {
        checks: [], maxRelDiff: Infinity, maxAbsDiff: Infinity,
        areaScale: residueScaleFloor(hData, phi.q),
        testPoints: [], maxDeg: maxOrder, numSamples: N,
        unbounded: true, lqdUnboundedSingular: true,
        warning: "could not find test points inside K",
      };
    }

    return verifyIdentityGeneric(phi, hData, options, {
      phiTaylorFn: phiTaylorAt_UQDLS,
      boundaryKernel(w) {
        const absW2 = Complex.abs2(w);
        if (absW2 < 1e-30) return null;
        return Complex.scale(Complex.inv(w), Math.log(absW2));
      },
      buildTestFunctions(phi, hData) {
        const tests = [];
        const finitePoles = (hData.poles || []).filter(
          p => Complex.abs2(p.a) > QD.ZERO_THRESHOLD
        );
        const a0Pole = _findA0Pole(hData);
        for (const b of testPoints) {
          for (let k = minOrder; k <= maxOrder; k++) {
            // f(w) = w/(w−b)^k = 1/(w−b)^{k−1} + b/(w−b)^k.
            // Residue of f·h at finite pole a_j (h ~ Σ_s C/(w−a_j)^s):
            //   C·(−1)^{s−1}·[ binom(k+s−3, s−1)/(a_j−b)^{k+s−2}
            //                + b·binom(k+s−2, s−1)/(a_j−b)^{k+s−1} ]
            // ∮_∂Ω f·h dw = -2πi · Σ_j Res_{a_j}(f·h) by sphere-residue
            // balance (Res_0 + Res_∞ = -Σ_j Res_{a_j}, and ∂Ω oriented CCW
            // around Ω picks up only the in-Ω residues at 0, ∞).
            //
            // The q/w simple pole at 0 contributes 0 to RHS (f·q/w analytic
            // at 0 since b ≠ 0). For higher-order q_l/w^l (l ≥ 2, HANDOFF
            // #24 case (b)): by sphere balance, the Res_0(f·q_l/w^l) is
            // automatically reflected via Res_∞ and Σ_j Res_{a_j}; since
            // f·q_l/w^l has no pole at any a_j ≠ 0, this works out as
            // Res_0(f·q_l/w^l) = -Res_∞(f·q_l/w^l), neither of which
            // appears in the rhsSum_finite formula. The identity therefore
            // holds without an explicit a=0 RHS contribution. The
            // q_l-dependent boundary shape enters LHS automatically via
            // the (●₀) γ-correction + (★)_Γ-pinned φ.
            let rhsSum = { re: 0, im: 0 };
            for (const pole of finitePoles) {
              const aMinusB = Complex.sub(pole.a, b);
              for (let sIdx = 0; sIdx < pole.principal.length; sIdx++) {
                const s = sIdx + 1;
                const C = pole.principal[sIdx];
                const sign = (s % 2 === 0) ? -1 : 1;
                // Term 1: 1/(w−b)^{k−1}
                const coef1 = QD.binomialCoeff(k + s - 3, s - 1);
                if (coef1 !== 0) {
                  const denom1 = Complex.pow(aMinusB, k + s - 2);
                  const term1 = Complex.div(C, denom1);
                  rhsSum = Complex.add(rhsSum, Complex.scale(term1, sign * coef1));
                }
                // Term 2: b/(w−b)^k
                const coef2 = QD.binomialCoeff(k + s - 2, s - 1);
                if (coef2 !== 0) {
                  const denom2 = Complex.pow(aMinusB, k + s - 1);
                  const term2 = Complex.mul(C, Complex.div(b, denom2));
                  rhsSum = Complex.add(rhsSum, Complex.scale(term2, sign * coef2));
                }
              }
            }
            // Higher-order pole at origin (HANDOFF #24 case (b)): the
            // Res_0(f · q_l/w^l) terms enter RHS with a sign / scaling
            // determined by the LQD-singular boundary identity. The
            // formula below was derived in HANDOFF #24 (closed-form from
            // expanding 1/(w−b)^k about w = 0 and reading off the 1/w
            // coefficient of f · q_l/w^l).
            if (a0Pole) {
              const signK = (k % 2 === 0) ? 1 : -1;
              for (let lIdx = 0; lIdx < a0Pole.principal.length; lIdx++) {
                const l = lIdx + 2;                          // q_l = principal[l−2]
                const ql = a0Pole.principal[lIdx];
                const coef = QD.binomialCoeff(k + l - 3, l - 2);
                if (coef === 0) continue;
                const denom = Complex.pow(b, k + l - 2);
                const term = Complex.div(ql, denom);
                rhsSum = Complex.add(rhsSum, Complex.scale(term, signK * coef));
              }
            }
            // polyPart contribution via Res_∞(f · h_polyPart) (HANDOFF #25):
            //   f · polyPart[i] · w^i = polyPart[i] · w^{i+1} / (w−b)^k
            // Laurent at w = ∞: w^{i+1}/(w−b)^k = Σ_n binom(k+n−1, n)·b^n·w^{i+1−k−n}.
            // 1/w coefficient ⇒ i+1−k−n = −1 ⇒ n = i+2−k. So:
            //   Res_∞(f · h_polyPart) = -Σ_{i≥k−2} polyPart[i]·binom(i+1, i+2−k)·b^{i+2−k}
            // The negative sign comes from Res_∞ := -[1/w coefficient].
            // Verified empirically (probe): when at least one finite pole
            // is present, this contribution closes the identity-verifier
            // gap to ~1e-12 (β-γ interaction case (b) and existing case
            // (a) polyPart+finite cases). The polyPart-only-no-finite-
            // poles edge case has a numerical-conditioning issue (very
            // large K geometry; b samples extend to ~|polyPart|/c^k) that
            // is NOT addressed by this contribution alone — see
            // TODO.md → LQD-identity-polyPart and HANDOFF.md §10.
            const polyPart = hData.polyPart || [];
            for (let i = 0; i < polyPart.length; i++) {
              const need = i + 2 - k;
              if (need < 0) continue;
              const coef = QD.binomialCoeff(i + 1, need);
              if (coef === 0) continue;
              const bPow = Complex.pow(b, need);
              const term = Complex.mul(polyPart[i], bPow);
              // rhsSum += -polyPart[i] · binom · b^need  (the Res_∞ value).
              rhsSum = Complex.add(rhsSum, Complex.scale(term, -coef));
            }
            const rhs = Complex.mul(minusTwoPiI, rhsSum);
            tests.push({
              label: 'w/(w-' + Complex.toString(b, 2) + ')^' + k,
              f: (w) => Complex.div(w, Complex.pow(Complex.sub(w, b), k)),
              residueRhs: rhs,
              tag: { k, b },
            });
          }
        }
        return tests;
      },
      resultFlags: {
        unbounded: true, lqdUnboundedSingular: true,
        testPoints, maxDeg: maxOrder,
      },
    });
  }

  // ===========================================================================
  // 8. Register Family.unboundedLQD_singular
  // ===========================================================================
  QD.Family.unboundedLQD_singular = {
    name: 'unboundedLQD_singular',
    enforceInDisk:  false,
    enforceOutDisk: true,
    matches(opts) { return !!(opts && opts.lqd && opts.unbounded && opts.singular); },

    normalizeOpts(opts, hData) {
      const c = opts.c;
      if (typeof c !== 'number' || !(c > 0)) {
        throw new Error("Family.unboundedLQD_singular: opts.c must be a positive number");
      }
      const q = opts.q || { re: 0, im: 0 };
      // Higher-order pole at origin (HANDOFF #24, case (b)): hData may
      // contain an a = 0 entry with principal = [q_2, …, q_{m₀+1}]
      // (length m₀). The simple-pole residue q_1 stays in opts.q.
      // Callers (parser/UI) are responsible for splitting q_1 off into
      // opts.q BEFORE constructing the principal list. We validate here.
      const hasPolyPart = !!(hData.polyPart && hData.polyPart.length > 0);
      const hasFinitePoles = (hData.poles || []).some(p => Complex.abs2(p.a) > QD.ZERO_THRESHOLD);
      const a0WithGamma = (hData.poles || []).some(p =>
        Complex.abs2(p.a) < QD.ZERO_THRESHOLD && p.principal && p.principal.length > 0
      );
      // h = q/w only (no finite poles, no polyPart, no a=0 γ-entry) has
      // no solution (Theorem 5.5.2-style). Otherwise the polyPart / finite
      // poles / γ provide the structure that pins φ.
      if (!hasFinitePoles && !hasPolyPart && !a0WithGamma && Complex.abs2(q) > QD.ZERO_THRESHOLD) {
        throw new Error(
          "Family.unboundedLQD_singular: no unbounded singular LQD exists for h = q/w " +
          "with no finite poles, no polynomial part, and no higher-order pole at 0 " +
          "(add a finite pole, a polyPart term, a higher-order residue at 0, or set q = 0)."
        );
      }
      return { lqd: true, unbounded: true, singular: true, c, q: Complex.clone(q) };
    },

    evalPhi: evalPhi_UQDLS,
    phiTaylorAt: phiTaylorAt_UQDLS,
    computeTargets(phi, hData) {
      return {
        A: computeTargetA_UQDLS(phi, hData),
        F: computeTargetF_UQDLS(phi, hData),
        G: computeTargetGamma_UQDLS(phi, hData),
      };
    },
    residual: residual_UQDLS,
    packPhi:  packPhi_UQDLS,
    unpackPhi: unpackPhi_UQDLS,
    canonicalizePhi: canonicalizePhi_UQDLS,
    initialGuess:          initialGuess_UQDLS,
    perturbedInitialGuess: perturbedInitialGuess_UQDLS,
    diverseInitialGuess:   diverseInitialGuess_UQDLS,
    continuationSolve:     continuationSolve_UQDLS,
    verifyQuadratureIdentity: verifyQuadratureIdentity_UQDLS,
  };
  QD.registerFamily('unboundedLQD_singular');

})();
