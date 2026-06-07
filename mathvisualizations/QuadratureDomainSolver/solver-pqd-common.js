// =============================================================================
// solver-pqd-common.js -- Shared machinery for the four power-weighted QD
// (PQD) families. Mirrors solver-lqd-common.js (QD.LqdCommon).
//
// The four PQD families
//   • powerQD              (bounded,   0∉Ω)   φ = (R#)^{1/α}
//   • powerQD_singular     (bounded,   0∈Ω)   φ = b_{z₀}·(R#)^{1/α}
//   • unboundedPQD         (unbounded, 0∉Ω)   φ = z·(r#)^{1/α}
//   • unboundedPQD_singular(unbounded, 0∈Ω)   φ = z·b_{z₀}·(r#)^{1/α}
// differ ONLY in (a) the rational kernel r#/R# and its Taylor (genuinely
// per-family — const/Laurent/Blaschke terms diverge, NOT shared here) and
// (b) a small prefactor + chain-rule term in the boundary sweep. Everything
// else — the continuous-arg phase unwrap, the αth-root, the deviation-refined
// boundary sampler, the R#-vanishing guard, and the identity-verifier's
// per-sample weighted-LHS accumulation + RHS residue sums — is identical and
// lives here. Each family supplies the divergent bits as small closures.
//
// API (QD.PqdCommon):
//   rHashVanishingGuard(samples)                       — winding/zero detection
//   sweepUnitCircle(phi, N, rHashTaylorFn, combine)    — continuous-arg ∂𝔻 sweep
//   boundaryMid(thMid, leftPt, phi, evalRHashFn, prefactorFn)  — refiner midpoint
//   sampleBoundaryViaSweep(phi, base, maxExtra, sweepFn, midFn, extraFallback)
//   accumulateWeightedLHS(samples, alpha, fOfW, skipNearZeroW2)  — Σ (unscaled)
//   boundedMonomialRHS(hData, k)                       — Σ residues for f=w^k
//   unboundedTestPointRHS(hData, k, b)                 — Σ residues for f=1/(w−b)^k
//
// `combine` and `prefactorFn` carry each family's exact prefactor/chain-rule
// expressions verbatim, so the migration is bit-identical to the prior inline
// copies (machine-precision identity unchanged).
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family) {
    throw new Error("solver-pqd-common.js: solver.js must be loaded first");
  }

  // ===========================================================================
  // R# non-vanishing guard (Q1.4). φ = (…)·(R#)^{1/α} is a single-valued
  // univalent map only if R# is non-vanishing on 𝔻̄: by the argument principle
  // the winding number of R#(e^{iθ}) about 0 equals the count of R# zeros
  // inside 𝔻, so winding ≠ 0 ⇒ a zero inside, and |R#| → 0 on a sample ⇒ a
  // zero touching ∂𝔻. Either case makes the αth root multi-valued and the
  // sampled "boundary" invalid; the verifier forces the identity to fail so
  // the spurious candidate is rejected. Used by the two BOUNDED PQD verifiers
  // (the unbounded r# is non-vanishing on 𝔻̄* by construction). One cheap pass
  // over the existing samples' stored `rH` values.
  // ===========================================================================
  function rHashVanishingGuard(samples) {
    const N = samples.length;
    if (N < 2) return { vanishes: false, winding: 0 };
    let netDArg = 0;
    let minMag2 = Infinity;
    let prevRaw = Math.atan2(samples[0].rH.im, samples[0].rH.re);
    for (let n = 1; n <= N; n++) {
      const rH = samples[n % N].rH;        // n === N closes the loop back to θ=0
      const m2 = rH.re * rH.re + rH.im * rH.im;
      if (m2 < minMag2) minMag2 = m2;
      const raw = Math.atan2(rH.im, rH.re);
      let d = raw - prevRaw;
      while (d >  Math.PI) d -= 2 * Math.PI;
      while (d < -Math.PI) d += 2 * Math.PI;
      netDArg += d;
      prevRaw = raw;
    }
    const winding = Math.round(netDArg / (2 * Math.PI));
    // |R#| ≈ 0 on a sample ⇒ R# touches the origin on ∂𝔻 (1e-12 in magnitude,
    // i.e. 1e-24 in |R#|²) — essentially a boundary zero.
    const vanishes = (winding !== 0) || (minMag2 < 1e-24);
    return { vanishes, winding };
  }

  // ===========================================================================
  // Anchored continuous αth-root branch.
  // ---------------------------------------------------------------------------
  // R# is analytic and non-vanishing on the (closed) disk where it is evaluated
  // (the rHashVanishingGuard above rejects any Newton iterate that violates
  // this), so a single-valued continuous arg(R#) — hence a single-valued
  // (R#)^{1/α} — exists there REGARDLESS of how far arg(R#) ranges. The
  // principal αth-root (Complex.cpow(·,1/α)) instead clamps arg to (−π,π] and
  // therefore lands on the WRONG sheet whenever α·arg(φ) leaves that interval
  // (e.g. a pole a with |arg a| > π/α): cprincipalRoot(aᵅ) is a rotated copy of
  // a, not a. That is the off-axis-pole bug. The fix is to ANCHOR the branch at
  // a point z* where φ is pinned and continue arg(R#) by unwrapping along a
  // straight segment anchorPt→z.
  //
  //   argContAt(phi, z, evalRHashFn, anchorArg0, anchorPt, K)
  //     → continuous arg(R#(z)), starting from argCont(R#(anchorPt)) = anchorArg0
  //       and unwrapping |Δarg| < π over K samples along anchorPt→z.
  //
  // WHO USES THIS: only the two BOUNDED PQD families (powerQD, powerQD_singular),
  // which anchor at z*=0 with anchorArg0 = α·arg(w0) (so the anchored root at
  // z=0 equals φ(0)=w0). The segment 0→z stays inside 𝔻̄, where R# is guaranteed
  // non-vanishing. The two UNBOUNDED families (φ = z·(r#)^{1/α}) deliberately do
  // NOT use this — the z prefactor carries the boundary winding while r# stays
  // near its positive-real value c^α (arg ≈ 0, inside the principal sector) for
  // ANY pole angle, so the principal root is already correct there (verified: they
  // solve univalently at every angle). Keeping them principal also avoids a
  // per-call continuation walk in the Schwarz fractal hot loop they share.
  //
  // K≈24 is safe in the realizable regime; a mis-track (|Δarg|>π between samples)
  // only happens for invalid candidates, which the identity verifier then rejects.
  // ===========================================================================
  function argContAt(phi, z, evalRHashFn, anchorArg0, anchorPt, K) {
    K = K || 24;
    let cont = anchorArg0;
    const rH0 = evalRHashFn(anchorPt, phi);
    // prevRaw tracks the PRINCIPAL arg for unwrapping; it may differ from
    // `cont` by a multiple of 2π (only the per-step increments matter).
    let prevRaw = Math.atan2(rH0.im, rH0.re);
    for (let i = 1; i <= K; i++) {
      const t = i / K;
      const zeta = {
        re: anchorPt.re + t * (z.re - anchorPt.re),
        im: anchorPt.im + t * (z.im - anchorPt.im),
      };
      const rH = evalRHashFn(zeta, phi);
      const raw = Math.atan2(rH.im, rH.re);
      let d = raw - prevRaw;
      while (d >  Math.PI) d -= 2 * Math.PI;
      while (d < -Math.PI) d += 2 * Math.PI;
      cont += d;
      prevRaw = raw;
    }
    return cont;
  }

  // (R#(z))^{1/α} on the anchored sheet. Returns the αth-root only — callers
  // multiply by any family prefactor (Blaschke b_{z₀}, the z factor, …).
  function phiAnchored(phi, z, evalRHashFn, anchorArg0, anchorPt, K) {
    const rH = evalRHashFn(z, phi);
    const invAlpha = 1 / phi.alpha;
    const mag = Math.pow(rH.re * rH.re + rH.im * rH.im, 0.5 * invAlpha);
    const a = argContAt(phi, z, evalRHashFn, anchorArg0, anchorPt, K) * invAlpha;
    return { re: mag * Math.cos(a), im: mag * Math.sin(a) };
  }

  // ===========================================================================
  // Continuous-arg unit-circle sweep (Q1.3). Walks z = e^{iθ} for N uniform
  // thetas, tracking the UNWRAPPED arg(R#) so the αth root stays on a single
  // continuous sheet (avoids the atan2 branch-cut artifact). The shared driver
  // computes z, R# + R#′ (via rHashTaylorFn order-1), the unwrapped contArg,
  // and root = (R#)^{1/α} on that sheet; then defers to `combine` for the
  // family-specific w and φ′:
  //
  //   combine(ctx) → { w, phiPrime }
  //   ctx = { z, theta, rH, rHp, root, rootMag, rootArg, mag2, contArg,
  //           alpha, invAlpha, phi }
  //
  // The unwrap is correct provided |Δarg(R#)| < π between consecutive samples
  // — true for a uniform grid of N ≥ ~50 on any well-behaved PQD. Every sample
  // stores `rH` for rHashVanishingGuard (harmless extra field for the unbounded
  // families, which don't run the guard).
  //
  // `anchorSpec` (optional) = { evalRHashFn, anchorPt, anchorArg0 }: when given,
  // the θ=0 sample's continuous arg is anchored to the interior value (walking
  // anchorPt→z=1 via argContAt) instead of the principal arg at θ=0. This puts
  // the WHOLE sweep on the same single global sheet as the interior φ(0)=w0
  // normalization, so the boundary is correct for off-axis poles. Omitted ⇒ the
  // legacy principal anchor (callers not yet migrated are unaffected; real-axis
  // poles are bit-identical, since the 0→1 walk reproduces the principal value
  // when arg(R#)=0 along the +real ray).
  // ===========================================================================
  function sweepUnitCircle(phi, N, rHashTaylorFn, combine, anchorSpec) {
    const alpha = phi.alpha;
    const invAlpha = 1 / alpha;
    const samples = new Array(N);
    let contArg = 0;
    let prevRawArg = 0;
    for (let i = 0; i < N; i++) {
      const theta = (2 * Math.PI * i) / N;
      const z = { re: Math.cos(theta), im: Math.sin(theta) };

      // R# and R#′ via the order-1 Taylor at z (closed-form).
      const rT = rHashTaylorFn(z, phi, 1);
      const rH = rT[0];
      const rHp = rT[1];

      const rawArg = Math.atan2(rH.im, rH.re);
      if (i === 0) {
        contArg = anchorSpec
          ? argContAt(phi, z, anchorSpec.evalRHashFn, anchorSpec.anchorArg0, anchorSpec.anchorPt)
          : rawArg;                      // anchored interior sheet, else principal at θ=0
      } else {
        let dArg = rawArg - prevRawArg;
        while (dArg >  Math.PI) dArg -= 2 * Math.PI;
        while (dArg < -Math.PI) dArg += 2 * Math.PI;
        contArg += dArg;
      }
      prevRawArg = rawArg;

      const mag2 = rH.re * rH.re + rH.im * rH.im;
      const rootMag = Math.pow(mag2, 0.5 * invAlpha);
      const rootArg = contArg * invAlpha;
      const root = { re: rootMag * Math.cos(rootArg), im: rootMag * Math.sin(rootArg) };

      const wp = combine({ z, theta, rH, rHp, root, rootMag, rootArg, mag2, contArg, alpha, invAlpha, phi });
      samples[i] = { theta, z, rH, contArg, w: wp.w, phiPrime: wp.phiPrime };
    }
    return samples;
  }

  // ===========================================================================
  // Continuous-arg midpoint evaluator for the curvature-aware refiner (§22).
  // Unwraps arg(R#) from the LEFT neighbour's contArg (|Δarg| < π for small
  // |Δθ|), forms root = (R#)^{1/α} on that sheet, and applies the family
  // prefactor: w = prefactorFn(z, phi) · root, where prefactorFn returns the
  // multiplier (b_{z₀}, z, z·b_{z₀}) or a falsy value for the identity
  // prefactor (powerQD, where w = root).
  // ===========================================================================
  function boundaryMid(thMid, leftPt, phi, evalRHashFn, prefactorFn) {
    const invAlpha = 1 / phi.alpha;
    const z = { re: Math.cos(thMid), im: Math.sin(thMid) };
    const rH = evalRHashFn(z, phi);
    const rawArg = Math.atan2(rH.im, rH.re);
    const refArg = leftPt.contArg;
    const refMod = refArg - 2 * Math.PI * Math.round(refArg / (2 * Math.PI));
    let dArg = rawArg - refMod;
    while (dArg >  Math.PI) dArg -= 2 * Math.PI;
    while (dArg < -Math.PI) dArg += 2 * Math.PI;
    const contArg = refArg + dArg;
    const wMag = Math.pow(rH.re * rH.re + rH.im * rH.im, 0.5 * invAlpha);
    const a = contArg * invAlpha;
    const root = { re: wMag * Math.cos(a), im: wMag * Math.sin(a) };
    const pre = prefactorFn ? prefactorFn(z, phi) : null;
    const w = pre ? Complex.mul(pre, root) : root;
    return { theta: thMid, w, contArg };
  }

  // ===========================================================================
  // Family.X.sampleBoundary hook body. Uniform continuous-arg sweep → coarse
  // polyline (carrying contArg) → curvature-aware deviation refinement
  // (QD.refineBoundaryByDeviation). Returns the `[{theta, w}, …]` shape the
  // generic sampler returns. `extraFallback`: when truthy, the refinement
  // budget defaults to N0·extraFallback if maxExtra is 0 (the singular /
  // unbounded families densify by 3×); when falsy, the budget is exactly
  // maxExtra (powerQD).
  // ===========================================================================
  function sampleBoundaryViaSweep(phi, baseSamples, maxExtra, sweepFn, midFn, extraFallback) {
    const N0 = Math.max(8, baseSamples | 0);
    const swept = sweepFn(phi, N0);
    const pts = swept.map(s => ({ theta: s.theta, w: s.w, contArg: s.contArg }));
    // θ=2π closure: R# winds 0× over ∂𝔻 ⇒ arg (hence contArg) returns to start.
    pts.push({ theta: 2 * Math.PI, w: { ...pts[0].w }, contArg: pts[0].contArg });
    const extra = extraFallback
      ? ((maxExtra | 0) || N0 * extraFallback)
      : (maxExtra | 0);
    return QD.refineBoundaryByDeviation(
      pts, (thMid, leftPt) => midFn(thMid, leftPt, phi),
      { maxPoints: N0 + extra });
  }

  // ===========================================================================
  // Identity-verifier LHS accumulator. Σ_n f(w_n)·(|w_n|²)^{α−1}·conj(w_n)·
  // φ′(z_n)·z_n over the swept boundary samples. The caller multiplies the
  // returned (unscaled) sum by its own ±1/(αN) factor (sign differs bounded vs
  // unbounded). `fOfW(w)` is the test function (w^k bounded, 1/(w−b)^k
  // unbounded). `skipNearZeroW2 > 0` drops samples with |w|² below it (the
  // singular families, where 0 ∈ Ω and the boundary can approach 0); pass 0 to
  // disable (the non-singular families).
  //
  // The weight is the SINGLE-VALUED form (|w|²)^{α−1}·conj(w) of the factor
  // conj(w)^α·w^{α−1} — a real power of the positive real |w|² (no branch cut),
  // valid for arbitrary real α (Complex.pow(·, α) would be multi-valued for
  // non-integer α). On ∂Ω of a non-singular PQD, 0 ∉ Ω̄ so |w| > 0.
  // ===========================================================================
  function accumulateWeightedLHS(samples, alpha, fOfW, skipNearZeroW2) {
    const skip = skipNearZeroW2 || 0;
    let acc = { re: 0, im: 0 };
    for (let n = 0; n < samples.length; n++) {
      const s = samples[n];
      const fVal = fOfW(s.w, n);
      const w2 = s.w.re * s.w.re + s.w.im * s.w.im;
      if (skip && w2 < skip) continue;
      const weight = Math.pow(w2, alpha - 1);
      let term = Complex.mul(fVal, Complex.conj(s.w));
      term = Complex.scale(term, weight);
      term = Complex.mul(term, s.phiPrime);
      term = Complex.mul(term, s.z);
      acc = Complex.add(acc, term);
    }
    return acc;
  }

  // ===========================================================================
  // RHS residue sums.
  // ---------------------------------------------------------------------------
  // Bounded PQDs (test monomials f = w^k): the quadrature RHS is
  //   Σ_j Σ_{s≥1} C_{j,s} · C(k, s−1) · a_j^{k−s+1}   (terms with s−1 > k drop).
  // ===========================================================================
  function boundedMonomialRHS(hData, k) {
    let rhs = { re: 0, im: 0 };
    for (const pole of hData.poles) {
      for (let sIdx = 0; sIdx < pole.principal.length; sIdx++) {
        const s = sIdx + 1;
        if (s - 1 > k) continue;
        const C = pole.principal[sIdx];
        const exp = k - s + 1;
        const aPow = exp === 0 ? { re: 1, im: 0 } : Complex.pow(pole.a, exp);
        const coef = QD.binomialCoeff(k, s - 1);
        rhs = Complex.add(rhs, Complex.scale(Complex.mul(C, aPow), coef));
      }
    }
    return rhs;
  }

  // Unbounded PQDs (test functions f = 1/(w−b)^k for b in the bounded
  // complement K): RHS = residues of f·h at the finite poles a_j plus the
  // pole-at-∞ (polynomial) contribution.
  //   finite: Σ_j Σ_{s≥1} C_{j,s} · (−1)^{s+1} · C(k+s−2, s−1) / (a_j − b)^{k+s−1}
  //   ∞     : −Σ_{l≥k−1} polyPart[l] · C(l, l−k+1) · b^{l−k+1}
  function unboundedTestPointRHS(hData, k, b) {
    let rhs = { re: 0, im: 0 };
    for (const pole of hData.poles) {
      const aMinusB = Complex.sub(pole.a, b);
      for (let sIdx = 0; sIdx < pole.principal.length; sIdx++) {
        const s = sIdx + 1;
        const C = pole.principal[sIdx];
        const sign = (s % 2 === 0) ? -1 : 1;
        const coef = QD.binomialCoeff(k + s - 2, s - 1);
        const denom = Complex.pow(aMinusB, k + s - 1);
        rhs = Complex.add(rhs, Complex.scale(Complex.div(C, denom), sign * coef));
      }
    }
    const polyPart = hData.polyPart || [];
    const m_inf = polyPart.length - 1;
    if (m_inf >= 0) {
      const lMin = Math.max(k - 1, 0);
      for (let l = lMin; l <= m_inf; l++) {
        const expo = l - k + 1;
        const bPow = expo === 0 ? { re: 1, im: 0 } : Complex.pow(b, expo);
        const coef = -QD.binomialCoeff(l, expo);
        rhs = Complex.add(rhs, Complex.scale(Complex.mul(polyPart[l], bPow), coef));
      }
    }
    return rhs;
  }

  // ===========================================================================
  // Realizability diagnostic (α-homotopy fold tracer). DIAGNOSTIC ONLY — it runs
  // many nested Newton solves, so it must NEVER be called on the interactive /
  // hot path; invoke it only AFTER a bounded-PQD solve has failed, to explain
  // WHY: for many quadrature configurations the requested-α PQD simply does not
  // exist. Given fixed quadrature data the bounded-PQD solution is unique (the
  // disk-automorphism gauge is fixed by canonicalizePhi), so
  //   "realizable at α=target"  ⟺  reachable by continuation in α from the
  // classical limit (α=1) up to a FOLD where the univalent branch turns back.
  // We seed at α≈1 (reliably in-basin via the standard pipeline), warm-start-
  // march α toward the target, and report where the branch folds.
  //
  // Returns { realizable, alphaMax, phi, reason } with reason ∈
  //   'realizable'             — reached target α with a univalent map (phi set)
  //   'fold-below-target'      — branch folds at alphaMax strictly between 1 and
  //                              target (no univalent PQD at α=target)
  //   'invalid-even-classical' — no valid QD even at α≈1 (the data is not a
  //                              realizable quadrature domain at all)
  //   'non-univalent'          — reached target α but the map is not univalent
  // alphaMax is the last α on the branch (≈ the fold when < target); it is
  // continuation-resolution-limited, so report it as approximate ("≈").
  //
  // opts: { alpha (target, required), w0?, singular? }. Self-contained — uses
  // QD.solveInverseQD for the α≈1 seed and QD.newtonSolve to march; both are on
  // QD by the time this is ever called.
  // ===========================================================================
  function diagnosePQDRealizability(hData, opts) {
    opts = opts || {};
    const target = opts.alpha;
    if (!(target > 0) || Math.abs(target - 1) < 1e-6) {
      // α ≈ 1 is the classical case — nothing to diagnose here.
      return { realizable: true, alphaMax: target || 1, phi: null, reason: 'realizable' };
    }
    const dir = target > 1 ? 1 : -1;     // march away from α=1 toward the target
    const alpha0 = 1 + dir * 2e-3;

    // (1) Seed near the classical limit via the standard pipeline. If even this
    //     fails, the data is not a realizable QD at all.
    const seedOpts = { alpha: alpha0, findAlternates: false };
    if (opts.w0) seedOpts.w0 = opts.w0;
    if (opts.singular) seedOpts.singular = true;
    let seed = null;
    try { seed = QD.solveInverseQD(hData, seedOpts); } catch (e) { seed = null; }
    if (!seed || !seed.success || !seed.primary || !seed.primary.phi) {
      return { realizable: false, alphaMax: 1, phi: null, reason: 'invalid-even-classical' };
    }

    // (2) March α → target, warm-starting from the previous solution.
    let phi = seed.primary.phi;
    let alpha = alpha0;
    let step = 0.03;
    const minStep = 1e-5, growth = 1.25, capStep = 0.12, maxSteps = 400;
    for (let i = 0; i < maxSteps && Math.abs(alpha - target) > 1e-9; i++) {
      const nextA = dir > 0 ? Math.min(target, alpha + step) : Math.max(target, alpha - step);
      const trial = QD.clonePhi(phi); trial.alpha = nextA;
      const res = QD.newtonSolve(trial, hData, { maxIter: 160 });
      if (res.success) { phi = res.phi; alpha = nextA; step = Math.min(step * growth, capStep); }
      else { step *= 0.5; if (step < minStep) break; }
    }

    // (3a) Reached the target α — confirm univalence (a valid PQD must be univalent).
    if (Math.abs(alpha - target) < 1e-9) {
      let univ = false;
      try { univ = QD.isBoundaryUnivalent(phi); } catch (e) { univ = false; }
      return univ
        ? { realizable: true, alphaMax: target, phi, reason: 'realizable' }
        : { realizable: false, alphaMax: target, phi, reason: 'non-univalent' };
    }

    // (3b) Folded before the target. Tighten the reported fold by bisection
    //      between the last good α and a known-failing α toward the target.
    let aLo = alpha;
    let aHi = dir > 0 ? Math.min(target, alpha + 0.1) : Math.max(target, alpha - 0.1);
    if (Math.abs(aHi - aLo) > 1e-6) {
      for (let b = 0; b < 12; b++) {
        const aMid = 0.5 * (aLo + aHi);
        const trial = QD.clonePhi(phi); trial.alpha = aMid;
        const res = QD.newtonSolve(trial, hData, { maxIter: 200 });
        if (res.success) { phi = res.phi; aLo = aMid; } else { aHi = aMid; }
      }
    }
    return { realizable: false, alphaMax: aLo, phi, reason: 'fold-below-target' };
  }

  // ===========================================================================
  // Continuation in α — the PQD continuation homotopy (shared by the three
  // families whose `continuationSolve` was previously a stub).
  //
  // WHY α and not residue-strength / continuation-in-c: those degenerate here.
  // Shrinking the residues toward 0 pulls Ω back to the disk around w₀, which
  // pushes the ORIGIN out of Ω — breaking the 0∈Ω ansatz of the SINGULAR
  // families; and a small conformal radius makes the unbounded seed blow up
  // (z_j = a_j/c → ∞). Continuing in the weight exponent α from the classical
  // limit avoids both: solve the near-classical α≈1 problem with the standard
  // pipeline (reliably in-basin there — same premise as diagnosePQDRealizability)
  // then warm-start-march α to the target. Returns the standard
  // { success, phi, residual, trace, method, error? } envelope and does NOT
  // canonicalize (the pipeline's evalCandidate does). On any failure the
  // pipeline falls through to multistart — today's behavior — so this can only
  // help or no-op.
  //
  // RECURSION GUARD: the α≈1 seed solve runs with usePhases.continuation:false,
  // so it can never re-enter this function.
  function continuationInAlpha(hData, norm, options = {}) {
    const target = norm.alpha;
    if (!(target > 0)) return { success: false, error: "continuation-in-α: no target α", trace: [] };
    const { growFactor = 1.6, shrinkFactor = 0.5, minStep = 1e-3, maxSteps = 80, newton = {} } = options;
    const dir = target >= 1 ? 1 : -1;
    // Start just off the classical limit, on the same side as the target.
    const alphaStart = 1 + dir * Math.min(0.05, Math.max(1e-3, Math.abs(target - 1) * 0.5));
    // Near-classical seed via the standard pipeline, with continuation DISABLED
    // (recursion guard). diverse on / deflation off mirrors a normal cold solve.
    const seedOpts = {
      alpha: alphaStart,
      unbounded: !!norm.unbounded, singular: !!norm.singular,
      c: norm.c, w0: norm.w0, q: norm.q,
      usePhases: { direct: true, continuation: false, multistart: true, diverse: true, deflation: false },
    };
    const seedR = QD.solveInverseQD(hData, seedOpts);
    if (!seedR.success || !seedR.primary) {
      return { success: false, error: "continuation-in-α: seed solve at α=" + alphaStart.toFixed(3) + " failed", trace: [] };
    }
    let phi = seedR.primary.phi;
    const trace = [{ alpha: alphaStart, ok: true, residual: seedR.primary.residual }];
    let last = alphaStart;
    let step = Math.max(Math.abs(target - alphaStart) * 0.25, minStep);
    for (let i = 0; i < maxSteps; i++) {
      if (Math.abs(last - target) < 1e-12) break;
      const nextA = dir > 0 ? Math.min(target, last + step) : Math.max(target, last - step);
      const trial = QD.clonePhi(phi); trial.alpha = nextA;
      const ns = QD.newtonSolve(trial, hData, newton);
      if (ns.success) {
        phi = ns.phi; last = nextA;
        trace.push({ alpha: nextA, ok: true, residual: ns.residual });
        step *= growFactor;
      } else {
        step *= shrinkFactor;
        trace.push({ alpha: nextA, ok: false, residual: ns.residual ?? null });
        if (step < minStep) {
          return { success: false, error: "continuation-in-α: step underflow at α=" + last.toFixed(4), phi, trace, lastAlpha: last };
        }
      }
    }
    if (Math.abs(last - target) > 1e-9) {
      return { success: false, error: "continuation-in-α: max steps reached at α=" + last.toFixed(4), phi, trace, lastAlpha: last };
    }
    // The α-march can cross a FOLD onto a branch that satisfies the algebraic
    // residual but is NOT a univalent QD (a spurious map). Verify the endpoint
    // is a genuine univalent QD; if not, report failure so the pipeline falls
    // through to multistart instead of returning — or letting the caller pick —
    // a spurious solution. This is what keeps continuation strictly "help or
    // no-op" (it never preempts multistart with a worse candidate).
    const fam = QD.selectFamily(norm);
    // Match the main pipeline's univalence resolution (default 500) so this gate
    // agrees with evalCandidate's verdict — a lower N can false-reject a
    // genuinely univalent unbounded boundary.
    const uN = options.univalenceSamples || 500;
    const univalent = QD.isBoundaryUnivalent(phi, uN);
    const idOK = fam.verifyQuadratureIdentity(phi, hData, { numSamples: uN }).maxRelDiff < (options.identityTol || 1e-6);
    if (!univalent || !idOK) {
      return { success: false, error: "continuation-in-α: endpoint not a univalent QD at α=" + target, phi, trace, lastAlpha: last };
    }
    return { success: true, phi, iterations: 0, residual: trace[trace.length - 1].residual, trace, method: "continuation-in-alpha" };
  }

  // ===========================================================================
  // Expose under QD.PqdCommon (+ back-compat alias QD._rHashVanishingGuard,
  // referenced by solver-pqd-singular.js and node-test.js).
  // ===========================================================================
  QD.PqdCommon = {
    rHashVanishingGuard,
    continuationInAlpha,
    argContAt,
    phiAnchored,
    sweepUnitCircle,
    boundaryMid,
    sampleBoundaryViaSweep,
    accumulateWeightedLHS,
    boundedMonomialRHS,
    unboundedTestPointRHS,
    diagnosePQDRealizability,
  };
  QD._rHashVanishingGuard = rHashVanishingGuard;

})();
