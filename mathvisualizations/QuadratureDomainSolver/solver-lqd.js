// =============================================================================
// solver-lqd.js -- Bounded non-singular log-weighted quadrature domains (LQDs)
//
// An LQD is a domain Ω satisfying
//     ∫_Ω f(w) / |w|² dA = ∮_∂Ω f(w) h(w) dw,    f ∈ L¹_a(Ω; ρ₀)
// with respect to the weight ρ₀(w) = |w|⁻², where h is rational.
//
// This file covers the BOUNDED NON-SINGULAR case (0 ∉ Ω̄, Ω bounded). The
// singular case (0 ∈ Ω, charge q parameter) and the unbounded case will be
// added in separate passes — they each get their own Family entry.
//
// Math summary (Chapter V of Graven 2026):
//
// • Riemann-map characterization (Theorem 5.4.1):
//     Ω ∈ QD₀ ⇔ φ_out extends to exp(rational)
//
// • Normalized parametrization (Eq. 5.5, non-singular):
//     φ(z) = w₀ · exp(r#(z)),   r ∈ A₀(𝔻*) rational,  r#(0) = 0
//
//   Matching the existing bounded-QD conjugation convention (poles at
//   1/conj(z_j), residues conjugated under the # reflection):
//     r#(z) = Σⱼ Σₖ conj(A_{j,k}) · z^k / (1 − conj(z_j) z)^k
//   so the unknowns are EXACTLY {z_j ∈ 𝔻, A_{j,k} ∈ ℂ}, same as bounded QD.
//
// • Locator constraints (●):
//     φ(z_j) = a_j ⇔ r#(z_j) = ln(a_j / w₀)         (principal branch)
//
// • Target equations (★), derived from Theorem 5.4.2:
//     A_{j,k}^target = Σ_{s=k..m_j} (s/k) · D_{j,s} · [t^s] ψ̃_j^k(t)
//   where
//     D_{j,s} = a_j · C_{j,s} + C_{j,s+1}    (C_{j, m_j + 1} ≡ 0)
//   are the residues of (w·h(w) + Res_∞ h) at a_j, and ψ̃_j(t) is the local
//   Taylor inverse of φ(z_j + ·) − a_j. STRUCTURALLY IDENTICAL to the bounded
//   QD formula in `computeTargetA` — the only difference is the modified
//   residues D_{j,s} replacing C_{j,s}.
//
// • Identity test (Green's theorem applied to f/|w|² with f ∈ L¹_a(Ω; ρ₀)):
//     ∫_Ω f/|w|² dA = ∮_∂Ω f · ln|w|²/w · dw      (no leading 1/2 with our
//                                                  contour-integral convention)
//
// All extension points hook into the Family interface defined in solver.js.
// This file consumes Complex, Taylor, and selected helpers from solver.js
// (which must be loaded first), and registers Family.boundedLQD.
// =============================================================================

(function () {
  'use strict';

  // Pull in the existing solver helpers we delegate to. In browser context
  // these are on `window.QD`; in Node, they're on the shared vm context.
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family) {
    throw new Error("solver-lqd.js: solver.js must be loaded first (no QD.Family found)");
  }

  // ===========================================================================
  // 1. Phi evaluation
  // ===========================================================================
  // r# evaluation (a rational kernel = "evalPhi with w₀ = 0") and its Taylor
  // expansion are shared with solver-lqd-singular.js via QD.LqdCommon.
  const evalRHash = QD.LqdCommon.evalRHash;
  const rHashTaylorAt = QD.LqdCommon.rHashTaylorAt;

  function evalPhi_LQD(z, phi) {
    const r = evalRHash(z, phi);
    const ea = Math.exp(r.re);
    const expR = { re: ea * Math.cos(r.im), im: ea * Math.sin(r.im) };
    return Complex.mul(phi.w0, expR);
  }

  // Taylor coefficients of φ(z) = w₀ · exp(r#(z)) at z = z0, up to t^L.
  //
  //   φ(z0 + t) = w₀ · exp(r#(z0)) · exp(r#(z0+t) − r#(z0))
  //
  // So compute the Taylor of r# at z0, strip its constant, exp-compose, and
  // multiply by the constant scalar w₀ · exp(r#(z0)).
  function phiTaylorAt_LQD(z0, phi, L) {
    const rT = rHashTaylorAt(z0, phi, L);       // Taylor of r# at z0
    const r0 = rT[0];                           // r#(z0)

    // Strip constant: rTilde = rT − r0
    const rTilde = rT.slice();
    rTilde[0] = { re: 0, im: 0 };

    // exp(rTilde) up to t^L
    const expRTilde = Taylor.exp(rTilde, L);

    // Multiplicative prefactor w₀ · exp(r#(z0))
    const ea = Math.exp(r0.re);
    const expR0 = { re: ea * Math.cos(r0.im), im: ea * Math.sin(r0.im) };
    const scale = Complex.mul(phi.w0, expR0);

    const out = new Array(L + 1);
    for (let l = 0; l <= L; l++) {
      out[l] = Complex.mul(scale, expRTilde[l]);
    }
    return out;
  }

  // ===========================================================================
  // 2. Target coefficients (★) and residual
  // ===========================================================================

  // Modified residues D_{j,s} = a_j · C_{j,s} + C_{j,s+1} and the (★) Faber
  // target loop are family-generic given a phiTaylorAt; both live in
  // QD.LqdCommon and are reused by solver-lqd-singular.js and (future)
  // unbounded LQD families.
  function computeTargetA_LQD(phi, hData) {
    return QD.LqdCommon.computeFaberTargetA(phi, hData, phiTaylorAt_LQD);
  }

  // Residual vector. Same shape as the bounded QD residual:
  //   for j: Re/Im(φ(z_j) − a_j)                  — locator (n complex)
  //   for j,k: Re/Im(A_{j,k} − target_{j,k})      — coefficient (d complex)
  //   + 1 gauge equation (Σ Im(conj(A_{j,1})) = 0; resolved sign in canonicalize)
  function residual_LQD(phi, hData, options) {
    options = options || {};
    const enforceGauge = options.enforceGauge !== false;
    const out = [];

    // (●) locator
    for (let j = 0; j < hData.poles.length; j++) {
      const zj = phi.branches[j].z;
      const phiZj = evalPhi_LQD(zj, phi);
      const diff = Complex.sub(phiZj, hData.poles[j].a);
      out.push(diff.re, diff.im);
    }

    // (★) coefficient
    const target = computeTargetA_LQD(phi, hData);
    for (let j = 0; j < hData.poles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], target[j][k]);
        out.push(diff.re, diff.im);
      }
    }

    // Gauge: fix the disk rotation z → e^{iθ} z by Σ Im(conj(A_{j,1})) = 0.
    // Note we use conj(A) because the parametrization uses conj(A_{j,k}) in
    // the (●) form — so the geometric "first coefficient" is conj(A_{j,1}).
    if (enforceGauge) {
      let imSum = 0;
      for (const br of phi.branches) {
        if (br.A.length > 0) imSum += -br.A[0].im;     // Im(conj(a)) = -Im(a)
      }
      out.push(imSum);
    }

    return out;
  }

  // Canonicalize the gauge sign. The Σ Im(conj(A_{j,1})) = 0 equation in
  // residual_LQD leaves a Z/2 sign ambiguity: rotation by π (z → −z) flips
  // z_j → −z_j and conj(A_{j,k}) → (−1)^k conj(A_{j,k}). We resolve it by
  // requiring Σ Re(conj(A_{j,1})) > 0 (i.e. φ'(0) ∝ Σ conj(A_{j,1}) points
  // along positive real, matching the standard Riemann-map normalization).
  function canonicalizePhi_LQD(phi) {
    let reSum = 0;
    for (const br of phi.branches) {
      if (br.A.length > 0) reSum += br.A[0].re;        // Re(conj(a)) = Re(a)
    }
    if (reSum >= 0) return phi;
    // 180° rotation: z_j → −z_j ; A_{j,k} → (−1)^k A_{j,k}
    //   (we apply the rotation in A-space; with k 0-indexed, k=0 ↔ A_{j,1}
    //    so sign flips; k=1 ↔ A_{j,2} so sign stays; etc.)
    return {
      family: 'boundedLQD',
      w0: Complex.clone(phi.w0),
      branches: phi.branches.map(br => ({
        z: Complex.neg(br.z),
        A: br.A.map((a, k) => (k % 2 === 0 ? Complex.neg(a) : Complex.clone(a))),
      })),
    };
  }

  // ===========================================================================
  // 3. Pack / unpack
  // ===========================================================================
  // Layout (matches bounded QD):
  //   for j = 1..n:                   [Re(z_j),     Im(z_j)]
  //   for j = 1..n: for k = 1..m_j:   [Re(A_{j,k}), Im(A_{j,k})]
  function packPhi_LQD(phi) {
    const v = [];
    for (const br of phi.branches) v.push(br.z.re, br.z.im);
    for (const br of phi.branches) for (const a of br.A) v.push(a.re, a.im);
    return v;
  }

  function unpackPhi_LQD(v, template) {
    const phi = {
      family: 'boundedLQD',
      w0: template.w0 ? Complex.clone(template.w0) : { re: 0, im: 0 },
      branches: [],
    };
    let idx = 0;
    for (let j = 0; j < template.branches.length; j++) {
      const z = { re: v[idx], im: v[idx + 1] };
      idx += 2;
      phi.branches.push({ z, A: [] });
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

  // ===========================================================================
  // 4. Initial guess via Theorem 5.1.1 (exp-of-QD bootstrap)
  // ===========================================================================
  //
  // For any bounded QD K with 0 ∉ K̄, e^K is an LQD. Conversely, for our
  // target LQD problem with poles {a_j} and residues {C_{j,s}}, an excellent
  // initial guess is to:
  //   (1) Translate the problem into log-space: solve a QD with "poles" at
  //       ln(a_j) and a (rough) corresponding h. We don't need to be exact —
  //       any near-disk that maps into a domain near Ω under exp will give
  //       Newton a good warm-start.
  //   (2) Exponentiate: read off the LQD parameters from φ_QD pushed through
  //       exp.
  //
  // For simplicity (and robustness when ln(a_j) has multi-valued issues),
  // we use a DIRECT geometric guess that works well in practice:
  //
  //   • z_j ≈ 0 (so φ(z_j) ≈ w₀; locator residual ≈ |a_j − w₀|, small for
  //     Ω near a disk around w₀). Actually we want z_j to be the preimage of
  //     a_j, so a better guess is z_j ≈ (a_j − w₀)/R where R is a guessed
  //     "radius" of Ω around w₀.
  //   • A_{j,k} chosen so that the local Taylor of φ at z_j roughly matches
  //     the residue structure. For the canonical one-point LQD (Theorem 5.3.2),
  //     φ(z) = w₀ exp(z√α), so A_{1,1} = √α (after matching conventions).
  //     For general h, A_{j,1} ≈ √(C_{j,1}/w₀) gives a sensible scale.
  //
  // This is a "good enough" first cut. If it fails, the multistart / diverse
  // phases of solveInverseQD pick up the slack.
  function initialGuess_LQD(hData, norm) {
    const w0 = norm.w0;
    const n = hData.poles.length;

    // Guess a "radius" of Ω: R such that |ln(a_j/w₀)| ~ R for all j.
    // (For a one-point LQD with charge α, R ~ √α.)
    let totalLog = 0;
    for (const p of hData.poles) {
      // ln(a_j/w₀) — principal branch.
      const ratio = Complex.div(p.a, w0);
      const mag = Math.hypot(Math.log(Complex.abs(ratio)), Math.atan2(ratio.im, ratio.re));
      if (mag > totalLog) totalLog = mag;
    }
    const R = Math.max(totalLog, 0.3);              // floor to keep z_j inside 𝔻
    const cap = 0.85;                               // keep |z_j| safely < 1

    const phi = { family: 'boundedLQD', w0: Complex.clone(w0), branches: [] };

    for (let j = 0; j < n; j++) {
      const p = hData.poles[j];
      const logRatio = {
        re: Math.log(Complex.abs(Complex.div(p.a, w0))),
        im: Math.atan2(p.a.im * w0.re - p.a.re * w0.im,
                       p.a.re * w0.re + p.a.im * w0.im),
      };
      // z_j ≈ logRatio / R, capped at cap to stay in 𝔻
      let z = Complex.scale(logRatio, 1 / R);
      const zr = Complex.abs(z);
      if (zr > cap) z = Complex.scale(z, cap / zr);

      // A_{j,k}: scale C_{j,k} by a factor that makes the (★) targets self-
      // consistent for a near-disk. For a disk of "radius" R, ψ̃(t) ≈ t/R,
      // so [t^s] ψ̃^k = R^{-k} δ_{sk}, and A_{j,k} ≈ D_{j,k} · R^{-k}.
      const A = [];
      const D = (() => {
        const out = new Array(p.principal.length);
        for (let s = 0; s < p.principal.length; s++) {
          const aC = Complex.mul(p.a, p.principal[s]);
          const next = (s + 1 < p.principal.length) ? p.principal[s + 1] : { re: 0, im: 0 };
          out[s] = Complex.add(aC, next);
        }
        return out;
      })();
      let Rk = 1;
      for (let k = 1; k <= p.principal.length; k++) {
        Rk *= R;
        A.push(Complex.scale(D[k - 1], 1 / Rk));
      }
      phi.branches.push({ z, A });
    }

    return phi;
  }

  function perturbedInitialGuess_LQD(hData, norm, rng, r) {
    const base = initialGuess_LQD(hData, norm);
    QD.LqdCommon.perturbBranchesInPlace(base.branches, rng, r || 0);
    return base;
  }

  function diverseInitialGuess_LQD(hData, norm, rng, r) {
    return {
      family: 'boundedLQD',
      w0: Complex.clone(norm.w0),
      branches: QD.LqdCommon.diverseSeedBranches(hData, rng),
    };
  }

  // No continuation strategy in the first cut — exp-of-QD bootstrap + the
  // shared multistart / diverse / deflation pipeline handle everything.
  // If a future case needs continuation (e.g. ramping a residue strength),
  // it slots in here.
  function continuationSolve_LQD(hData, norm, opts) {
    // No-op: signal failure so the solver pipeline falls through to multistart.
    return { success: false, error: "continuation not implemented for LQD", trace: [] };
  }

  // ===========================================================================
  // 5. Identity verification (modified Cauchy-kernel test class)
  // ===========================================================================
  //
  // Test functions f(w) = 1/(w − b)^k with b ∈ Ω* and b ≠ 0. These are
  // analytic in Ω (since b ∉ Ω), L¹ against ρ₀ = |w|⁻² near 0 (since
  // 0 ∈ Ω* by assumption, so f stays bounded inside Ω near 0), and
  // bounded near ∂Ω.
  //
  // LHS (Green's theorem):
  //   ∫_Ω f/|w|² dA = ∮_∂Ω f(w) · ln|w|²/w · dw
  // Numerically, sampling z_n = e^{iθ_n}, w_n = φ(z_n):
  //   LHS = (1/N) Σ_n (1/(w_n − b)^k) · ln|w_n|²/w_n · φ'(z_n) · z_n
  //
  // RHS (residue theorem at poles a_j of h inside Ω):
  //   RHS = Σ_j Σ_s C_{j,s} · (−1)^{s−1} · binom(k+s−2, s−1) / (a_j − b)^{k+s−1}
  function verifyQuadratureIdentity_LQD(phi, hData, options) {
    options = options || {};
    const N             = options.numSamples ?? 500;
    const maxOrder      = options.maxDegree ?? 3;
    const numTestPoints = options.numTestPoints ?? 3;

    // Sample boundary
    const samples = new Array(N);
    for (let n = 0; n < N; n++) {
      const theta = 2 * Math.PI * n / N;
      const z = { re: Math.cos(theta), im: Math.sin(theta) };
      const taylor = phiTaylorAt_LQD(z, phi, 1);
      samples[n] = { z, w: taylor[0], phiPrime: taylor[1] };
    }

    // Test points: b = 0 is always in Ω* by non-singular assumption; pair
    // with a couple of points well outside the bounding box of ∂Ω.
    let minRe = Infinity, maxRe = -Infinity, minIm = Infinity, maxIm = -Infinity;
    for (const s of samples) {
      if (s.w.re < minRe) minRe = s.w.re;
      if (s.w.re > maxRe) maxRe = s.w.re;
      if (s.w.im < minIm) minIm = s.w.im;
      if (s.w.im > maxIm) maxIm = s.w.im;
    }
    const spanRe = maxRe - minRe, spanIm = maxIm - minIm;
    const testPoints = [{ re: 0, im: 0 }];
    for (let i = 1; i < numTestPoints; i++) {
      // Points "far" outside — at the corners of an expanded bounding box.
      const ang = 2 * Math.PI * (i - 1) / Math.max(numTestPoints - 1, 1);
      const r = Math.max(spanRe, spanIm) * 1.5 + 1;
      testPoints.push({
        re: 0.5 * (minRe + maxRe) + r * Math.cos(ang),
        im: 0.5 * (minIm + maxIm) + r * Math.sin(ang),
      });
    }

    // Identity-scale (avoids 0/0 trouble for LHS≈RHS≈0): rough |area|.
    // For LQDs the natural scale is Σ |C_{j,1}| (sum of first-order residues).
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
        // --- LHS via trapezoidal rule ---
        let lhs = { re: 0, im: 0 };
        for (let n = 0; n < N; n++) {
          const s = samples[n];
          const diff = Complex.sub(s.w, b);
          const dPow = Complex.pow(diff, k);
          const fVal = Complex.inv(dPow);                 // 1/(w-b)^k

          // ln|w|²/w
          const absW2 = Complex.abs2(s.w);
          if (absW2 < 1e-30) continue;                    // skip if boundary touches 0
          const lnAbsW2 = Math.log(absW2);
          const lnTerm = Complex.scale(Complex.inv(s.w), lnAbsW2);

          let term = Complex.mul(fVal, lnTerm);
          term = Complex.mul(term, s.phiPrime);
          term = Complex.mul(term, s.z);
          lhs = Complex.add(lhs, term);
        }
        lhs = Complex.scale(lhs, 1 / N);

        // --- RHS via residues at finite poles ---
        let rhs = { re: 0, im: 0 };
        for (const pole of hData.poles) {
          const aMinusB = Complex.sub(pole.a, b);
          for (let sIdx = 0; sIdx < pole.principal.length; sIdx++) {
            const s = sIdx + 1;
            const C = pole.principal[sIdx];
            const sign = (s % 2 === 0) ? -1 : 1;          // (-1)^{s-1}
            const coef = QD.binomialCoeff
              ? QD.binomialCoeff(k + s - 2, s - 1)
              : binom(k + s - 2, s - 1);
            const denom = Complex.pow(aMinusB, k + s - 1);
            const term = Complex.div(C, denom);
            rhs = Complex.add(rhs, Complex.scale(term, sign * coef));
          }
        }

        const diff = Complex.sub(lhs, rhs);
        const absDiff = Complex.abs(diff);
        const scale = Math.max(Complex.abs(lhs), Complex.abs(rhs), scaleRef);
        const relDiff = absDiff / scale;
        if (relDiff > maxRelDiff) maxRelDiff = relDiff;
        if (absDiff > maxAbsDiff) maxAbsDiff = absDiff;
        checks.push({ bIdx: pIdx, k, lhs, rhs, absDiff, relDiff });
      }
    }

    return {
      checks,
      maxRelDiff,
      maxAbsDiff,
      areaScale: scaleRef,
      testPoints,
      maxDeg: maxOrder,
      numSamples: N,
      lqd: true,
    };
  }

  // Local fallback for the binomial coefficient (in case QD.binomialCoeff
  // isn't exported — it's an internal helper in solver.js).
  function binom(n, k) {
    if (k < 0 || k > n) return 0;
    let r = 1;
    for (let i = 0; i < k; i++) r = r * (n - i) / (i + 1);
    return r;
  }

  // ===========================================================================
  // 6. Register Family.boundedLQD
  // ===========================================================================
  QD.Family.boundedLQD = {
    name: 'boundedLQD',
    enforceInDisk:  true,
    enforceOutDisk: false,
    // Dispatch criterion: any opts bag with `lqd: true` and no `unbounded`.
    // Walked by selectFamily in registration order (boundedLQD before
    // boundedQD, since boundedQD is the catch-all default).
    matches(opts) { return !!(opts && opts.lqd && !opts.unbounded); },

    normalizeOpts(opts, hData) {
      let w0 = opts.w0;
      if (!w0) {
        // Default: centroid of poles. Warn (via thrown error) if centroid is at 0.
        let sumRe = 0, sumIm = 0;
        for (const p of hData.poles) { sumRe += p.a.re; sumIm += p.a.im; }
        const n = hData.poles.length;
        w0 = n > 0 ? { re: sumRe / n, im: sumIm / n } : { re: 1, im: 0 };
      }
      if (Complex.abs2(w0) < QD.ZERO_THRESHOLD) {
        throw new Error("Family.boundedLQD: w₀ = φ(0) must be nonzero (0 ∉ Ω̄ required)");
      }
      return { lqd: true, w0 };
    },

    evalPhi: evalPhi_LQD,
    phiTaylorAt: phiTaylorAt_LQD,

    computeTargets(phi, hData) {
      return { A: computeTargetA_LQD(phi, hData), F: null };
    },

    residual: residual_LQD,
    packPhi:  packPhi_LQD,
    unpackPhi: unpackPhi_LQD,
    canonicalizePhi: canonicalizePhi_LQD,

    initialGuess: initialGuess_LQD,
    perturbedInitialGuess: perturbedInitialGuess_LQD,
    diverseInitialGuess: diverseInitialGuess_LQD,
    continuationSolve: continuationSolve_LQD,
    verifyQuadratureIdentity: verifyQuadratureIdentity_LQD,
  };

  // Register with the dispatch order. This makes the internal selectFamily
  // (inside solver.js's closure) aware of boundedLQD — the previous monkey-
  // patch only updated QD.selectFamily and missed the internal call site.
  if (QD.registerFamily) {
    QD.registerFamily('boundedLQD');
  } else {
    throw new Error("solver-lqd.js: QD.registerFamily not found (solver.js too old?)");
  }

})();
