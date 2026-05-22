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
// LIMITATIONS / TODO:
//   • Higher-order pole at 0 in h (a_j = 0 with m_j ≥ 2). Per the thesis the
//     parametrization extends to add c_l/z^l terms to r#; we don't yet handle
//     this. Solver throws a clear "deferred" error in normalizeOpts.
//   • Polynomial part of h not yet matched by the solver (same TODO as in
//     solver-uqd-lqd.js).
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
  // 1. φ evaluation: c · |z₀| · z · b_{z₀}(z) · exp(r#(z) − r#(∞))
  // ===========================================================================
  function evalPhi_UQDLS(z, phi) {
    const r = evalRHash(z, phi);
    const rInf = rHashAtInfinity(phi);
    const rEff = Complex.sub(r, rInf);
    const ea = Math.exp(rEff.re);
    const expR = { re: ea * Math.cos(rEff.im), im: ea * Math.sin(rEff.im) };
    const b = blaschkeEval(z, phi.z0);
    const absZ0 = Complex.abs(phi.z0);
    // c · |z_0| · z · b · expR
    const scale = Complex.scale(z, phi.c * absZ0);
    return Complex.mul(Complex.mul(scale, b), expR);
  }

  // ===========================================================================
  // 2. φ-Taylor at z = z_c
  // ===========================================================================
  // φ(z_c + t) = c·|z₀|·(z_c+t)·b_{z₀}(z_c+t)·exp(r#(z_c+t) − r#(∞))
  //            = K · lin(t) · bT(t) · expRTilde(t)
  // with
  //   K          = c·|z₀|·exp(r#(z_c) − r#(∞))
  //   lin(t)     = z_c + t                                (Taylor [z_c, 1, 0, ...])
  //   bT(t)      = b_{z₀}(z_c + t)                        (blaschkeTaylor)
  //   expRTilde  = exp(r#(z_c+t) − r#(z_c))               (Taylor.exp with zero constant)
  function phiTaylorAt_UQDLS(zc, phi, L) {
    const rT = rHashTaylorAt(zc, phi, L);
    const rInf = rHashAtInfinity(phi);
    const r0minusInf = Complex.sub(rT[0], rInf);

    const rTilde = rT.slice();
    rTilde[0] = { re: 0, im: 0 };
    const expRTilde = Taylor.exp(rTilde, L);

    const ea = Math.exp(r0minusInf.re);
    const expR0 = { re: ea * Math.cos(r0minusInf.im), im: ea * Math.sin(r0minusInf.im) };
    const absZ0 = Complex.abs(phi.z0);
    const K = Complex.scale(expR0, phi.c * absZ0);

    const lin = Taylor.zero(L + 1);
    lin[0] = Complex.clone(zc);
    if (L >= 1) lin[1] = { re: 1, im: 0 };

    const bT = blaschkeTaylor(zc, phi.z0, L);

    const t1 = Taylor.mul(lin, bT, L);
    const t2 = Taylor.mul(t1, expRTilde, L);
    const out = new Array(L + 1);
    for (let l = 0; l <= L; l++) out[l] = Complex.mul(K, t2[l]);
    return out;
  }

  function computeTargetA_UQDLS(phi, hData) {
    return computeFaberTargetA(phi, hData, phiTaylorAt_UQDLS);
  }

  // ===========================================================================
  // 3. Residual: (●) 2n + (★) 2M + (●₀ q-eq) 2
  // ===========================================================================
  function residual_UQDLS(phi, hData, options) {
    options = options || {};
    const out = [];

    for (let j = 0; j < hData.poles.length; j++) {
      const phiZj = evalPhi_UQDLS(phi.branches[j].z, phi);
      const diff = Complex.sub(phiZj, hData.poles[j].a);
      out.push(diff.re, diff.im);
    }

    const target = computeTargetA_UQDLS(phi, hData);
    for (let j = 0; j < hData.poles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], target[j][k]);
        out.push(diff.re, diff.im);
      }
    }

    // (●₀) q-equation:
    //   q = ln(c²|z₀|²) + r̃#(z₀) + conj(r̃#(1/conj(z₀)))
    // with r̃#(z) := r#(z) − r#(∞).
    const rInf = rHashAtInfinity(phi);
    const rZ0 = evalRHash(phi.z0, phi);
    const absZ02 = Complex.abs2(phi.z0);
    const oneOverConjZ0 = Complex.scale(phi.z0, 1 / absZ02);     // 1/conj(z₀)
    const rInvZ0 = evalRHash(oneOverConjZ0, phi);
    const rTildeZ0  = Complex.sub(rZ0, rInf);
    const rTildeInv = Complex.sub(rInvZ0, rInf);
    const sum = Complex.add(rTildeZ0, Complex.conj(rTildeInv));
    const lnFactor = Math.log(phi.c * phi.c * absZ02);
    const lhs = { re: sum.re + lnFactor, im: sum.im };
    const diff = Complex.sub(phi.q, lhs);
    out.push(diff.re, diff.im);

    return out;
  }

  // c > 0 pins the rotation gauge; no Z/2 ambiguity.
  function canonicalizePhi_UQDLS(phi) { return phi; }

  // ===========================================================================
  // 4. Pack / unpack — schema-driven
  // ===========================================================================
  // Layout: [{z_j}_{j=1..n}, z₀, {A_{j,k}}].  z_j and z₀ both clamped to 𝔻*;
  // z₀ additionally bounded above (maxR=1000) to keep Newton out of the
  // deferred z₀ → ∞ degeneracy.
  const SCHEMA_UQDLS = [
    { kind: 'branchesZ', clamp: { side: 'out', cap: 1.0001 } },
    { kind: 'complex',   name: 'z0', clamp: { side: 'out', cap: 1.0001, maxR: 1000 } },
    { kind: 'branchesA' },
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
  function initialGuess_UQDLS(hData, norm) {
    const c = norm.c;
    const q = norm.q;

    let zj_guess = null, A_guess = null, z0_guess = null;

    // Try companion bootstrap (only if at least one finite pole).
    if (hData.poles.length > 0) {
      try {
        const companion = QD.solveInverseQD(hData, {
          lqd: true, unbounded: true, c,
          identityTol: 1e-3, autoEscalate: false, findAlternates: false,
        });
        if (companion.success && companion.primary && companion.primary.phi) {
          const phiUQDL = companion.primary.phi;
          zj_guess = phiUQDL.branches.map(br => Complex.clone(br.z));
          A_guess  = phiUQDL.branches.map(br => br.A.map(Complex.clone));

          // z₀ = argmin |φ_UQDL(z)| on |z| = 1.01, pushed slightly outward.
          // (Non-singular companion has 0 ∉ Ω̄, so |φ| > 0 everywhere; we
          // use the closest-approach point as a starting guess for z₀.)
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
      zj_guess = hData.poles.map(p => {
        if (Complex.abs2(p.a) < 1e-30) return { re: 2, im: 0 };
        let z = Complex.scale(p.a, 1 / c);
        const r = Complex.abs(z);
        if (r < 1.05) z = Complex.scale(z, 1.05 / Math.max(r, 1e-15));
        return z;
      });
      A_guess = hData.poles.map(p => {
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

    return {
      family: 'unboundedLQD_singular',
      unbounded: true,
      c, q: Complex.clone(q),
      z0: z0_guess,
      w0: undefined,
      branches: zj_guess.map((z, j) => ({ z, A: A_guess[j].map(Complex.clone) })),
    };
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
    return base;
  }

  function diverseInitialGuess_UQDLS(hData, norm, rng, r) {
    const c = norm.c, q = norm.q;
    const mz0 = Math.exp(Math.log(1.05) + rng() * Math.log(30 / 1.05));
    const pz0 = 2 * Math.PI * rng();
    return {
      family: 'unboundedLQD_singular',
      unbounded: true,
      c, q: Complex.clone(q),
      z0: { re: mz0 * Math.cos(pz0), im: mz0 * Math.sin(pz0) },
      w0: undefined,
      branches: diverseSeedBranches(hData, rng, { zMin: 1.05, zMax: 30 }),
    };
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
        for (const b of testPoints) {
          for (let k = minOrder; k <= maxOrder; k++) {
            // f(w) = w/(w−b)^k = 1/(w−b)^{k−1} + b/(w−b)^k
            // Residue at a_j:
            //   C·(−1)^{s−1}·[ binom(k+s−3, s−1)/(a_j−b)^{k+s−2}
            //                + b·binom(k+s−2, s−1)/(a_j−b)^{k+s−1} ]
            // q/w pole at 0 contributes 0 (f·q/w analytic at 0).
            let rhsSum = { re: 0, im: 0 };
            for (const pole of hData.poles) {
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
      // Any pole at a = 0 in hData implies higher-order singularity at 0
      // (since the simple residue belongs to opts.q). Deferred to follow-up.
      for (const p of hData.poles) {
        if (Complex.abs2(p.a) < 1e-20) {
          throw new Error(
            "Family.unboundedLQD_singular: pole at a = 0 in hData detected — " +
            "higher-order pole at origin not yet implemented (deferred to follow-up). " +
            "For order-1 pole at 0, use opts.q instead of an hData entry."
          );
        }
      }
      // h = q/w only (no finite poles, nonzero q) has no solution (Theorem 5.5.2-style).
      if (hData.poles.length === 0 && Complex.abs2(q) > 1e-20) {
        throw new Error(
          "Family.unboundedLQD_singular: no unbounded singular LQD exists for h = q/w " +
          "with no finite poles (you can add a finite pole, or set q = 0)."
        );
      }
      return { lqd: true, unbounded: true, singular: true, c, q: Complex.clone(q) };
    },

    evalPhi: evalPhi_UQDLS,
    phiTaylorAt: phiTaylorAt_UQDLS,
    computeTargets(phi, hData) {
      return { A: computeTargetA_UQDLS(phi, hData), F: null };
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
