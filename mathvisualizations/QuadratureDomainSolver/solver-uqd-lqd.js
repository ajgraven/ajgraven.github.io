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
// IMPORTANT — polynomial part of h:  Per the user's answer (Q2), polynomial
// terms in h affect r# (not separate F_l unknowns). The full inverse-Faber
// machinery accounting for a polynomial-part h is deferred to a follow-up
// (this file currently handles only finite-pole h; the polyPart contribution
// to the IDENTITY VERIFIER's RHS is included so the verifier doesn't lie
// when the user does feed in a polyPart, but the SOLVER will not yet match
// it correctly — flagged as TODO).
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
  // 1. Phi evaluation: φ(z) = c · z · exp(r#(z) − r#(∞))
  // ===========================================================================
  function evalPhi_UQDL(z, phi) {
    const r = QD.LqdCommon.evalRHash(z, phi);
    const rInf = rHashAtInfinity(phi);
    const rEff = Complex.sub(r, rInf);          // r#(z) − r#(∞)
    const ea = Math.exp(rEff.re);
    const expR = { re: ea * Math.cos(rEff.im), im: ea * Math.sin(rEff.im) };
    const cz = Complex.scale(z, phi.c);
    return Complex.mul(cz, expR);
  }

  // ===========================================================================
  // 2. Taylor of φ at z = z_c
  // ===========================================================================
  // φ(z_c + t) = c (z_c + t) · exp(r#(z_c) − r#(∞)) · exp(r#(z_c+t) − r#(z_c))
  //            = K · (z_c + t) · exp(rTilde(t))
  // with K = c · exp(r#(z_c) − r#(∞)) and rTilde(t) = r#(z_c+t) − r#(z_c).
  function phiTaylorAt_UQDL(zc, phi, L) {
    const rT = QD.LqdCommon.rHashTaylorAt(zc, phi, L);   // Taylor of r# at z_c
    const rInf = rHashAtInfinity(phi);
    const r0minusInf = Complex.sub(rT[0], rInf);         // r#(z_c) − r#(∞)

    const rTilde = rT.slice();
    rTilde[0] = { re: 0, im: 0 };
    const expRTilde = Taylor.exp(rTilde, L);

    const ea = Math.exp(r0minusInf.re);
    const expR0 = { re: ea * Math.cos(r0minusInf.im), im: ea * Math.sin(r0minusInf.im) };
    const K = Complex.scale(expR0, phi.c);

    const lin = Taylor.zero(L + 1);
    lin[0] = Complex.clone(zc);
    if (L >= 1) lin[1] = { re: 1, im: 0 };

    const linTimesExp = Taylor.mul(lin, expRTilde, L);
    const out = new Array(L + 1);
    for (let l = 0; l <= L; l++) out[l] = Complex.mul(K, linTimesExp[l]);
    return out;
  }

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
  // 5. Residual
  // ===========================================================================
  // Layout (square system: 2n + 2M):
  //   (●)   2n   φ(z_j) = a_j
  //   (★)   2M   A_{j,k} = target via modified-residue Faber
  // The ∞-gauge r#(∞) = 0 is built into evalPhi / phiTaylorAt via the
  // r#(z) − r#(∞) subtraction, so no explicit equation here.
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
    return out;
  }

  // No Z/2 ambiguity in the unbounded case — c > 0 already pins everything.
  function canonicalizePhi_UQDL(phi) { return phi; }

  // ===========================================================================
  // 6. Pack / unpack — schema-driven
  // ===========================================================================
  // Same layout as bounded LQD non-singular: [{z_j}, {A_{j,k}}]. z_j ∈ 𝔻*
  // so clamp to keep |z| ≥ 1.0001.
  const SCHEMA_UQDL = [
    { kind: 'branchesZ', clamp: { side: 'out', cap: 1.0001 } },
    { kind: 'branchesA' },
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

    return {
      family: 'unboundedLQD',
      unbounded: true,
      c, w0: undefined,
      branches,
    };
  }

  function perturbedInitialGuess_UQDL(hData, norm, rng, r) {
    const base = initialGuess_UQDL(hData, norm);
    QD.LqdCommon.perturbBranchesInPlace(base.branches, rng, r || 0,
      { side: 'out', zCap: 1.05, zScale: 1.10 });
    return base;
  }

  function diverseInitialGuess_UQDL(hData, norm, rng, r) {
    return {
      family: 'unboundedLQD',
      unbounded: true,
      c: norm.c, w0: undefined,
      branches: QD.LqdCommon.diverseSeedBranches(hData, rng, { zMin: 1.05, zMax: 30 }),
    };
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
      return { A: computeTargetA_UQDL(phi, hData), F: null };
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
