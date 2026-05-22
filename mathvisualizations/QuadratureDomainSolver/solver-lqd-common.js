// =============================================================================
// solver-lqd-common.js -- Shared machinery for LQD families
//
// Both bounded non-singular LQD (solver-lqd.js) and bounded singular LQD
// (solver-lqd-singular.js) share:
//
//   • modified-residue construction  D_{j,s} = a_j · C_{j,s} + C_{j,s+1}
//   • a "rational kernel" evaluator  r#(z), r#-Taylor-at-z_0
//       both reuse the bounded-QD primitive with w_0 = 0
//   • the (★) inverse-Faber target loop, parameterized by a
//     family-specific phiTaylorAt
//   • multistart helpers:
//       perturbedFromBase(phi, rng, r)   — Gaussian perturbation of any LQD phi
//       diverseSeedBranches(hData, rng)  — log-uniform {z_j, A_{j,k}} sampling
//   • identity-verifier scaffolding:
//       sampleBoundaryWithDerivative(phi, N, phiTaylorAt)  — uniform ∂Ω samples
//       compareLhsRhs(lhs, rhs, scaleRef) → {absDiff, relDiff}
//
// Each LQD family supplies the *family-specific* bits: gauge equation,
// q-equation if applicable, test-function class, canonicalization. The
// common machinery lives here.
//
// The upcoming unbounded LQD variants will reuse all of this (the Faber loop
// is family-agnostic given phiTaylorAt and modified residues).
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family) {
    throw new Error("solver-lqd-common.js: solver.js must be loaded first");
  }

  // ===========================================================================
  // Blaschke factor b_{z_0} (thesis form) — shared between bounded singular
  // and unbounded singular LQDs.
  //
  //   b_{z_0}(z) = -(conj(z_0)/|z_0|) · (z - z_0) / (1 - conj(z_0) z)
  //
  // Möbius-like map with a simple zero at z = z_0 and a simple pole at
  // z = 1/conj(z_0). For z_0 ∈ 𝔻 (bounded singular LQD): zero in 𝔻, pole in
  // 𝔻*. For z_0 ∈ 𝔻* (unbounded singular LQD): zero in 𝔻*, pole in 𝔻. Same
  // formula either way; the (1 - conj(z_0)z) denominator vanishes at z =
  // 1/conj(z_0) on the OPPOSITE side of the unit circle from z_0, so b_{z_0}
  // is locally smooth on the side containing z_0.
  //
  // Identity used by the q-equation derivations: b_{z_0} · b_{z_0}^# ≡ 1
  // (direct algebra; holds for either |z_0| < 1 or |z_0| > 1).
  // ===========================================================================
  function blaschkeEval(z, z0) {
    const absZ0 = Complex.abs(z0);
    if (absZ0 < 1e-14) throw new Error("blaschkeEval: z_0 = 0 not supported");
    const z0C = Complex.conj(z0);
    const phaseFactor = Complex.scale(z0C, -1 / absZ0);         // -conj(z_0)/|z_0|
    const num = Complex.sub(z, z0);
    const denom = Complex.sub(Complex.ONE(), Complex.mul(z0C, z));
    return Complex.mul(phaseFactor, Complex.div(num, denom));
  }

  // Taylor of b_{z_0}(z) at z = z_c, up to t^L. Closed-form expansion of
  // the two factors. Requires α := 1 - conj(z_0) z_c ≠ 0 (i.e., z_c ≠
  // 1/conj(z_0), the location of b's pole).
  function blaschkeTaylor(zc, z0, L) {
    const absZ0 = Complex.abs(z0);
    if (absZ0 < 1e-14) throw new Error("blaschkeTaylor: z_0 = 0 not supported");
    const z0C = Complex.conj(z0);
    const phaseFactor = Complex.scale(z0C, -1 / absZ0);

    const num = Taylor.zero(L + 1);
    num[0] = Complex.sub(zc, z0);
    if (L >= 1) num[1] = { re: 1, im: 0 };

    const alpha = Complex.sub(Complex.ONE(), Complex.mul(z0C, zc));
    const alphaInv = Complex.inv(alpha);
    const beta = Complex.mul(z0C, alphaInv);
    const invDenom = Taylor.zero(L + 1);
    let betaPow = Complex.ONE();
    for (let l = 0; l <= L; l++) {
      invDenom[l] = Complex.mul(betaPow, alphaInv);
      betaPow = Complex.mul(betaPow, beta);
    }

    const product = Taylor.mul(num, invDenom, L);
    const out = new Array(L + 1);
    for (let l = 0; l <= L; l++) out[l] = Complex.mul(phaseFactor, product[l]);
    return out;
  }

  // ===========================================================================
  // r#(∞) — closed-form sum for the unbounded LQD parametrization
  // ===========================================================================
  // r#(z) = Σ_j Σ_k conj(A_{j,k}) z^k / (1 - conj(z_j) z)^k
  // r#(∞) = Σ_j Σ_k conj(A_{j,k}) · (-1)^k / conj(z_j)^k
  //
  // Used by the unbounded LQD families to absorb the ∞-gauge into the
  // parametrization (φ uses exp(r#(z) - r#(∞)) so that φ'(∞) = c exactly).
  function rHashAtInfinity(phi) {
    let result = { re: 0, im: 0 };
    for (const br of phi.branches) {
      if (br.A.length === 0) continue;
      const zjC = Complex.conj(br.z);
      const zjCinv = Complex.inv(zjC);
      let zjCinvPow = { re: 1, im: 0 };                          // (1/conj(z_j))^0
      for (let k = 1; k <= br.A.length; k++) {
        zjCinvPow = Complex.mul(zjCinvPow, zjCinv);
        const sign = (k % 2 === 0) ? 1 : -1;                     // (-1)^k
        const AkC = Complex.conj(br.A[k - 1]);
        result = Complex.add(result, Complex.scale(Complex.mul(AkC, zjCinvPow), sign));
      }
    }
    return result;
  }

  // ===========================================================================
  // Modified residues — D_{j,s} = a_j · C_{j,s} + C_{j,s+1}
  // (with C_{j, m_j+1} ≡ 0). The bounded-QD inverse Faber formula carries
  // through unchanged for LQDs with C → D (Theorem 5.4.2).
  // ===========================================================================
  function modifiedResidues(hData) {
    return hData.poles.map(pole => {
      const C = pole.principal;
      const m = C.length;
      const D = new Array(m);
      for (let s = 0; s < m; s++) {
        const aC = Complex.mul(pole.a, C[s]);
        const next = (s + 1 < m) ? C[s + 1] : { re: 0, im: 0 };
        D[s] = Complex.add(aC, next);
      }
      return D;
    });
  }

  // ===========================================================================
  // r# evaluator -- reuse the bounded-QD primitive with w_0 = 0.
  //
  // r#(z) = Σ_j Σ_k conj(A_{j,k}) z^k / (1 − conj(z_j) z)^k
  //
  // We need a Family.boundedQD-tagged stub so the dispatcher routes to the
  // right evalPhi. (We can't just call QD.Family.boundedQD.evalPhi directly
  // because evalPhi/phiTaylorAt are top-level dispatchers — but the dispatch
  // path is identical and cheap.)
  // ===========================================================================
  function evalRHash(z, phi) {
    return QD.Family.boundedQD.evalPhi(z, {
      w0: { re: 0, im: 0 },
      branches: phi.branches,
    });
  }

  function rHashTaylorAt(z0, phi, L) {
    return QD.Family.boundedQD.phiTaylorAt(z0, {
      w0: { re: 0, im: 0 },
      branches: phi.branches,
    }, L);
  }

  // ===========================================================================
  // (★) Faber target loop -- family-generic given:
  //   • phiTaylorAt(z_j, phi, m_j)     — local Taylor of the family's full φ
  //   • D[j] = [D_{j,1}, …, D_{j,m_j}] — modified residues at pole j
  //
  // Returns target[j][k-1] = A_{j,k}^target.
  // ===========================================================================
  function computeFaberTargetA(phi, hData, phiTaylorFn) {
    // LQD inverse Faber at each finite pole, using MODIFIED residues
    // D_{j,s} = a_j · C_{j,s} + C_{j,s+1}. Underlying primitive is shared
    // with classical QDs via QD.Faber.
    const D_all = modifiedResidues(hData);
    const target = [];
    for (let j = 0; j < hData.poles.length; j++) {
      const D = D_all[j];
      const mj = D.length;
      const zj = phi.branches[j].z;

      const phiT = phiTaylorFn(zj, phi, mj);
      // φ̃(t) = φ(z_j + t) − a_j: drop constant (locator handles its value).
      const phiTilde = Taylor.zero(mj + 1);
      for (let i = 1; i <= mj; i++) phiTilde[i] = Complex.clone(phiT[i]);

      target.push(QD.Faber.inverseFaberAtPole(D, phiTilde));
    }
    return target;
  }

  // ===========================================================================
  // Multistart helpers
  // ===========================================================================
  // Gaussian-perturb in place: z_j gets re/im noise then clamp; A_{j,k} gets
  // multiplicative-on-re + additive-on-im noise.
  //
  // opts.side controls the disk-side constraint applied AFTER perturbation:
  //   'in'  (default) — keep |z| < zCap     (bounded LQD: z_j ∈ 𝔻)
  //   'out'           — keep |z| > zCap     (unbounded LQD: z_j ∈ 𝔻*)
  function perturbBranchesInPlace(branches, rng, r = 0, opts = {}) {
    const sigma = (opts.sigmaBase ?? 0.15) + (opts.sigmaSlope ?? 0.25) * r;
    const side  = opts.side ?? 'in';
    const cap   = opts.zCap   ?? (side === 'in' ? 0.9  : 1.05);
    const scale = opts.zScale ?? (side === 'in' ? 0.85 : 1.10);
    for (const br of branches) {
      br.z = {
        re: br.z.re + sigma * (rng() - 0.5),
        im: br.z.im + sigma * (rng() - 0.5),
      };
      const rr = Math.hypot(br.z.re, br.z.im);
      if (side === 'in' && rr > cap) {
        br.z.re *= scale / rr; br.z.im *= scale / rr;
      } else if (side === 'out' && rr < cap) {
        const s = rr > 1e-12 ? scale / rr : scale;
        br.z.re *= s; br.z.im *= s;
      }
      for (let k = 0; k < br.A.length; k++) {
        br.A[k] = {
          re: br.A[k].re * (1 + sigma * (rng() - 0.5)),
          im: br.A[k].im + sigma * (rng() - 0.5),
        };
      }
    }
  }

  // Log-uniform branch seed (z_j ∈ 𝔻 ring, A_{j,k} log-uniform magnitudes).
  function diverseSeedBranches(hData, rng, opts = {}) {
    const zMin = opts.zMin ?? 0.05, zMax = opts.zMax ?? 0.95;
    const aMin = opts.aMin ?? 0.1,  aMax = opts.aMax ?? 3.0;
    const branches = [];
    for (const p of hData.poles) {
      const mz = Math.exp(Math.log(zMin) + rng() * Math.log(zMax / zMin));
      const pz = 2 * Math.PI * rng();
      const z = { re: mz * Math.cos(pz), im: mz * Math.sin(pz) };
      const A = [];
      for (let k = 0; k < p.principal.length; k++) {
        const ma = Math.exp(Math.log(aMin) + rng() * Math.log(aMax / aMin));
        const pa = 2 * Math.PI * rng();
        A.push({ re: ma * Math.cos(pa), im: ma * Math.sin(pa) });
      }
      branches.push({ z, A });
    }
    return branches;
  }

  // ===========================================================================
  // Identity-verifier scaffolding
  // ===========================================================================
  // Sample {z_n, w_n, φ'(z_n)} on z = e^{iθ} for N uniform thetas, using the
  // family's phiTaylorAt to get both φ and φ' at once.
  function sampleBoundaryWithDerivative(phi, N, phiTaylorFn) {
    const samples = new Array(N);
    for (let n = 0; n < N; n++) {
      const theta = 2 * Math.PI * n / N;
      const z = { re: Math.cos(theta), im: Math.sin(theta) };
      const taylor = phiTaylorFn(z, phi, 1);
      samples[n] = { z, w: taylor[0], phiPrime: taylor[1] };
    }
    return samples;
  }

  // Standard (LHS, RHS) → diff metric used by all LQD verifiers. scaleRef is
  // the "natural" magnitude floor (typically Σ|C_{j,1}| + |q|) to keep the
  // relative diff meaningful when both sides are near zero.
  function compareLhsRhs(lhs, rhs, scaleRef) {
    const diff = Complex.sub(lhs, rhs);
    const absDiff = Complex.abs(diff);
    const scale = Math.max(Complex.abs(lhs), Complex.abs(rhs), scaleRef);
    return { absDiff, relDiff: absDiff / scale };
  }

  // Standard scale-floor reference: Σ |C_{j,1}| + |q| (q optional).
  function residueScaleFloor(hData, qOpt) {
    let s = 0;
    for (const pole of hData.poles) {
      if (pole.principal.length > 0) s += Complex.abs(pole.principal[0]);
    }
    if (qOpt) s += Complex.abs(qOpt);
    return s || 1;
  }

  // ===========================================================================
  // Generic identity-verifier skeleton (R4)
  // ---------------------------------------------------------------------------
  // All LQD identity verifiers compute LHS via trapezoidal rule on ∂Ω and
  // RHS via a closed-form residue sum. The differences are:
  //   • which test functions to use (monomials, 1/(w-b)^k, w/(w-b)^k, …)
  //   • how to integrate f against the boundary kernel
  //   • how to compute the RHS residue sum
  //
  // The skeleton handles boundary sampling, scale-floor heuristics, the per-
  // test bookkeeping, and the result struct. Each family supplies:
  //
  //   spec = {
  //     phiTaylorFn(z, phi, L)               — family's phiTaylorAt
  //     buildTestFunctions(phi, hData, opts) — returns [{ label, f, residueRhs }]
  //       where f(w) is analytic in Ω and residueRhs is a closed-form complex
  //       value for ∮f h dw (typically computed from h's residues at the a_j
  //       plus, for unbounded, residue at ∞).
  //     boundaryKernel(w)                    — returns the LHS integrand
  //       multiplier (e.g. ln|w|²/w for bounded LQDs).
  //     scaleRefExtra(hData, phi)            — extra contribution to the
  //       relative-error floor; defaults to 0.
  //     resultFlags                          — e.g. { lqdSingular: true }
  //   }
  //
  // Returns the canonical verifier result struct:
  //   { checks: [{label, lhs, rhs, absDiff, relDiff}], maxRelDiff, maxAbsDiff,
  //     areaScale, numSamples, maxDeg, ...resultFlags }
  // ===========================================================================
  function verifyIdentityGeneric(phi, hData, options, spec) {
    const N = options.numSamples ?? 500;
    const samples = sampleBoundaryWithDerivative(phi, N, spec.phiTaylorFn);

    // Boundary kernel pre-applied per-sample for efficiency.
    const samplesK = samples.map(s => {
      const kernel = spec.boundaryKernel(s.w);
      return { ...s, kernel };
    });

    const scaleRef = residueScaleFloor(hData, phi.q) + (spec.scaleRefExtra?.(hData, phi) ?? 0);
    const tests = spec.buildTestFunctions(phi, hData, options);

    const checks = [];
    let maxRelDiff = 0;
    let maxAbsDiff = 0;

    for (const test of tests) {
      // LHS = (2π/N) · i · Σ_n  f(w_n) · kernel(w_n) · z_n · φ'(z_n)
      //   (the boundary integral ∮_∂Ω f · kernel · dw with dw = i z φ'(z) dθ
      //    and the 2π/N trapezoid weight, matching the existing verifiers'
      //    conventions for ∮ … dw with no leading factor.)
      let lhs = { re: 0, im: 0 };
      for (const s of samplesK) {
        if (!s.kernel) continue;
        const fv = test.f(s.w);
        if (!fv) continue;
        let term = Complex.mul(fv, s.kernel);
        term = Complex.mul(term, s.z);
        term = Complex.mul(term, s.phiPrime);
        lhs = Complex.add(lhs, term);
      }
      lhs = Complex.scale(lhs, 2 * Math.PI / N);
      lhs = { re: -lhs.im, im: lhs.re };           // multiply by i

      const rhs = test.residueRhs;
      const { absDiff, relDiff } = compareLhsRhs(lhs, rhs, scaleRef);
      if (relDiff > maxRelDiff) maxRelDiff = relDiff;
      if (absDiff > maxAbsDiff) maxAbsDiff = absDiff;
      checks.push({ label: test.label, lhs, rhs, absDiff, relDiff, ...(test.tag || {}) });
    }

    return {
      checks, maxRelDiff, maxAbsDiff,
      areaScale: scaleRef,
      numSamples: N,
      maxDeg: options.maxDegree ?? 3,
      ...(spec.resultFlags || {}),
    };
  }

  // ===========================================================================
  // Expose under QD.LqdCommon
  // ===========================================================================
  QD.LqdCommon = {
    blaschkeEval,
    blaschkeTaylor,
    rHashAtInfinity,
    modifiedResidues,
    evalRHash,
    rHashTaylorAt,
    computeFaberTargetA,
    perturbBranchesInPlace,
    diverseSeedBranches,
    sampleBoundaryWithDerivative,
    compareLhsRhs,
    residueScaleFloor,
    verifyIdentityGeneric,
  };

})();
