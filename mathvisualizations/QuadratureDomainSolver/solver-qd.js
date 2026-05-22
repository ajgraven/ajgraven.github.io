// =============================================================================
// solver-qd.js -- Bounded classical quadrature domains (Family.boundedQD)
//
// A classical bounded QD satisfies
//     ∫_Ω f(w) dA = ∮_∂Ω f(w) h(w) dw,   f ∈ A(Ω) ∩ L¹(Ω)
// with h rational and poles in Ω (the "quadrature nodes" at finite poles a_j).
//
// Riemann-map parametrization:
//     φ(z) = w_0 + Σ_j Σ_k conj(A_{j,k}) · z^k / (1 - conj(z_j) z)^k,  z ∈ 𝔻,
//     w_0 := φ(0) (user input).
//
// Locator (●): φ(z_j) = a_j  (one complex equation per pole).
// Faber-target (★): A_{j,k} = Σ_{s ≥ k} (s/k) · C_{j,s} · [t^s] ψ̃_j^k(t)
//   where ψ̃_j is the local inverse of φ(z_j + ·) − a_j.
// Gauge: Σ_j Im(A_{j,1}) = 0 + canonicalize Re sign.
//
// Companion file to solver-uqd.js (unbounded variant). Shared infrastructure
// (Newton driver, family registry, boundary sampler, top-level orchestrator)
// lives in solver.js.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family) {
    throw new Error("solver-qd.js: solver.js must be loaded first");
  }

  // ===========================================================================
  // 1. φ evaluation
  // ===========================================================================
  //   φ(z) = w_0 + Σ_j Σ_k conj(A_{j,k}) · z^k / (1 - conj(z_j) z)^k
  function evalPhi_QD(z, phi) {
    let result = Complex.clone(phi.w0);
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

  // ===========================================================================
  // 2. Taylor expansion of φ at z = z_0 up to order L
  // ===========================================================================
  // Closed form for u_j(z) = z/(1 − conj(z_j) z) at z = z_0:
  //   u_j(z_0) = z_0 / α,   u_j^{(l)}(z_0)/l! = conj(z_j)^{l-1} / α^{l+1}   (l ≥ 1)
  // with α = 1 − conj(z_j) z_0.
  function phiTaylorAt_QD(z0, phi, L) {
    const result = Taylor.zero(L + 1);
    result[0] = Complex.clone(phi.w0);

    for (const br of phi.branches) {
      const zjC = Complex.conj(br.z);
      const alpha = Complex.sub(Complex.ONE(), Complex.mul(zjC, z0));
      const alphaInv = Complex.inv(alpha);

      const uT = Taylor.zero(L + 1);
      uT[0] = Complex.mul(z0, alphaInv);
      if (L >= 1) {
        let zjcPow = { re: 1, im: 0 };                              // conj(z_j)^0
        let alphaInvPow = Complex.mul(alphaInv, alphaInv);          // 1/α^2
        for (let l = 1; l <= L; l++) {
          uT[l] = Complex.mul(zjcPow, alphaInvPow);
          zjcPow = Complex.mul(zjcPow, zjC);
          alphaInvPow = Complex.mul(alphaInvPow, alphaInv);
        }
      }

      let uPow = Taylor.truncate(uT, L);                            // u^1
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
  // 3. Target A_{j,k} from (★)
  // ===========================================================================
  function computeTargetA_QD(phi, hData) {
    // Inverse Faber transform at each finite pole. Shared with classical
    // unbounded QDs and (with modified residues) all LQDs via QD.Faber.
    const target = [];
    for (let j = 0; j < hData.poles.length; j++) {
      const principal = hData.poles[j].principal;
      const mj = principal.length;
      const zj = phi.branches[j].z;

      const phiT = phiTaylorAt_QD(zj, phi, mj);
      const phiTilde = Taylor.zero(mj + 1);
      for (let i = 1; i <= mj; i++) phiTilde[i] = Complex.clone(phiT[i]);

      target.push(QD.Faber.inverseFaberAtPole(principal, phiTilde));
    }
    return target;
  }

  // ===========================================================================
  // 4. Residual, pack, unpack, canonicalize
  // ===========================================================================
  function residual_QD(phi, hData, options = {}) {
    const enforceGauge = options.enforceGauge !== false;
    const out = [];

    // (●) φ(z_j) = a_j
    for (let j = 0; j < hData.poles.length; j++) {
      const phiZj = evalPhi_QD(phi.branches[j].z, phi);
      const diff = Complex.sub(phiZj, hData.poles[j].a);
      out.push(diff.re, diff.im);
    }
    // (★) A_{j,k} = target
    const target = computeTargetA_QD(phi, hData);
    for (let j = 0; j < hData.poles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], target[j][k]);
        out.push(diff.re, diff.im);
      }
    }
    // Gauge: Σ Im(A_{j,1}) = 0
    if (enforceGauge) {
      let imSum = 0;
      for (const br of phi.branches) if (br.A.length > 0) imSum += br.A[0].im;
      out.push(imSum);
    }
    return out;
  }

  function packPhi_QD(phi) {
    const v = [];
    for (const br of phi.branches) v.push(br.z.re, br.z.im);
    for (const br of phi.branches) for (const a of br.A) v.push(a.re, a.im);
    return v;
  }

  function unpackPhi_QD(v, template) {
    const phi = {
      unbounded: false,
      c: undefined,
      w0: Complex.clone(template.w0),
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
    return phi;
  }

  // Z/2 sign: Σ_j Re(A_{j,1}) > 0. If negative, apply z → -z.
  function canonicalizePhi_QD(phi) {
    let reSum = 0;
    for (const br of phi.branches) if (br.A.length > 0) reSum += br.A[0].re;
    if (reSum >= 0) return phi;
    return {
      w0: Complex.clone(phi.w0),
      branches: phi.branches.map(br => ({
        z: Complex.neg(br.z),
        A: br.A.map((a, k) => (k % 2 === 0 ? Complex.neg(a) : Complex.clone(a))),
      })),
    };
  }

  // ===========================================================================
  // 5. Initial guesses
  // ===========================================================================
  function diskInitialGuess_QD(hData, w0, scale = null) {
    const n = hData.poles.length;

    let totalC = 0;
    for (const p of hData.poles) {
      if (p.principal.length > 0) totalC += Complex.abs(p.principal[0]);
    }
    if (totalC === 0) totalC = 1;
    let R = scale !== null ? scale : Math.sqrt(totalC);
    let maxR = 0;
    for (const p of hData.poles) {
      const d = Complex.abs(Complex.sub(p.a, w0));
      if (d > maxR) maxR = d;
    }
    if (R < 1.5 * maxR) R = 1.5 * maxR;
    if (R === 0) R = 1;

    const phi = { unbounded: false, w0: Complex.clone(w0), c: undefined, branches: [] };
    for (let j = 0; j < n; j++) {
      const z = Complex.scale(Complex.sub(hData.poles[j].a, w0), 1 / R);
      const A = [];
      const mj = hData.poles[j].principal.length;
      let Rk = 1;
      for (let k = 1; k <= mj; k++) {
        Rk *= R;
        A.push(Complex.scale(hData.poles[j].principal[k - 1], 1 / Rk));
      }
      phi.branches.push({ z, A });
    }
    return phi;
  }

  function perturbedInitialGuess_QD(hData, w0, rng, r = 0) {
    const base = diskInitialGuess_QD(hData, w0);
    const sigma = 0.15 + 0.25 * r;
    for (const br of base.branches) {
      br.z = {
        re: br.z.re + sigma * (rng() - 0.5),
        im: br.z.im + sigma * (rng() - 0.5),
      };
      const rr = Math.hypot(br.z.re, br.z.im);
      if (rr > 0.9) { br.z.re *= 0.85 / rr; br.z.im *= 0.85 / rr; }
      for (let k = 0; k < br.A.length; k++) {
        br.A[k] = {
          re: br.A[k].re * (1 + sigma * (rng() - 0.5)),
          im: br.A[k].im + sigma * (rng() - 0.5),
        };
      }
    }
    return base;
  }

  // ===========================================================================
  // 6. Continuation (homotopy in pole distance)
  // ===========================================================================
  function continuationSolve_QD(hData, w0, options = {}) {
    const {
      tStart       = 0.10,
      tStartMin    = 1e-3,
      growFactor   = 1.6,
      shrinkFactor = 0.5,
      minStep      = 5e-4,
      maxSteps     = 80,
      newton       = {},
    } = options;

    const trace = [];
    let t = tStart;
    let phi = diskInitialGuess_QD(QD.scaleHDataPoles(hData, t, w0), w0);

    let warmupResult;
    while (true) {
      warmupResult = QD.newtonSolve(phi, QD.scaleHDataPoles(hData, t, w0), newton);
      if (warmupResult.success) { phi = warmupResult.phi; break; }
      t *= shrinkFactor;
      if (t < tStartMin) {
        return {
          success: false,
          error: "continuation: warmup failed even at t=" + t.toExponential(2),
          phi: warmupResult.phi,
          trace,
        };
      }
      phi = diskInitialGuess_QD(QD.scaleHDataPoles(hData, t, w0), w0);
    }
    trace.push({ t, ok: true, residual: warmupResult.residual });

    let lastSuccessT = t;
    let stepSize = Math.max(t, 0.1);
    for (let step = 0; step < maxSteps; step++) {
      if (lastSuccessT >= 1.0 - 1e-12) break;
      const nextT = Math.min(1.0, lastSuccessT + stepSize);
      const hData_t = QD.scaleHDataPoles(hData, nextT, w0);
      const result = QD.newtonSolve(phi, hData_t, newton);
      if (result.success) {
        phi = result.phi;
        lastSuccessT = nextT;
        trace.push({ t: nextT, ok: true, residual: result.residual });
        stepSize *= growFactor;
      } else {
        stepSize *= shrinkFactor;
        trace.push({ t: nextT, ok: false, residual: result.residual ?? null });
        if (stepSize < minStep) {
          return {
            success: false,
            error: "continuation: step size underflow at t=" + lastSuccessT.toFixed(4),
            phi, trace, lastT: lastSuccessT,
          };
        }
      }
    }
    if (lastSuccessT < 1.0 - 1e-9) {
      return {
        success: false,
        error: "continuation: max steps reached at t=" + lastSuccessT.toFixed(4),
        phi, trace, lastT: lastSuccessT,
      };
    }
    return {
      success: true, phi, iterations: 0,
      residual: trace[trace.length - 1].residual,
      trace, method: "continuation",
    };
  }

  // ===========================================================================
  // 7. Identity verification — monomials w^k (analytic in Ω, k = 0..maxDegree)
  //
  // ∫_Ω w^k dA = ∮_∂Ω w^k h(w) dw written as
  //   LHS  (1/N) Σ w_n^k · conj(w_n) · φ'(z_n) · z_n  (Green / Stokes)
  //   RHS  Σ_j Σ_s C_{j,s} · binom(k, s-1) · a_j^{k-s+1}  (residue sum)
  // Matches the original solver.js byte-for-byte.
  // ===========================================================================
  function verifyQuadratureIdentity_QD(phi, hData, options = {}) {
    const N = options.numSamples ?? 500;
    const totalDeg = hData.poles.reduce((s, p) => s + p.principal.length, 0);
    const K = options.maxDegree ?? Math.max(totalDeg, 4);

    const samples = new Array(N);
    for (let n = 0; n < N; n++) {
      const theta = 2 * Math.PI * n / N;
      const z = { re: Math.cos(theta), im: Math.sin(theta) };
      const taylor = phiTaylorAt_QD(z, phi, 1);
      samples[n] = { z, w: taylor[0], phiPrime: taylor[1] };
    }

    let areaScale = 0;
    for (const pole of hData.poles) {
      if (pole.principal.length > 0) areaScale += Complex.abs(pole.principal[0]);
    }
    if (areaScale === 0) areaScale = 1;

    const checks = [];
    let maxRelDiff = 0;
    let maxAbsDiff = 0;

    for (let k = 0; k <= K; k++) {
      let lhs = { re: 0, im: 0 };
      for (let n = 0; n < N; n++) {
        const s = samples[n];
        const wPow = Complex.pow(s.w, k);
        let term = Complex.mul(wPow, Complex.conj(s.w));
        term = Complex.mul(term, s.phiPrime);
        term = Complex.mul(term, s.z);
        lhs = Complex.add(lhs, term);
      }
      lhs = Complex.scale(lhs, 1 / N);

      let rhs = { re: 0, im: 0 };
      for (const pole of hData.poles) {
        for (let sIdx = 0; sIdx < pole.principal.length; sIdx++) {
          const s = sIdx + 1;
          if (s - 1 > k) continue;
          const C = pole.principal[sIdx];
          const exp = k - s + 1;
          const aPow = exp === 0 ? { re: 1, im: 0 } : Complex.pow(pole.a, exp);
          const coef = QD.binomialCoeff(k, s - 1);
          const term = Complex.scale(Complex.mul(C, aPow), coef);
          rhs = Complex.add(rhs, term);
        }
      }

      const diff = Complex.sub(lhs, rhs);
      const absDiff = Complex.abs(diff);
      const scale = Math.max(Complex.abs(lhs), Complex.abs(rhs), areaScale);
      const relDiff = absDiff / scale;
      if (relDiff > maxRelDiff) maxRelDiff = relDiff;
      if (absDiff > maxAbsDiff) maxAbsDiff = absDiff;
      checks.push({ k, lhs, rhs, absDiff, relDiff });
    }

    return {
      checks, maxRelDiff, maxAbsDiff,
      areaScale, maxDeg: K, numSamples: N,
    };
  }

  // ===========================================================================
  // 8. Register Family.boundedQD
  // ===========================================================================
  QD.Family.boundedQD = {
    name: 'boundedQD',
    enforceInDisk:  true,
    enforceOutDisk: false,
    matches(opts) { return true; },           // catch-all (checked last)

    normalizeOpts(opts, hData) {
      let w0 = opts.w0;
      if (!w0) {
        let sumRe = 0, sumIm = 0;
        for (const p of hData.poles) { sumRe += p.a.re; sumIm += p.a.im; }
        const n = hData.poles.length;
        w0 = n > 0 ? { re: sumRe / n, im: sumIm / n } : { re: 0, im: 0 };
      }
      return { w0 };
    },

    evalPhi: evalPhi_QD,
    phiTaylorAt: phiTaylorAt_QD,
    computeTargets(phi, hData) {
      return { A: computeTargetA_QD(phi, hData), F: null };
    },
    residual: residual_QD,
    packPhi: packPhi_QD,
    unpackPhi: unpackPhi_QD,
    canonicalizePhi: canonicalizePhi_QD,
    initialGuess(hData, norm) { return diskInitialGuess_QD(hData, norm.w0); },
    perturbedInitialGuess(hData, norm, rng, r) {
      return perturbedInitialGuess_QD(hData, norm.w0, rng, r);
    },
    diverseInitialGuess(hData, norm, rng, r) {
      return QD.diverseInitialGuess(hData, norm, rng, r);
    },
    continuationSolve(hData, norm, opts) {
      return continuationSolve_QD(hData, norm.w0, opts);
    },
    verifyQuadratureIdentity: verifyQuadratureIdentity_QD,
  };
  QD.registerFamily('boundedQD');

  // Export helpers under QD for other modules that need them (notably the
  // companion-QD bootstrap used by Family.boundedLQD_singular).
  QD.diskInitialGuess        = diskInitialGuess_QD;
  QD.perturbedInitialGuess   = perturbedInitialGuess_QD;
  QD.computeTargetA          = computeTargetA_QD;
  QD.continuationSolve       = continuationSolve_QD;
  QD.verifyQuadratureIdentity = verifyQuadratureIdentity_QD;

})();
