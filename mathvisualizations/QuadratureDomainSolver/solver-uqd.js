// =============================================================================
// solver-uqd.js -- Unbounded classical quadrature domains (Family.unboundedQD)
//
// An unbounded simply-connected Ω with bounded complement K, satisfying
//     ∫_Ω f(w) dA = ∮_∂Ω f(w) h(w) dw,   f ∈ A_0(Ω)
// (the algebra A_0(Ω) of analytic functions on Ω vanishing at ∞).
//
// Riemann-map parametrization (φ: 𝔻* → Ω, φ(∞) = ∞, φ'(∞) = c > 0):
//     φ(z) = c·z + Σ_{l=0..m_∞} F_l / z^l
//          + Σ_j Σ_k conj(A_{j,k}) · z^k / (1 - conj(z_j) z)^k.
// z_j ∈ 𝔻* (poles of φ at 1/conj(z_j) ∈ 𝔻).
//
// Locator (●): φ(z_j) = a_j  (finite-pole locations).
// Faber (★)_A: A_{j,k} = Σ_{s≥k} (s/k) C_{j,s} [t^s] ψ̃_j^k(t)   (finite poles)
// Faber (★)_F: F_l = Σ_{l'≥l} conj(C_{∞,l'} · [u^{l'-l}] g(u)^{l'})   (∞ poles)
// Gauge: implicit — c > 0 pins the disk rotation; no Σ-Im constraint.
//
// Companion file to solver-qd.js (bounded variant). Shared infrastructure
// lives in solver.js.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Family) {
    throw new Error("solver-uqd.js: solver.js must be loaded first");
  }

  // ===========================================================================
  // 1. φ evaluation
  // ===========================================================================
  function evalPhi_UQD(z, phi) {
    let result = Complex.scale(z, phi.c);

    if (phi.polyA && phi.polyA.length > 0) {
      result = Complex.add(result, phi.polyA[0]);
      let zPow = Complex.clone(z);                      // z^1
      for (let l = 1; l < phi.polyA.length; l++) {
        result = Complex.add(result, Complex.div(phi.polyA[l], zPow));
        if (l + 1 < phi.polyA.length) zPow = Complex.mul(zPow, z);
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

  // ===========================================================================
  // 2. Taylor of φ at z = z_0
  // ===========================================================================
  function phiTaylorAt_UQD(z0, phi, L) {
    const result = Taylor.zero(L + 1);

    // c·z at z = z_0:  [c·z_0, c, 0, ...]
    result[0] = Complex.scale(z0, phi.c);
    if (L >= 1) result[1] = { re: phi.c, im: 0 };

    // Polynomial-part F_l / z^l Taylor at z = z_0.
    if (phi.polyA && phi.polyA.length > 0) {
      result[0] = Complex.add(result[0], phi.polyA[0]);
      if (phi.polyA.length >= 2) {
        const z0inv = Complex.inv(z0);
        let z0invPowL = Complex.clone(z0inv);
        for (let l = 1; l < phi.polyA.length; l++) {
          let z0invPowLM = Complex.clone(z0invPowL);
          let binom = 1;
          for (let m = 0; m <= L; m++) {
            const sign = (m % 2 === 0) ? 1 : -1;
            const coef = Complex.scale(
              Complex.mul(phi.polyA[l], z0invPowLM), sign * binom);
            result[m] = Complex.add(result[m], coef);
            if (m < L) {
              z0invPowLM = Complex.mul(z0invPowLM, z0inv);
              binom = binom * (l + m) / (m + 1);
            }
          }
          if (l + 1 < phi.polyA.length) z0invPowL = Complex.mul(z0invPowL, z0inv);
        }
      }
    }

    // Finite-pole contributions (same shape as bounded case).
    for (const br of phi.branches) {
      const zjC = Complex.conj(br.z);
      const alpha = Complex.sub(Complex.ONE(), Complex.mul(zjC, z0));
      const alphaInv = Complex.inv(alpha);
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
  // 3. Laurent at infinity
  // ===========================================================================
  //   φ(z) = c·z + f_0 + f_1/z + f_2/z² + …
  //   f_l = F_l + Σ_{j,k} conj(A_{j,k}) · (-1)^k · binom(k+l-1, l) / conj(z_j)^{k+l}
  function phiLaurentAtInfinity_UQD(phi, L) {
    if (L <= 0) return [];
    const f = new Array(L);
    for (let l = 0; l < L; l++) f[l] = { re: 0, im: 0 };

    if (phi.polyA) {
      for (let l = 0; l < L && l < phi.polyA.length; l++) {
        f[l] = Complex.add(f[l], phi.polyA[l]);
      }
    }
    for (const br of phi.branches) {
      const zjC = Complex.conj(br.z);
      const zjCinv = Complex.inv(zjC);
      for (let k = 1; k <= br.A.length; k++) {
        const AkC = Complex.conj(br.A[k - 1]);
        const sign = (k % 2 === 0) ? 1 : -1;
        let zjCinvPow = Complex.pow(zjCinv, k);
        let binom = 1;
        for (let l = 0; l < L; l++) {
          const contrib = Complex.scale(Complex.mul(AkC, zjCinvPow), sign * binom);
          f[l] = Complex.add(f[l], contrib);
          if (l + 1 < L) {
            zjCinvPow = Complex.mul(zjCinvPow, zjCinv);
            binom = binom * (k + l) / (l + 1);
          }
        }
      }
    }
    return f;
  }

  // ===========================================================================
  // 4. Targets — A_{j,k} (finite poles) and F_l (∞ pole)
  // ===========================================================================
  function computeTargetA_UQD(phi, hData) {
    // Same Faber primitive as bounded QD; only difference is z_j ∈ 𝔻*.
    const target = [];
    for (let j = 0; j < hData.poles.length; j++) {
      const principal = hData.poles[j].principal;
      const mj = principal.length;
      const zj = phi.branches[j].z;

      const phiT = phiTaylorAt_UQD(zj, phi, mj);
      const phiTilde = Taylor.zero(mj + 1);
      for (let i = 1; i <= mj; i++) phiTilde[i] = Complex.clone(phiT[i]);

      target.push(QD.Faber.inverseFaberAtPole(principal, phiTilde));
    }
    return target;
  }

  function computeTargetF_UQD(phi, hData) {
    // Inverse Faber transform at ∞ for the polynomial part of h.
    const polyPart = hData.polyPart || [];
    const m_inf = polyPart.length - 1;
    if (m_inf < 0) return [];
    const f = phiLaurentAtInfinity_UQD(phi, m_inf);
    return QD.Faber.inverseFaberAtInfinity(polyPart, f, phi.c);
  }

  // ===========================================================================
  // 5. Residual, pack, unpack, canonicalize
  // ===========================================================================
  function residual_UQD(phi, hData, options = {}) {
    const out = [];

    // (●)
    for (let j = 0; j < hData.poles.length; j++) {
      const phiZj = evalPhi_UQD(phi.branches[j].z, phi);
      const diff = Complex.sub(phiZj, hData.poles[j].a);
      out.push(diff.re, diff.im);
    }
    // (★)_A
    const targetA = computeTargetA_UQD(phi, hData);
    for (let j = 0; j < hData.poles.length; j++) {
      const A = phi.branches[j].A;
      for (let k = 0; k < A.length; k++) {
        const diff = Complex.sub(A[k], targetA[j][k]);
        out.push(diff.re, diff.im);
      }
    }
    // (★)_F (polynomial part)
    if (phi.polyA && phi.polyA.length > 0) {
      const targetF = computeTargetF_UQD(phi, hData);
      for (let l = 0; l < phi.polyA.length; l++) {
        const diff = Complex.sub(phi.polyA[l], targetF[l]);
        out.push(diff.re, diff.im);
      }
    }
    // No gauge constraint for unbounded (c > 0 fixes it).
    return out;
  }

  function packPhi_UQD(phi) {
    const v = [];
    for (const br of phi.branches) v.push(br.z.re, br.z.im);
    for (const br of phi.branches) for (const a of br.A) v.push(a.re, a.im);
    if (phi.polyA) for (const F of phi.polyA) v.push(F.re, F.im);
    return v;
  }

  function unpackPhi_UQD(v, template) {
    const phi = {
      unbounded: true,
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

  // No gauge ambiguity in the unbounded case (c > 0 pins the rotation).
  function canonicalizePhi_UQD(phi) { return phi; }

  // ===========================================================================
  // 6. Initial guesses
  // ===========================================================================
  function unboundedInitialGuess_UQD(hData, cUser) {
    let minA = Infinity;
    for (const p of hData.poles) {
      const m = Complex.abs(p.a);
      if (m < minA) minA = m;
    }
    const cap = isFinite(minA) && minA > 0 ? 0.5 * minA : Math.min(1, cUser);
    const effC = Math.min(cUser, cap);
    const phi = {
      unbounded: true, c: cUser, w0: undefined, polyA: [], branches: [],
    };
    for (const p of hData.poles) {
      let z;
      if (Complex.abs2(p.a) < 1e-30) {
        z = { re: 2, im: 0 };
      } else {
        z = Complex.scale(p.a, 1 / effC);
      }
      const A = [];
      let cPow = 1;
      for (let k = 1; k <= p.principal.length; k++) {
        cPow *= effC;
        A.push(Complex.scale(p.principal[k - 1], 1 / cPow));
      }
      phi.branches.push({ z, A });
    }
    if (hData.polyPart && hData.polyPart.length > 0) {
      let cPowL = 1;
      for (let l = 0; l < hData.polyPart.length; l++) {
        if (l > 0) cPowL *= cUser;
        phi.polyA.push(Complex.scale(hData.polyPart[l], cPowL));
      }
    }
    return phi;
  }

  function perturbedUnboundedInitialGuess_UQD(hData, c, rng, r = 0) {
    const base = unboundedInitialGuess_UQD(hData, c);
    const sigma = 0.15 + 0.25 * r;
    for (const br of base.branches) {
      br.z = {
        re: br.z.re + sigma * (rng() - 0.5),
        im: br.z.im + sigma * (rng() - 0.5)
      };
      const rr = Math.hypot(br.z.re, br.z.im);
      if (rr < 1.05) {
        const scale = 1.05 / Math.max(rr, 1e-9);
        br.z.re *= scale; br.z.im *= scale;
      }
      for (let k = 0; k < br.A.length; k++) {
        br.A[k] = {
          re: br.A[k].re * (1 + sigma * (rng() - 0.5)),
          im: br.A[k].im + sigma * (rng() - 0.5)
        };
      }
    }
    for (let l = 0; l < base.polyA.length; l++) {
      base.polyA[l] = {
        re: base.polyA[l].re * (1 + sigma * (rng() - 0.5))
              + (Math.abs(base.polyA[l].re) < 1e-9 ? sigma * (rng() - 0.5) : 0),
        im: base.polyA[l].im + sigma * (rng() - 0.5),
      };
    }
    return base;
  }

  // ===========================================================================
  // 7. Continuation in c
  // ===========================================================================
  function continuationInC_UQD(hData, cTarget, options = {}) {
    const {
      cStart       = null,
      growFactor   = 1.6,
      shrinkFactor = 0.5,
      minStep      = 1e-4,
      maxSteps     = 80,
      newton       = {},
    } = options;

    let minA = Infinity;
    for (const p of hData.poles) {
      const m = Complex.abs(p.a);
      if (m > 0 && m < minA) minA = m;
    }
    const startGuess = cStart ?? Math.min(cTarget, isFinite(minA) ? 0.25 * minA : 0.25);
    if (startGuess <= 0) {
      return { success: false, error: "continuationInC: invalid starting c", trace: [] };
    }

    const trace = [];
    let c = startGuess;
    let phi = unboundedInitialGuess_UQD(hData, c);

    let warmup;
    while (true) {
      warmup = QD.newtonSolve(phi, hData, newton);
      if (warmup.success) { phi = warmup.phi; break; }
      c *= shrinkFactor;
      if (c < minStep) {
        return {
          success: false,
          error: "continuationInC: warmup failed even at c=" + c.toExponential(2),
          phi: warmup.phi, trace,
        };
      }
      phi = unboundedInitialGuess_UQD(hData, c);
    }
    trace.push({ c, ok: true, residual: warmup.residual });

    if (c >= cTarget - 1e-12) {
      return { success: true, phi, iterations: 0, residual: warmup.residual, trace, method: "continuation-in-c" };
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
            error: "continuationInC: step underflow at c=" + lastSuccessC.toFixed(4) +
                   " (target c=" + cTarget.toFixed(4) + ")",
            phi, trace, lastC: lastSuccessC,
          };
        }
      }
    }
    if (lastSuccessC < cTarget - 1e-9) {
      return {
        success: false,
        error: "continuationInC: max steps reached at c=" + lastSuccessC.toFixed(4),
        phi, trace, lastC: lastSuccessC,
      };
    }
    return {
      success: true, phi, iterations: 0,
      residual: trace[trace.length - 1].residual,
      trace, method: "continuation-in-c",
    };
  }

  // ===========================================================================
  // 8. Identity verifier — test functions f(w) = 1/(w − b)^k for b ∈ K
  // ===========================================================================
  function verifyQuadratureIdentity_UQD(phi, hData, options = {}) {
    const N             = options.numSamples ?? 500;
    const maxOrder      = options.maxDegree ?? 3;
    const numTestPoints = options.numTestPoints ?? 3;

    const samples = new Array(N);
    for (let n = 0; n < N; n++) {
      const theta = 2 * Math.PI * n / N;
      const z = { re: Math.cos(theta), im: Math.sin(theta) };
      const taylor = phiTaylorAt_UQD(z, phi, 1);
      samples[n] = { z, w: taylor[0], phiPrime: taylor[1] };
    }

    let cx = 0, cy = 0;
    for (const s of samples) { cx += s.w.re; cy += s.w.im; }
    cx /= N; cy /= N;
    const centroid = { re: cx, im: cy };
    let maxDev = 0;
    for (const s of samples) {
      const d = Math.hypot(s.w.re - cx, s.w.im - cy);
      if (d > maxDev) maxDev = d;
    }
    const testPoints = [centroid];
    for (let i = 1; i < numTestPoints; i++) {
      const ang = 2 * Math.PI * (i - 1) / Math.max(numTestPoints - 1, 1);
      const r = 0.18 * maxDev;
      testPoints.push({ re: cx + r * Math.cos(ang), im: cy + r * Math.sin(ang) });
    }

    let areaScale = 0;
    for (const pole of hData.poles) {
      if (pole.principal.length > 0) areaScale += Complex.abs(pole.principal[0]);
    }
    if (areaScale === 0) areaScale = 1;

    const checks = [];
    let maxRelDiff = 0;
    let maxAbsDiff = 0;

    for (let pIdx = 0; pIdx < testPoints.length; pIdx++) {
      const b = testPoints[pIdx];
      for (let k = 1; k <= maxOrder; k++) {
        let lhs = { re: 0, im: 0 };
        for (let n = 0; n < N; n++) {
          const s = samples[n];
          const diff = Complex.sub(s.w, b);
          const dPow = Complex.pow(diff, k);
          const fVal = Complex.inv(dPow);
          let term = Complex.mul(fVal, Complex.conj(s.w));
          term = Complex.mul(term, s.phiPrime);
          term = Complex.mul(term, s.z);
          lhs = Complex.add(lhs, term);
        }
        lhs = Complex.scale(lhs, -1 / N);

        let rhs = { re: 0, im: 0 };
        for (const pole of hData.poles) {
          const aMinusB = Complex.sub(pole.a, b);
          for (let sIdx = 0; sIdx < pole.principal.length; sIdx++) {
            const s = sIdx + 1;
            const C = pole.principal[sIdx];
            const sign = (s % 2 === 0) ? -1 : 1;
            const coef = QD.binomialCoeff(k + s - 2, s - 1);
            const expon = k + s - 1;
            const denom = Complex.pow(aMinusB, expon);
            const term = Complex.div(C, denom);
            rhs = Complex.add(rhs, Complex.scale(term, sign * coef));
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

        const diff = Complex.sub(lhs, rhs);
        const absDiff = Complex.abs(diff);
        const scale = Math.max(Complex.abs(lhs), Complex.abs(rhs), areaScale);
        const relDiff = absDiff / scale;
        if (relDiff > maxRelDiff) maxRelDiff = relDiff;
        if (absDiff > maxAbsDiff) maxAbsDiff = absDiff;
        checks.push({ bIdx: pIdx, k, lhs, rhs, absDiff, relDiff });
      }
    }

    return {
      checks, maxRelDiff, maxAbsDiff, areaScale, testPoints,
      maxDeg: maxOrder, numSamples: N, unbounded: true,
    };
  }

  // ===========================================================================
  // 9. Register Family.unboundedQD
  // ===========================================================================
  QD.Family.unboundedQD = {
    name: 'unboundedQD',
    enforceInDisk:  false,
    enforceOutDisk: true,
    matches(opts) { return !!(opts && opts.unbounded); },

    normalizeOpts(opts, hData) {
      const c = opts.c;
      if (typeof c !== 'number' || !(c > 0)) {
        throw new Error("Family.unboundedQD: opts.c must be a positive number");
      }
      return { unbounded: true, c };
    },

    evalPhi: evalPhi_UQD,
    phiTaylorAt: phiTaylorAt_UQD,
    computeTargets(phi, hData) {
      return { A: computeTargetA_UQD(phi, hData), F: computeTargetF_UQD(phi, hData) };
    },
    residual: residual_UQD,
    packPhi: packPhi_UQD,
    unpackPhi: unpackPhi_UQD,
    canonicalizePhi: canonicalizePhi_UQD,
    initialGuess(hData, norm) { return unboundedInitialGuess_UQD(hData, norm.c); },
    perturbedInitialGuess(hData, norm, rng, r) {
      return perturbedUnboundedInitialGuess_UQD(hData, norm.c, rng, r);
    },
    diverseInitialGuess(hData, norm, rng, r) {
      return QD.diverseInitialGuess(hData, norm, rng, r);
    },
    continuationSolve(hData, norm, opts) {
      return continuationInC_UQD(hData, norm.c, opts);
    },
    verifyQuadratureIdentity: verifyQuadratureIdentity_UQD,
  };
  QD.registerFamily('unboundedQD');

  // Exports
  QD.unboundedInitialGuess          = unboundedInitialGuess_UQD;
  QD.perturbedUnboundedInitialGuess = perturbedUnboundedInitialGuess_UQD;
  QD.continuationInC                = continuationInC_UQD;
  QD.phiLaurentAtInfinity           = phiLaurentAtInfinity_UQD;
  QD.computeTargetF                 = computeTargetF_UQD;
  QD.verifyQuadratureIdentityUnbounded = verifyQuadratureIdentity_UQD;

})();
