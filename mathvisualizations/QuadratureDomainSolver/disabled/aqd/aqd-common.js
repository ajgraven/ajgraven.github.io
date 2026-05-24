// =============================================================================
// aqd-common.js -- Shared infrastructure for Algebraic Quadrature Domains.
//
// AQDs are domains Ω with weight ρ = |α|² where α = R′ for some rational
// primitive R. (LQDs, with α = w⁻¹, are NOT AQDs because R = log w is not
// rational; they live in their own solver tab.)
//
// Data shape for R (mirrors hData exactly so the same UI widget can drive
// both inputs):
//
//   R = {
//     polyPart: [c_0, c_1, ..., c_{m_∞}],            // R_poly(w) = Σ c_k w^k
//     poles:    [ { a: Complex,
//                   principal: [D_1, D_2, ..., D_m]   // Σ D_s / (w − a)^s
//                 }, ...
//                ],
//   }
//
// Note: with R rational of this shape, α = R′ AUTOMATICALLY has zero residue
// at every finite pole (R′ contributes only orders ≥ 2 to each pole), which
// is precisely the AQD admissibility condition. No validator needed.
//
// The unifying inverse-problem fact (Thm 6.4.1, bounded case, no R-poles in Ω̄):
//
//     (R∘φ)(z) = R(w₀) + r#(z)         on |z| = 1,
//
// where r(z) = Φφ⁻¹(AnalyticIn[h/R′, Ω⊥])(z). The parametric unknown is the
// rational function R∘φ : 𝔻̄ → ℂ. φ itself is recovered numerically by
// inverting R at each evaluation point.
//
// Stage 1 status: math kernel for R (eval, derivative, Schwarz reflection,
// Newton inversion). Composition with φ and principal-part extraction land
// in Stage 2 once we have the R∘φ ansatz nailed down.
// =============================================================================
'use strict';

(function (global) {

  const QD = global.QD || (global.QD = {});
  const Aqd = QD.Aqd || (QD.Aqd = {});
  Aqd.version = '0.1.0-stage1';

  // We reuse the existing Complex namespace. In the browser it's a global;
  // in node-test.js it's pulled from the vm context. Either way: just refer
  // to `Complex` and `Taylor` as if they were globals.
  const C = (typeof Complex !== 'undefined') ? Complex
          : (typeof global.Complex !== 'undefined' ? global.Complex : null);
  if (!C) throw new Error("aqd-common.js: Complex (complex.js) must be loaded first");

  // ===========================================================================
  // R evaluation
  // ===========================================================================
  //   R(w) = Σ_k polyPart[k] · w^k  +  Σ_j Σ_s principal_j[s-1] / (w − a_j)^s
  // Polynomial part is Horner from highest down.
  function evalR(R, w) {
    let result = { re: 0, im: 0 };

    // Polynomial part by Horner.
    const poly = R.polyPart || [];
    if (poly.length > 0) {
      result = C.clone(poly[poly.length - 1]);
      for (let k = poly.length - 2; k >= 0; k--) {
        result = C.add(C.mul(result, w), poly[k]);
      }
    }

    // Pole contributions.
    for (const pole of (R.poles || [])) {
      const diff = C.sub(w, pole.a);
      const d2 = diff.re * diff.re + diff.im * diff.im;
      if (d2 < 1e-300) return { re: NaN, im: NaN };   // exactly at a pole
      // 1/diff, 1/diff^2, ...
      let invPow = C.inv(diff);
      const invStep = C.clone(invPow);
      for (let s = 1; s <= pole.principal.length; s++) {
        result = C.add(result, C.mul(pole.principal[s - 1], invPow));
        if (s < pole.principal.length) invPow = C.mul(invPow, invStep);
      }
    }
    return result;
  }

  // ===========================================================================
  // α = R'(w)
  // ===========================================================================
  //   R'(w) = Σ_k k · polyPart[k] · w^{k-1}
  //         + Σ_j Σ_s (-s) · principal_j[s-1] / (w − a_j)^{s+1}
  function evalAlpha(R, w) {
    let result = { re: 0, im: 0 };

    const poly = R.polyPart || [];
    if (poly.length > 1) {
      // Build the derivative polynomial in place via Horner.
      result = C.scale(poly[poly.length - 1], poly.length - 1);
      for (let k = poly.length - 2; k >= 1; k--) {
        result = C.add(C.mul(result, w), C.scale(poly[k], k));
      }
    }

    for (const pole of (R.poles || [])) {
      const diff = C.sub(w, pole.a);
      const d2 = diff.re * diff.re + diff.im * diff.im;
      if (d2 < 1e-300) return { re: NaN, im: NaN };
      // 1/diff^2, 1/diff^3, ...
      const invDiff = C.inv(diff);
      let invPow = C.mul(invDiff, invDiff);             // start at (w-a)^{-2}
      for (let s = 1; s <= pole.principal.length; s++) {
        result = C.add(result, C.scale(C.mul(pole.principal[s - 1], invPow), -s));
        if (s < pole.principal.length) invPow = C.mul(invPow, invDiff);
      }
    }
    return result;
  }

  // ===========================================================================
  // Schwarz reflection R#(z) := conj(R(1/conj(z)))
  // ===========================================================================
  // Useful for assembling boundary identities. For z on |z|=1 (z̄ = 1/z), this
  // equals conj(R(z)).
  function evalRSharp(R, z) {
    const d2 = z.re * z.re + z.im * z.im;
    if (d2 < 1e-300) {
      // Limit at z=0 corresponds to R(∞); use the polynomial leading behavior.
      const poly = R.polyPart || [];
      if (poly.length <= 1) {
        // R(∞) = polyPart[0] (constant); the principal parts go to 0 at ∞.
        return C.conj(poly[0] || { re: 0, im: 0 });
      }
      // R has positive-degree polynomial part → R(∞) = ∞.
      return { re: Infinity, im: 0 };
    }
    // w = 1/conj(z) = z / (z·conj(z)) = z / |z|²  -- but actually
    //   1/conj(z) = (z) / |z|^2  ... wait: conj(z) = (re, -im), so
    //   1/conj(z) = conj(z̄)^{-1}  = ... cleanest: w_re + i w_im = 1/(re - i·im)
    //   = (re + i·im) / (re² + im²)  = z / |z|².
    const w = { re: z.re / d2, im: z.im / d2 };
    const Rw = evalR(R, w);
    return C.conj(Rw);
  }

  // ===========================================================================
  // Newton inversion of R: given target T, solve R(w) = T from seed.
  // ===========================================================================
  // Returns { w, converged: bool, iters: int }. Converged when both
  //   |R(w) − T| < absTol  AND  |Δw| < relTol · max(|w|, 1).
  // Diverges-out (|w| > maxAbs) or |R'(w)| < jacFloor aborts.
  function invertR(R, target, seed, opts = {}) {
    const maxIter = opts.maxIter ?? 32;
    const absTol  = opts.absTol  ?? 1e-13;
    const relTol  = opts.relTol  ?? 1e-13;
    const jacFloor = opts.jacFloor ?? 1e-14;
    const maxAbs  = opts.maxAbs  ?? 1e8;

    let w = C.clone(seed);
    for (let iter = 0; iter < maxIter; iter++) {
      const Rw = evalR(R, w);
      const diff = C.sub(Rw, target);
      const absDiff = C.abs(diff);

      const Rp = evalAlpha(R, w);
      const absRp = C.abs(Rp);
      if (!isFinite(absRp) || absRp < jacFloor) {
        // Saddle: try a tiny perturbation; if it persists, bail.
        w = { re: w.re + 1e-8, im: w.im - 1e-8 };
        continue;
      }
      const step = C.div(diff, Rp);
      const absStep = C.abs(step);
      const scale = Math.max(C.abs(w), 1);

      // Update before convergence check (so converged step also gets applied).
      w = C.sub(w, step);
      if (!isFinite(w.re) || !isFinite(w.im) || C.abs(w) > maxAbs) {
        return { w, converged: false, iters: iter + 1, reason: 'diverged' };
      }
      if (absDiff < absTol && absStep < relTol * scale) {
        return { w, converged: true, iters: iter + 1 };
      }
    }
    return { w, converged: false, iters: maxIter, reason: 'maxIter' };
  }

  // ===========================================================================
  // Sanity helpers (used by tests and the UI to surface bad input).
  // ===========================================================================

  // Sum of pole orders (= number of finite pole DOF) plus polynomial-part DOF.
  function rationalDegreeOfFreedom(R) {
    let dof = (R.polyPart?.length ?? 0);
    for (const p of (R.poles || [])) dof += p.principal.length;
    return dof;
  }

  // Is R the identity R(w) = w?  This corresponds to a classical QD; the user
  // is better served by the QD/LQD tab.
  function isIdentityR(R) {
    const poly = R.polyPart || [];
    const polesOk = (R.poles?.length ?? 0) === 0;
    const polyOk = poly.length === 2
                && Math.hypot(poly[0].re, poly[0].im) < 1e-12
                && Math.abs(poly[1].re - 1) < 1e-12
                && Math.abs(poly[1].im)     < 1e-12;
    return polesOk && polyOk;
  }

  // ===========================================================================
  // Expose
  // ===========================================================================
  Aqd.evalR        = evalR;
  Aqd.evalAlpha    = evalAlpha;
  Aqd.evalRSharp   = evalRSharp;
  Aqd.invertR      = invertR;
  Aqd.rationalDegreeOfFreedom = rationalDegreeOfFreedom;
  Aqd.isIdentityR  = isIdentityR;

  // ===========================================================================
  // Taylor of R around a point w₀: [R(w₀), R'(w₀)/1!, ..., R^(L)(w₀)/L!]
  // ===========================================================================
  // Polynomial part: contribute via binomial expansion of (w₀ + h)^k.
  // Pole contribution D/(w − a)^s = D · (w − a)^{-s}: expand around w₀ with
  //   1/(w₀ − a + h)^s = (w₀ − a)^{-s} · Σ_l C(-s, l) (h/(w₀-a))^l
  //                   = (w₀ − a)^{-s} · Σ_l (-1)^l · binom(s+l-1, l) (h/(w₀-a))^l.
  // (Requires w₀ to not be a pole of R.)
  function rTaylorAt(R, w0, L) {
    const out = Taylor.zero(L + 1);

    // Polynomial part
    const poly = R.polyPart || [];
    if (poly.length > 0) {
      for (let k = 0; k < poly.length; k++) {
        const ck = poly[k];
        if (ck.re === 0 && ck.im === 0) continue;
        // Contribution: c_k · Σ_{l=0..min(k,L)} binom(k,l) · w0^{k-l} · h^l
        const lMax = Math.min(k, L);
        let wPow = C.pow(w0, k - lMax);             // start at w0^{k-lMax}
        // Roll: at l, factor is binom(k,l) · w0^{k-l}. We start from l=lMax
        // and work down so we can multiply by w0 each step.
        let binom = (function () {
          let b = 1;
          for (let i = 0; i < lMax; i++) b = b * (k - i) / (i + 1);
          return b;
        })();
        for (let l = lMax; l >= 0; l--) {
          const term = C.scale(C.mul(ck, wPow), binom);
          out[l] = C.add(out[l], term);
          if (l > 0) {
            wPow = C.mul(wPow, w0);                 // w0^{k-(l-1)} = w0 · prev
            binom = binom * l / (k - l + 1);
          }
        }
      }
    }

    // Pole contributions
    for (const pole of (R.poles || [])) {
      const wma = C.sub(w0, pole.a);
      if (C.abs2(wma) < 1e-300) {
        return null;                                // w0 is at this R-pole
      }
      const wmaInv = C.inv(wma);
      for (let s = 1; s <= pole.principal.length; s++) {
        const D = pole.principal[s - 1];
        if (D.re === 0 && D.im === 0) continue;
        // base = D · (w0 - a)^{-s}
        let base = C.mul(D, C.pow(wmaInv, s));
        let binom = 1;
        for (let l = 0; l <= L; l++) {
          const sign = (l % 2 === 0) ? 1 : -1;
          out[l] = C.add(out[l], C.scale(base, sign * binom));
          if (l < L) {
            base = C.mul(base, wmaInv);             // (w0-a)^{-(s+l+1)}
            binom = binom * (s + l) / (l + 1);      // binom(s+l, l+1) = ...
          }
        }
      }
    }

    return out;
  }

  // ===========================================================================
  // Taylor of 1/R'(w) at a point a (Stage 2 requires this at each h-pole).
  // Built from rTaylorAt of length L+2 (to get R^(1)..R^(L+1)) and Taylor
  // reciprocal.
  // ===========================================================================
  function oneOverRprimeTaylorAt(R, a, L) {
    // Need R'(a + h) Taylor up to h^L. R'(a+h)'s [h^l] coeff = R^(l+1)(a)/l!.
    // From rTaylorAt(R, a, L+1): [h^{l+1}] coeff = R^(l+1)(a)/(l+1)!.
    // So R'(a+h)[h^l] = (l+1) · rTaylorAt(R, a, L+1)[l+1].
    const rT = rTaylorAt(R, a, L + 1);
    if (!rT) return null;                           // a is an R-pole
    const Rp = Taylor.zero(L + 1);
    for (let l = 0; l <= L; l++) {
      Rp[l] = C.scale(rT[l + 1], l + 1);
    }
    if (C.abs2(Rp[0]) < 1e-30) return null;         // R'(a) = 0 (R'-zero)
    return Taylor.reciprocal(Rp, L);
  }

  // ===========================================================================
  // Principal-part preprocessing for AQDs: residues of h/R' at h-pole a_j.
  // ===========================================================================
  //   h(w) = Σ_{s=1..m} C_{j,s} / (w − a_j)^s + (regular at a_j)
  //   1/R'(w) = Σ_{k≥0} g_{j,k} (w − a_j)^k                     [Taylor at a_j]
  //   (h/R')(w) = Σ_{s,k} C_{j,s} · g_{j,k} · (w − a_j)^{k − s} + (regular)
  //   Principal part: coefficient of (w − a_j)^{-t} for t ≥ 1:
  //       P_{j,t} = Σ_{s=t..m} C_{j,s} · g_{j, s − t}.
  // Returns [P_{j,1}, ..., P_{j,m}] complex array, or null if R has a pole or
  // R' has a zero AT a_j (which would violate Stage 2 admissibility).
  function principalPartsHOverRprime(hPole, R) {
    const a = hPole.a;
    const m = hPole.principal.length;
    if (m === 0) return [];
    // Need g_0..g_{m-1}, i.e., 1/R' Taylor up to order m-1.
    const g = oneOverRprimeTaylorAt(R, a, m - 1);
    if (!g) return null;
    const P = new Array(m);
    for (let t = 1; t <= m; t++) {
      let acc = { re: 0, im: 0 };
      for (let s = t; s <= m; s++) {
        const Cs = hPole.principal[s - 1];
        const gk = g[s - t];
        acc = C.add(acc, C.mul(Cs, gk));
      }
      P[t - 1] = acc;
    }
    return P;
  }

  Aqd.rTaylorAt                = rTaylorAt;
  Aqd.oneOverRprimeTaylorAt    = oneOverRprimeTaylorAt;
  Aqd.principalPartsHOverRprime = principalPartsHOverRprime;

  // ---------------------------------------------------------------------------
  // Stage 2 primitives (to come):
  //   Aqd.composeRphi(R, phiData, z)        -> Complex   (analytic kernel)
  //   Aqd.principalPartsInOmegaClosure(...) -> data for inverse Faber transform
  //   Aqd.boundarySampleAQD(...)            -> per-point R-inversion
  // ---------------------------------------------------------------------------

  // Node-test interop: expose on module.exports if running under Node.
  if (typeof module !== 'undefined' && module.exports) {
    module.exports.Aqd = Aqd;
  }

}(typeof window !== 'undefined' ? window
   : typeof global   !== 'undefined' ? global
   : globalThis));
