// =============================================================================
// schwarz-common.js -- Schwarz-reflection dynamics core.
//
// Given a Riemann map φ (from the inverse-tab solver), this module builds the
// Schwarz function S(w) = F(ψ(w)) and the anti-meromorphic reflection
// σ(w) = conj(S(w)), where ψ = φ⁻¹ and F is the "Schwarz extension" of φ
// across the unit circle.
//
// All six inverse families are supported. Each adapter constructs φ, φ', and
// F in closed form from the inverse-Faber data on `phi`:
//
//   boundedQD             φ = w₀ + Σⱼ Σₖ conj(A_{j,k})·u_j(z)^k
//                         F = conj(w₀) + R##(z)
//                         where u_j(z) = z/(1-conj(z_j)z) and
//                         R##(z) := Σⱼ Σₖ A_{j,k}/(z-z_j)^k.
//
//   unboundedQD           φ = c·z + Σ_l F_l/z^l + Σ branches
//                         F = (c/z) + Σ_l conj(F_l)·z^l + R##(z)
//                         (Laurent at ∞; phi.polyA holds [F_0..F_{m-1}].)
//
//   boundedQDRational     φ = P(z)/Q(z)  (Direct-tab entry only)
//                         F = z^{q-p} · P̃(z)/Q̃(z)   (reverse-conjugate trick)
//
//   boundedLQD            φ = w₀ · exp(r#(z))            r# = Σ branches
//                         F = conj(w₀) · exp(R##(z))
//
//   boundedLQD_singular   φ = γ · b_{z₀}(z) · exp(r#(z))      z₀ ∈ 𝔻
//                         F = conj(γ) · b#_{z₀}(z) · exp(R##(z))
//
//   unboundedLQD          φ = c·z · exp(r#(z) − r#(∞))
//                         F = (c/z) · exp(R##(z) − conj(r#(∞)))
//
//   unboundedLQD_singular φ = c·|z₀|·z·b_{z₀}(z)·exp(r#(z)−r#(∞))   z₀ ∈ 𝔻*
//                         F = (c·|z₀|/z)·b#_{z₀}(z)·exp(R##(z)−conj r#(∞))
//
// Iteration:  starting at w₀ ∈ Ω, repeatedly apply σ. Track the smallest n for
// which σⁿ(w₀) ∉ Ω (lands in the "fundamental tile" Ω^c, escape time = n);
// for unbounded Ω also track diverging orbits |σⁿ(w₀)| → ∞ ("escaping set").
//
// "In Ω" is tested by point-in-polygon against the sampled ∂Ω of the inverse
// solution. Cheap and unambiguous given the curve we already have.
//
// API (on QD.Schwarz):
//   buildSchwarzFromPhi(phi, hData, boundaryPts)
//     → { sigma(w, seedHint?) → {re, im, _z} | null,
//         psi(w, seedHint?)   → {re, im} | null,
//         evalPhi, evalF, isInOmega, escapeR, family, unbounded, adapter }
//   buildSchwarzFromRational(phi, boundaryPts) → same shape (rational P/Q path)
//   pointInPolygon(pt, polyPts) → bool   (even-odd; polyPts may be CCW or CW)
//   polygonBounds(polyPts) → { center, radius }
//   escapeTime(w0, schwarz, opts) → { kind, n, lastW }
//       kind: 'fundamental' (σⁿ left Ω), 'escaped' (|σⁿ| > escapeR; unbounded
//              only), 'interior' (still in Ω after maxIter), 'invalid'
//              (Newton failed; ψ undefined)
//   makeOrbit(w0, schwarz, opts) → [{re,im}, ...]   for click-to-orbit
// =============================================================================

(function (global) {
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' && module.exports ? module.exports : (global.QD || (global.QD = {})));

  const C =
    (typeof Complex !== 'undefined' ? Complex : null)
    || (QD && QD.Complex)
    || (typeof window !== 'undefined' && window.Complex)
    || (global && global.Complex);
  if (!C) throw new Error("schwarz-common.js: Complex namespace not found");

  const Schwarz = QD.Schwarz || (QD.Schwarz = {});

  // ===========================================================================
  // Polynomial helpers (ascending-power Complex[]).
  // ===========================================================================
  function polyEval(coeffs, z) {
    // Horner from the top.
    if (!coeffs.length) return { re: 0, im: 0 };
    let v = { re: coeffs[coeffs.length - 1].re, im: coeffs[coeffs.length - 1].im };
    for (let k = coeffs.length - 2; k >= 0; k--) {
      v = C.add(C.mul(v, z), coeffs[k]);
    }
    return v;
  }
  function polyDerivEval(coeffs, z) {
    // Evaluate Σ k·c_k·z^{k-1} via Horner on the derivative coefficients.
    if (coeffs.length < 2) return { re: 0, im: 0 };
    const d = coeffs.length - 1;
    let v = C.scale(coeffs[d], d);
    for (let k = d - 1; k >= 1; k--) {
      v = C.add(C.mul(v, z), C.scale(coeffs[k], k));
    }
    return v;
  }

  // ===========================================================================
  // Newton in z for φ(z) = w, where φ is given by (evalPhi, derivPhi).
  // Returns null if convergence fails.
  // ===========================================================================
  function newtonInvert(evalPhi, derivPhi, w, zSeed, opts) {
    opts = opts || {};
    const maxIter = opts.maxIter || 40;
    const tol     = opts.tol     || 1e-12;
    let z = { re: zSeed.re, im: zSeed.im };
    for (let it = 0; it < maxIter; it++) {
      const fz = C.sub(evalPhi(z), w);
      if (Math.hypot(fz.re, fz.im) < tol) return z;
      const dfz = derivPhi(z);
      const denom = Math.hypot(dfz.re, dfz.im);
      if (denom < 1e-300) return null;
      const step = C.div(fz, dfz);
      z = C.sub(z, step);
      // Diverging guard
      if (!isFinite(z.re) || !isFinite(z.im) || Math.hypot(z.re, z.im) > 1e8) return null;
    }
    // Final residual check
    const fz = C.sub(evalPhi(z), w);
    if (Math.hypot(fz.re, fz.im) < tol * 100) return z;
    return null;
  }

  // ===========================================================================
  // Point-in-polygon (even-odd rule), tolerant of CCW or CW orientation.
  //
  // Two flavors:
  //   • pointInPolygon(pt, polyPts)              — naïve O(N) per query.
  //   • pointInPolygonIndexed(pt, idx)           — O(N / nBins) via a
  //     y-binned edge list built once by buildPolygonIndex.
  //
  // The CPU Schwarz renderer calls the test once per pixel × maxIter
  // (~256² × 64 ≈ 4M calls per frame); with ~200-vertex boundaries the
  // indexed path is ~6–10× faster end-to-end. The naïve path is kept for
  // one-off tests and backward compat.
  // ===========================================================================
  function pointInPolygon(pt, polyPts) {
    let inside = false;
    const x = pt.re, y = pt.im;
    for (let i = 0, j = polyPts.length - 1; i < polyPts.length; j = i++) {
      const yi = polyPts[i].im, xi = polyPts[i].re;
      const yj = polyPts[j].im, xj = polyPts[j].re;
      const crosses = ((yi > y) !== (yj > y)) &&
        (x < (xj - xi) * (y - yi) / (yj - yi || 1e-300) + xi);
      if (crosses) inside = !inside;
    }
    return inside;
  }

  // Build a y-binned edge index over a closed polygon. Each bin lists the
  // edges whose y-extent crosses that bin; a ray-cast query only needs to
  // walk the bin containing the test point's y.
  //
  // For the polygons we typically deal with (~50–500 boundary samples) a
  // 32-bin grid keeps each bin's edge count to ~3–30 — small enough that
  // the inner loop fits in L1 cache.
  function buildPolygonIndex(polyPts, nBins) {
    nBins = nBins || 32;
    const N = polyPts.length;
    let yMin = Infinity, yMax = -Infinity, xMin = Infinity, xMax = -Infinity;
    for (let i = 0; i < N; i++) {
      const p = polyPts[i];
      if (p.im < yMin) yMin = p.im;
      if (p.im > yMax) yMax = p.im;
      if (p.re < xMin) xMin = p.re;
      if (p.re > xMax) xMax = p.re;
    }
    const span = (yMax - yMin) || 1;
    const dy = span / nBins;
    // Flat-packed bins: edges_i = polygon-edge-vertex index pair stored as
    // adjacent ints (i, j). This avoids per-edge object allocation.
    const binEdges = new Array(nBins);
    for (let b = 0; b < nBins; b++) binEdges[b] = [];
    for (let i = 0, j = N - 1; i < N; j = i++) {
      const yi = polyPts[i].im, yj = polyPts[j].im;
      const lo = (yi < yj) ? yi : yj;
      const hi = (yi < yj) ? yj : yi;
      let b0 = Math.floor((lo - yMin) / dy);
      let b1 = Math.floor((hi - yMin) / dy);
      if (b0 < 0) b0 = 0;
      if (b1 >= nBins) b1 = nBins - 1;
      for (let b = b0; b <= b1; b++) { binEdges[b].push(i); binEdges[b].push(j); }
    }
    return { polyPts, binEdges, nBins, yMin, yMax, xMin, xMax, dy };
  }

  function pointInPolygonIndexed(pt, idx) {
    const y = pt.im;
    if (y < idx.yMin || y > idx.yMax) return false;
    // x-extent rejection (rare hit for typical convex-ish boundaries, but free).
    if (pt.re < idx.xMin || pt.re > idx.xMax) return false;
    let bin = Math.floor((y - idx.yMin) / idx.dy);
    if (bin < 0) bin = 0;
    if (bin >= idx.nBins) bin = idx.nBins - 1;
    const edges = idx.binEdges[bin];
    const pts = idx.polyPts;
    const x = pt.re;
    let inside = false;
    for (let k = 0; k < edges.length; k += 2) {
      const i = edges[k], j = edges[k + 1];
      const pi = pts[i], pj = pts[j];
      const yi = pi.im, yj = pj.im;
      if ((yi > y) !== (yj > y)) {
        const xi = pi.re, xj = pj.re;
        if (x < (xj - xi) * (y - yi) / (yj - yi || 1e-300) + xi) inside = !inside;
      }
    }
    return inside;
  }

  // Polygon centroid + radius (rough), used for escape-radius defaults.
  function polygonBounds(polyPts) {
    let cx = 0, cy = 0;
    for (const p of polyPts) { cx += p.re; cy += p.im; }
    cx /= polyPts.length; cy /= polyPts.length;
    let r = 0;
    for (const p of polyPts) {
      const d = Math.hypot(p.re - cx, p.im - cy);
      if (d > r) r = d;
    }
    return { center: { re: cx, im: cy }, radius: r };
  }

  // For unbounded φ (Laurent), Newton seed: if |w| is large, z ≈ w/c is good.
  // For w closer to K, fall back to last-seed or boundary heuristic.
  function bestSeedUnbounded(phi, w, lastSeed) {
    if (lastSeed && Math.hypot(lastSeed.re, lastSeed.im) > QD.DISK_CLAMP_OUT) return lastSeed;
    const c = phi.c || 1;
    const cand = { re: w.re / c, im: w.im / c };
    if (Math.hypot(cand.re, cand.im) > 1.05) return cand;
    // Pull outside the disk along the same ray.
    const m = Math.hypot(cand.re, cand.im);
    if (m < 1e-12) return { re: 1.1, im: 0 };
    return { re: cand.re * 1.1 / m, im: cand.im * 1.1 / m };
  }

  // ===========================================================================
  // Family adapters: each returns { evalPhi, derivPhi, F, isDiskInside, seed }.
  // ===========================================================================
  // -- Shared inverse-Faber-branch sums --------------------------------------
  // For each branch j with z_j inside 𝔻 and principal-part coefficients
  // A_{j,1..m_j}, the contribution to φ is
  //   Σ_k conj(A_{j,k}) · u_j(z)^k     where  u_j(z) = z / (1 − conj(z_j)·z).
  // u_j'(z) = 1 / (1 − conj(z_j)·z)^2  (a clean chain-rule identity).
  // On |z|=1, conj(u_j(z))^k = 1 / (z − z_j)^k, which extends as the principal
  // part of the Schwarz function at the conjugate-reflected pole z_j ∈ 𝔻.
  function branchPhiContribution(phi, z) {
    let acc = { re: 0, im: 0 };
    const branches = phi.branches || [];
    for (const br of branches) {
      const zjC = C.conj(br.z);
      const denom = C.sub({ re: 1, im: 0 }, C.mul(zjC, z));
      if (Math.hypot(denom.re, denom.im) < 1e-300) continue;
      const u = C.div(z, denom);
      let uPow = { re: 1, im: 0 };
      for (let k = 0; k < br.A.length; k++) {
        uPow = C.mul(uPow, u);
        acc = C.add(acc, C.mul(C.conj(br.A[k]), uPow));
      }
    }
    return acc;
  }
  function branchPhiDeriv(phi, z) {
    // d/dz Σ_k conj(A_{j,k}) · u^k  =  Σ_k k · conj(A_{j,k}) · u^{k-1} · u'
    //                              =  (1/denom²) · Σ_k k · conj(A_{j,k}) · u^{k-1}.
    let acc = { re: 0, im: 0 };
    const branches = phi.branches || [];
    for (const br of branches) {
      const zjC = C.conj(br.z);
      const denom = C.sub({ re: 1, im: 0 }, C.mul(zjC, z));
      if (Math.hypot(denom.re, denom.im) < 1e-300) continue;
      const u = C.div(z, denom);
      const denom2 = C.mul(denom, denom);
      let uPowKm1 = { re: 1, im: 0 };               // u^{k-1}, start at k=1
      let inner = { re: 0, im: 0 };
      for (let k = 1; k <= br.A.length; k++) {
        inner = C.add(inner, C.mul(C.scale(C.conj(br.A[k - 1]), k), uPowKm1));
        uPowKm1 = C.mul(uPowKm1, u);
      }
      acc = C.add(acc, C.div(inner, denom2));
    }
    return acc;
  }
  // Schwarz-extension contribution of the branches on the OPPOSITE side of |z|=1.
  // For bounded (z ∈ 𝔻) the contribution is Σ_j Σ_k A_{j,k}/(z − z_j)^k.
  // For unbounded (z ∈ 𝔻*) the same formula is used; z_j ∈ 𝔻 stays in 𝔻,
  // so 1/(z − z_j)^k is regular on 𝔻*.
  function branchSchwarzContribution(phi, z) {
    let acc = { re: 0, im: 0 };
    const branches = phi.branches || [];
    for (const br of branches) {
      const d = C.sub(z, br.z);
      if (Math.hypot(d.re, d.im) < 1e-300) continue;
      const dInv = C.inv(d);
      let dInvPow = { re: 1, im: 0 };               // 1/(z−z_j)^k starts at k=0
      for (let k = 0; k < br.A.length; k++) {
        dInvPow = C.mul(dInvPow, dInv);             // → 1/(z−z_j)^{k+1}
        acc = C.add(acc, C.mul(br.A[k], dInvPow));
      }
    }
    return acc;
  }

  // -- Polynomial-h β-extension helpers (UQDL / UQDLS, HANDOFF #22) ----------
  // The unbounded LQD parametrization is
  //   φ(z) = c · z · exp(r̃#(z) + B(1/z))
  // where B(1/z) = Σ_{l=1..N} β_l / z^l encodes h's polynomial part. These
  // helpers mirror evalB_OverZ / friends in solver-uqd-lqd.js so the Schwarz
  // module can evaluate the full φ.
  // evalBOverZ is shared with LqdCommon (HANDOFF #27 dedupe); the derivative
  // and Schwarz-reflection variants are Schwarz-specific.
  const evalBOverZ = QD.LqdCommon.evalB_OverZ;
  // d/dz B(1/z) = -Σ_l l · β_l / z^{l+1}.
  function evalBOverZDeriv(phi, z) {
    const beta = phi.lqdBeta || [];
    if (beta.length === 0) return { re: 0, im: 0 };
    const zInv = C.inv(z);
    let pow = C.mul(zInv, zInv);                    // 1/z^2 (l=1 term)
    let acc = { re: 0, im: 0 };
    for (let l = 1; l <= beta.length; l++) {
      acc = C.add(acc, C.mul(C.scale(beta[l - 1], -l), pow));
      if (l < beta.length) pow = C.mul(pow, zInv);
    }
    return acc;
  }
  // Schwarz reflection of B(1/z): conj(B(1/conj(z))) = Σ_l conj(β_l) · z^l.
  // Polynomial in z (no constant term), evaluated via Horner from highest
  // degree down.
  function evalBConjOfZ(phi, z) {
    const beta = phi.lqdBeta || [];
    if (beta.length === 0) return { re: 0, im: 0 };
    let acc = { re: 0, im: 0 };
    for (let l = beta.length; l >= 1; l--) {
      acc = C.add(C.mul(acc, z), C.conj(beta[l - 1]));
    }
    return C.mul(acc, z);                            // multiply by z once (no β_0 constant)
  }

  // γ-merge: prepend a synthetic branch at z=z₀ holding phi.lqdGamma. Mirrors
  // _phiWithSyntheticBranch in solver-uqd-lqd-singular.js. With this merge,
  // the existing branchPhiContribution / branchPhiDeriv /
  // branchSchwarzContribution helpers automatically pick up the synth
  // contribution required for higher-order-pole-at-origin support (HANDOFF #24).
  function withSyntheticBranch(phi) {
    const gamma = phi.lqdGamma || [];
    if (gamma.length === 0 || !phi.z0) return phi;
    return {
      ...phi,
      branches: (phi.branches || []).concat([{
        z: { re: phi.z0.re, im: phi.z0.im },
        A: gamma.map(c => ({ re: c.re, im: c.im })),
      }]),
    };
  }

  function adaptBounded(phi) {
    // φ(z) = w_0 + Σ branches contribution.
    const w0 = phi.w0 || { re: 0, im: 0 };
    function evalPhi(z) {
      return C.add(w0, branchPhiContribution(phi, z));
    }
    function derivPhi(z) {
      return branchPhiDeriv(phi, z);
    }
    // F(z) = conj(w_0) + Σ_j Σ_k A_{j,k} / (z − z_j)^k.
    // Meromorphic on 𝔻 with poles at the z_j (all inside 𝔻).
    function evalF(z) {
      return C.add(C.conj(w0), branchSchwarzContribution(phi, z));
    }
    // Estimate "linearization at 0" for seeding: φ'(0) = Σ k=1 conj(A_{j,1}).
    // (Each branch contributes conj(A_{j,1}) at first order.)
    let dphi0 = { re: 0, im: 0 };
    for (const br of phi.branches || []) {
      if (br.A.length > 0) dphi0 = C.add(dphi0, C.conj(br.A[0]));
    }
    function seedFor(w, last) {
      if (last && Math.hypot(last.re, last.im) < 0.9999) return last;
      const denomMag = Math.hypot(dphi0.re, dphi0.im);
      if (denomMag > 1e-12) {
        const cand = C.div(C.sub(w, w0), dphi0);
        const m = Math.hypot(cand.re, cand.im);
        if (m < 0.95) return cand;
        return { re: cand.re * 0.9 / m, im: cand.im * 0.9 / m };
      }
      return { re: 0, im: 0 };
    }
    return {
      family: 'boundedQD',
      w0,
      evalPhi, derivPhi, evalF,
      seedFor,
      acceptZ: (z) => Math.hypot(z.re, z.im) < 1 - 1e-9,
    };
  }

  function adaptUnbounded(phi) {
    const c = phi.c;
    // Inverse-solver stores Laurent F_l as `phi.polyA`; some Direct-tab paths
    // use `F`. Accept either.
    const F = phi.polyA || phi.F || [];
    const m = F.length;
    // φ(z) = c·z + Σ_{l=0..m-1} F_l / z^l + Σ branches contribution.
    function evalPhi(z) {
      let acc = C.scale(z, c);
      if (m > 0) {
        const zInv = C.inv(z);
        let zInvPow = { re: 1, im: 0 };
        for (let l = 0; l < m; l++) {
          acc = C.add(acc, C.mul(F[l], zInvPow));
          zInvPow = C.mul(zInvPow, zInv);
        }
      }
      acc = C.add(acc, branchPhiContribution(phi, z));
      return acc;
    }
    function derivPhi(z) {
      // φ'(z) = c − Σ_{l=1..m-1} l·F_l / z^{l+1} + branch deriv.
      let acc = { re: c, im: 0 };
      if (m > 1) {
        const zInv = C.inv(z);
        let zInvPow = C.mul(zInv, zInv);            // z^{-2}
        for (let l = 1; l < m; l++) {
          acc = C.sub(acc, C.mul(C.scale(F[l], l), zInvPow));
          zInvPow = C.mul(zInvPow, zInv);
        }
      }
      acc = C.add(acc, branchPhiDeriv(phi, z));
      return acc;
    }
    // G(z) (Schwarz extension to 𝔻*):
    //   G(z) = c/z + Σ_{l=0..m-1} conj(F_l) · z^l + Σ branches contribution.
    function evalF(z) {
      const zInv = C.inv(z);
      let acc = C.scale(zInv, c);
      let zPow = { re: 1, im: 0 };
      for (let l = 0; l < m; l++) {
        acc = C.add(acc, C.mul(C.conj(F[l]), zPow));
        zPow = C.mul(zPow, z);
      }
      acc = C.add(acc, branchSchwarzContribution(phi, z));
      return acc;
    }
    return {
      family: 'unboundedQD',
      c, F,
      evalPhi, derivPhi, evalF,
      seedFor: (w, last) => bestSeedUnbounded(phi, w, last),
      acceptZ: (z) => Math.hypot(z.re, z.im) > 1 + 1e-9,
    };
  }

  // ===========================================================================
  // LQD support — bounded, unbounded, ± singular. The four LQD families share
  // an inverse-Faber rational kernel r#(z) (mathematically identical to the
  // bounded-classical-QD φ with w₀ = 0). Each family wraps r# differently:
  //
  //   boundedLQD            φ(z) = w₀ · exp(r#(z))
  //   boundedLQD_singular   φ(z) = γ · b_{z₀}(z) · exp(r#(z))           z₀ ∈ 𝔻
  //   unboundedLQD          φ(z) = c·z · exp(r#(z) − r#(∞))
  //   unboundedLQD_singular φ(z) = c·|z₀|·z·b_{z₀}(z)·exp(r#(z) − r#(∞)) z₀ ∈ 𝔻*
  //
  // The Schwarz extension F is derived from the identity
  //   conj(u_j(z))|_{|z|=1} = 1/(z − z_j),
  // which turns r#(z) into R##(z) := Σⱼ Σₖ A_{j,k}/(z − zⱼ)ᵏ on the boundary.
  // Combined with conj(z) = 1/z and the Blaschke #-identity b·b# = 1, we get:
  //
  //   F (bounded)             = conj(w₀) · exp(R##(z))
  //   F (bounded singular)    = conj(γ) · b#_{z₀}(z) · exp(R##(z))
  //   F (unbounded)           = (c/z) · exp(R##(z) − conj(r#(∞)))
  //   F (unbounded singular)  = (c·|z₀|/z) · b#_{z₀}(z) · exp(R##(z) − conj(r#(∞)))
  //
  // Derivatives use the log-derivative trick:
  //   (uv)' / (uv) = u'/u + v'/v.
  // φ'/φ is therefore a sum of cheap pieces:
  //    bounded LQD                : r#'(z)
  //    bounded singular LQD       : b'(z)/b(z) + r#'(z)
  //    unbounded LQD              : 1/z + r#'(z)
  //    unbounded singular LQD     : 1/z + b'(z)/b(z) + r#'(z)
  // and r#'(z) is the already-existing branchPhiDeriv. b'/b = 1/(z−z₀) +
  // conj(z₀)/(1 − conj(z₀)z) is a 4-line closed-form.
  // ===========================================================================
  function cexp(a) {
    const ea = Math.exp(a.re);
    return { re: ea * Math.cos(a.im), im: ea * Math.sin(a.im) };
  }
  function blaschkeEval(z, z0) {
    const absZ0 = Math.hypot(z0.re, z0.im);
    if (absZ0 < 1e-14) throw new Error("blaschkeEval: z₀ = 0 not supported");
    const z0C = C.conj(z0);
    const phase = { re: -z0C.re / absZ0, im: -z0C.im / absZ0 };
    const num = C.sub(z, z0);
    const den = C.sub({ re: 1, im: 0 }, C.mul(z0C, z));
    if (Math.hypot(den.re, den.im) < 1e-300) return { re: Infinity, im: 0 };
    return C.mul(phase, C.div(num, den));
  }
  function blaschkeSchwarz(z, z0) {
    // b#_{z₀}(z) = -(z₀/|z₀|) · (1 - conj(z₀)z) / (z - z₀)
    const absZ0 = Math.hypot(z0.re, z0.im);
    if (absZ0 < 1e-14) throw new Error("blaschkeSchwarz: z₀ = 0 not supported");
    const phase = { re: -z0.re / absZ0, im: -z0.im / absZ0 };
    const num = C.sub({ re: 1, im: 0 }, C.mul(C.conj(z0), z));
    const den = C.sub(z, z0);
    if (Math.hypot(den.re, den.im) < 1e-300) return { re: Infinity, im: 0 };
    return C.mul(phase, C.div(num, den));
  }
  function blaschkeLogDeriv(z, z0) {
    // d/dz log b_{z₀}(z) = 1/(z−z₀) + conj(z₀)/(1 − conj(z₀)z).
    const z0C = C.conj(z0);
    const t1 = C.inv(C.sub(z, z0));
    const t2 = C.div(z0C, C.sub({ re: 1, im: 0 }, C.mul(z0C, z)));
    return C.add(t1, t2);
  }
  function rHashAtInfinity(phi) {
    // r#(∞) = Σⱼ Σₖ conj(A_{j,k}) · (-1)^k / conj(z_j)^k.
    // Mirrors the closed-form in QD.LqdCommon.rHashAtInfinity; reproduced
    // here so schwarz-common.js stays self-contained.
    let acc = { re: 0, im: 0 };
    for (const br of phi.branches || []) {
      if (!br.A.length) continue;
      const zjCinv = C.inv(C.conj(br.z));
      let zjCinvPow = { re: 1, im: 0 };
      for (let k = 1; k <= br.A.length; k++) {
        zjCinvPow = C.mul(zjCinvPow, zjCinv);
        const sign = (k % 2 === 0) ? 1 : -1;
        acc = C.add(acc, C.scale(C.mul(C.conj(br.A[k - 1]), zjCinvPow), sign));
      }
    }
    return acc;
  }

  // Common bounded-LQD seed: φ(0) ≈ w₀ (since r#(0) = 0 → exp = 1),
  // φ'(0) ≈ w₀ · Σⱼ conj(A_{j,1}).
  function seedBoundedLQD(phi, w0, w, last) {
    if (last && Math.hypot(last.re, last.im) < 0.9999) return last;
    let dr0 = { re: 0, im: 0 };
    for (const br of phi.branches || []) {
      if (br.A.length > 0) dr0 = C.add(dr0, C.conj(br.A[0]));
    }
    const dphi0 = C.mul(w0, dr0);
    if (Math.hypot(dphi0.re, dphi0.im) < 1e-12) return { re: 0, im: 0 };
    const cand = C.div(C.sub(w, w0), dphi0);
    const m = Math.hypot(cand.re, cand.im);
    if (m < 0.95) return cand;
    return { re: cand.re * 0.9 / m, im: cand.im * 0.9 / m };
  }

  function adaptBoundedLQD(phi) {
    const w0 = phi.w0 || { re: 0, im: 0 };
    function evalPhi(z) {
      return C.mul(w0, cexp(branchPhiContribution(phi, z)));
    }
    function derivPhi(z) {
      // φ' = φ · r#'(z).
      const phiV = evalPhi(z);
      return C.mul(phiV, branchPhiDeriv(phi, z));
    }
    function evalF(z) {
      // F(z) = conj(w₀) · exp(R##(z)).
      return C.mul(C.conj(w0), cexp(branchSchwarzContribution(phi, z)));
    }
    return {
      family: 'boundedLQD',
      evalPhi, derivPhi, evalF,
      seedFor: (w, last) => seedBoundedLQD(phi, w0, w, last),
      acceptZ: (z) => Math.hypot(z.re, z.im) < 1 - 1e-9,
    };
  }

  function adaptBoundedLQD_singular(phi) {
    const gamma = phi.gamma;
    const z0    = phi.z0;
    function evalPhi(z) {
      const b = blaschkeEval(z, z0);
      return C.mul(C.mul(gamma, b), cexp(branchPhiContribution(phi, z)));
    }
    function derivPhi(z) {
      // φ' = φ · (b'/b + r#').
      const phiV = evalPhi(z);
      return C.mul(phiV, C.add(blaschkeLogDeriv(z, z0), branchPhiDeriv(phi, z)));
    }
    function evalF(z) {
      // F(z) = conj(γ) · b#_{z₀}(z) · exp(R##(z)).
      const bs = blaschkeSchwarz(z, z0);
      return C.mul(C.mul(C.conj(gamma), bs), cexp(branchSchwarzContribution(phi, z)));
    }
    // Seed: φ(0) = γ · b(0) · 1 = γ · (z₀/|z₀|) · (-z₀)/(1) · (-1)
    //              = γ · (z₀/|z₀|) · z₀ = γ · |z₀|·(z₀/|z₀|²) ... too fragile.
    // Use a robust generic seed inside 𝔻 and let the 4-seed retry ladder
    // handle landing in the right basin.
    function seedFor(w, last) {
      if (last && Math.hypot(last.re, last.im) < 0.9999) return last;
      return { re: 0.3, im: 0.3 };
    }
    return {
      family: 'boundedLQD_singular',
      evalPhi, derivPhi, evalF,
      seedFor,
      acceptZ: (z) => Math.hypot(z.re, z.im) < 1 - 1e-9,
    };
  }

  function adaptUnboundedLQD(phi) {
    // φ(z) = c·z·exp(r̃#(z) + B(1/z)), where B(1/z) = Σ_l β_l/z^l comes
    // from polynomial-h support (HANDOFF #22). β is `phi.lqdBeta`; the
    // evalBOverZ / evalBOverZDeriv / evalBConjOfZ helpers above mirror
    // solver-uqd-lqd.js's evaluator.
    const c    = phi.c;
    const rInf = rHashAtInfinity(phi);              // r#(∞)
    const rInfConj = C.conj(rInf);                  // for F
    function evalPhi(z) {
      const r = branchPhiContribution(phi, z);
      const b = evalBOverZ(phi, z);
      const eff = C.add(C.sub(r, rInf), b);
      return C.mul(C.scale(z, c), cexp(eff));
    }
    function derivPhi(z) {
      // φ' = φ · (1/z + r#'(z) + B(1/z)').
      const phiV = evalPhi(z);
      const tot = C.add(C.add(C.inv(z), branchPhiDeriv(phi, z)),
                        evalBOverZDeriv(phi, z));
      return C.mul(phiV, tot);
    }
    function evalF(z) {
      // F(z) = (c/z) · exp(R##(z) − conj(r#(∞)) + conj(B(z))) where
      // conj(B(z)) = Σ_l conj(β_l)·z^l is the Schwarz reflection of B(1/z).
      const R = branchSchwarzContribution(phi, z);
      const bConj = evalBConjOfZ(phi, z);
      const eff = C.add(C.sub(R, rInfConj), bConj);
      const cOverZ = C.scale(C.inv(z), c);
      return C.mul(cOverZ, cexp(eff));
    }
    function seedFor(w, last) {
      if (last && Math.hypot(last.re, last.im) > QD.DISK_CLAMP_OUT) return last;
      const cand = { re: w.re / c, im: w.im / c };
      const m = Math.hypot(cand.re, cand.im);
      if (m > 1.05) return cand;
      if (m < 1e-12) return { re: 1.1, im: 0 };
      return { re: cand.re * 1.1 / m, im: cand.im * 1.1 / m };
    }
    return {
      family: 'unboundedLQD',
      c, rInf,
      evalPhi, derivPhi, evalF,
      seedFor,
      acceptZ: (z) => Math.hypot(z.re, z.im) > 1 + 1e-9,
    };
  }

  function adaptUnboundedLQD_singular(phi) {
    // φ(z) = c·|z₀|·z·b_{z₀}(z)·exp(r̃#(z) + B(1/z)).
    // β extends r̃# at ∞ (HANDOFF #22). γ (HANDOFF #24) extends r̃# with a
    // synthetic Möbius branch anchored at z=z₀. We merge γ into the branches
    // up-front via `withSyntheticBranch` so the existing
    // branchPhiContribution / branchPhiDeriv / branchSchwarzContribution
    // pickup it up transparently. β remains a separate exp-argument term.
    const phiX  = withSyntheticBranch(phi);
    const c     = phi.c;
    const z0    = phi.z0;
    const absZ0 = Math.hypot(z0.re, z0.im);
    const rInf  = rHashAtInfinity(phiX);
    const rInfConj = C.conj(rInf);
    function evalPhi(z) {
      const r = branchPhiContribution(phiX, z);
      const b = evalBOverZ(phi, z);
      const eff = C.add(C.sub(r, rInf), b);
      const bl = blaschkeEval(z, z0);
      const scale = C.scale(z, c * absZ0);
      return C.mul(C.mul(scale, bl), cexp(eff));
    }
    function derivPhi(z) {
      // φ' = φ · (1/z + b'/b + r#'(z) + B(1/z)').
      const phiV = evalPhi(z);
      const tot = C.add(
        C.add(C.inv(z), blaschkeLogDeriv(z, z0)),
        C.add(branchPhiDeriv(phiX, z), evalBOverZDeriv(phi, z))
      );
      return C.mul(phiV, tot);
    }
    function evalF(z) {
      // F(z) = (c·|z₀|/z) · b#_{z₀}(z) · exp(R##(z) − conj(r#(∞)) + conj(B(z))).
      const R = branchSchwarzContribution(phiX, z);
      const bConj = evalBConjOfZ(phi, z);
      const eff = C.add(C.sub(R, rInfConj), bConj);
      const bs = blaschkeSchwarz(z, z0);
      const cOverZ = C.scale(C.inv(z), c * absZ0);
      return C.mul(C.mul(cOverZ, bs), cexp(eff));
    }
    function seedFor(w, last) {
      if (last && Math.hypot(last.re, last.im) > QD.DISK_CLAMP_OUT) return last;
      const denom = c * absZ0;
      if (denom < 1e-12) return { re: 1.1, im: 0 };
      const cand = { re: w.re / denom, im: w.im / denom };
      const m = Math.hypot(cand.re, cand.im);
      if (m > 1.05) return cand;
      if (m < 1e-12) return { re: 1.1, im: 0 };
      return { re: cand.re * 1.1 / m, im: cand.im * 1.1 / m };
    }
    return {
      family: 'unboundedLQD_singular',
      c, z0,
      evalPhi, derivPhi, evalF,
      seedFor,
      acceptZ: (z) => Math.hypot(z.re, z.im) > 1 + 1e-9,
    };
  }

  function adaptBoundedRational(phi) {
    // phi = { rational: true, P: Complex[], Q: Complex[], w0?, ... }.
    // φ(z) = P(z) / Q(z), assumed analytic on 𝔻̄ (Q has no roots there).
    const P = phi.P, Q = phi.Q;
    const p = P.length - 1, q = Q.length - 1;
    function evalPhiRat(z) {
      return C.div(polyEval(P, z), polyEval(Q, z));
    }
    function derivPhiRat(z) {
      // (P/Q)' = (P'Q − PQ')/Q²
      const Pv = polyEval(P, z), Qv = polyEval(Q, z);
      const Pp = polyDerivEval(P, z), Qp = polyDerivEval(Q, z);
      return C.div(C.sub(C.mul(Pp, Qv), C.mul(Pv, Qp)), C.mul(Qv, Qv));
    }
    // Schwarz extension F: φ#(z) = z^{q−p} · P̃(z) / Q̃(z), where X̃ is the
    // reverse-conjugate of X. This is the same construction the direct
    // rational kernel uses; it gives σ on Ω.
    const Ptil = reverseConjugate(P);
    const Qtil = reverseConjugate(Q);
    function evalF(z) {
      let zPow = { re: 1, im: 0 };
      const e = q - p;
      if (e > 0) for (let k = 0; k < e; k++) zPow = C.mul(zPow, z);
      else if (e < 0) {
        const zInv = C.inv(z);
        for (let k = 0; k < -e; k++) zPow = C.mul(zPow, zInv);
      }
      return C.mul(zPow, C.div(polyEval(Ptil, z), polyEval(Qtil, z)));
    }
    // Seed: linearization around z=0 (φ(0) = P[0]/Q[0]; φ'(0) = (P[1]Q[0] − P[0]Q[1])/Q[0]²).
    const coeff0 = C.div(P[0] || { re: 0, im: 0 }, Q[0] || { re: 1, im: 0 });
    const coeff1Num = C.sub(
      C.mul(P[1] || { re: 0, im: 0 }, Q[0]),
      C.mul(P[0] || { re: 0, im: 0 }, Q[1] || { re: 0, im: 0 }),
    );
    const coeff1 = C.div(coeff1Num, C.mul(Q[0], Q[0]));
    function seedFor(w, last) {
      if (last && Math.hypot(last.re, last.im) < 0.9999) return last;
      const m1 = Math.hypot(coeff1.re, coeff1.im);
      if (m1 < 1e-12) return { re: 0, im: 0 };
      const cand = C.div(C.sub(w, coeff0), coeff1);
      const m = Math.hypot(cand.re, cand.im);
      if (m < 0.95) return cand;
      return { re: cand.re * 0.9 / m, im: cand.im * 0.9 / m };
    }
    return {
      family: 'boundedQDRational',
      P, Q,
      evalPhi: evalPhiRat,
      derivPhi: derivPhiRat,
      evalF,
      seedFor,
      acceptZ: (z) => Math.hypot(z.re, z.im) < 1 - 1e-9,
    };
  }

  function reverseConjugate(P) {
    // X̃(z) = Σ_k conj(X_{deg X − k}) z^k
    const n = P.length;
    const out = new Array(n);
    for (let k = 0; k < n; k++) out[k] = C.conj(P[n - 1 - k]);
    return out;
  }

  // ===========================================================================
  // Public builder: phi (+ optional hData) → { sigma, isInOmega, w0, ... } or null.
  // ===========================================================================
  function buildSchwarzFromPhi(phi, hData, boundaryPts) {
    if (!phi) return null;
    // Family-driven dispatch. Classical bounded/unbounded QDs don't set
    // phi.family, so we detect them by the `unbounded` flag (see HANDOFF
    // gotcha #1 — classical phis intentionally leave the family tag unset).
    let adapter;
    switch (phi.family) {
      case 'boundedLQD':            adapter = adaptBoundedLQD(phi); break;
      case 'boundedLQD_singular':   adapter = adaptBoundedLQD_singular(phi); break;
      case 'unboundedLQD':          adapter = adaptUnboundedLQD(phi); break;
      case 'unboundedLQD_singular': adapter = adaptUnboundedLQD_singular(phi); break;
      default:
        adapter = phi.unbounded ? adaptUnbounded(phi) : adaptBounded(phi, hData);
    }
    return buildFromAdapter(adapter, phi, boundaryPts);
  }
  function buildSchwarzFromRational(phi, boundaryPts) {
    if (!phi || !phi.P || !phi.Q) return null;
    const adapter = adaptBoundedRational(phi);
    return buildFromAdapter(adapter, phi, boundaryPts);
  }

  function buildFromAdapter(adapter, phi, boundaryPts) {
    // Cached "last successful seed" — naive single-cell cache (renderer will
    // own per-pixel seeds when iterating across the grid).
    let lastSeed = null;

    function psi(w, seedHint) {
      const seed = seedHint || lastSeed || adapter.seedFor(w, null);
      const z = newtonInvert(adapter.evalPhi, adapter.derivPhi, w, seed, { maxIter: 40, tol: 1e-12 });
      if (z && adapter.acceptZ(z)) { lastSeed = z; return z; }
      // Fallback: try one more time with the "fresh" linearization seed.
      const seed2 = adapter.seedFor(w, null);
      const z2 = newtonInvert(adapter.evalPhi, adapter.derivPhi, w, seed2, { maxIter: 60, tol: 1e-11 });
      if (z2 && adapter.acceptZ(z2)) { lastSeed = z2; return z2; }
      return null;
    }
    function sigma(w, seedHint) {
      const z = psi(w, seedHint);
      if (!z) return null;
      // Pole-near-origin sentinel. F doesn't generally have a pole at z=0
      // (classical bounded F has poles at the z_j ∈ 𝔻 instead; LQD-bounded F
      // has essential singularities at z_j; for unbounded families ψ never
      // lands near 0). The check survives mostly as a safety: if ψ converges
      // to a tiny z by accident, evalF likely produces a non-finite value
      // that the next guard catches anyway.
      if (Math.hypot(z.re, z.im) < 1e-14) return null;
      let Sv;
      try { Sv = adapter.evalF(z); }
      catch (e) { return null; }
      if (!Sv || !isFinite(Sv.re) || !isFinite(Sv.im)) return null;
      return { re: Sv.re, im: -Sv.im, _z: z };       // conj
    }

    // Polygon-based in-Ω test (works for bounded and unbounded). For unbounded
    // Ω the user passes the inner boundary of the bounded complement K; "in Ω"
    // means OUTSIDE the polygon. We expose both predicates and let the caller
    // pick the right one.
    //
    // Build the y-binned polygon index once per Schwarz handle; each
    // per-pixel test then walks only the bin's edges (~6–10× speedup on
    // typical boundaries with the CPU pyramid renderer).
    const polyIdx = (boundaryPts && boundaryPts.length >= 4)
      ? buildPolygonIndex(boundaryPts) : null;
    function inOmegaBounded(w) {
      return polyIdx ? pointInPolygonIndexed(w, polyIdx)
                     : pointInPolygon(w, boundaryPts);
    }
    function inOmegaUnbounded(w) {
      return !(polyIdx ? pointInPolygonIndexed(w, polyIdx)
                       : pointInPolygon(w, boundaryPts));
    }
    const isInOmega = phi.unbounded ? inOmegaUnbounded : inOmegaBounded;

    // Escape-radius default for unbounded: a few times the K-diameter.
    let escapeR = Infinity;
    if (phi.unbounded && boundaryPts && boundaryPts.length) {
      const b = polygonBounds(boundaryPts);
      escapeR = b.radius * 30 + 10;
    }

    return {
      family:    adapter.family,
      sigma,
      psi,
      evalPhi:   adapter.evalPhi,
      evalF:     adapter.evalF,
      isInOmega,
      escapeR,
      unbounded: !!phi.unbounded,
      // For UI use:
      w0:        phi.w0 || (phi.unbounded ? { re: Infinity, im: Infinity } : null),
      adapter,
    };
  }

  // ===========================================================================
  // Escape-time computation for a single starting point.
  //   Returns { kind, n, lastW }:
  //     kind = 'fundamental' (σⁿ left Ω into Ω^c)
  //          | 'escaped'     (|σⁿ| > escapeR; unbounded only)
  //          | 'interior'    (still in Ω after maxIter)
  //          | 'invalid'     (Newton failed; ψ undefined)
  // ===========================================================================
  function escapeTime(w0, schwarz, opts) {
    opts = opts || {};
    const maxIter  = opts.maxIter  || 64;
    const escapeR  = opts.escapeR  || schwarz.escapeR || Infinity;
    // initialSeedHint: spatial warm-start for the first ψ call. The caller
    // (a raster renderer) passes the previous pixel's converged z so the
    // Newton inverse starts from a nearby point in 𝔻 — typically converges
    // in 1–3 iterations instead of 5–10 from a cold seed. Also returned in
    // the result as `firstZ` so the caller can chain across pixels.
    let seedHint = opts.initialSeedHint || null;

    let w = { re: w0.re, im: w0.im };
    if (!schwarz.isInOmega(w)) return { kind: 'fundamental', n: 0, lastW: w, firstZ: null };

    let firstZ = null;
    for (let n = 1; n <= maxIter; n++) {
      const next = schwarz.sigma(w, seedHint);
      if (!next) return { kind: 'invalid', n: n - 1, lastW: w, firstZ };
      if (firstZ === null) firstZ = next._z;
      seedHint = next._z;
      w = { re: next.re, im: next.im };
      if (!isFinite(w.re) || !isFinite(w.im)) return { kind: 'escaped', n, lastW: w, firstZ };
      if (Math.hypot(w.re, w.im) > escapeR)   return { kind: 'escaped', n, lastW: w, firstZ };
      if (!schwarz.isInOmega(w))              return { kind: 'fundamental', n, lastW: w, firstZ };
    }
    return { kind: 'interior', n: maxIter, lastW: w, firstZ };
  }

  // Make an explicit orbit polyline {w0, σ(w0), σ²(w0), ...} stopping at the
  // same conditions as escapeTime (no further iteration once we leave Ω).
  function makeOrbit(w0, schwarz, opts) {
    opts = opts || {};
    const maxIter = opts.maxIter || 64;
    const escapeR = opts.escapeR || schwarz.escapeR || Infinity;
    const pts = [{ re: w0.re, im: w0.im }];
    let w = pts[0];
    let seedHint = null;
    for (let n = 0; n < maxIter; n++) {
      if (!schwarz.isInOmega(w)) break;
      const next = schwarz.sigma(w, seedHint);
      if (!next) break;
      seedHint = next._z;
      w = { re: next.re, im: next.im };
      pts.push(w);
      if (!isFinite(w.re) || !isFinite(w.im)) break;
      if (Math.hypot(w.re, w.im) > escapeR)   break;
    }
    return pts;
  }

  // ===========================================================================
  // Export.
  // ===========================================================================
  Schwarz.buildSchwarzFromPhi      = buildSchwarzFromPhi;
  Schwarz.buildSchwarzFromRational = buildSchwarzFromRational;
  Schwarz.pointInPolygon           = pointInPolygon;
  Schwarz.buildPolygonIndex        = buildPolygonIndex;
  Schwarz.pointInPolygonIndexed    = pointInPolygonIndexed;
  Schwarz.escapeTime               = escapeTime;
  Schwarz.makeOrbit                = makeOrbit;
  Schwarz.polygonBounds            = polygonBounds;
  // For tests:
  Schwarz._adaptBounded            = adaptBounded;
  Schwarz._adaptUnbounded          = adaptUnbounded;
  Schwarz._adaptBoundedRational    = adaptBoundedRational;
  Schwarz._reverseConjugate        = reverseConjugate;

})(typeof globalThis !== 'undefined' ? globalThis : this);
