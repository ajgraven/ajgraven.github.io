// =============================================================================
// critical-set.js  —  Critical points of φ and their w-plane images.
//
// A Riemann map φ: 𝔻 → Ω (or 𝔻* → Ω for unbounded families) is univalent
// iff φ'(z) ≠ 0 throughout the relevant disk.  Zeros of φ' that lie INSIDE
// 𝔻 (or 𝔻*) signal univalence loss — the map folds onto itself.  Zeros
// near |z| = 1 predict imminent degeneracy as parameters vary.
//
// This module finds those zeros numerically and reports each one as
//   {z, w = φ(z), |z|, in-the-relevant-disk?, severity}.
//
// API:
//   QD.findCriticalPoints(phi, opts) → {
//     points: [
//       { z: {re,im},
//         w: {re,im},
//         absZ: Number,
//         inDomain: Boolean,     // |z|<1 for bounded, |z|>1 for unbounded
//         severity: 'critical' | 'near' | 'safe' },
//       ...
//     ],
//     stats: { nSeeds, nConverged, nUnique }
//   }
//
// Algorithm:
//   1. Seed a polar grid covering |z| ∈ [innerR, outerR] (where the bounds
//      depend on bounded vs unbounded; the relevant disk is always one of
//      𝔻 or 𝔻*, the other side is sampled too so we can show "safe" roots).
//   2. For each seed, run complex Newton on f(z) = φ'(z) using
//      phiTaylorAt(z, phi, 2) — a₁ = f(z), a₂ = f'(z)/2.  Newton update
//      z ← z − a₁ / (2·a₂).
//   3. Stop on convergence (|a₁| < TOL), divergence (|z| > MAX_R), or
//      iteration limit.
//   4. Deduplicate accepted roots by snapping to a 1e-5 grid.
//   5. For each unique root, evaluate w = evalPhi(z, phi).  Compute
//      severity:
//        bounded family:
//          inDomain = |z| < 1
//          severity = 'critical' if |z| < 1 − BAND
//                     'near'     if  |1 − |z|| ≤ BAND
//                     'safe'     otherwise
//        unbounded family: symmetric (relevant disk is |z| > 1).
//
// Note on dedup tolerance: 1e-5 in z catches genuinely-equal roots while
// preserving distinct close roots (typical critical points of a high-order
// φ are separated by at least 1e-3).
//
// This module is pure (no DOM, no WebGL) so it loads cleanly in node-test.js.
// =============================================================================

(function (global) {
  'use strict';

  // QD-namespace resolution — same idiom every solver file uses.
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD) return;

  const TOL_CONVERGE = 1e-10;   // |φ'(z)| accepted as zero
  const TOL_DERIV    = 1e-30;   // avoid division by ≈0 (Newton step explosion)
  const MAX_NEWTON   = 40;
  const MAX_Z        = 100;     // declare "escaped" beyond this
  const DEDUP_TOL    = 1e-5;    // snap-grid spacing for dedup
  const SEVERITY_BAND = 0.05;   // |1 − |z|| ≤ this counts as "near"

  // ---------------------------------------------------------------------------
  // Default seed grid: covers both the bounded (|z|<1) and unbounded (|z|>1)
  // sides so we can also report safely-outside critical points.  For bounded
  // families the relevant disk is 𝔻; for unbounded it's 𝔻*.  Seeds outside
  // |z| ≤ 3 are unnecessary in practice — any further roots are pushed back
  // by Newton anyway, and almost no QD has critical points far from the unit
  // circle.
  // ---------------------------------------------------------------------------
  function _defaultSeeds(unbounded) {
    const seeds = [];
    // Radii: dense around |z|=1, coarser further out.
    const radii = [0.10, 0.25, 0.45, 0.65, 0.82, 0.92, 0.98,
                   1.02, 1.08, 1.20, 1.45, 1.80, 2.40];
    // Angles: every 30°.
    const nAng = 12;
    for (const r of radii) {
      for (let k = 0; k < nAng; k++) {
        const th = 2 * Math.PI * k / nAng;
        seeds.push({ re: r * Math.cos(th), im: r * Math.sin(th) });
      }
    }
    // Add the origin and ∞-side anchor for symmetry.
    seeds.push({ re: 0, im: 0 });
    // (We don't add an explicit ∞ seed — Newton would diverge; the |z|=2.4
    //  ring is already in the unbounded-family regime.)
    return seeds;
  }

  // ---------------------------------------------------------------------------
  // Complex Newton on φ'(z) = 0.  Returns {z, ok} where ok=true means the
  // iteration converged to within TOL_CONVERGE.  Caller deduplicates and
  // classifies.
  // ---------------------------------------------------------------------------
  function _newtonStep(z, phi) {
    // phiTaylorAt(z0, phi, L) returns [a_0, a_1, ..., a_L] where
    //   φ(z0 + t) = a_0 + a_1·t + a_2·t² + ...
    // So a_1 = φ'(z0), a_2 = φ''(z0)/2.  Newton on f = φ' uses
    //   z ← z − f/f' = z − a_1 / (2·a_2).
    let aT;
    try {
      aT = QD.phiTaylorAt(z, phi, 2);
    } catch (e) {
      return { ok: false, reason: 'taylor-throw' };
    }
    if (!aT || aT.length < 3) return { ok: false, reason: 'taylor-short' };
    const f  = aT[1];                       // φ'(z)
    const f2 = aT[2];                       // φ''(z)/2  (== d/dz φ' / 2)
    const fpAbs2 = f.re * f.re + f.im * f.im;
    if (fpAbs2 < TOL_CONVERGE * TOL_CONVERGE) {
      return { ok: true, z, fpAbs2 };
    }
    const dpAbs2 = f2.re * f2.re + f2.im * f2.im;
    if (4 * dpAbs2 < TOL_DERIV) {
      return { ok: false, reason: 'flat-derivative', f, f2 };
    }
    // step = f / (2·f2)
    const denom = 2 * dpAbs2;
    // Complex division (a + bi)/(c + di) = (ac + bd + i(bc − ad)) / (c²+d²)
    // Here numerator is f, denominator is 2·f2 = (2 f2.re) + i (2 f2.im).
    const cr = 2 * f2.re, ci = 2 * f2.im;
    const dr2 = cr * cr + ci * ci;
    const stepRe = (f.re * cr + f.im * ci) / dr2;
    const stepIm = (f.im * cr - f.re * ci) / dr2;
    return {
      ok: false,
      next: { re: z.re - stepRe, im: z.im - stepIm },
      fpAbs2,
    };
  }

  function _newton(seed, phi) {
    let z = { re: seed.re, im: seed.im };
    for (let iter = 0; iter < MAX_NEWTON; iter++) {
      const r = _newtonStep(z, phi);
      if (r.ok) {
        return { ok: true, z: r.z, fpAbs2: r.fpAbs2, iters: iter };
      }
      if (!r.next) {
        return { ok: false, reason: r.reason || 'no-step' };
      }
      z = r.next;
      const az2 = z.re * z.re + z.im * z.im;
      if (!isFinite(az2) || az2 > MAX_Z * MAX_Z) {
        return { ok: false, reason: 'diverged' };
      }
    }
    return { ok: false, reason: 'max-iter' };
  }

  // ---------------------------------------------------------------------------
  // Snap a complex number to a coarse grid for dedup.
  // ---------------------------------------------------------------------------
  function _snapKey(z) {
    const inv = 1 / DEDUP_TOL;
    const kr = Math.round(z.re * inv);
    const ki = Math.round(z.im * inv);
    return kr + ',' + ki;
  }

  // ---------------------------------------------------------------------------
  // Severity classifier.
  //   bounded   family: relevant disk = 𝔻 (|z| < 1)
  //   unbounded family: relevant disk = 𝔻* (|z| > 1)
  // ---------------------------------------------------------------------------
  function _classify(absZ, unbounded) {
    const dist = absZ - 1;   // signed distance from the unit circle (radial)
    if (!unbounded) {
      // Bounded: inDomain ⇔ |z| < 1.
      const inDomain = absZ < 1;
      let severity;
      if (Math.abs(dist) <= SEVERITY_BAND) severity = 'near';
      else if (absZ < 1 - SEVERITY_BAND)   severity = 'critical';
      else                                  severity = 'safe';
      return { inDomain, severity };
    }
    // Unbounded: inDomain ⇔ |z| > 1.
    const inDomain = absZ > 1;
    let severity;
    if (Math.abs(dist) <= SEVERITY_BAND) severity = 'near';
    else if (absZ > 1 + SEVERITY_BAND)   severity = 'critical';
    else                                  severity = 'safe';
    return { inDomain, severity };
  }

  // ---------------------------------------------------------------------------
  // findCriticalPoints  — public entry point
  // ---------------------------------------------------------------------------
  function findCriticalPoints(phi, opts) {
    opts = opts || {};
    if (!phi) return { points: [], stats: { nSeeds: 0, nConverged: 0, nUnique: 0 } };

    const unbounded = !!phi.unbounded || /unbounded/i.test(phi.family || '');
    const seeds     = opts.seeds || _defaultSeeds(unbounded);

    const uniqueMap = new Map();  // snap-key → best (smallest |φ'|) {z, fpAbs2}
    let nConverged = 0;

    for (const seed of seeds) {
      const r = _newton(seed, phi);
      if (!r.ok) continue;
      nConverged++;
      const key = _snapKey(r.z);
      const cur = uniqueMap.get(key);
      if (!cur || r.fpAbs2 < cur.fpAbs2) {
        uniqueMap.set(key, { z: r.z, fpAbs2: r.fpAbs2 });
      }
    }

    const points = [];
    for (const { z } of uniqueMap.values()) {
      const absZ = Math.hypot(z.re, z.im);
      const { inDomain, severity } = _classify(absZ, unbounded);
      let w;
      try {
        w = QD.evalPhi(z, phi);
      } catch (e) {
        // Skip points where φ itself blows up (shouldn't happen at a true
        // critical point of φ', but a numerically-stalled Newton landing
        // exactly on a branch pole could).
        continue;
      }
      if (!w || !isFinite(w.re) || !isFinite(w.im)) continue;
      points.push({
        z: { re: z.re, im: z.im },
        w: { re: w.re, im: w.im },
        absZ,
        inDomain,
        severity,
      });
    }

    // Sort: critical first, then near, then safe; within a group, by |z|.
    const sevRank = { critical: 0, near: 1, safe: 2 };
    points.sort((a, b) => {
      const da = sevRank[a.severity] - sevRank[b.severity];
      if (da !== 0) return da;
      return a.absZ - b.absZ;
    });

    return {
      points,
      stats: { nSeeds: seeds.length, nConverged, nUnique: points.length },
    };
  }

  // Expose under QD.CriticalSet for namespace hygiene, and a top-level
  // alias QD.findCriticalPoints for convenience.
  QD.CriticalSet = {
    findCriticalPoints,
    _defaultSeeds, _newton, _classify, _snapKey,   // exposed for tests
  };
  QD.findCriticalPoints = findCriticalPoints;

})(typeof globalThis !== 'undefined' ? globalThis : this);
