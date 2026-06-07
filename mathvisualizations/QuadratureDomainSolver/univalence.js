// =============================================================================
// univalence.js  —  Special univalence criteria for the solved Riemann map φ.
//
// Beyond plain univalence (boundary self-intersection, isBoundaryUnivalent),
// classical geometric-function theory classifies the shape of Ω = φ(𝔻) (or
// φ(𝔻*) for unbounded families):
//
//   • star-like   — every ray from the center c stays in Ω.
//   • spiral-like — a logarithmic-spiral generalization of star-like.
//   • convex      — Ω is convex (bounded families only).
//
// All three reduce to sign/argument tests of quantities built from φ, φ′, φ″
// on the boundary z = e^{iθ}.  By the minimum principle for the harmonic Re of
// the relevant analytic function, the boundary test is equivalent to the
// interior condition WHEN φ is univalent — so we evaluate on ∂𝔻 only.
//
// Star-function g(z):
//   bounded   (center c = φ(0) = phi.w0):   g = z·φ′/(φ − c)
//   unbounded (star-like w.r.t. ∞):         g = z·φ′/φ
//   star-like   ⇔ min_θ Re(g) > 0
//   spiral-like ⇔ {arg g(e^{iθ})} fit in an open arc of width < π
//                 (optimal λ = −center of the covering arc)
//
// Convex (bounded only):  q = 1 + z·φ″/φ′ ;  convex ⇔ min_θ Re(q) > 0.
//   Unbounded Ω (exterior of a bounded complement) is essentially never
//   convex, so convexity is reported N/A there.
//
// Hierarchy:  convex ⟹ star-like ⟹ spiral-like, and star-like/convex ⟹
//   univalent.  A criterion is only asserted "yes" when φ is univalent.
//
// API:
//   QD.classifyUnivalence(phi, opts) → {
//     bounded, univalent, center: {re,im} | 'infinity',
//     starLike:   { is, margin },                  // margin = min Re(g)
//     spiralLike: { is, arcWidth, angleDeg },      // angleDeg = optimal λ (deg)
//     convex:     { is, margin } | { na: true },   // na for unbounded
//     notes: [ ... ],
//   }
//   opts: { samples = 360, univalent }   (univalent only affects `notes`)
//
// Reuses QD.phiTaylorAt(z, phi, 2) → [a0, a1, a2] with φ = a0, φ′ = a1,
// φ″ = 2·a2 (the standard Taylor convention; see critical-set.js).  No new
// numeric primitives.  Pure + DOM-free; loaded page-side only (like
// critical-set.js), not bundled into the solver Workers.
// =============================================================================

(function (global) {
  'use strict';

  // QD-namespace resolution — same idiom every solver file uses.
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD) return;

  const C = QD.Complex;

  const DEFAULT_SAMPLES = 360;
  const MARGIN_TOL      = 1e-6;   // strict-inequality tolerance band
  const PHIP_FLOOR2     = 1e-24;  // |φ′|² below this ⇒ boundary critical point

  // Normalize an angle to (−π, π].
  function _wrapPi(a) {
    while (a >  Math.PI) a -= 2 * Math.PI;
    while (a <= -Math.PI) a += 2 * Math.PI;
    return a;
  }

  // Largest angular gap among a set of angles on the circle. Returns
  // { arcWidth, centerAngle } where arcWidth = 2π − largestGap is the width of
  // the minimal arc covering all the angles, and centerAngle is that arc's
  // midpoint. If all angles coincide, arcWidth ≈ 0.
  function _coveringArc(args) {
    const n = args.length;
    if (n === 0) return { arcWidth: 0, centerAngle: 0 };
    const sorted = args.slice().sort((a, b) => a - b);
    let maxGap = -Infinity, gapEnd = sorted[0];
    for (let i = 0; i < n; i++) {
      const lo = sorted[i];
      const hi = (i + 1 < n) ? sorted[i + 1] : sorted[0] + 2 * Math.PI;
      const gap = hi - lo;
      if (gap > maxGap) { maxGap = gap; gapEnd = hi; }
    }
    const arcWidth = 2 * Math.PI - maxGap;
    // Covered arc runs from gapEnd around to gapStart+2π; its midpoint:
    const centerAngle = _wrapPi(gapEnd + arcWidth / 2);
    return { arcWidth, centerAngle };
  }

  // ---------------------------------------------------------------------------
  // classifyUnivalence — public entry point.
  // ---------------------------------------------------------------------------
  function classifyUnivalence(phi, opts) {
    opts = opts || {};
    const N = Math.max(16, opts.samples | 0 || DEFAULT_SAMPLES);
    const bounded = !phi.unbounded;
    const notes = [];

    // Center c = φ(0) for bounded families (where the star-function subtracts
    // it). Unbounded families are star-like w.r.t. ∞ ⇒ no subtraction.
    let center = 'infinity';
    if (bounded) {
      center = (phi.w0 && Number.isFinite(phi.w0.re)) ? phi.w0 : QD.evalPhi({ re: 0, im: 0 }, phi);
    }

    let starMargin  = Infinity;     // min Re(g)
    let convexMargin = Infinity;    // min Re(q), bounded only
    let convexIndeterminate = false;
    const argsG = [];

    for (let i = 0; i < N; i++) {
      const theta = (2 * Math.PI * i) / N;
      const z = { re: Math.cos(theta), im: Math.sin(theta) };
      const t = QD.phiTaylorAt(z, phi, 2);     // [a0, a1, a2]
      const w  = t[0];                          // φ
      const dp = t[1];                          // φ′
      const ddp = { re: 2 * t[2].re, im: 2 * t[2].im };   // φ″ = 2·a2

      // Star-function g = z·φ′ / (φ − c)   (bounded) or z·φ′/φ (unbounded).
      const num = C.mul(z, dp);
      const den = bounded ? C.sub(w, center) : w;
      const g = C.div(num, den);
      if (g.re < starMargin) starMargin = g.re;
      argsG.push(Math.atan2(g.im, g.re));

      // Convexity q = 1 + z·φ″/φ′  (bounded only).
      if (bounded) {
        if (C.abs2(dp) < PHIP_FLOOR2) {
          convexIndeterminate = true;          // φ′ ≈ 0 on ∂𝔻 (cusp / fold)
        } else {
          const q = C.add({ re: 1, im: 0 }, C.div(C.mul(z, ddp), dp));
          if (q.re < convexMargin) convexMargin = q.re;
        }
      }
    }

    const arc = _coveringArc(argsG);

    const starLike = { is: starMargin > MARGIN_TOL, margin: starMargin };
    const spiralLike = {
      is: arc.arcWidth < Math.PI - MARGIN_TOL,
      arcWidth: arc.arcWidth,
      angleDeg: _wrapPi(-arc.centerAngle) * 180 / Math.PI,   // optimal λ
    };

    let convex;
    if (!bounded) {
      convex = { na: true };
      notes.push('Convexity is reported only for bounded Ω.');
    } else if (convexIndeterminate) {
      convex = { is: false, margin: 0, indeterminate: true };
      notes.push('φ′ vanishes on ∂𝔻 (cusp / fold) — convexity indeterminate.');
    } else {
      convex = { is: convexMargin > MARGIN_TOL, margin: convexMargin };
    }

    const univalent = (opts.univalent === undefined) ? null : !!opts.univalent;
    if (univalent === false) {
      notes.push('Map is not univalent — these sufficient conditions presuppose a univalent Ω.');
    }

    return {
      bounded,
      univalent,
      center,
      starLike,
      spiralLike,
      convex,
      samples: N,
      notes,
    };
  }

  QD.Univalence = { classifyUnivalence };
  QD.classifyUnivalence = classifyUnivalence;

})(typeof globalThis !== 'undefined' ? globalThis : this);
