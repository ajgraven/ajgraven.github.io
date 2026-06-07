// =============================================================================
// schwarz-forward.js — Phase S5: forward-dynamics kernels.
//
// Companion to schwarz-inverse.js (S1-S3) and schwarz-analysis.js (S4).
// Pure math; UI integration in schwarz-ui.js.
//
// Exports on QD.Schwarz:
//   canonicalSeeds(schwarz)              → seed points worth tracing per family
//   iterateCurveForward(pts, sw, k)      → σ-images of a polyline, k iterations
//   findCycles(sw, n, opts)              → period-n cycles via Newton
//   sampleSweepSeeds(kind, params)       → evenly-spaced seeds on a line/circle
//   domainColoringField(sw, viewport, …) → per-pixel HSL field of σ(w)
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Schwarz) {
    throw new Error("schwarz-forward.js: schwarz-common.js must be loaded first");
  }

  // ---------------------------------------------------------------------------
  // 1. canonicalSeeds — per-family natural orbit-tracker seed points (H7).
  //
  //   boundedQD            : φ(0) = w₀                  (the "center" of Ω)
  //   unboundedQD          : φ(c⁻¹·w₀) — not meaningful; fall back to centroid
  //                          of branches. We use 0 if Ω contains it, else the
  //                          centroid of φ-images of z_j.
  //   boundedLQD           : φ(0) = w₀
  //   boundedLQD_singular  : φ(z₀)  (Blaschke center)
  //   unboundedLQD         : centroid of branches' a_j (the user-supplied
  //                          quadrature node mean).
  //   unboundedLQD_singular: same + φ(z₀) when defined
  //
  // Each entry: { w, label }. Caller iterates σ from each to get the orbit.
  // ---------------------------------------------------------------------------
  function canonicalSeeds(schwarz) {
    if (!schwarz || !schwarz._phi) return [];
    const phi      = schwarz._phi;
    const family   = schwarz.family;
    const out      = [];

    const addIfInOmega = (w, label) => {
      if (!w || !isFinite(w.re) || !isFinite(w.im)) return;
      try { if (!schwarz.isInOmega(w)) return; } catch (_) { return; }
      out.push({ w, label });
    };

    // φ(0) = w₀ for bounded families. PQDs share the same convention
    // (φ(0) = (R#(0))^{1/α} = (w₀^α)^{1/α} = w₀).
    if (family === 'boundedQD' || family === 'boundedLQD' || family === 'powerQD'
        || family === 'powerQD_singular') {
      const w0 = phi.w0 || { re: 0, im: 0 };
      addIfInOmega(w0, 'φ(0)');
    }
    // Blaschke center for singular families: z = z₀ (for singular PQDs
    // φ(z₀) = 0, the origin ∈ Ω).
    if (family === 'boundedLQD_singular' || family === 'unboundedLQD_singular'
        || family === 'powerQD_singular' || family === 'unboundedPQD_singular') {
      if (phi.z0) {
        let w;
        try { w = schwarz.evalPhi(phi.z0); } catch (_) { w = null; }
        if (w) addIfInOmega(w, 'φ(z₀)');
      }
    }
    // For unbounded families: try the centroid of pole positions a_j.
    if (family === 'unboundedQD' || family === 'unboundedLQD'
        || family === 'unboundedLQD_singular' || family === 'unboundedPQD'
        || family === 'unboundedPQD_singular') {
      const branches = phi.branches || [];
      if (branches.length > 0) {
        let cx = 0, cy = 0, n = 0;
        for (const br of branches) {
          // Branch z_j is in the disk side; the corresponding w-space "pole" is
          // φ(1/conj(z_j)) (the Schwarz reflection). Use that.
          const absZj = Math.hypot(br.z.re, br.z.im);
          if (absZj < 1e-12) continue;
          const zR = { re: br.z.re / (absZj * absZj), im: br.z.im / (absZj * absZj) };
          let w; try { w = schwarz.evalPhi(zR); } catch (_) { continue; }
          if (!w || !isFinite(w.re) || !isFinite(w.im)) continue;
          cx += w.re; cy += w.im; n++;
        }
        if (n > 0) addIfInOmega({ re: cx / n, im: cy / n }, 'centroid');
      }
    }
    // Always include the origin if it's in Ω (the "natural" point for many
    // analytic-function visualisations).
    addIfInOmega({ re: 0, im: 0 }, 'w=0');
    return out;
  }

  // ---------------------------------------------------------------------------
  // 2. iterateCurveForward — apply σ to every vertex of a polyline, k times.
  //
  // Returns an array of length (k+1): the original polyline, then σ(pts), then
  // σ²(pts), … At each iteration step, vertices that have left Ω or
  // produced an invalid σ result are filtered out. If all vertices leave, the
  // remaining iterations are empty arrays.
  //
  // Used by E11 (forward-image of user-drawn curves) and H8 (orbit-family
  // sweep — apply iterateCurveForward to evenly-spaced seeds).
  // ---------------------------------------------------------------------------
  function iterateCurveForward(pts, schwarz, k) {
    k = (k != null) ? Math.max(0, k | 0) : 5;
    if (!pts || pts.length === 0 || !schwarz || !schwarz.sigma) return [[]];
    const out = [pts.map(p => ({ re: p.re, im: p.im }))];
    let current = out[0];
    for (let it = 0; it < k; it++) {
      const next = [];
      for (const w of current) {
        try {
          if (!schwarz.isInOmega(w)) continue;
        } catch (_) { continue; }
        let sv;
        try { sv = schwarz.sigma(w); } catch (_) { continue; }
        if (!sv || !isFinite(sv.re) || !isFinite(sv.im)) continue;
        next.push({ re: sv.re, im: sv.im });
      }
      out.push(next);
      if (next.length === 0) {
        // Pad with empty arrays so output length is always k+1.
        for (let pad = it + 1; pad < k; pad++) out.push([]);
        break;
      }
      current = next;
    }
    return out;
  }

  // ---------------------------------------------------------------------------
  // 3. findCycles — period-n cycles via Newton on σⁿ(w) = w  (E10).
  //
  // Approach: grid-search over Ω with M seeds (default 24 × 24 = 576). At
  // each seed, do Newton on G(w) = σⁿ(w) − w. F'(w) = (σⁿ)'(w) = Π_{k=0..n-1}
  // σ'(σ^k(w)) — too expensive to compute exactly per step; use forward
  // finite differences. Dedup converged roots in w-space; trace the cycle by
  // iterating σ from each found point.
  //
  // Note: this is intentionally a coarse global search, not a guaranteed
  // exhaustive enumeration. n=1 (fixed points) is typically reliable; n=2,3,…
  // become progressively harder due to spurious convergence and basin
  // narrowing. Caller should treat results as advisory.
  // ---------------------------------------------------------------------------
  function findCycles(schwarz, n, opts) {
    opts = opts || {};
    n = Math.max(1, n | 0);
    const M       = opts.gridSize || 18;
    const maxIter = opts.maxIter  || 30;
    const tol     = opts.tol      || 1e-8;
    const h       = 1e-6;

    if (!schwarz || !schwarz.sigma) return [];

    // Bounding box of Ω from the boundary polygon (or a default).
    const bdy = schwarz.boundaryPts || schwarz._boundaryPts || [];
    let minRe = -2, maxRe = 2, minIm = -2, maxIm = 2;
    if (bdy.length > 0) {
      minRe = Infinity; maxRe = -Infinity;
      minIm = Infinity; maxIm = -Infinity;
      for (const p of bdy) {
        if (p.re < minRe) minRe = p.re;
        if (p.re > maxRe) maxRe = p.re;
        if (p.im < minIm) minIm = p.im;
        if (p.im > maxIm) maxIm = p.im;
      }
      // Inset slightly so we don't seed on the boundary.
      const dx = maxRe - minRe, dy = maxIm - minIm;
      minRe += 0.05 * dx; maxRe -= 0.05 * dx;
      minIm += 0.05 * dy; maxIm -= 0.05 * dy;
    }

    // sigma^n(w) — accumulate n forward iterations.
    function sigmaN(w) {
      let cur = w;
      for (let it = 0; it < n; it++) {
        if (!schwarz.isInOmega(cur)) return null;
        let sv;
        try { sv = schwarz.sigma(cur); } catch (_) { return null; }
        if (!sv || !isFinite(sv.re) || !isFinite(sv.im)) return null;
        cur = sv;
      }
      return cur;
    }

    const roots = [];
    for (let iy = 0; iy < M; iy++) {
      for (let ix = 0; ix < M; ix++) {
        let w = {
          re: minRe + (ix + 0.5) * (maxRe - minRe) / M,
          im: minIm + (iy + 0.5) * (maxIm - minIm) / M,
        };
        if (!schwarz.isInOmega(w)) continue;
        let converged = false;
        for (let it = 0; it < maxIter; it++) {
          const sN = sigmaN(w);
          if (!sN) break;
          const diffR = sN.re - w.re, diffI = sN.im - w.im;
          if (Math.hypot(diffR, diffI) < tol) { converged = true; break; }
          // G'(w) = (σⁿ)'(w) − 1.
          const sNh = sigmaN({ re: w.re + h, im: w.im });
          if (!sNh) break;
          const fpR = (sNh.re - sN.re) / h - 1;
          const fpI = (sNh.im - sN.im) / h - 0;       // ∂/∂x of (im − im) = 0
          const denom = fpR * fpR + fpI * fpI;
          if (denom < 1e-30) break;
          const stepR = -(diffR * fpR + diffI * fpI) / denom;
          const stepI = -(diffI * fpR - diffR * fpI) / denom;
          const nw = { re: w.re + stepR, im: w.im + stepI };
          if (!schwarz.isInOmega(nw)) break;
          w = nw;
        }
        if (!converged) continue;
        // Dedup against existing roots.
        let isDup = false;
        for (const r of roots) {
          if (Math.hypot(r.re - w.re, r.im - w.im) < 1e-4) { isDup = true; break; }
        }
        if (!isDup) roots.push(w);
      }
    }

    // Convert each root into a cycle by iterating σ until we return.
    const cycles = [];
    for (const r of roots) {
      const pts = [r];
      let cur = r;
      for (let it = 0; it < n + 1; it++) {
        if (!schwarz.isInOmega(cur)) break;
        let sv; try { sv = schwarz.sigma(cur); } catch (_) { break; }
        if (!sv || !isFinite(sv.re) || !isFinite(sv.im)) break;
        if (Math.hypot(sv.re - r.re, sv.im - r.im) < 1e-4 && it > 0) break;
        pts.push(sv);
        cur = sv;
      }
      // Detect period: smallest k such that pts[k] ≈ pts[0].
      let period = pts.length - 1;
      for (let k = 1; k < pts.length; k++) {
        if (Math.hypot(pts[k].re - r.re, pts[k].im - r.im) < 1e-4) {
          period = k; break;
        }
      }
      // Skip cycle if it's a sub-period of a shorter cycle already in `cycles`.
      let isSubperiod = false;
      for (const c of cycles) {
        if (period % c.period === 0) {
          // Does any point in pts match any point in c.points?
          for (const p of pts) {
            for (const q of c.points) {
              if (Math.hypot(p.re - q.re, p.im - q.im) < 1e-4) {
                isSubperiod = true; break;
              }
            }
            if (isSubperiod) break;
          }
        }
        if (isSubperiod) break;
      }
      if (!isSubperiod) {
        cycles.push({ period, points: pts.slice(0, period) });
      }
    }
    return cycles;
  }

  // ---------------------------------------------------------------------------
  // 4. sampleSweepSeeds — evenly-spaced seeds on a line or circle (H8).
  //
  //   kind   : 'line' | 'circle'
  //   params : line   → { from: {re,im}, to: {re,im}, n }
  //            circle → { center: {re,im}, radius, n }
  // ---------------------------------------------------------------------------
  function sampleSweepSeeds(kind, params) {
    if (!params || !params.n || params.n < 1) return [];
    const n = params.n | 0;
    const out = [];
    if (kind === 'line') {
      const a = params.from, b = params.to;
      for (let i = 0; i < n; i++) {
        const t = (n === 1) ? 0.5 : i / (n - 1);
        out.push({ re: a.re + t * (b.re - a.re), im: a.im + t * (b.im - a.im) });
      }
    } else if (kind === 'circle') {
      const c = params.center, r = params.radius;
      for (let i = 0; i < n; i++) {
        const th = 2 * Math.PI * i / n;
        out.push({ re: c.re + r * Math.cos(th), im: c.im + r * Math.sin(th) });
      }
    }
    return out;
  }

  // ---------------------------------------------------------------------------
  // 5. domainColoringField — per-pixel HSL of σ(w) over a viewport (F6).
  //
  // Output: Uint8ClampedArray of size 4·W·H (RGBA, premultiplied). Pixels
  // outside Ω get fully transparent.
  //   hue        = arg(σ)/τ                  (0–1)
  //   saturation = 1
  //   lightness  = sigmoid(log|σ|)           (0.15–0.85)
  //
  // 256×256 takes ~150–500 ms depending on family. Caller may want to debounce.
  // ---------------------------------------------------------------------------
  function domainColoringField(schwarz, viewport, opts) {
    opts = opts || {};
    const W = opts.W || 256, H = opts.H || 256;
    if (!schwarz || !schwarz.sigma) return new Uint8ClampedArray(4 * W * H);

    const buf = new Uint8ClampedArray(4 * W * H);
    const dx = (viewport.reMax - viewport.reMin) / (W - 1);
    const dy = (viewport.imMax - viewport.imMin) / (H - 1);
    for (let iy = 0; iy < H; iy++) {
      const im = viewport.imMin + iy * dy;
      for (let ix = 0; ix < W; ix++) {
        const re = viewport.reMin + ix * dx;
        const w = { re, im };
        const p = 4 * (iy * W + ix);
        if (!schwarz.isInOmega(w)) {
          buf[p] = 0; buf[p + 1] = 0; buf[p + 2] = 0; buf[p + 3] = 0;
          continue;
        }
        let sv;
        try { sv = schwarz.sigma(w); } catch (_) { sv = null; }
        if (!sv || !isFinite(sv.re) || !isFinite(sv.im)) {
          buf[p] = 0; buf[p + 1] = 0; buf[p + 2] = 0; buf[p + 3] = 0;
          continue;
        }
        const h = (Math.atan2(sv.im, sv.re) / (2 * Math.PI) + 1) % 1;     // 0–1
        const mag = Math.hypot(sv.re, sv.im);
        const l = 0.5 + 0.35 * Math.tanh(0.5 * Math.log(Math.max(1e-12, mag)));
        const rgb = _hslToRgb(h, 1, l);
        buf[p]     = rgb[0];
        buf[p + 1] = rgb[1];
        buf[p + 2] = rgb[2];
        buf[p + 3] = 200;                                                  // semi-opaque
      }
    }
    return buf;
  }

  function _hslToRgb(h, s, l) {
    function hueRgb(p, q, t) {
      if (t < 0) t += 1; if (t > 1) t -= 1;
      if (t < 1/6) return p + (q - p) * 6 * t;
      if (t < 1/2) return q;
      if (t < 2/3) return p + (q - p) * (2/3 - t) * 6;
      return p;
    }
    if (s === 0) return [l * 255, l * 255, l * 255];
    const q = l < 0.5 ? l * (1 + s) : l + s - l * s;
    const p = 2 * l - q;
    return [
      Math.round(hueRgb(p, q, h + 1/3) * 255),
      Math.round(hueRgb(p, q, h)       * 255),
      Math.round(hueRgb(p, q, h - 1/3) * 255),
    ];
  }

  // ---------------------------------------------------------------------------
  // Wire onto QD.Schwarz.
  // ---------------------------------------------------------------------------
  QD.Schwarz.canonicalSeeds       = canonicalSeeds;
  QD.Schwarz.iterateCurveForward  = iterateCurveForward;
  QD.Schwarz.findCycles           = findCycles;
  QD.Schwarz.sampleSweepSeeds     = sampleSweepSeeds;
  QD.Schwarz.domainColoringField  = domainColoringField;
})();
