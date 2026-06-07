// =============================================================================
// symmetry.js  —  Detect the symmetry group of a quadrature domain Ω from its
// conformal map φ (a lightweight detector; partially delivers roadmap #11).
//
//   QD.detectSymmetry(phi, opts) → {
//     rotationalOrder,   // n for a Z_n / D_n domain (Infinity ⇒ a disk);  1 ⇒ none
//     reflectionAxes,    // [angleRad, …] in [0, π) of the mirror lines (through center)
//     center,            // {re, im} — the symmetry center (boundary centroid)
//     confidence,        // 0..1
//     continuous,        // true for a rotationally-continuous Ω (disk)
//   }
//
// HOW IT WORKS — the conformal map intertwines disk and plane symmetries. If Ω
// has Z_n symmetry about c then φ(e^{2πi/n}·z) = c + e^{2πi/n}·(φ(z) − c); on the
// boundary z = e^{iθ} this is a pure index shift by M/n samples, so the test is
// EXACT (no interpolation) when the sample count M is divisible by n. Likewise a
// reflection across the axis through c at angle α corresponds to the disk
// reflection θ ↦ 2β − θ (β the on-axis sample), again an exact index map. We
// therefore sample ∂Ω at M = 2520 points (divisible by every order 2..12 except
// 11) and test the shift / reflection index identities directly. Robust for
// non-star-shaped boundaries (deltoid cusps included) since it never resamples by
// polar angle. Pure + DOM-free — loads in node-test like critical-set.js / cusps.js.
// =============================================================================

(function (global) {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD) return;

  const DEFAULT_M        = 2520;   // divisible by 2,3,4,5,6,7,8,9,10,12
  const DEFAULT_MAXORDER = 12;
  const DEFAULT_TOL      = 1e-4;   // relative to domain scale; tolerant of solved-φ noise
  const CIRCLE_TOL       = 1e-3;   // radius flatness ⇒ continuous (a disk)

  const TRIVIAL = () => ({
    rotationalOrder: 1, reflectionAxes: [], center: { re: 0, im: 0 },
    confidence: 0, continuous: false,
  });

  function detectSymmetry(phi, opts) {
    opts = opts || {};
    if (!phi || typeof QD.evalPhi !== 'function') return TRIVIAL();

    const M        = opts.samples  || DEFAULT_M;
    const maxOrder = opts.maxOrder || DEFAULT_MAXORDER;
    const rotTol   = (opts.rotTol  != null ? opts.rotTol  : DEFAULT_TOL);
    const reflTol  = (opts.reflTol != null ? opts.reflTol : DEFAULT_TOL);

    // Sample ∂Ω and its centroid (the only sensible symmetry center).
    let cx = 0, cy = 0;
    const w = new Array(M);
    for (let k = 0; k < M; k++) {
      const th = 2 * Math.PI * k / M;
      let p;
      try { p = QD.evalPhi({ re: Math.cos(th), im: Math.sin(th) }, phi); }
      catch (e) { return TRIVIAL(); }
      if (!p || !isFinite(p.re) || !isFinite(p.im)) return TRIVIAL();
      w[k] = p; cx += p.re; cy += p.im;
    }
    cx /= M; cy /= M;
    const center = { re: cx, im: cy };

    // Recenter; domain scale = max boundary radius; radius range ⇒ circle test.
    const pc = new Array(M);
    let scale = 0, rMin = Infinity, rMax = 0, kFar = 0, kNear = 0;
    for (let k = 0; k < M; k++) {
      const dr = w[k].re - cx, di = w[k].im - cy;
      pc[k] = { re: dr, im: di };
      const r = Math.hypot(dr, di);
      if (r > scale) scale = r;
      if (r > rMax) { rMax = r; kFar = k; }
      if (r < rMin) { rMin = r; kNear = k; }
    }
    if (!(scale > 0)) return TRIVIAL();

    const continuous = rMax > 0 && (rMax - rMin) / rMax < CIRCLE_TOL;

    // ---- rotational order: largest n | M with the exact shift identity --------
    let order = 1, rotMargin = 1;
    const rotEps = rotTol * scale;
    for (let n = maxOrder; n >= 2; n--) {
      if (M % n !== 0) continue;
      const s = M / n;
      const ang = 2 * Math.PI / n, cs = Math.cos(ang), sn = Math.sin(ang);
      let maxErr = 0;
      for (let k = 0; k < M; k++) {
        const p = pc[k];
        const rotRe = p.re * cs - p.im * sn;
        const rotIm = p.re * sn + p.im * cs;
        const t = pc[(k + s) % M];
        const e = Math.hypot(rotRe - t.re, rotIm - t.im);
        if (e > maxErr) { maxErr = e; if (maxErr > rotEps) break; }
      }
      if (maxErr <= rotEps) { order = n; rotMargin = maxErr / scale; break; }
    }

    // ---- reflection: test the axis through the farthest (then nearest) point --
    // Refl across the line through c at angle α maps p ↦ e^{2iα}·conj(p), which on
    // the boundary is the index map k ↦ 2·kStar − k (the on-axis sample is fixed).
    const reflEps = reflTol * scale;
    function reflAxisAngle(kStar) {
      const alpha = Math.atan2(pc[kStar].im, pc[kStar].re);
      const c2 = Math.cos(2 * alpha), s2 = Math.sin(2 * alpha);
      let maxErr = 0;
      for (let k = 0; k < M; k++) {
        const cr = pc[k].re, ci = -pc[k].im;            // conj(p)
        const reRe = cr * c2 - ci * s2;
        const reIm = cr * s2 + ci * c2;
        const idx = ((2 * kStar - k) % M + M) % M;
        const t = pc[idx];
        const e = Math.hypot(reRe - t.re, reIm - t.im);
        if (e > maxErr) { maxErr = e; if (maxErr > reflEps) break; }
      }
      return maxErr <= reflEps ? alpha : null;
    }

    let axes = [];
    if (!continuous) {
      let a0 = reflAxisAngle(kFar);
      if (a0 === null) a0 = reflAxisAngle(kNear);
      if (a0 !== null) {
        // A D_n group has n mirror lines spaced π/n; for n=1 a single axis.
        for (let j = 0; j < order; j++) {
          let ax = a0 + j * Math.PI / order;
          ax = ((ax % Math.PI) + Math.PI) % Math.PI;     // canonical axis angle in [0, π)
          axes.push(ax);
        }
        axes.sort((p, q) => p - q);
      }
    }

    const confidence = continuous ? 1
      : (order > 1 || axes.length > 0) ? Math.max(0.3, 1 - rotMargin * 200) : 0;

    return {
      rotationalOrder: continuous ? Infinity : order,
      reflectionAxes: axes,
      center,
      confidence,
      continuous,
    };
  }

  QD.Symmetry = { detectSymmetry };
  QD.detectSymmetry = detectSymmetry;

})(typeof globalThis !== 'undefined' ? globalThis : this);
