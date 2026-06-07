// =============================================================================
// cusps.js  —  Boundary cusp detection + (p,q)-type classification for φ.
//
// The boundary ∂Ω = φ(∂𝔻) develops a CUSP exactly where the conformal map's
// derivative vanishes on the unit circle: φ′(e^{iθ₀}) = 0. The order m of that
// zero fixes the local shape. Writing the boundary curve γ(s) = φ(ẑ·e^{is})
// near the singular point ẑ = e^{iθ₀} as a power series in the arclength
// parameter s,
//        γ(s) − γ(0) = b_p s^p + b_q s^q + …      (p < q, b_p,b_q complex)
// the leading exponent is p = m+1 and the next REAL-independent exponent q is
// generically m+2. We report the pair (p,q): m = 1 ⇒ (2,3), the ordinary
// 3⁄2-power cusp seen on the cardioid/deltoid; higher m ⇒ (m+1, m+2).
//
// A strictly univalent domain has φ′ ≠ 0 on the closed disk, so an EXACT cusp
// is non-generic — it appears precisely as a zero of φ′ crosses ∂𝔻 (the
// realizability boundary). We therefore also report PROXIMITY: each boundary-
// adjacent zero of φ′ is listed with its distance d = ‖z‖ − 1 from ∂𝔻 and an
// `isCusp` flag (true once d is within tolerance); a near-but-not-on-boundary
// zero is an "incipient" cusp whose type is what it WILL be at the bifurcation.
//
// Two independent estimators (per the design): the (p,q) type is computed
// EXACTLY from φ's Taylor coefficients (we have closed-form φ), and a NUMERICAL
// log–log fit of ‖γ(θ₀+δ)−γ(θ₀)‖ vs δ cross-checks it (and distinguishes a
// genuine cusp, leading slope > 1, from a smooth near-cusp, slope ≈ 1).
//
// Interior-angle / corner estimation is deliberately deferred: for these
// analytic rational/algebraic φ a φ′-zero yields interior angle 0 or 2π, so a
// general turning-angle estimator is a follow-up.
//
// API:
//   QD.classifyCusps(phi, opts) → {
//     cusps: [ {
//       thetaDeg,            // boundary angle θ₀ of the (incipient) cusp, deg
//       w: {re,im},          // φ(e^{iθ₀}) — the cusp tip in the w-plane
//       zZero: {re,im},      // the actual zero of φ′ (may be slightly off ∂𝔻)
//       dist,                // d = ‖zZero‖ − 1 (signed); |d| ≈ 0 ⇒ real cusp
//       orderM,              // order m of the zero of φ′
//       type: [p, q],        // leading boundary exponents, p = m+1
//       typeLabel,           // "ordinary 3⁄2-cusp" (2,3) | "(p,q) cusp"
//       isCusp,              // |dist| < cuspTol — an actual boundary cusp now
//       numeric: { pLeading },// numerical leading-exponent estimate (→ p for a
//                            //   real cusp, → 1 for a smooth near-cusp)
//       confidence,          // 0..1 agreement of pLeading with the exact p
//     }, ... ],              // sorted by |dist| (closest-to-boundary first)
//     samples, notes,
//   }
//   opts: { cuspTol=5e-3, taylorOrder=8, relTol=1e-2, indepTol=1e-3,
//           fitDeltas } — all tunable; defaults below.
//
// `relTol` is a RELATIVE-gap tolerance: a Taylor coefficient counts as zero
// when it is < relTol · (largest coefficient). This is deliberately loose
// (1e-2) because at a multiple root of φ′ the root is located only to ≈√(tol),
// so the would-be-zero coefficients carry that much error — a relative gap is
// robust where an absolute threshold would misread the order.
//
// Reuses QD.findCriticalPoints (zeros of φ′, app/critical-set.js), QD.phiTaylorAt
// (exact Taylor coeffs, family-agnostic), QD.evalPhi, and QD.Complex. Pure +
// DOM-free; loaded page-side only (like critical-set.js / univalence.js), not
// bundled into the solver Workers.
// =============================================================================

(function (global) {
  'use strict';

  // QD-namespace resolution — same idiom every solver file uses.
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD) return;

  const C = QD.Complex;

  const DEFAULT_CUSP_TOL   = 5e-3;   // |‖z‖−1| below this ⇒ an actual cusp
  const DEFAULT_TAYLOR_L   = 8;      // Taylor order for φ at the singular point
  const DEFAULT_REL_TOL    = 1e-2;   // |coeff| < relTol·scale ⇒ treated as zero
  const DEFAULT_INDEP_TOL  = 1e-3;   // |Im(b_j/b_p)| above this ⇒ independent
  const THETA_WINDOW       = 0.05;   // half-width (rad) of the θ-refine scan
  const THETA_WINDOW_N     = 41;     // samples across the refine window
  const FIT_DELTA0         = 4e-3;   // largest δ for the numerical log–log fit
  const FIT_RATIO          = 0.5;    // geometric δ ratio
  const FIT_N              = 7;      // number of δ samples

  // --- tiny local truncated-series helpers (avoid a hard QD.Taylor dep) ------
  // Series are Complex[] with index = power of s; entry j is the coeff of s^j.
  function _seriesZero(L) {
    const r = new Array(L + 1);
    for (let i = 0; i <= L; i++) r[i] = { re: 0, im: 0 };
    return r;
  }
  // (p · q) truncated to degree L.
  function _seriesMul(p, q, L) {
    const r = _seriesZero(L);
    for (let i = 0; i <= L && i < p.length; i++) {
      if (p[i].re === 0 && p[i].im === 0) continue;
      for (let j = 0; j <= L - i && j < q.length; j++) {
        // r[i+j] += p[i]·q[j]
        r[i + j].re += p[i].re * q[j].re - p[i].im * q[j].im;
        r[i + j].im += p[i].re * q[j].im + p[i].im * q[j].re;
      }
    }
    return r;
  }

  // ---------------------------------------------------------------------------
  // Order of the zero of φ′ at z, read from φ's Taylor coeffs a = [a0,a1,…].
  // φ′(z+t) = Σ_{k≥0} (k+1)·a_{k+1}·t^k, so φ′(z)=a1. At a zero a1≈0 and the
  // order m = (first index j≥2 with a_j ≠ 0) − 1. Returns m≥1, or null if φ′
  // doesn't actually vanish here (a1 not small) or no nonzero term is found.
  // ---------------------------------------------------------------------------
  function _orderOfZero(a, relTol) {
    if (!a || a.length < 3) return null;
    // Scale for the relative zero-test: the largest coefficient (k ≥ 1).
    let scale = 0;
    for (let k = 1; k < a.length; k++) scale = Math.max(scale, Math.hypot(a[k].re, a[k].im));
    if (scale === 0) return null;
    const tol = relTol * scale;
    if (Math.hypot(a[1].re, a[1].im) > tol) return null;   // φ′ not vanishing
    for (let j = 2; j < a.length; j++) {
      if (Math.hypot(a[j].re, a[j].im) > tol) return j - 1;
    }
    return null;
  }

  // ---------------------------------------------------------------------------
  // Boundary expansion coefficients b_j (coeff of s^j) of
  //   γ(s) − γ(0) = Σ_{k ≥ m+1} a_k · t(s)^k,   t(s) = ẑ·(e^{is} − 1).
  // We sum only k ≥ m+1 so the result reflects the CUSP SHAPE, independent of
  // the small a_1..a_m terms that, for a near-cusp (ẑ slightly off the true
  // zero), merely encode the displacement d rather than the singularity.
  // t(s) = Σ_{j≥1} ẑ·(i^j/j!) s^j  (no constant term).
  // ---------------------------------------------------------------------------
  function _boundaryCoeffs(a, zHat, mPlus1, L) {
    // Build T (the t(s) series) up to degree L.
    const T = _seriesZero(L);
    let ip = { re: 1, im: 0 };   // i^0
    let fact = 1;
    for (let j = 1; j <= L; j++) {
      ip = { re: -ip.im, im: ip.re };   // ip *= i  (so ip = i^j)
      fact *= j;                        // j!
      // T[j] = ẑ · (i^j / j!)
      const cr = ip.re / fact, ci = ip.im / fact;
      T[j] = { re: zHat.re * cr - zHat.im * ci, im: zHat.re * ci + zHat.im * cr };
    }
    // Accumulate b = Σ_{k=mPlus1..L} a_k · T^k. Build T^k incrementally.
    const b = _seriesZero(L);
    let Tk = T;                          // T^1
    for (let k = 1; k <= L; k++) {
      if (k >= mPlus1 && k < a.length) {
        const ak = a[k];
        for (let j = 0; j <= L; j++) {
          // b[j] += a_k · Tk[j]
          b[j].re += ak.re * Tk[j].re - ak.im * Tk[j].im;
          b[j].im += ak.re * Tk[j].im + ak.im * Tk[j].re;
        }
      }
      if (k < L) Tk = _seriesMul(Tk, T, L);
    }
    return b;
  }

  // From the boundary coeffs b, find the two leading exponents (p,q): p = first
  // j with |b_j| above tol, q = first j>p with b_j NOT a real multiple of b_p.
  function _leadingExponents(b, relTol, indepTol) {
    let scale = 0;
    for (let j = 1; j < b.length; j++) scale = Math.max(scale, Math.hypot(b[j].re, b[j].im));
    if (scale === 0) return null;
    const tol = relTol * scale;
    let p = -1;
    for (let j = 1; j < b.length; j++) {
      if (Math.hypot(b[j].re, b[j].im) > tol) { p = j; break; }
    }
    if (p < 0) return null;
    const bp = b[p];
    for (let j = p + 1; j < b.length; j++) {
      if (Math.hypot(b[j].re, b[j].im) <= tol) continue;
      // ratio = b_j / b_p; independent ⇔ |Im(ratio)| (relative) > indepTol.
      const ratio = C.div(b[j], bp);
      const rmag = Math.hypot(ratio.re, ratio.im) || 1;
      if (Math.abs(ratio.im) / rmag > indepTol) return { p, q: j };
    }
    return { p, q: null };   // all higher terms collinear with the leading one
  }

  // Least-squares slope of y vs x.
  function _slope(xs, ys) {
    const n = xs.length;
    let sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; }
    const den = n * sxx - sx * sx;
    if (Math.abs(den) < 1e-300) return NaN;
    return (n * sxy - sx * sy) / den;
  }

  // Numerical estimate of the LEADING exponent at boundary angle θ₀: the slope
  // of log‖γ(θ₀+δ)−γ(θ₀)‖ vs log δ → p for a real cusp, → 1 for a smooth
  // near-cusp. Independent of the Taylor path — a robust cross-check on p and a
  // confidence source. (The next exponent q is NOT estimated numerically: the
  // perpendicular-component slope is dominated by floating-point noise at these
  // δ scales, so q is taken from the exact Taylor path only.)
  function _numericalLeadingExponent(phi, theta0, opts) {
    const evalAt = (th) => QD.evalPhi({ re: Math.cos(th), im: Math.sin(th) }, phi);
    let g0;
    try { g0 = evalAt(theta0); } catch (e) { return null; }
    const logD = [], logR = [];
    let d = opts.fitDelta0;
    for (let k = 0; k < opts.fitN; k++) {
      let g;
      try { g = evalAt(theta0 + d); } catch (e) { return null; }
      const r = Math.hypot(g.re - g0.re, g.im - g0.im);
      if (r > 0) { logD.push(Math.log(d)); logR.push(Math.log(r)); }
      d *= opts.fitRatio;
    }
    if (logD.length < 2) return null;
    return _slope(logD, logR);
  }

  // ---------------------------------------------------------------------------
  // Refine the boundary angle θ₀ to the local minimum of |φ′(e^{iθ})|² in a
  // small window around the seed angle (the φ′-zero, projected to ∂𝔻). Cheap
  // dense scan — robust and good enough to anchor the Taylor + numerical passes.
  // ---------------------------------------------------------------------------
  function _refineTheta(phi, thetaSeed, halfWidth, nSamp) {
    let best = thetaSeed, bestVal = Infinity;
    for (let i = 0; i < nSamp; i++) {
      const th = thetaSeed - halfWidth + (2 * halfWidth) * i / (nSamp - 1);
      const z = { re: Math.cos(th), im: Math.sin(th) };
      let t;
      try { t = QD.phiTaylorAt(z, phi, 1); } catch (e) { continue; }
      const v = C.abs2(t[1]);
      if (v < bestVal) { bestVal = v; best = th; }
    }
    return best;
  }

  // ---------------------------------------------------------------------------
  // classifyCusps — public entry point.
  // ---------------------------------------------------------------------------
  function classifyCusps(phi, opts) {
    opts = opts || {};
    const cuspTol  = opts.cuspTol   != null ? opts.cuspTol   : DEFAULT_CUSP_TOL;
    const L        = opts.taylorOrder|| DEFAULT_TAYLOR_L;
    const relTol   = opts.relTol    != null ? opts.relTol    : DEFAULT_REL_TOL;
    const indepTol = opts.indepTol  != null ? opts.indepTol  : DEFAULT_INDEP_TOL;
    const fitOpts  = {
      fitDelta0: opts.fitDelta0 || FIT_DELTA0,
      fitRatio:  opts.fitRatio  || FIT_RATIO,
      fitN:      opts.fitN      || FIT_N,
    };
    const notes = [];

    if (!phi || typeof QD.findCriticalPoints !== 'function') {
      return { cusps: [], samples: 0, notes: ['cusp analysis unavailable'] };
    }

    // 1. Candidate cusps = boundary-adjacent zeros of φ′ (severity 'near').
    let crit;
    try { crit = QD.findCriticalPoints(phi); } catch (e) {
      return { cusps: [], samples: 0, notes: ['critical-point search failed'] };
    }
    const near = (crit.points || []).filter(p => p.severity === 'near');

    const cusps = [];
    for (const cp of near) {
      const z = cp.z;
      const absZ = cp.absZ;
      const dist = absZ - 1;                       // signed distance from ∂𝔻

      // Order m of the zero of φ′, read AT the zero (well-defined off-circle).
      let aZero;
      try { aZero = QD.phiTaylorAt(z, phi, L); } catch (e) { continue; }
      const m = _orderOfZero(aZero, relTol);
      if (m == null) continue;                     // not a genuine φ′-zero here

      // Boundary point ẑ = e^{iθ₀}; refine θ₀ to the |φ′| min along ∂𝔻.
      const theta0 = _refineTheta(phi, Math.atan2(z.im, z.re), THETA_WINDOW, THETA_WINDOW_N);
      const zHat = { re: Math.cos(theta0), im: Math.sin(theta0) };
      let w;
      try { w = QD.evalPhi(zHat, phi); } catch (e) { continue; }

      // Exact (p,q) from the boundary expansion at ẑ (summing k ≥ m+1).
      let aHat;
      try { aHat = QD.phiTaylorAt(zHat, phi, L); } catch (e) { aHat = aZero; }
      const b = _boundaryCoeffs(aHat, zHat, m + 1, L);
      const exps = _leadingExponents(b, relTol, indepTol);
      const p = exps ? exps.p : (m + 1);
      const q = (exps && exps.q != null) ? exps.q : null;
      const type = [p, q != null ? q : (p + 1)];   // fall back to generic q=p+1

      // Numerical cross-check on the leading exponent + confidence.
      const pNum = _numericalLeadingExponent(phi, theta0, fitOpts);
      const numeric = (pNum != null && isFinite(pNum)) ? { pLeading: pNum } : null;
      const confidence = numeric ? Math.max(0, 1 - Math.abs(pNum - type[0])) : 0;

      const isCusp = Math.abs(dist) < cuspTol;
      const ordinary = (type[0] === 2 && type[1] === 3);
      const typeLabel = ordinary ? 'ordinary 3⁄2-cusp' : `(${type[0]},${type[1]}) cusp`;

      cusps.push({
        thetaDeg: theta0 * 180 / Math.PI,
        w: { re: w.re, im: w.im },
        zZero: { re: z.re, im: z.im },
        dist,
        orderM: m,
        type,
        typeLabel,
        isCusp,
        numeric,                 // { pLeading } — numerical leading-exponent estimate
        confidence,
      });
    }

    // Dedupe by boundary angle: a MULTIPLE zero of φ′ converges (linearly) to a
    // cluster of slightly-separated roots that all map to the same boundary
    // point, so findCriticalPoints can return several of them. Merge entries
    // whose θ₀ agree (mod 2π) within ANG_DEDUP, keeping the one closest to ∂𝔻.
    const ANG_DEDUP = 1e-2;   // rad (~0.57°)
    cusps.sort((u, v) => Math.abs(u.dist) - Math.abs(v.dist));
    const merged = [];
    for (const c of cusps) {
      const dup = merged.find(e => {
        let da = Math.abs(e.thetaDeg - c.thetaDeg) * Math.PI / 180;
        da = Math.min(da, 2 * Math.PI - da);     // wraparound (e.g. 180° vs −180°)
        return da < ANG_DEDUP;
      });
      if (!dup) merged.push(c);                  // first (smallest |dist|) wins
    }
    merged.sort((u, v) => Math.abs(u.dist) - Math.abs(v.dist));

    if (merged.length === 0) {
      // Surface the nearest φ′-zero distance as a "how close to a cusp" gauge.
      let nearest = Infinity;
      for (const cp of (crit.points || [])) nearest = Math.min(nearest, Math.abs(cp.absZ - 1));
      if (isFinite(nearest)) notes.push('No boundary cusps; nearest φ′-zero at d = ' + nearest.toFixed(3));
      else notes.push('No critical points of φ′ found.');
    }

    return { cusps: merged, samples: near.length, notes };
  }

  QD.Cusps = { classifyCusps };
  QD.classifyCusps = classifyCusps;

})(typeof globalThis !== 'undefined' ? globalThis : this);
