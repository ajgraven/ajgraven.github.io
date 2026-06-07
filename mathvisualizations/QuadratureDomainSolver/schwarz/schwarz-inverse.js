// =============================================================================
// schwarz-inverse.js — σ⁻¹ kernel + preimage / tile-tree builder (Phase S1).
//
// Math:
//   σ(w)     = conj(F(ψ(w)))    where F is the per-family Schwarz-pullback
//                                of φ and ψ = φ⁻¹.
//   σ⁻¹(w′)  = φ(F⁻¹(conj(w′)))
//
// "Inverting F" means solving F(z) = target for z in the relevant disk
// (𝔻 for bounded, 𝔻* for unbounded). For boundedQD (Phase S1 scope), F has
// the rational form
//
//   F(z) = conj(w₀) + Σⱼ Σₖ A_{j,k} / (z − z_j)^k
//
// Clearing the denominator L(z) = Πⱼ (z − z_j)^{m_j} gives a polynomial
// equation of degree d_L = Σⱼ m_j in z (when conj(w₀) ≠ target). We use
// QD.Direct.polynomialRoots (Durand–Kerner) to get all roots, then filter to
// |z| < 1 − ε to keep only the in-disk preimages. Each surviving z is mapped
// through φ to give a preimage in Ω.
//
// API:
//   QD.Schwarz.sigmaInverse(w, schwarz)
//       → [{re, im}, ...]  preimages of w under σ (length ≤ deg(σ))
//   QD.Schwarz.buildPreimageTree(seed, schwarz, opts)
//       → { generations, edges, truncatedByBudget }
//
// Phase S2 will extend the family-dispatch in `sigmaInverse` to the other
// five families. Phase S3 will reuse `sigmaInverse` for the limit-set chaos
// game + box-counting dimension.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Complex || !QD.Schwarz) {
    throw new Error("schwarz-inverse.js: schwarz-common.js must be loaded first");
  }

  const C = QD.Complex;

  // ---------------------------------------------------------------------------
  // Polynomial helpers (ascending-power Complex[]) — shared QD.Poly
  // (poly-helpers.js, code-review CR3). These are the non-trimming primitives,
  // matching this module's historical semantics (degree feeds the σ⁻¹ root
  // count, so trailing coefficients must be preserved). polyScaleC is the
  // complex-scalar scale.
  // ---------------------------------------------------------------------------
  if (!QD.Poly) throw new Error("schwarz-inverse.js: QD.Poly not found — poly-helpers.js must load first");
  const polyZero        = QD.Poly.zero;
  const polyOne         = QD.Poly.one;
  const polyAdd         = QD.Poly.add;
  const polyMul         = QD.Poly.mul;
  const polyScaleC      = QD.Poly.scale;
  const polyLinearPower = QD.Poly.linearPower;

  // ---------------------------------------------------------------------------
  // boundedQD: invert F(z) = conj(w₀) + Σⱼ Σₖ A_{j,k} / (z − z_j)^k to find
  // all z with F(z) = target.
  //
  // Algorithm:
  //   For each pole j of multiplicity m_j:
  //     L_j(z)   = (z − z_j)^{m_j}                  (pole's denominator)
  //     L⁻ⱼ(z)  = Πⱼ'≠j L_{j'}(z)                  (the "other poles" product)
  //     P_j(z)   = Σₖ A_{j,k} · (z − z_j)^{m_j − k}  (numerator after clearing)
  //   F(z) − target = 0 ⇔ Σⱼ P_j(z)·L⁻ⱼ(z) + (conj(w₀) − target)·L(z) = 0.
  // ---------------------------------------------------------------------------
  function _boundedQDInvertF(phi, target) {
    const branches = phi.branches || [];
    if (branches.length === 0) {
      // F(z) = conj(w₀); inversion is either everything or nothing.
      return [];
    }

    // L_j and L⁻ⱼ for every j. Precompute L_j once; L⁻ⱼ = L / L_j is built
    // by multiplying every L_{j'} for j' ≠ j.
    const m  = branches.map(br => br.A.length);
    const Lj = branches.map((br, j) => polyLinearPower(br.z, m[j]));
    // Build L = Πⱼ L_j once.
    let L = polyOne();
    for (const lj of Lj) L = polyMul(L, lj);

    // For each j, build L⁻ⱼ via polyMul of all Lⱼ' for j'≠j. For small n
    // (typically ≤ 4) this is fine; for larger n we could compute via prefix
    // products, but this isn't a hot path (called once per σ⁻¹ call, not
    // millions).
    const Lmj = branches.map((_, j) => {
      let acc = polyOne();
      for (let jp = 0; jp < branches.length; jp++) {
        if (jp === j) continue;
        acc = polyMul(acc, Lj[jp]);
      }
      return acc;
    });

    // P_j(z) = Σₖ A_{j,k} · (z − z_j)^{m_j − k}, ascending in k from 1 to m_j.
    // For k=m_j the term is A_{j,m_j}·(z−z_j)^0 = constant; for k=1 it's
    // A_{j,1}·(z−z_j)^{m_j−1} (highest-degree term).
    const Pj = branches.map((br, j) => {
      let acc = polyZero();
      for (let k = 1; k <= m[j]; k++) {
        const Ajk = br.A[k - 1];
        const power = polyLinearPower(br.z, m[j] - k);
        acc = polyAdd(acc, polyScaleC(power, Ajk));
      }
      return acc;
    });

    // LHS = Σⱼ P_j · L⁻ⱼ + (conj(w₀) − target) · L
    const w0 = phi.w0 || { re: 0, im: 0 };
    const c0 = { re: w0.re - target.re, im: -w0.im - target.im };  // conj(w₀) − target

    let LHS = polyScaleC(L, c0);
    for (let j = 0; j < branches.length; j++) {
      LHS = polyAdd(LHS, polyMul(Pj[j], Lmj[j]));
    }
    return LHS;
  }

  // ---------------------------------------------------------------------------
  // Generic Newton-based F-inverter (S2). Used for the 5 non-boundedQD
  // families. Distributes nSeeds across the relevant disk (𝔻 for bounded
  // families, 𝔻* for unbounded), runs Newton on F(z) − target = 0 from
  // each, validates by F-residual, dedupes, and filters by |z|.
  //
  // Holomorphic F: F'(z) = ∂F/∂x (Cauchy-Riemann). We approximate via
  // forward finite difference of step h = 1e-7 — same epsilon used by
  // solver.js numericalJacobian.
  //
  // This sacrifices the perfect closed-form root multiplicity of S1's
  // polynomial method (we may miss roots if seeding doesn't cover their
  // basins of attraction) in exchange for unified family support. For
  // typical presets the principal-branch roots are reliably found; if a
  // future application needs richer multi-branch coverage, the
  // boundedQD-style closed-form polynomial method generalises to LQDs via
  // the log-step + principal-branch polynomial approach (sketched in the
  // §13 S2 design).
  // ---------------------------------------------------------------------------
  function _sigmaInverseViaNewton(target, schwarz, opts) {
    opts = opts || {};
    const side    = opts.side    || 'in';      // 'in' = |z|<1, 'out' = |z|>1
    const nSeeds  = opts.nSeeds  || 32;
    const maxIter = opts.maxIter || 40;
    const tol     = opts.tol     || 1e-10;
    const evalF   = schwarz.evalF;
    if (!evalF) return [];

    // Geometric radius grid; angles offset between odd/even seeds for
    // coverage without too much rotational symmetry collisions.
    function seedAt(i) {
      const angleStep = (2 * Math.PI) / Math.max(8, Math.ceil(nSeeds / 4));
      const angle = (i * angleStep) + ((i & 1) ? 0.27 : 0);
      const ringIdx = Math.floor(i / 8);
      if (side === 'in') {
        const r = 0.15 + 0.75 * (ringIdx / Math.max(1, Math.ceil(nSeeds / 8) - 1));
        return { re: r * Math.cos(angle), im: r * Math.sin(angle) };
      } else {
        const r = 1.20 + 2.50 * (ringIdx / Math.max(1, Math.ceil(nSeeds / 8) - 1));
        return { re: r * Math.cos(angle), im: r * Math.sin(angle) };
      }
    }

    const roots = [];
    for (let s = 0; s < nSeeds; s++) {
      let z = seedAt(s);
      let converged = false;
      for (let it = 0; it < maxIter; it++) {
        let Fz;
        try { Fz = evalF(z); } catch (_) { break; }
        const diffR = Fz.re - target.re;
        const diffI = Fz.im - target.im;
        if (Math.hypot(diffR, diffI) < tol) { converged = true; break; }
        // Numerical F'(z) by forward difference along the real axis.
        const h = 1e-7;
        let Fzp;
        try { Fzp = evalF({ re: z.re + h, im: z.im }); }
        catch (_) { break; }
        const fpR = (Fzp.re - Fz.re) / h;
        const fpI = (Fzp.im - Fz.im) / h;
        const denom = fpR * fpR + fpI * fpI;
        if (denom < 1e-30) break;
        // Complex Newton step: z ← z − (F(z) − target) / F'(z).
        const stepR = -(diffR * fpR + diffI * fpI) / denom;
        const stepI = -(diffI * fpR - diffR * fpI) / denom;
        const nz = { re: z.re + stepR, im: z.im + stepI };
        // Crude bounding to keep the iterate from running off to ∞.
        const nr = Math.hypot(nz.re, nz.im);
        if (!isFinite(nr) || nr > 1e4) break;
        z = nz;
      }
      if (!converged) continue;

      // Validate by residual.
      let Fz;
      try { Fz = evalF(z); } catch (_) { continue; }
      const err = Math.hypot(Fz.re - target.re, Fz.im - target.im);
      if (err > 1e-6) continue;

      // Side filter.
      const absZ = Math.hypot(z.re, z.im);
      if (side === 'in'  && absZ >= 1 - 1e-6) continue;
      if (side === 'out' && absZ <= 1 + 1e-6) continue;

      // Dedupe against already-found roots.
      let isDup = false;
      for (const r of roots) {
        if (Math.hypot(r.re - z.re, r.im - z.im) < 1e-5) { isDup = true; break; }
      }
      if (!isDup) roots.push(z);
    }
    return roots;
  }

  // Validate a candidate σ-preimage by the round-trip σ(σ⁻¹(w)) ≈ w
  // test. Used after both polynomial root-finding (boundedQD) and the
  // Newton-based generic inverter (other 5 families).
  function _validatePreimage(wPre, w, schwarz) {
    try {
      const sw = schwarz.sigma(wPre);
      if (!sw) return false;
      const err = Math.hypot(sw.re - w.re, sw.im - w.im);
      return err < 1e-3;
    } catch (_) { return false; }
  }

  // ---------------------------------------------------------------------------
  // sigmaInverse(w, schwarz) — family-dispatched.
  //   boundedQD       : closed-form polynomial (S1)
  //   unboundedQD     : closed-form polynomial in 𝔻* (S2)
  //   boundedLQD*     : numerical Newton in 𝔻  (S2)
  //   unboundedLQD*   : numerical Newton in 𝔻* (S2)
  // ---------------------------------------------------------------------------
  function sigmaInverse(w, schwarz) {
    if (!schwarz || !schwarz._phi) {
      throw new Error("sigmaInverse: schwarz handle missing _phi (was it built via buildSchwarzFromPhi?)");
    }
    const phi    = schwarz._phi;
    const family = schwarz.family;
    const target = C.conj(w);                                     // σ⁻¹(w) = φ(F⁻¹(conj(w)))

    // ----- boundedQD: closed-form polynomial F-inversion (S1). -----
    if (family === 'boundedQD') {
      const lhs = _boundedQDInvertF(phi, target);
      while (lhs.length > 1 && Math.hypot(lhs[lhs.length - 1].re, lhs[lhs.length - 1].im) < 1e-300) {
        lhs.pop();
      }
      let rootsZ;
      try { rootsZ = QD.Direct.polynomialRoots(lhs); }
      catch (_) { return []; }
      const preimages = [];
      for (const z of rootsZ) {
        if (Math.hypot(z.re, z.im) >= 1 - 1e-6) continue;
        const wPre = schwarz.evalPhi(z);
        if (!_validatePreimage(wPre, w, schwarz)) continue;
        preimages.push(wPre);
      }
      return preimages;
    }

    // ----- powerQD with INTEGER α: same polynomial structure as boundedQD
    //       with the αth lifting applied: F(z) = (R(z))^{1/α} = target ⇔
    //       R(z) = target^α. Use the existing inverter with a synthetic
    //       w₀ = w₀^α (so `conj(synth_w₀) = conj(w₀)^α = r(∞)`) and
    //       target' = target^α. For NON-integer α, target^α is not a
    //       polynomial in z after clearing denominators, so we fall through
    //       to the generic Newton-based F-inverter below (which uses the
    //       adapter's evalF = (R(z))^{1/α} closure, valid for any α).
    if (family === 'powerQD' && Number.isInteger(phi.alpha)) {
      const alpha = phi.alpha;
      let targetPow = { re: 1, im: 0 };
      for (let k = 0; k < alpha; k++) targetPow = C.mul(targetPow, target);
      let w0Pow = { re: 1, im: 0 };
      for (let k = 0; k < alpha; k++) w0Pow = C.mul(w0Pow, phi.w0 || { re: 0, im: 0 });
      const phiSynth = Object.assign({}, phi, { w0: w0Pow });
      const lhs = _boundedQDInvertF(phiSynth, targetPow);
      while (lhs.length > 1 && Math.hypot(lhs[lhs.length - 1].re, lhs[lhs.length - 1].im) < 1e-300) {
        lhs.pop();
      }
      let rootsZ;
      try { rootsZ = QD.Direct.polynomialRoots(lhs); }
      catch (_) { return []; }
      const preimages = [];
      for (const z of rootsZ) {
        if (Math.hypot(z.re, z.im) >= 1 - 1e-6) continue;
        let wPre;
        try { wPre = schwarz.evalPhi(z); } catch (_) { continue; }
        if (!_validatePreimage(wPre, w, schwarz)) continue;
        // Dedupe in w-space (the αth-root branching can produce equal w's
        // from distinct z's).
        let isDup = false;
        for (const wp of preimages) {
          if (Math.hypot(wp.re - wPre.re, wp.im - wPre.im) < 1e-5) { isDup = true; break; }
        }
        if (!isDup) preimages.push(wPre);
      }
      return preimages;
    }

    // ----- All other families (the 5 LQD/unbounded variants AND non-integer
    //       powerQD): numerical Newton via the F closure on the schwarz
    //       handle. The disk-side filter depends on whether φ maps 𝔻
    //       (bounded) or 𝔻* (unbounded) onto Ω. powerQD is bounded ⇒
    //       side='in', non-lqd ⇒ 24 seeds.
    const isUnbounded = (family === 'unboundedQD'
                      || family === 'unboundedLQD'
                      || family === 'unboundedLQD_singular'
                      || family === 'unboundedPQD'
                      || family === 'unboundedPQD_singular');
    const side    = isUnbounded ? 'out' : 'in';
    // More seeds for LQD families because the principal-branch landscape
    // is more complex (exp factor); fewer for the pure-rational cases.
    const lqd     = family.endsWith('LQD') || family.endsWith('LQD_singular');
    const nSeeds  = lqd ? 64 : 24;
    const rootsZ  = _sigmaInverseViaNewton(target, schwarz, { side, nSeeds });
    const preimages = [];
    for (const z of rootsZ) {
      let wPre;
      try { wPre = schwarz.evalPhi(z); } catch (_) { continue; }
      if (!_validatePreimage(wPre, w, schwarz)) continue;
      // Dedupe in w-space (different z's can converge to numerically-equal w).
      let isDup = false;
      for (const wp of preimages) {
        if (Math.hypot(wp.re - wPre.re, wp.im - wPre.im) < 1e-5) { isDup = true; break; }
      }
      if (!isDup) preimages.push(wPre);
    }
    return preimages;
  }

  // ---------------------------------------------------------------------------
  // Preimage tree builder (Phase S1 visualization core).
  //
  // Tree nodes by generation; generation k has up to deg(σ)^k nodes capped by
  // visualBudget. Each child carries an `edge` reference (index of parent
  // within previous generation) so the renderer can draw parent→child lines.
  // ---------------------------------------------------------------------------
  function buildPreimageTree(seed, schwarz, opts) {
    opts = opts || {};
    const depth         = (opts.depth != null) ? Math.max(0, opts.depth | 0) : 4;
    const visualBudget  = (opts.visualBudget != null) ? Math.max(1, opts.visualBudget | 0) : 4096;

    const generations = [[{ re: seed.re, im: seed.im }]];
    const edges = [];                                       // [{fromGen, fromIdx, toGen, toIdx}, ...]
    let truncatedByBudget = false;
    let totalNodes = 1;

    for (let g = 0; g < depth; g++) {
      const parents = generations[g];
      const next = [];
      for (let p = 0; p < parents.length; p++) {
        if (totalNodes >= visualBudget) { truncatedByBudget = true; break; }
        let preimages;
        try { preimages = sigmaInverse(parents[p], schwarz); }
        catch (_) { preimages = []; }
        for (const wPre of preimages) {
          if (totalNodes >= visualBudget) { truncatedByBudget = true; break; }
          edges.push({ fromGen: g, fromIdx: p, toGen: g + 1, toIdx: next.length });
          next.push(wPre);
          totalNodes++;
        }
        if (truncatedByBudget) break;
      }
      generations.push(next);
      if (truncatedByBudget) break;
      if (next.length === 0) break;                         // no further preimages
    }

    return { generations, edges, truncatedByBudget };
  }

  // ---------------------------------------------------------------------------
  // Phase S3: Limit-set chaos game.
  //
  // Random walk through σ⁻¹: at each step, pick one of the preimages
  // uniformly at random and continue. After a burnIn warmup, the trajectory
  // densely approximates the limit set — the boundary of the tiling set.
  //
  // Returns a Float64Array(2·n) of interleaved [re, im] coords.
  //
  // Restart logic: if at any step σ⁻¹ yields zero preimages (the current
  // point has no σ-preimage in Ω), we restart from a random Ω^c seed by
  // sampling within a coarse bounding box of the boundary polygon. After
  // 100 consecutive failed restarts we give up and return what we have.
  // ---------------------------------------------------------------------------
  function sampleLimitSet(schwarz, opts) {
    opts = opts || {};
    const n      = (opts.n != null) ? Math.max(0, opts.n | 0) : 50000;
    const burnIn = (opts.burnIn != null) ? Math.max(0, opts.burnIn | 0) : 500;
    const rng    = opts.rng || Math.random;

    const out = new Float64Array(2 * n);
    if (n === 0) return out;

    // Build a bounding box of the boundary polygon (lazily inferred from the
    // schwarz handle). We sample restart seeds from a bbox padded 1.5×.
    const bdy = schwarz.boundaryPts || schwarz._boundaryPts || [];
    let bbMinRe = -2, bbMaxRe = 2, bbMinIm = -2, bbMaxIm = 2;
    if (bdy.length > 0) {
      bbMinRe = Infinity; bbMaxRe = -Infinity;
      bbMinIm = Infinity; bbMaxIm = -Infinity;
      for (const p of bdy) {
        if (p.re < bbMinRe) bbMinRe = p.re;
        if (p.re > bbMaxRe) bbMaxRe = p.re;
        if (p.im < bbMinIm) bbMinIm = p.im;
        if (p.im > bbMaxIm) bbMaxIm = p.im;
      }
      const dx = bbMaxRe - bbMinRe;
      const dy = bbMaxIm - bbMinIm;
      bbMinRe -= 0.25 * dx; bbMaxRe += 0.25 * dx;
      bbMinIm -= 0.25 * dy; bbMaxIm += 0.25 * dy;
    }

    function randomOmegaCSeed() {
      // Reject-sample inside the padded bbox for a point in Ω^c.
      for (let tries = 0; tries < 200; tries++) {
        const w = {
          re: bbMinRe + rng() * (bbMaxRe - bbMinRe),
          im: bbMinIm + rng() * (bbMaxIm - bbMinIm),
        };
        if (!schwarz.isInOmega(w)) return w;
      }
      // No luck — fall back to a corner of the bbox (often in Ω^c for typical Ω).
      return { re: bbMaxRe, im: bbMaxIm };
    }

    let cur = opts.seed ? { re: opts.seed.re, im: opts.seed.im } : randomOmegaCSeed();
    let written = 0;
    let consecRestarts = 0;
    const total = burnIn + n;
    let stepIdx = 0;

    while (written < n && consecRestarts < 100) {
      let preimages;
      try { preimages = sigmaInverse(cur, schwarz); }
      catch (_) { preimages = []; }
      if (preimages.length === 0) {
        cur = randomOmegaCSeed();
        consecRestarts++;
        continue;
      }
      consecRestarts = 0;
      const pick = preimages[Math.floor(rng() * preimages.length)] || preimages[0];
      cur = pick;
      stepIdx++;
      if (stepIdx > burnIn) {
        out[2 * written]     = cur.re;
        out[2 * written + 1] = cur.im;
        written++;
      }
      // Failsafe budget — give up after 20·total steps to avoid infinite loops
      // on pathological φ.
      if (stepIdx > 20 * total) break;
    }

    // Trim if we exited early.
    if (written < n) return out.subarray(0, 2 * written);
    return out;
  }

  // ---------------------------------------------------------------------------
  // Box-counting dimension (S3, paired with sampleLimitSet).
  //
  // For a sequence of box sizes ε, count the number of grid cells of size ε
  // that contain at least one point. log(count) ≈ −d · log(ε) + const.
  // Least-squares slope estimate gives d (box-counting = Hausdorff for nice sets).
  //
  // Inputs:
  //   points    Float64Array of interleaved [re, im, re, im, ...] (from
  //             sampleLimitSet), OR an array of {re, im} (also accepted).
  //   opts.boxSizes  optional explicit ε array (geometric; default 2^{-3..-12}).
  //
  // Returns { boxSizes: number[], counts: number[], slope: number,
  //           intercept: number, dim: number }.
  //
  // Uses a Map<string, true> for occupancy — fine up to ~10^6 boxes.
  // ---------------------------------------------------------------------------
  function boxCountingDimension(points, opts) {
    opts = opts || {};
    const sizes = opts.boxSizes || [
      0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625,
      0.001953125, 0.0009765625, 0.00048828125, 0.000244140625,
    ];

    // Normalise input into Float64Array of interleaved coords.
    let pts;
    if (ArrayBuffer.isView(points)) {
      pts = points;
    } else if (Array.isArray(points)) {
      pts = new Float64Array(points.length * 2);
      for (let i = 0; i < points.length; i++) {
        pts[2 * i]     = points[i].re;
        pts[2 * i + 1] = points[i].im;
      }
    } else {
      throw new TypeError('boxCountingDimension: points must be Float64Array or {re,im}[]');
    }

    const nPts = pts.length / 2;
    const counts = new Array(sizes.length).fill(0);
    for (let s = 0; s < sizes.length; s++) {
      const eps = sizes[s];
      const occ = new Map();
      for (let i = 0; i < nPts; i++) {
        const bx = Math.floor(pts[2 * i]     / eps);
        const by = Math.floor(pts[2 * i + 1] / eps);
        const key = bx + ',' + by;
        if (!occ.has(key)) occ.set(key, true);
      }
      counts[s] = occ.size;
    }

    // Linear regression on (log ε, log N(ε)). Slope = −d.
    let sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
    let nValid = 0;
    for (let s = 0; s < sizes.length; s++) {
      if (counts[s] < 2) continue;                               // skip degenerate
      const lx = Math.log(sizes[s]);
      const ly = Math.log(counts[s]);
      sumX += lx; sumY += ly; sumXY += lx * ly; sumXX += lx * lx;
      nValid++;
    }
    let slope = NaN, intercept = NaN, dim = NaN;
    if (nValid >= 2) {
      const denom = nValid * sumXX - sumX * sumX;
      if (denom !== 0) {
        slope     = (nValid * sumXY - sumX * sumY) / denom;
        intercept = (sumY - slope * sumX) / nValid;
        dim       = -slope;
      }
    }
    return { boxSizes: sizes, counts, slope, intercept, dim };
  }

  // ---------------------------------------------------------------------------
  // Wire onto the QD.Schwarz namespace.
  // ---------------------------------------------------------------------------
  QD.Schwarz.sigmaInverse           = sigmaInverse;
  QD.Schwarz.buildPreimageTree      = buildPreimageTree;
  QD.Schwarz.sampleLimitSet         = sampleLimitSet;
  QD.Schwarz.boxCountingDimension   = boxCountingDimension;

  // Lightweight handle augmenter: schwarz-common.js doesn't know about this
  // file, so we monkey-patch buildSchwarzFromPhi to stash the phi reference on
  // the returned handle (used by sigmaInverse). This is a safe no-op for older
  // handles that still carry _phi via the existing buildFromAdapter path.
  const origBuild = QD.Schwarz.buildSchwarzFromPhi;
  if (origBuild && !origBuild._sigmaInversePatched) {
    QD.Schwarz.buildSchwarzFromPhi = function patchedBuildSchwarzFromPhi(phi, hData, boundaryPts) {
      const handle = origBuild(phi, hData, boundaryPts);
      if (handle) {
        handle._phi          = phi;
        handle._boundaryPts  = boundaryPts;
      }
      return handle;
    };
    QD.Schwarz.buildSchwarzFromPhi._sigmaInversePatched = true;
  }
  const origBuildRat = QD.Schwarz.buildSchwarzFromRational;
  if (origBuildRat && !origBuildRat._sigmaInversePatched) {
    QD.Schwarz.buildSchwarzFromRational = function patchedBuildSchwarzFromRational(phi, boundaryPts) {
      const handle = origBuildRat(phi, boundaryPts);
      if (handle) {
        handle._phi          = phi;
        handle._boundaryPts  = boundaryPts;
      }
      return handle;
    };
    QD.Schwarz.buildSchwarzFromRational._sigmaInversePatched = true;
  }
})();
