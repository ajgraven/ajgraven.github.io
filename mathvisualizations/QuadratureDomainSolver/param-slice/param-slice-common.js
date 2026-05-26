// =============================================================================
// param-slice-common.js -- Pure math kernel for the Parameter-slice tab.
//
// Provides:
//   • ParamRef descriptors      — what "a sweepable parameter" is
//   • listAvailableParams(...)  — enumerate all sweepable params for a scenario
//   • applyParam(scenario,...)  — clone+mutate a scenario at one parameter
//   • formatParamLabel(ref)     — human-readable axis label
//   • classifyResult(result)    — map a solver result to a pixel classification
//   • CLASS_COLORS, CLASS_ORDER — categorical color LUT
//
// No DOM dependencies — same source runs in the main thread and in the
// per-tile Web Worker (see param-slice-pool.js).
//
// Scenario shape (the unit of work for one pixel):
//   {
//     hData: { poles: [{a, principal:[...]}, ...], polyPart: [...] },
//     norm:  { w0?, c?, q?, lqd?, unbounded?, singular? },
//     opts:  any extra QD.solveInverseQD options (newton/continuation/etc.)
//   }
// =============================================================================

(function (global) {
  'use strict';
  const QD = global.QD;
  if (!QD || !QD.Complex) {
    // Loaded before solver — tolerate this; consumers re-resolve QD lazily.
  }
  const C = (QD && QD.Complex) || null;

  // ---------------------------------------------------------------------------
  // Parameter descriptors
  // ---------------------------------------------------------------------------
  //   { kind: 'residueRe', poleIdx, residueIdx }
  //   { kind: 'residueIm', poleIdx, residueIdx }
  //   { kind: 'poleRe',    poleIdx }
  //   { kind: 'poleIm',    poleIdx }
  //   { kind: 'polyRe',    degree }
  //   { kind: 'polyIm',    degree }
  //   { kind: 'cReal' }
  //   { kind: 'qRe' } | { kind: 'qIm' }
  //   { kind: 'w0Re' } | { kind: 'w0Im' }
  // ---------------------------------------------------------------------------

  // Modes (mirror MODES in ui.js)
  const MODE_FAMILY_TAG = {
    'bounded':                undefined,
    'unbounded':              undefined,
    'lqd-bounded':            'boundedLQD',
    'lqd-bounded-singular':   'boundedLQD_singular',
    'lqd-unbounded':          'unboundedLQD',
    'lqd-unbounded-singular': 'unboundedLQD_singular',
  };
  function modeHasC(mode)        { return mode === 'unbounded' || mode === 'lqd-unbounded' || mode === 'lqd-unbounded-singular'; }
  function modeHasQ(mode)        { return mode === 'lqd-bounded-singular' || mode === 'lqd-unbounded-singular'; }
  function modeHasW0Manual(mode) { return mode === 'bounded' || mode === 'lqd-bounded' || mode === 'lqd-bounded-singular'; }
  function modeAllowsPoly(mode)  { return mode === 'unbounded' || mode === 'lqd-unbounded' || mode === 'lqd-unbounded-singular'; }

  // Subscript helper
  const SUBS = ['₀','₁','₂','₃','₄','₅','₆','₇','₈','₉'];
  function sub(n) { return String(n).split('').map(d => SUBS[+d] || d).join(''); }

  // Format a ParamRef as a human-readable label.
  function formatParamLabel(ref) {
    switch (ref.kind) {
      case 'residueRe': return `Re(C${sub(ref.poleIdx+1)},${sub(ref.residueIdx+1)})`;
      case 'residueIm': return `Im(C${sub(ref.poleIdx+1)},${sub(ref.residueIdx+1)})`;
      case 'poleRe':    return `Re(a${sub(ref.poleIdx+1)})`;
      case 'poleIm':    return `Im(a${sub(ref.poleIdx+1)})`;
      case 'polyRe':    return `Re(C∞,${sub(ref.degree)})`;
      case 'polyIm':    return `Im(C∞,${sub(ref.degree)})`;
      case 'cReal':     return `c (conformal radius)`;
      case 'qRe':       return `Re(q)`;
      case 'qIm':       return `Im(q)`;
      case 'w0Re':      return `Re(w₀)`;
      case 'w0Im':      return `Im(w₀)`;
    }
    return '?';
  }

  // Read the current value of a parameter from a scenario. Used to default
  // the [min, max] range to (cur−1, cur+1) in the UI.
  function readParam(scenario, ref) {
    const { hData, norm } = scenario;
    switch (ref.kind) {
      case 'residueRe': return hData.poles[ref.poleIdx].principal[ref.residueIdx].re;
      case 'residueIm': return hData.poles[ref.poleIdx].principal[ref.residueIdx].im;
      case 'poleRe':    return hData.poles[ref.poleIdx].a.re;
      case 'poleIm':    return hData.poles[ref.poleIdx].a.im;
      case 'polyRe':    return (hData.polyPart[ref.degree] || {re:0,im:0}).re;
      case 'polyIm':    return (hData.polyPart[ref.degree] || {re:0,im:0}).im;
      case 'cReal':     return norm.c;
      case 'qRe':       return norm.q ? norm.q.re : 0;
      case 'qIm':       return norm.q ? norm.q.im : 0;
      case 'w0Re':      return norm.w0 ? norm.w0.re : 0;
      case 'w0Im':      return norm.w0 ? norm.w0.im : 0;
    }
    return 0;
  }

  // Enumerate all ParamRefs that make sense for the current scenario+mode.
  // Returns [{ ref, label }, ...].
  function listAvailableParams(scenario, mode) {
    const out = [];
    const push = (ref) => out.push({ ref, label: formatParamLabel(ref) });
    const { hData } = scenario;

    // Pole positions and residues
    for (let j = 0; j < hData.poles.length; j++) {
      push({ kind: 'poleRe', poleIdx: j });
      push({ kind: 'poleIm', poleIdx: j });
      const ord = hData.poles[j].principal.length;
      for (let s = 0; s < ord; s++) {
        push({ kind: 'residueRe', poleIdx: j, residueIdx: s });
        push({ kind: 'residueIm', poleIdx: j, residueIdx: s });
      }
    }
    // Polynomial-at-∞ coefficients
    if (modeAllowsPoly(mode)) {
      for (let l = 0; l < hData.polyPart.length; l++) {
        push({ kind: 'polyRe', degree: l });
        push({ kind: 'polyIm', degree: l });
      }
    }
    // c, q, w0
    if (modeHasC(mode))        push({ kind: 'cReal' });
    if (modeHasQ(mode))      { push({ kind: 'qRe' }); push({ kind: 'qIm' }); }
    if (modeHasW0Manual(mode)) { push({ kind: 'w0Re' }); push({ kind: 'w0Im' }); }
    return out;
  }

  // ---------------------------------------------------------------------------
  // Scenario clone + applyParam
  // ---------------------------------------------------------------------------
  // Deep-clone the small shapes (hData + norm) so per-pixel mutations don't
  // bleed across tiles. The opts bag is shared (treated as immutable).
  function _cloneComplex(z) { return z ? { re: z.re, im: z.im } : z; }
  function _clonePole(p)    { return { a: _cloneComplex(p.a), principal: p.principal.map(_cloneComplex) }; }
  function _cloneHData(h)   { return { poles: h.poles.map(_clonePole), polyPart: (h.polyPart || []).map(_cloneComplex) }; }
  function _cloneNorm(n) {
    const out = {};
    for (const k of Object.keys(n)) {
      const v = n[k];
      if (v && typeof v === 'object' && 're' in v && 'im' in v) out[k] = _cloneComplex(v);
      else out[k] = v;
    }
    return out;
  }
  function cloneScenario(s) {
    return { hData: _cloneHData(s.hData), norm: _cloneNorm(s.norm), opts: s.opts };
  }

  // Apply one parameter value to a scenario in place. Returns the scenario.
  function applyParamInPlace(scenario, ref, value) {
    const { hData, norm } = scenario;
    switch (ref.kind) {
      case 'residueRe': hData.poles[ref.poleIdx].principal[ref.residueIdx].re = value; break;
      case 'residueIm': hData.poles[ref.poleIdx].principal[ref.residueIdx].im = value; break;
      case 'poleRe':    hData.poles[ref.poleIdx].a.re = value; break;
      case 'poleIm':    hData.poles[ref.poleIdx].a.im = value; break;
      case 'polyRe': {
        while (hData.polyPart.length <= ref.degree) hData.polyPart.push({re:0,im:0});
        hData.polyPart[ref.degree].re = value;
        break;
      }
      case 'polyIm': {
        while (hData.polyPart.length <= ref.degree) hData.polyPart.push({re:0,im:0});
        hData.polyPart[ref.degree].im = value;
        break;
      }
      case 'cReal':  norm.c = value; break;
      case 'qRe':    if (!norm.q) norm.q = {re:0,im:0}; norm.q.re = value; break;
      case 'qIm':    if (!norm.q) norm.q = {re:0,im:0}; norm.q.im = value; break;
      case 'w0Re':   if (!norm.w0) norm.w0 = {re:0,im:0}; norm.w0.re = value; break;
      case 'w0Im':   if (!norm.w0) norm.w0 = {re:0,im:0}; norm.w0.im = value; break;
    }
    return scenario;
  }

  // Convenience: clone first, then mutate. Used in tests; the worker uses
  // applyParamInPlace on its own per-tile clone for speed.
  function applyParam(scenario, ref, value) {
    return applyParamInPlace(cloneScenario(scenario), ref, value);
  }

  // ---------------------------------------------------------------------------
  // Classification
  // ---------------------------------------------------------------------------
  // Maps a solver result to one of a small set of categorical classes
  // (mirroring the existing solver's failure taxonomy + the
  // checkLqdPolynomialGap capability refusal).
  //
  // `result` is whatever the worker returns: either a {success,...} bag from
  // QD.solveInverseQD or a {success,phi,residual,iterations,...} bag from a
  // warm-start QD.newtonSolve. For warm-start successes we still need
  // univalent/identityOK to classify "valid" — the worker attaches those.

  const CLASS_VALID            = 'valid';
  const CLASS_IDENTITY_FAIL    = 'identity-fail';
  const CLASS_UNIVALENCE_FAIL  = 'univalence-fail';
  const CLASS_NEWTON_DIVERGED  = 'newton-diverged';
  const CLASS_NO_ROOT          = 'no-root';
  const CLASS_CAPABILITY       = 'capability-refused';
  const CLASS_UNCLASSIFIED     = 'unclassified';

  const CLASS_ORDER = [
    CLASS_VALID,
    CLASS_IDENTITY_FAIL,
    CLASS_UNIVALENCE_FAIL,
    CLASS_NEWTON_DIVERGED,
    CLASS_NO_ROOT,
    CLASS_CAPABILITY,
    CLASS_UNCLASSIFIED,
  ];

  // RGB triples (0..255). Brightness for VALID is modulated by iter count.
  const CLASS_COLORS = {
    [CLASS_VALID]:           [ 40, 170,  80],   // green
    [CLASS_IDENTITY_FAIL]:   [220, 200,  60],   // yellow
    [CLASS_UNIVALENCE_FAIL]: [230, 140,  40],   // orange
    [CLASS_NEWTON_DIVERGED]: [200,  60,  60],   // red
    [CLASS_NO_ROOT]:         [140, 140, 140],   // gray
    [CLASS_CAPABILITY]:      [ 60,  70,  90],   // dark slate
    [CLASS_UNCLASSIFIED]:    [220,  80, 220],   // magenta (debug)
  };

  const CLASS_LABELS = {
    [CLASS_VALID]:           'Valid QD',
    [CLASS_IDENTITY_FAIL]:   'Identity fails',
    [CLASS_UNIVALENCE_FAIL]: 'Boundary self-intersects',
    [CLASS_NEWTON_DIVERGED]: 'Newton diverged',
    [CLASS_NO_ROOT]:         'No algebraic root',
    [CLASS_CAPABILITY]:      'Capability refused',
    [CLASS_UNCLASSIFIED]:    'Unclassified',
  };

  function classifyResult(result) {
    if (!result) return { cls: CLASS_UNCLASSIFIED, iterations: 0 };
    if (!result.success) {
      const err = (result.error || '').toString();
      if (/no algebraic root/i.test(err)) return { cls: CLASS_NO_ROOT, iterations: 0, err };
      // Capability-refused: only genuine "not yet implemented" / "deferred"
      // gates qualify, not every normalizeOpts throw (which prefixes
      // "solveInverseQD: " indiscriminately for any input-shape error).
      if (/not yet implemented|deferred|higher-order pole/i.test(err)) {
        return { cls: CLASS_CAPABILITY, iterations: 0, err };
      }
      if (/iter|line search|jacobian|singular/i.test(err)) {
        return { cls: CLASS_NEWTON_DIVERGED, iterations: result.iterations || 0, err };
      }
      return { cls: CLASS_UNCLASSIFIED, iterations: 0, err };
    }
    const univ = result.univalent !== false;
    const idOK = result.identityOK !== false;
    if (univ && idOK) return { cls: CLASS_VALID, iterations: result.iterations || 0, residual: result.residual };
    if (!univ)        return { cls: CLASS_UNIVALENCE_FAIL, iterations: result.iterations || 0 };
    return                  { cls: CLASS_IDENTITY_FAIL,    iterations: result.iterations || 0 };
  }

  // Map a classification to an RGBA byte tuple.  Iter-count-based intensity
  // applied only to the VALID class.
  function colorFor(classification) {
    const base = CLASS_COLORS[classification.cls] || CLASS_COLORS[CLASS_UNCLASSIFIED];
    if (classification.cls === CLASS_VALID) {
      // Brightness ∝ 1 / (1 + iter/20). Few-iter solves → bright; many-iter → dim.
      const iter = classification.iterations || 0;
      const k = 1 / (1 + iter / 20);
      const bg = 32;   // dark background for the dim end
      return [
        Math.round(bg + (base[0] - bg) * k),
        Math.round(bg + (base[1] - bg) * k),
        Math.round(bg + (base[2] - bg) * k),
        255,
      ];
    }
    return [base[0], base[1], base[2], 255];
  }

  // ---------------------------------------------------------------------------
  // solveOnePoint -- The actual per-pixel solve. Shared by the Web Worker
  // handler and the main-thread fallback so they stay in lock-step.
  //
  // Returns a "result tile" entry:
  //   { cls, iterations, residual, warmUsed, phiSerialized, errSample }
  //
  // Reads the QD namespace from the global (`self.QD` in workers,
  // `window.QD` on the main thread).
  // ---------------------------------------------------------------------------
  function solveOnePoint(scenario, point, warmHint, expectedFamilyTag) {
    const QD = global.QD;
    if (!QD) {
      return { cls: CLASS_UNCLASSIFIED, iterations: 0, residual: NaN,
               warmUsed: false, phiSerialized: null,
               errSample: '(QD namespace not loaded)' };
    }
    const s = cloneScenario(scenario);
    for (let k = 0; k < point.length; k++) applyParamInPlace(s, point[k].ref, point[k].value);
    return _solveScenarioBody(QD, s, warmHint, expectedFamilyTag);
  }
  // Variant of solveOnePoint that mutates a caller-supplied `scratch`
  // scenario in place instead of cloning the base scenario per pixel.
  //
  // SAFETY INVARIANT: every point in a sweep is assumed to assign values
  // for the SAME set of ParamRefs (just different values). The worker /
  // pool dispatches always pass `point` with the same ref shape derived
  // from the sweep's axes, so this holds in practice. Cross-sweep contamination
  // is impossible because the scratch is freshly cloned per `tile` message.
  //
  // Saves the deep-clone of hData (poles + Complex objects) per pixel —
  // ~30 small object allocations × #pixels per tile.
  function solveOnePointWithScratch(scratch, point, warmHint, expectedFamilyTag) {
    const QD = global.QD;
    if (!QD) {
      return { cls: CLASS_UNCLASSIFIED, iterations: 0, residual: NaN,
               warmUsed: false, phiSerialized: null,
               errSample: '(QD namespace not loaded)' };
    }
    for (let k = 0; k < point.length; k++) {
      applyParamInPlace(scratch, point[k].ref, point[k].value);
    }
    return _solveScenarioBody(QD, scratch, warmHint, expectedFamilyTag);
  }

  // Shared body of solveOnePoint / solveOnePointWithScratch: the actual
  // solver call + classification.
  function _solveScenarioBody(QD, s, warmHint, expectedFamilyTag) {
    const opts = Object.assign({}, s.opts || {});
    if (s.norm.w0)        opts.w0 = s.norm.w0;
    if (s.norm.c != null) opts.c  = s.norm.c;
    if (s.norm.q)         opts.q  = s.norm.q;
    if (s.norm.lqd)       opts.lqd = true;
    if (s.norm.unbounded) opts.unbounded = true;
    if (s.norm.singular)  opts.singular = true;

    const canWarm = warmHint &&
      warmHint.family === expectedFamilyTag &&
      warmHint.branches && warmHint.branches.length === s.hData.poles.length &&
      warmHint.branches.every((br, j) =>
        s.hData.poles[j] && br.A.length === s.hData.poles[j].principal.length);

    let resultBag;
    try {
      if (canWarm) {
        const init = QD.clonePhi(warmHint);
        if (init.w0 && s.norm.w0) { init.w0.re = s.norm.w0.re; init.w0.im = s.norm.w0.im; }
        if (s.norm.c != null)     init.c = s.norm.c;
        if (init.q && s.norm.q)   { init.q.re = s.norm.q.re; init.q.im = s.norm.q.im; }
        // Speculative tighter maxIter when the warm hint carries a
        // coarse-pass iteration count (`_coarseIter`, set by the
        // adaptive renderer's nearestPhi). A refined sub-pixel whose
        // neighbour converged in N iters virtually always converges in
        // ~N too; cap at max(12, 2N) to bail early on the few-percent
        // of pixels that land in a different basin. The retry below
        // guarantees we never misclassify those as `newton-diverged`.
        const baseNewton  = opts.newton || { maxIter: 40, tolerance: 1e-9 };
        const baseMaxIter = baseNewton.maxIter || 40;
        const hintIter    = warmHint && warmHint._coarseIter;
        const tightMax    = hintIter
          ? Math.min(baseMaxIter, Math.max(12, hintIter * 2))
          : baseMaxIter;
        let ns = QD.newtonSolve(init, s.hData,
          Object.assign({}, baseNewton, { maxIter: tightMax }));
        if (!ns.success && tightMax < baseMaxIter) {
          // Speculative cap missed — retry once with the full budget.
          // (For non-cap failures this just costs one extra short run.)
          ns = QD.newtonSolve(init, s.hData,
            Object.assign({}, baseNewton, { maxIter: baseMaxIter }));
        }
        if (ns.success) {
          const family = QD.selectFamily(s.norm);
          const phi = family.canonicalizePhi(ns.phi);
          const univalent = QD.isBoundaryUnivalent(phi, opts.univalenceSamples || 64);
          const id = family.verifyQuadratureIdentity(phi, s.hData,
            { numSamples: opts.univalenceSamples || 64 });
          resultBag = {
            success: true, phi, residual: ns.residual,
            iterations: ns.iterations, univalent,
            identityOK: id.maxRelDiff < (opts.identityTol || 1e-6),
            warmUsed: true,
          };
        } else {
          resultBag = _wrapFullSolve(QD.solveInverseQD(s.hData, opts));
        }
      } else {
        resultBag = _wrapFullSolve(QD.solveInverseQD(s.hData, opts));
      }
    } catch (err) {
      resultBag = { success: false, error: String(err && err.message || err) };
    }
    const cls = classifyResult(resultBag);
    return {
      cls: cls.cls,
      iterations: cls.iterations,
      residual: resultBag.residual,
      warmUsed: !!resultBag.warmUsed,
      phiSerialized: (cls.cls === CLASS_VALID && resultBag.phi)
        ? QD.clonePhi(resultBag.phi) : null,
      errSample: resultBag.success ? null : (resultBag.error || '(unknown error)'),
    };
  }

  function _wrapFullSolve(full) {
    if (!full.success) return { success: false, error: full.error };
    const p = full.primary || {};
    return {
      success: true, phi: p.phi, residual: p.residual, iterations: p.iterations,
      univalent: p.univalent, identityOK: p.identityOK,
    };
  }

  // ---------------------------------------------------------------------------
  // Adaptive-mesh helpers
  // ---------------------------------------------------------------------------
  const UNKNOWN_CLASS = 255;
  // Stable integer encoding for the classification names — fits in a Uint8.
  const CLASS_TO_IDX = (function () {
    const m = {};
    for (let i = 0; i < CLASS_ORDER.length; i++) m[CLASS_ORDER[i]] = i;
    return m;
  })();
  const IDX_TO_CLASS = CLASS_ORDER.slice();

  // True iff the 4 corners of the stride-sized cell anchored at (c, r) all
  // have a known classification that agrees. classGrid is a flat
  // Uint8Array of length n0*n1; UNKNOWN_CLASS encodes "not yet sampled".
  function cornersAgree(classGrid, n0, n1, c, r, stride) {
    const c1 = Math.min(c + stride, n0 - 1);
    const r1 = Math.min(r + stride, n1 - 1);
    const k = classGrid[r * n0 + c];
    if (k === UNKNOWN_CLASS) return false;
    return k === classGrid[r * n0 + c1]
        && k === classGrid[r1 * n0 + c]
        && k === classGrid[r1 * n0 + c1];
  }

  // Stronger predicate that also vetoes subdivision-skip when the four
  // corners span a too-large iter-count range. Iter count modulates the
  // VALID class's brightness (see `colorFor`) — without this check, a
  // uniformly-VALID region with a sharp iter gradient renders as visible
  // coarse blocks. For non-VALID classes iter is colour-irrelevant, so
  // we only apply the iter test when classes agree on VALID.
  //
  //   opts.iterDelta : max allowed (max-min) iter spread across corners
  //                    before refinement is triggered. Default 8 (~10–15%
  //                    perceived brightness shift). Pass Infinity to get
  //                    cornersAgree's behaviour exactly.
  //
  // Cheap: 4 grid lookups + at most 4 iter comparisons.
  function cellIsHomogeneous(classGrid, iterGrid, n0, n1, c, r, stride, opts) {
    if (!cornersAgree(classGrid, n0, n1, c, r, stride)) return false;
    const validIdx = CLASS_TO_IDX[CLASS_VALID];
    const k = classGrid[r * n0 + c];
    if (k !== validIdx) return true;
    const iterDelta = (opts && opts.iterDelta != null) ? opts.iterDelta : 8;
    if (!isFinite(iterDelta)) return true;
    const c1 = Math.min(c + stride, n0 - 1);
    const r1 = Math.min(r + stride, n1 - 1);
    const i00 = iterGrid[r * n0 + c];
    const i01 = iterGrid[r * n0 + c1];
    const i10 = iterGrid[r1 * n0 + c];
    const i11 = iterGrid[r1 * n0 + c1];
    const lo = Math.min(i00, i01, i10, i11);
    const hi = Math.max(i00, i01, i10, i11);
    return (hi - lo) <= iterDelta;
  }

  // The 5 new sample points generated by subdividing a stride-cell at (c,r):
  // the four edge-midpoints and the cell center. Returns [{c,r}, ...] within
  // the grid bounds; caller dedupes against previously-evaluated cells.
  function subdivisionPoints(c, r, stride, n0, n1) {
    const half = stride >> 1;
    if (half === 0) return [];
    const out = [];
    const candidates = [
      [c + half,   r],
      [c,          r + half],
      [c + half,   r + half],
      [c + stride, r + half],
      [c + half,   r + stride],
    ];
    for (const [cc, rr] of candidates) {
      if (cc < n0 && rr < n1) out.push({ c: cc, r: rr });
    }
    return out;
  }

  // ---------------------------------------------------------------------------
  // Exports
  // ---------------------------------------------------------------------------
  const exports = {
    // descriptors
    formatParamLabel, listAvailableParams, readParam,
    // mutation
    applyParam, applyParamInPlace, cloneScenario,
    // classification
    classifyResult, colorFor,
    CLASS_VALID, CLASS_IDENTITY_FAIL, CLASS_UNIVALENCE_FAIL,
    CLASS_NEWTON_DIVERGED, CLASS_NO_ROOT, CLASS_CAPABILITY, CLASS_UNCLASSIFIED,
    CLASS_ORDER, CLASS_COLORS, CLASS_LABELS,
    // mode helpers
    modeHasC, modeHasQ, modeHasW0Manual, modeAllowsPoly, MODE_FAMILY_TAG,
    // solver entry
    solveOnePoint, solveOnePointWithScratch,
    // adaptive-mesh helpers
    UNKNOWN_CLASS, CLASS_TO_IDX, IDX_TO_CLASS,
    cornersAgree, cellIsHomogeneous, subdivisionPoints,
  };

  global.ParamSlice = exports;
  if (typeof module !== 'undefined' && module.exports) module.exports = exports;
})(typeof self !== 'undefined' ? self : (typeof window !== 'undefined' ? window : globalThis));
