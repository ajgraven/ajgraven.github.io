// =============================================================================
// solver-cmax.js  —  Estimate the maximum conformal radius c* for which an
// UNBOUNDED quadrature domain still exists.
//
// For an unbounded family the conformal radius c = φ′(∞) is a FREE GAUGE the
// user picks; the solver enforces it and solves for the rest of the map
// φ: 𝔻* → Ω. As c grows the domain grows, and at a critical scale c* its
// boundary develops a cusp and then self-overlaps (the classic Hele-Shaw /
// Polubarinova–Galin blow-up). This module finds c* automatically.
//
//   QD.estimateMaxConformalRadius(hData, baseOpts, solveFn, ctl) → Promise<{
//     found,        // true ⇔ a finite c* was bracketed
//     cMax,         // the estimated largest valid c (= cLowValid when found)
//     cLowValid,    // largest c proven valid (also set when !found, if any)
//     phiAtMax,     // the solved φ at cMax (the largest clean domain)
//     trace,        // [{ c, valid, residual, phase }] — the probe history
//     reason,       // 'bracketed' | 'no-invalid-below-ceiling' | 'no-valid-at-start'
//     ceiling?,     // search ceiling (only on 'no-invalid-below-ceiling')
//     solves,       // number of solveFn invocations spent
//   }>
//
// DESIGN — dependency-injected + pure-ish orchestration:
//   • `solveFn(hData, opts) → result | Promise<result>` is INJECTED. In the
//     browser the UI passes the off-thread QD.PrimarySolverWorker.solve; Node
//     tests pass the synchronous QD.solveInverseQD. Both are await-safe. This
//     keeps the bracket/bisection logic unit-testable with the real solver and
//     with cheap synthetic stubs.
//   • TWO-REGIME VALIDITY GATE. A c is "valid" if the solve is on the physical
//     QD branch. We need TWO criteria because no single one is reliable over the
//     whole range up to c* (lesson from h = 1.5/w + 0.5/w², where c* ≈ 1.46):
//       (1) GENUINE-QD regime (away from the cusp): univalent AND satisfies the
//           quadrature identity (primary.identityOK). This anchors us on a real
//           QD branch and rejects spurious roots that solve (★) but aren't QDs.
//       (2) CUSP regime (near c*): as a φ′ zero migrates onto |z| = 1 the domain
//           Ω develops a cusp and its complement (the hole K) thins to nothing —
//           so the identity's interior test points can no longer clear ∂Ω and
//           the verifier turns UNRELIABLE (maxRelDiff drifts 1e-6 → O(1)) even
//           though the map is still a valid univalent QD. There the geometric
//           CUSP CRITERION is the only trustworthy signal: track g = max|z| over
//           φ′ zeros (QD.findCriticalPoints); the boundary exists while g < 1 and
//           c* is where g → 1. Once g is high (≥ CUSP_NEAR) we gate on
//           univalence + g < 1 and ignore the (unreliable) identity.
//     `success` alone is never the gate — solveInverseQD returns success:true for
//     converged-but-non-univalent fallback roots; we check the fields explicitly
//     and force identityCheck on. When φ′-zero info is unavailable (e.g. synthetic
//     test stubs) we fall back to criterion (1) alone — the classic gate.
//   • WARM-START GAUGE INJECTION (the c-slider-bug lesson, HANDOFF #61): the
//     solver TRUSTS the warm seed's gauge. When we reuse the last valid φ as a
//     seed for a new c we MUST set seed.c = c first, else Newton re-converges to
//     the stale-c solution. We clone + inject c on every warm probe.
//   • CONFIRM-INVALID GUARD: a warm-started solve can fail merely because the
//     seed was bad. Before accepting "invalid" at a given c we retry once
//     WITHOUT a seed (fresh multistart). A bad seed must never shrink c*.
//
// Algorithm: establish a GENUINE-QD lower bracket cLo → grow (warm-started, with
// finer steps once near the cusp) until an invalid cHi (or the ceiling) → bisect
// [cLo,cHi] on the two-regime gate to relTol. cLo tracks the largest known-valid
// c (+ its φ). `mechanism` reports whether c* was set by a cusp or a fold/identity
// boundary.
//
// Pure + DOM-free; loaded page-side only (like critical-set.js / cusps.js), not
// bundled into the solver Workers — it DRIVES a worker rather than running in one.
// It additionally reads QD.findCriticalPoints (also page-side) for the cusp gate.
// =============================================================================

(function () {
  'use strict';

  // QD-namespace resolution — same idiom every solver/analysis file uses.
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD) return;

  const C = QD.Complex;

  // Cusp-regime thresholds. Once a φ′ zero is within CUSP_NEAR of |z| = 1 the
  // domain is treated as "near the cusp": the identity check is no longer trusted
  // and the geometric criterion (g < 1) gates validity. CRIT_BACKSTEP allows a
  // little non-monotonicity in g (sampling noise) without falling off the branch.
  const CUSP_NEAR     = 0.95;
  const CRIT_BACKSTEP = 0.02;
  const CONVERGED_RESID = 1e-6;   // primary.residual above this ⇒ not converged

  // g(φ) = max |z| over φ′ zeros near/inside the unit circle (|z| < 1.05). As the
  // gauge c grows a zero migrates outward; g → 1 marks the cusp (univalence loss).
  // Returns NaN when φ′-zero info is unavailable (no QD.findCriticalPoints, or a
  // synthetic stub φ with no real critical points) — callers then use the genuine
  // identity gate alone.
  function critModulus(phi) {
    if (!phi || typeof QD.findCriticalPoints !== 'function') return NaN;
    try {
      const cs = QD.findCriticalPoints(phi);
      let mx = 0, any = false;
      for (const pt of (cs && cs.points) || []) {
        if (pt.absZ < 1.05) { any = true; if (pt.absZ > mx) mx = pt.absZ; }
      }
      return any ? mx : NaN;
    } catch (_) { return NaN; }
  }

  function residualOf(r) {
    return (r && r.primary && typeof r.primary.residual === 'number') ? r.primary.residual : null;
  }
  function isConverged(r) {
    const res = residualOf(r);
    return res == null || res <= CONVERGED_RESID;
  }

  // Classify a solve result against the two-regime gate (see header).
  //   prevCrit : g at the last accepted c, for branch-continuity in the cusp gate
  // Returns { valid, kind:'genuine'|'cusp'|'invalid', phi, residual, crit }.
  function classify(r, prevCrit) {
    const res = residualOf(r);
    if (!(r && r.success === true && r.primary && r.primary.univalent === true)) {
      return { valid: false, kind: 'invalid', phi: null, residual: res, crit: NaN };
    }
    const p = r.primary;
    const g = critModulus(p.phi);
    // (1) Genuine QD: univalent + quadrature identity holds. Anchors the branch.
    if (p.identityOK === true) {
      return { valid: true, kind: 'genuine', phi: p.phi, residual: res, crit: g };
    }
    // (2) Cusp regime: univalent + converged + a φ′ zero is near |z| = 1 (but still
    // inside), continuing the branch monotonically. The identity is unreliable
    // here (thin hole), so we accept on geometry alone. Requires real g info.
    if (g >= CUSP_NEAR && g < 1 && isConverged(r) &&
        (!(prevCrit > 0) || g >= prevCrit - CRIT_BACKSTEP)) {
      return { valid: true, kind: 'cusp', phi: p.phi, residual: res, crit: g };
    }
    return { valid: false, kind: 'invalid', phi: null, residual: res, crit: g };
  }

  async function estimateMaxConformalRadius(hData, baseOpts, solveFn, ctl) {
    baseOpts = baseOpts || {};
    ctl = ctl || {};
    if (typeof solveFn !== 'function') {
      throw new Error('estimateMaxConformalRadius: solveFn must be a function');
    }

    // ---- controls / defaults -------------------------------------------------
    // growFactor is modest (was 1.6): the physical QD branch can fold, and big
    // multiplicative jumps fall off it onto spurious roots. cuspGrow is finer
    // still, used once g ≥ CUSP_NEAR so we don't overshoot the cusp. maxSolves is
    // raised to fund the smaller steps. identitySamples lifts the verifier's
    // contour-integral resolution for every probe (decoupled from univalence).
    const growFactor  = (ctl.growFactor  > 1) ? ctl.growFactor  : 1.2;
    const cuspGrow    = (ctl.cuspGrow    > 1) ? ctl.cuspGrow    : 1.06;
    const relTol      = (ctl.relTol      > 0) ? ctl.relTol      : 1e-3;
    const maxSolves   = (ctl.maxSolves   > 0) ? ctl.maxSolves   : 80;
    const shrinkTries = (ctl.shrinkTries >= 0) ? (ctl.shrinkTries | 0) : 8;
    const idSamples   = (ctl.identitySamples > 0) ? (ctl.identitySamples | 0) : 3000;
    const progress    = (typeof ctl.progress === 'function') ? ctl.progress : null;

    // Starting c: the caller's hint (the live slider value), else a heuristic
    // mirroring continuationInC_UQD — 0.25 × smallest finite pole magnitude.
    let cStart = ctl.cStart;
    if (!(cStart > 0) || !isFinite(cStart)) {
      let minA = Infinity;
      for (const p of (hData.poles || [])) {
        const m = C.abs(p.a);
        if (m > 0 && m < minA) minA = m;
      }
      cStart = isFinite(minA) ? 0.25 * minA : 0.25;
    }
    // Search ceiling: generous relative cap, never below an absolute floor.
    const cCeiling = (ctl.cCeiling > 0) ? ctl.cCeiling : Math.max(100, cStart * 256);

    // ---- probing -------------------------------------------------------------
    let solves = 0;
    const trace = [];
    function record(c, valid, residual, phase) {
      trace.push({ c, valid, residual, phase });
      if (progress) { try { progress({ phase, c, valid, step: trace.length, solves }); } catch (_) {} }
    }

    // Build per-c options: force the identity gate on with a high-resolution
    // contour integral (idSamples), and inject c into a clone of the warm seed
    // (gauge injection) so Newton can't drift back to a stale c.
    //
    // adaptiveSamples:false (#11) — disable the verifier's near-cusp node
    // escalation here: this estimator already switches to the GEOMETRIC criterion
    // near the cusp (where escalation would otherwise fire), so escalating the
    // identity there is wasted work; away from the cusp the gate never trips. The
    // fixed high idSamples remains the genuine-QD-regime resolution.
    function mergeOpts(c, seedPhi) {
      const opts = Object.assign({}, baseOpts,
        { c: c, identityCheck: true, identitySamples: idSamples, adaptiveSamples: false });
      if (seedPhi) {
        const seed = QD.clonePhi(seedPhi);
        seed.c = c;
        opts.warmPhi = seed;
      } else {
        delete opts.warmPhi;
      }
      return opts;
    }

    async function callSolve(c, seedPhi) {
      solves++;
      try {
        return await solveFn(hData, mergeOpts(c, seedPhi));
      } catch (e) {
        // A solver throw at this c ⇒ no determinable valid QD here.
        return { success: false, error: String((e && e.message) || e) };
      }
    }

    // Classify c against the two-regime gate, with the confirm-invalid guard: a
    // warm failure is retried fresh (full multistart) before being believed, so a
    // bad seed never shrinks c*. `prevCrit` carries g forward for cusp-continuity.
    async function probe(c, seedPhi, prevCrit) {
      let cls = classify(await callSolve(c, seedPhi), prevCrit);
      if (cls.valid) return cls;
      if (seedPhi) {
        const cls2 = classify(await callSolve(c, null), prevCrit);
        if (cls2.valid) return cls2;
        cls = cls2 || cls;
      }
      return cls;
    }

    // ---- 1) lower bracket: a GENUINE QD (identity-confirmed) anchors the branch
    let cLo = cStart;
    let p = await probe(cLo, null, 0);
    record(cLo, p.valid, p.residual, 'bracket-lo');
    let shrinks = 0;
    while (!(p.valid && p.kind === 'genuine') && shrinks < shrinkTries) {
      cLo *= 0.5; shrinks++;
      p = await probe(cLo, null, 0);
      record(cLo, p.valid, p.residual, 'bracket-lo');
    }
    if (!(p.valid && p.kind === 'genuine')) {
      return { found: false, cMax: null, cLowValid: null, phiAtMax: null,
               trace, reason: 'no-valid-at-start', solves };
    }
    let loPhi = p.phi;
    let loCrit = (p.crit > 0) ? p.crit : 0;

    // ---- 2) grow (warm) to an invalid upper bracket cHi --------------------
    // Step finely once near the cusp (loCrit ≥ CUSP_NEAR) so we don't leap over c*
    // onto the non-univalent side or off the branch.
    const stepFactor = () => (loCrit >= CUSP_NEAR ? cuspGrow : growFactor);
    let cHi = cLo * stepFactor();
    let bracketed = false;
    while (cHi <= cCeiling && solves < maxSolves) {
      const pg = await probe(cHi, loPhi, loCrit);
      record(cHi, pg.valid, pg.residual, 'bracket-hi');
      if (pg.valid) {
        cLo = cHi; loPhi = pg.phi;
        if (pg.crit > 0) loCrit = pg.crit;
        cHi = cLo * stepFactor();
      } else {
        bracketed = true;
        break;
      }
    }
    if (!bracketed) {
      return { found: false, cMax: null, cLowValid: cLo, phiAtMax: loPhi,
               trace, reason: 'no-invalid-below-ceiling', ceiling: cCeiling,
               critAtMax: loCrit, solves };
    }

    // ---- 3) bisect [cLo, cHi] on the two-regime gate ----------------------
    while ((cHi - cLo) / cHi > relTol && solves < maxSolves) {
      const mid = 0.5 * (cLo + cHi);
      const pm = await probe(mid, loPhi, loCrit);
      record(mid, pm.valid, pm.residual, 'bisect');
      if (pm.valid) { cLo = mid; loPhi = pm.phi; if (pm.crit > 0) loCrit = pm.crit; }
      else { cHi = mid; }
    }

    // Mechanism that set c*: a forming cusp (a φ′ zero reached |z| = 1) vs a fold /
    // identity-existence boundary away from any cusp.
    const mechanism = (loCrit >= CUSP_NEAR) ? 'cusp' : 'fold';
    // Confidence (#11): how trustworthy the reported c* is, in [0, 1]. Two factors:
    //   • mechanism cleanliness — a cusp whose g climbed close to 1 is an
    //     unambiguous geometric signal (high); a fold / identity-boundary is
    //     inherently softer (≈0.5 floor).
    //   • bracket tightness — the final [cLo, cHi] relative width (≤ relTol).
    const clamp01 = (x) => Math.max(0, Math.min(1, x));
    const relWidth = (cHi > 0) ? (cHi - cLo) / cHi : 1;
    const bracketConf = clamp01(1 - relWidth);
    const mechConf = (mechanism === 'cusp')
      ? 0.6 + 0.4 * clamp01((loCrit - CUSP_NEAR) / (1 - CUSP_NEAR))
      : 0.5;
    const confidence = clamp01(mechConf * bracketConf);
    return { found: true, cMax: cLo, cLowValid: cLo, phiAtMax: loPhi,
             trace, reason: 'bracketed', mechanism, critAtMax: loCrit,
             confidence, solves };
  }

  QD.estimateMaxConformalRadius = estimateMaxConformalRadius;

})();
