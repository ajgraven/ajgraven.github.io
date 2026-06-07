// =============================================================================
// ui-modes.js -- Mode descriptors + aggressiveness presets for the Inverse tab.
//
// Extracted from ui.js by the Phase-3 UI modularization (item E). Exposes a
// QD_UI.installModes(uiCtx) factory (see ui-domain-plot.js for the template).
// ui.js installs it EARLY (right after buildW0 is available) and captures:
//   ({ MODES, PRESETS, modeDescriptor, currentPresetList } =
//        window.QD_UI.installModes(uiCtx));
// so every existing reference (MODES[...], modeDescriptor(), PRESETS[...]) is
// unchanged.
//
// The only ui.js closure these bodies need is buildW0 (called by five families'
// buildNorm). It's read as `ui.buildW0` AT CALL TIME (i.e. at solve time), so it
// only needs to be present on uiCtx before the first solve — not at install.
// `QD` is the global solver namespace (window.QD); the preset arrays
// (QD_PRESETS_BOUNDED, …) are top-level globals from ui-presets.js. `state`
// inside buildNorm/extraHContrib is the explicit PARAMETER passed by
// buildNormalization(hData) — unchanged from the original.
// =============================================================================

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  global.QD_UI.installModes = function installModes(ui) {
    const state = ui.state;   // for modeDescriptor()'s MODES[state.mode] lookup

    // =======================================================================
    // MODE DESCRIPTORS (R5)
    // -----------------------------------------------------------------------
    // Single source of truth for everything that varies between QD/LQD modes:
    //   • the family tag expected on phi
    //   • which UI cards are visible
    //   • which preset list to populate the dropdown with
    //   • how to build the `norm` and route into solver opts
    //   • the vector-field "external" label
    //   • whether auto-escalate runs on solve failure
    //
    // Adding a new mode (e.g. the upcoming unbounded LQDs) is one entry here +
    // one radio in index.html + per-family solver file. No more if/else chains
    // scattered across setMode / buildNormalization / applyNormToOpts /
    // quickSolveAndRender / currentPresetList.
    // =======================================================================
    const MODES = {
      'bounded': {
        label: 'Bounded QD',
        familyTag:        undefined,           // legacy: untagged phi (boundedQD)
        cards: { w0: true, c: false, poly: false, q: false, alpha: false },
        hint: null,
        presets:          () => QD_PRESETS_BOUNDED,
        externalFieldLabel: 'External field: w − h̄(w)',
        externalFieldKind:  'qd',              // 'qd' = w − h̄;  'lqd' = ln|w|²/w̄ − h̄
        vectorFieldOriginAbs2Floor: 1e-30,     // origin not in Ω, no special clip
        extraHContrib:    null,                // no extra terms beyond polyPart + finite poles
        autoEscalate:     true,
        buildNorm(hData, state) {
          const w0 = ui.buildW0(hData);
          if (w0.error) return w0;
          return { w0: w0.w0 };
        },
        applyNorm(opts, norm) { opts.w0 = norm.w0; },
        warmStartUpdate(initPhi, norm) { initPhi.w0 = { re: norm.w0.re, im: norm.w0.im }; },
      },
      'pqd-bounded': {
        label: 'Bounded PQD',
        familyTag:        'powerQD',
        cards: { w0: true, c: false, poly: false, q: false, alpha: true },
        hint:             'pqd-hint',
        presets:          () => QD_PRESETS_BOUNDED_PQD,
        externalFieldLabel: 'External field: (1/α) w |w|^(2(α−1)) − h̄(w)',
        externalFieldKind:  'pqd',
        vectorFieldOriginAbs2Floor: 1e-30,
        extraHContrib:    null,
        autoEscalate:     false,               // PQD existence has a realizability
                                               // floor; auto-escalate doesn't help.
        // In 'auto' w₀-mode w₀ = the live centroid of the poles, recomputed on every
        // solve (so it tracks continuously while a pole is dragged and stays
        // interior — avoids a stale w₀ drifting OUT of Ω on a long drag). The user
        // can still pin a manual w₀.
        buildNorm(hData, state) {
          // α from the dedicated PQD input: any real α > 0, α ≠ 1 (QA milestone).
          const alpha = +state.alpha;
          if (!(alpha > 0) || alpha === 1) {
            return { error: 'PQD power α must be a real number > 0 with α ≠ 1 (α = 1 is classical bounded QD).' };
          }
          const w0 = ui.buildW0(hData);                  // manual value, or the centroid
          if (w0.error) return w0;
          let center = w0.w0;
          // A bounded PQD needs w₀ ≠ 0 (0 ∉ Ω). In auto mode, if the centroid lands
          // ~0 (e.g. poles symmetric about the origin) fall back to the
          // dominant-|residue| pole (always interior, ≠ 0) — the same candidate
          // order as bootstrapW0_PQD, but without the per-frame nested classical
          // solve, so it stays cheap during a drag.
          if (state.w0Mode !== 'manual' && QD.Complex.abs(center) < 1e-9) {
            let best = null, bestMag = -1;
            for (const p of hData.poles) {
              const mag = p.principal.length ? QD.Complex.abs(p.principal[0]) : 0;
              if (mag > bestMag && QD.Complex.abs(p.a) > 1e-9) { bestMag = mag; best = p.a; }
            }
            if (best) center = best;
          }
          if (QD.Complex.abs(center) < 1e-12) {
            return { error: 'Bounded PQD requires w₀ ≠ 0 (0 ∉ Ω assumed).' };
          }
          return { w0: center, alpha };
        },
        applyNorm(opts, norm) {
          // buildNorm always supplies w0 now (centroid in auto, manual otherwise),
          // so opts.w0 is set every solve. The `if` only guards the degenerate
          // no-w0 shape; solver-side bootstrapW0_PQD remains the fallback when w0 is
          // genuinely absent (e.g. headless/direct callers).
          if (norm.w0) opts.w0 = norm.w0;
          opts.alpha = norm.alpha;
        },
        warmStartUpdate(initPhi, norm) {
          if (norm.w0) initPhi.w0 = { re: norm.w0.re, im: norm.w0.im };
          initPhi.alpha = norm.alpha;
        },
      },
      'pqd-bounded-singular': {
        label: 'Bounded singular PQD',
        familyTag:        'powerQD_singular',
        cards: { w0: true, c: false, poly: false, q: false, alpha: true },
        hint:             'pqd-singular-hint',
        presets:          () => QD_PRESETS_BOUNDED_PQD_SINGULAR,
        externalFieldLabel: 'External field: (1/α) w |w|^(2(α−1)) − h̄(w)',
        externalFieldKind:  'pqd',
        vectorFieldOriginAbs2Floor: 1e-30,
        extraHContrib:    null,
        autoEscalate:     false,
        buildNorm(hData, state) {
          // Singular PQD: 0 ∈ Ω (the origin is the Blaschke-zero image φ(z₀)=0),
          // but w₀ = φ(0) is a DIFFERENT interior point and must be nonzero (it
          // appears in the hardwired constant w₀^α/|z₀|^α).
          const w0 = ui.buildW0(hData);
          if (w0.error) return w0;
          if (QD.Complex.abs(w0.w0) < 1e-12) {
            return { error: 'Bounded singular PQD requires w₀ = φ(0) ≠ 0 (a non-origin interior point). Set a manual w₀.' };
          }
          const alpha = +state.alpha;
          if (!(alpha > 0) || alpha === 1) {
            return { error: 'PQD power α must be a real number > 0 with α ≠ 1 (α = 1 is classical bounded QD).' };
          }
          return { w0: w0.w0, alpha, singular: true };
        },
        applyNorm(opts, norm) {
          opts.w0 = norm.w0;
          opts.alpha = norm.alpha;
          opts.singular = true;
        },
        warmStartUpdate(initPhi, norm) {
          initPhi.w0 = { re: norm.w0.re, im: norm.w0.im };
          initPhi.alpha = norm.alpha;
        },
      },
      'pqd-unbounded': {
        label: 'Unbounded PQD',
        familyTag:        'unboundedPQD',
        cards: { w0: false, c: true, poly: true, q: false, alpha: true },
        hint:             'pqd-unbounded-hint',
        presets:          () => QD_PRESETS_UNBOUNDED_PQD,
        externalFieldLabel: 'External field: (1/α) w |w|^(2(α−1)) − h̄(w)',
        externalFieldKind:  'pqd',
        vectorFieldOriginAbs2Floor: 1e-30,
        extraHContrib:    null,
        autoEscalate:     false,
        buildNorm(hData, state) {
          // Unbounded PQD: φ(z)=z·(r#)^{1/α} on 𝔻*, r#(∞)=c^α. c is the conformal
          // radius (user input, as classical unbounded QD); α any real > 0, α ≠ 1.
          const c = +state.c;
          if (!(c > 0) || !isFinite(c)) return { error: 'c must be a positive number' };
          const alpha = +state.alpha;
          if (!(alpha > 0) || alpha === 1) {
            return { error: 'PQD power α must be a real number > 0 with α ≠ 1 (α = 1 is classical unbounded QD).' };
          }
          return { c, alpha, unbounded: true };
        },
        applyNorm(opts, norm) { opts.unbounded = true; opts.c = norm.c; opts.alpha = norm.alpha; },
        warmStartUpdate(initPhi, norm) { initPhi.c = norm.c; initPhi.alpha = norm.alpha; },
      },
      'pqd-unbounded-singular': {
        label: 'Unbounded singular PQD',
        familyTag:        'unboundedPQD_singular',
        cards: { w0: false, c: true, poly: true, q: false, alpha: true },
        hint:             'pqd-unbounded-singular-hint',
        presets:          () => QD_PRESETS_UNBOUNDED_PQD_SINGULAR,
        externalFieldLabel: 'External field: (1/α) w |w|^(2(α−1)) − h̄(w)',
        externalFieldKind:  'pqd',
        vectorFieldOriginAbs2Floor: 1e-30,
        extraHContrib:    null,
        autoEscalate:     false,
        buildNorm(hData, state) {
          // Unbounded singular PQD: 0 ∈ Ω (origin-preimage z₀ ∈ 𝔻*, φ(z₀)=0).
          // No q; the z₀-closure is r(z₀)=0 (Prop 4.6.3). c is the conformal radius.
          const c = +state.c;
          if (!(c > 0) || !isFinite(c)) return { error: 'c must be a positive number' };
          const alpha = +state.alpha;
          if (!(alpha > 0) || alpha === 1) {
            return { error: 'PQD power α must be a real number > 0 with α ≠ 1 (α = 1 is classical unbounded QD).' };
          }
          return { c, alpha, unbounded: true, singular: true };
        },
        applyNorm(opts, norm) { opts.unbounded = true; opts.singular = true; opts.c = norm.c; opts.alpha = norm.alpha; },
        warmStartUpdate(initPhi, norm) { initPhi.c = norm.c; initPhi.alpha = norm.alpha; },
      },
      'unbounded': {
        label: 'Unbounded QD',
        familyTag:        undefined,           // legacy: untagged phi (unboundedQD)
        cards: { w0: false, c: true, poly: true, q: false },
        hint: null,
        presets:          () => QD_PRESETS_UNBOUNDED,
        externalFieldLabel: 'External field: w − h̄(w)',
        externalFieldKind:  'qd',
        vectorFieldOriginAbs2Floor: 1e-30,
        extraHContrib:    null,
        autoEscalate:     true,
        buildNorm(hData, state) {
          const c = +state.c;
          if (!(c > 0) || !isFinite(c)) return { error: 'c must be a positive number' };
          return { c, unbounded: true };
        },
        applyNorm(opts, norm) { opts.unbounded = true; opts.c = norm.c; },
        warmStartUpdate(initPhi, norm) { initPhi.c = norm.c; },
      },
      'lqd-bounded': {
        label: 'Bounded LQD',
        familyTag:        'boundedLQD',
        cards: { w0: true, c: false, poly: false, q: false },
        hint:             'lqd-hint',
        presets:          () => LQD_PRESETS_BOUNDED,
        externalFieldLabel: 'External field: ln|w|²/w̄ − h̄(w)',
        externalFieldKind:  'lqd',
        vectorFieldOriginAbs2Floor: 1e-30,     // 0 ∉ Ω̄, no special clip
        extraHContrib:    null,
        autoEscalate:     false,                // existence is constrained (Thm 5.3.2)
        buildNorm(hData, state) {
          const w0 = ui.buildW0(hData);
          if (w0.error) return w0;
          if (QD.Complex.abs(w0.w0) < 1e-12) {
            return { error: 'LQD mode requires w₀ ≠ 0 (non-singular: 0 ∉ Ω̄). Set a manual w₀.' };
          }
          return { w0: w0.w0, lqd: true };
        },
        applyNorm(opts, norm) { opts.lqd = true; opts.w0 = norm.w0; },
        warmStartUpdate(initPhi, norm) { initPhi.w0 = { re: norm.w0.re, im: norm.w0.im }; },
      },
      'lqd-unbounded': {
        label: 'Unbounded LQD',
        familyTag:        'unboundedLQD',
        cards: { w0: false, c: true, poly: true, q: false },
        hint:             'lqd-unbounded-hint',
        presets:          () => LQD_PRESETS_UNBOUNDED,
        externalFieldLabel: 'External field: ln|w|²/w̄ − h̄(w)',
        externalFieldKind:  'lqd',
        vectorFieldOriginAbs2Floor: 1e-30,    // 0 ∈ K, no special clip
        extraHContrib:    null,
        autoEscalate:     false,
        buildNorm(hData, state) {
          const c = +state.c;
          if (!(c > 0) || !isFinite(c)) return { error: 'c must be a positive number' };
          return { c, lqd: true, unbounded: true };
        },
        applyNorm(opts, norm) { opts.unbounded = true; opts.lqd = true; opts.c = norm.c; },
        warmStartUpdate(initPhi, norm) { initPhi.c = norm.c; },
      },
      'lqd-unbounded-singular': {
        label: 'Unbounded singular LQD',
        familyTag:        'unboundedLQD_singular',
        cards: { w0: false, c: true, poly: true, q: true },
        hint:             'lqd-unbounded-singular-hint',
        presets:          () => LQD_PRESETS_UNBOUNDED_SINGULAR,
        externalFieldLabel: 'External field: ln|w|²/w̄ − h̄(w)',
        externalFieldKind:  'lqd',
        vectorFieldOriginAbs2Floor: 1e-4,      // 0 ∈ Ω; clip arrows near origin
        extraHContrib(w, hData, phi, state) {
          // Singular LQD: h has an extra q/w pole at the origin.
          const q = (phi && phi.q) ? phi.q : QD.Complex.parse(state.q) || { re: 0, im: 0 };
          const denQ = w.re * w.re + w.im * w.im;
          if (denQ < 1e-30) return { re: 0, im: 0 };
          return {
            re: (q.re * w.re + q.im * w.im) / denQ,
            im: (q.im * w.re - q.re * w.im) / denQ,
          };
        },
        autoEscalate:     false,
        buildNorm(hData, state) {
          const c = +state.c;
          if (!(c > 0) || !isFinite(c)) return { error: 'c must be a positive number' };
          const q = QD.Complex.parse(state.q);
          if (!q) return { error: 'Invalid value for q' };
          return { c, q, lqd: true, unbounded: true, singular: true };
        },
        applyNorm(opts, norm) {
          opts.unbounded = true; opts.lqd = true; opts.singular = true;
          opts.c = norm.c; opts.q = norm.q;
        },
        warmStartUpdate(initPhi, norm) {
          initPhi.c = norm.c;
          initPhi.q = { re: norm.q.re, im: norm.q.im };
        },
      },
      'lqd-bounded-singular': {
        label: 'Bounded singular LQD',
        familyTag:        'boundedLQD_singular',
        cards: { w0: true, c: false, poly: false, q: true },
        hint:             'lqd-singular-hint',
        presets:          () => LQD_PRESETS_BOUNDED_SINGULAR,
        externalFieldLabel: 'External field: ln|w|²/w̄ − h̄(w)',
        externalFieldKind:  'lqd',
        vectorFieldOriginAbs2Floor: 1e-4,      // 0 ∈ Ω; clip arrows near origin
        // Singular LQDs add a simple pole of h at w = 0 with residue q.
        extraHContrib(w, hData, phi, state) {
          const q = (phi && phi.q) ? phi.q : QD.Complex.parse(state.q) || { re: 0, im: 0 };
          const denQ = w.re * w.re + w.im * w.im;
          if (denQ < 1e-30) return { re: 0, im: 0 };
          return {
            re: (q.re * w.re + q.im * w.im) / denQ,
            im: (q.im * w.re - q.re * w.im) / denQ,
          };
        },
        autoEscalate:     false,
        buildNorm(hData, state) {
          const w0 = ui.buildW0(hData);
          if (w0.error) return w0;
          if (QD.Complex.abs(w0.w0) < 1e-12) {
            return { error: 'Singular LQD requires w₀ = φ(0) ≠ 0 (preimage 0 ↔ z_0 ≠ 0). Set a manual w₀.' };
          }
          const q = QD.Complex.parse(state.q);
          if (!q) return { error: 'Invalid value for q' };
          return { w0: w0.w0, q, lqd: true, singular: true };
        },
        applyNorm(opts, norm) {
          opts.lqd = true; opts.singular = true; opts.w0 = norm.w0; opts.q = norm.q;
        },
        warmStartUpdate(initPhi, norm) {
          initPhi.w0 = { re: norm.w0.re, im: norm.w0.im };
          initPhi.q  = { re: norm.q.re,  im: norm.q.im  };
        },
      },
    };

    function modeDescriptor() { return MODES[state.mode] || MODES['bounded']; }

    function currentPresetList() {
      return modeDescriptor().presets();
    }

    // =======================================================================
    // Aggressiveness presets
    // -----------------------------------------------------------------------
    // Each entry tunes the four cost knobs of the solver:
    //
    //   numRestarts         — multistart budget AND base for diverse/deflation
    //                         phases AND foreground alternates loop
    //   newton.maxIter      — per-Newton-attempt iteration cap
    //   newton.tolerance    — residual at which Newton declares success
    //   continuation.tStart — initial step in the pole-distance continuation
    //   continuation.growFactor — how aggressively to grow t each successful step
    //   bgAltChunks         — number of background search rounds after a solve
    //   bgAltChunkSize      — restarts per background round
    //
    // Total background alternate-search restarts = bgAltChunks × bgAltChunkSize.
    // To make presets more/less aggressive, just edit the numbers here.
    // "exhaustive" is also wired to the "Try harder" button in the UI.
    // =======================================================================
    const PRESETS = {

      //              | numRestarts |  Newton              |  Continuation              |  bgAltChunks × size
      //              | (a3 + alts) |  maxIter   tolerance |  tStart    growFactor      |  → total bg restarts
      quick: {
        numRestarts:    3,
        newton:       { maxIter:  40, tolerance: 1e-8  },
        continuation: { tStart: 0.20, growFactor: 2.0 },
        bgAltChunks:    8,
        bgAltChunkSize: 4,
      },

      standard: {
        numRestarts:    8,
        newton:       { maxIter:  80, tolerance: 1e-10 },
        continuation: { tStart: 0.10, growFactor: 1.6  },
        bgAltChunks:   20,
        bgAltChunkSize: 6,
      },

      thorough: {
        numRestarts:   20,
        newton:       { maxIter: 150, tolerance: 1e-12 },
        continuation: { tStart: 0.05, growFactor: 1.4  },
        bgAltChunks:   40,
        bgAltChunkSize: 8,
      },

      // Used by the "Try harder" button (and auto-escalation, when enabled in
      // the search-options panel). Much larger multistart budget; deflation is
      // implicit (always on in solveInverseQD once spurious roots appear).
      exhaustive: {
        numRestarts:   60,
        newton:       { maxIter: 200, tolerance: 1e-12 },
        continuation: { tStart: 0.03, growFactor: 1.3  },
        bgAltChunks:   60,
        bgAltChunkSize: 10,
      },

    };

    return { MODES, PRESETS, modeDescriptor, currentPresetList };
  };
})(typeof window !== 'undefined' ? window : globalThis);
