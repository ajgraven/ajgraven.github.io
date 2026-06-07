// =============================================================================
// ui-solve.js -- The solve -> render -> analyze pipeline for the Inverse tab.
//
// Extracted from ui.js by the Phase-3 UI modularization (item E). This is the
// solve subsystem as ONE cohesive module (the three clusters that the plan
// listed separately -- solve, output, analysis -- are mutually recursive and
// share the stale-result token _solveAndRenderToken, so keeping them together
// keeps every intra-pipeline call a bare same-scope call and the token fully
// module-internal). Exposes QD_UI.installSolve(uiCtx); ui.js captures the entry
// points (scheduleSolve / solveAndRender / showSolution / ...) into the
// forward-declared lets with their ORIGINAL names, so all call sites are
// unchanged.
//
// Every ui.js dependency is destructured below with its original name, so the
// function bodies are VERBATIM moves. They are all present on uiCtx by the time
// installSolve runs (the install is at the end of ui.js, after every dependency
// is defined). Small shared render helpers (escapeHTML / formatExp / setStatus)
// and the option/hData builders + publishPrimarySolution stay in ui.js (used by
// the retained Try-harder / cross-tab handlers too) and are read from here.
// =============================================================================

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  global.QD_UI.installSolve = function installSolve(ui) {
    // ---- Injected ui.js dependencies (original names; bodies verbatim) ----
    const state                 = ui.state;
    const $                     = ui.$;
    const MODES                 = ui.MODES;
    const PRESETS               = ui.PRESETS;
    const modeDescriptor        = ui.modeDescriptor;
    const debounce              = ui.debounce;
    const plot                  = ui.plot;
    const setStatus             = ui.setStatus;
    const escapeHTML            = ui.escapeHTML;
    const escapeAttr            = ui.escapeAttr;
    const formatExp             = ui.formatExp;
    const applyModeVisuals      = ui.applyModeVisuals;
    const buildHData            = ui.buildHData;
    const buildNormalization    = ui.buildNormalization;
    const buildSolverOptions    = ui.buildSolverOptions;
    const buildAltSearchOptions = ui.buildAltSearchOptions;
    const applyNormToOpts       = ui.applyNormToOpts;
    const publishPrimarySolution = ui.publishPrimarySolution;

// ---------- Solving (debounced) ------------------------------------------
// Debounce dropped from 250 ms → 60 ms after P0.2 moved solveAndRender off
// the main thread (the worker absorbs the cost; the debounce only needs to
// coalesce keystroke bursts).
const scheduleSolve = debounce(() => { solveAndRender(); }, 60);

// Snappier path used while a slider is being dragged: rAF-throttled, warm-
// starts from the previous solution, skips multistart / continuation /
// alternate-search, and renders without re-fitting the view.
let _quickSolveRaf = null;
function scheduleQuickSolve() {
  if (_quickSolveRaf) return;
  _quickSolveRaf = requestAnimationFrame(() => {
    _quickSolveRaf = null;
    quickSolveAndRender();
  });
}

// Reduced boundary-sample budget for the LIVE (drag) path's univalence +
// identity checks. The debounced full solve re-verifies at the full
// state.samples; during a drag we only need a cheap consistency signal, so we
// cap the per-frame sampling to keep each live frame snappy.
const LIVE_SAMPLES = 96;

// Token bumped on every live (drag) frame. A worker result whose token is no
// longer current — a newer frame superseded it — is dropped so a late paint
// can't clobber newer state. Mirrors _solveAndRenderToken for the full solve.
let _liveSolveToken = 0;

function quickSolveAndRender() {
  const built = buildHData();
  if (!built || built.error) return;
  const norm = buildNormalization(built);
  if (norm.error) return;

  const preset = PRESETS[state.aggressiveness];
  const unbounded = !!norm.unbounded;
  const desc = modeDescriptor();
  const expectedFamilyTag = desc.familyTag;

  // Pick the live Newton seed on the MAIN thread (cheap): warm-start from the
  // previous solution when family tag / branch structure / bounded-mode all
  // match, otherwise a fresh family.initialGuess. The expensive part — the
  // Newton solve plus univalence + identity checks — then runs OFF the main
  // thread via PSW.liveSolve (QD.liveSolveStep), so a slow frame never freezes
  // the tab. Only this seed selection and the final paint stay on the main
  // thread.
  const family = QD.selectFamily(norm);
  let initPhi = null;
  const prev = state.current && state.current.success && state.current.primary
             ? state.current.primary.phi : null;
  if (prev &&
      prev.family === expectedFamilyTag &&
      !!prev.unbounded === unbounded &&
      prev.branches.length === built.poles.length &&
      prev.branches.every((br, j) => br.A.length === built.poles[j].principal.length)) {
    initPhi = QD.clonePhi(prev);
    desc.warmStartUpdate(initPhi, norm);
  } else {
    initPhi = family.initialGuess(built, norm);
  }

  state.altSearchToken++;
  state.altSearchActive = false;
  $('#alt-search-indicator').classList.add('hidden');

  // §23 (mid-drag) regime switch — DETECTION only, never an inline full solve
  // (that synchronous solve every frame was what froze the tab). When the
  // warm-start can't represent the live pole config (Newton fails, or the
  // result is inconsistent with the current regime) we kick the debounced
  // worker solve — solveAndRender carries autoSwitchSingular for the four PQD
  // modes — and still paint the live current-family result so the boundary
  // tracks the drag until the worker catches up on the next pause.
  const isPqdAuto = state.autoSwitchSingular && /^pqd-/.test(state.mode);
  const isSingularMode = state.mode.endsWith('-singular');

  // Live Newton + verification budget. Reduced sample count keeps each frame
  // cheap; the debounced full solve re-verifies at full state.samples.
  const liveOpts = {
    newton: { ...preset.newton, maxIter: 30 },
    numSamples: Math.min(state.samples, LIVE_SAMPLES),
    wantOriginInside: isPqdAuto,
  };

  const myToken = ++_liveSolveToken;
  const PSW = QD.PrimarySolverWorker;
  // Run the live step on the dedicated live worker (cancel-and-replace per
  // frame); fall back to a synchronous main-thread liveSolveStep when no worker
  // is available (file:// origin, unit tests).
  const runLive = (PSW && typeof PSW.liveSolve === 'function')
    ? PSW.liveSolve(built, initPhi, liveOpts)
    : Promise.resolve().then(() => QD.liveSolveStep(built, initPhi, liveOpts));

  runLive.then((result) => {
    if (myToken !== _liveSolveToken) return;       // superseded by a newer frame
    if (!result || !result.success) {
      // Warm-start diverged (often a mid-drag regime crossing). Hand off to the
      // debounced worker solve, which can multistart and auto-switch the regime.
      scheduleSolve();
      return;
    }

    // Is the warm-start result consistent with the current regime? "consistent"
    // = univalent + identity + (origin Ω-membership matches singular/non-singular
    // mode). If not (a clean origin crossing, or the invalid-ansatz case where
    // the R# non-vanishing guard trips), kick the debounced worker solve to
    // switch regime when the drag pauses; we still paint the result below so
    // nothing appears stuck.
    if (isPqdAuto) {
      const consistent = result.univalent && result.identityOK &&
                         result.originInside === isSingularMode;
      if (!consistent) scheduleSolve();
    }

    const sol = {
      success: true,
      phi: result.phi,
      method: 'live',
      univalent: result.univalent,
      identity: result.identity,
      identityOK: result.identityOK,
      residual: result.residual,
      iterations: result.iterations,
    };

    state.current = {
      success: true,
      primary: sol,
      alternates: [],
      hData: built,
      w0Used: norm.w0 || (sol.phi && sol.phi.w0),
      cUsed:  norm.c,
      unbounded,
      attempts: [],
    };
    state.selectedSolutionIdx = 0;
    publishPrimarySolution();

    showSolution(sol, built, /*isPrimary=*/ false);
    refreshAlternatesPanel();
    // Live-refresh the status-panel cards (geometry / cusps / observables) so
    // they track the drag. Throttled + accuracy-deferred (see scheduleLiveAnalysis);
    // the drag-end full solve runs the authoritative pass.
    scheduleLiveAnalysis();
  }).catch(() => { /* superseded / aborted live job — ignore */ });
}

// Token that increments every time solveAndRender() is called. Used to
// discard stale worker results when the user edits faster than solves
// complete.
let _solveAndRenderToken = 0;

// The status-panel analyses (geometric properties, boundary singularities,
// geometry & accuracy) share their OWN token, separate from the solve token, so
// the cheap passes can refresh live during a pole/slider drag (the quick-solve
// path) without being tied to the full-solve cadence. A newer analysis request
// supersedes any older in-flight idle callback.
let _analysisToken = 0;
let _liveAnalysisLast = 0;          // last live-analysis timestamp (throttle)

// solveAndRender — the main (debounced) solve pipeline. Steps:
//   1. buildHData() + buildNormalization() from the DOM; bail on parse errors.
//   2. Bump _solveAndRenderToken (stale-result guard) and run the solve on the
//      warm PrimarySolverWorker (falls back to a sync QD.solveInverseQD in the
//      no-worker / unit-test path). The worker preempts any in-flight solve.
//   3. Auto-escalate: if the standard preset failed and the mode allows it
//      (MODES[x].autoEscalate — LQDs opt out since non-existence is genuine),
//      retry once with the exhaustive preset.
//   4. After EVERY await, re-check myToken !== _solveAndRenderToken and bail if
//      a newer call superseded this one (prevents stale paints).
//   5. Stash the result on state.current (+ hData/w0Used/cUsed/unbounded) and
//      publishPrimarySolution() so the other tabs see it.
//   6. §23: if the solver auto-switched the PQD singular⇄non-singular regime,
//      reflectFamilyMode() updates the compact domain-type control WITHOUT
//      re-solving, and the alternate search uses the switched norm (normFromPhi).
//   7. showSolution → refreshAlternatesPanel → startBackgroundAltSearch →
//      scheduleGeomClassification (the async geometric-univalence card).
// Map a failed solve into one line of plain-language, mode-aware guidance the
// user can act on, prepended to the raw solver `reason:` dump (B7). The intent
// is to answer "what do I change?" rather than just reporting the failure.
function failureGuidance(mode) {
  // Guidance prose lives in QD.Strings.guidance (ui-strings.js); fall back to the
  // inline defaults if the strings module isn't present.
  const G = (typeof window !== 'undefined' && window.QD && window.QD.Strings && window.QD.Strings.guidance) || {};
  const tips = [];
  if (state.aggressiveness !== 'exhaustive') {
    tips.push(G.tryHarder || "try the “Try harder (exhaustive search)” button or raise Aggressiveness");
  }
  if (/^lqd-/.test(mode)) {
    // LQD existence is genuinely bounded (Thm 5.3.2 / 5.6.2): not every h
    // admits a log-weighted QD. Smaller residues / different c can help.
    tips.push(G.lqd || "this h may have no log-weighted QD — try smaller residues, or adjust c");
  } else if (/^pqd-/.test(mode)) {
    // PQD realizability needs the residue magnitude above a threshold
    // (C > (pᵃ − w₀ᵃ)²/α²) and an interior w₀.
    tips.push(G.pqd || "PQDs need a large-enough residue and an interior w₀ — try a bigger |C| or move w₀");
  } else {
    tips.push(G.poles || "move poles away from each other and the boundary, or adjust residue magnitudes");
  }
  return (G.noSolutionPrefix || 'No quadrature domain found. Suggestions: ') + tips.join('; ') + '.';
}

function solveAndRender() {
  const built = buildHData();
  if (!built) {
    setStatus({ kind: 'err', text: 'No poles entered.' });
    return;
  }
  if (built.error) {
    setStatus({ kind: 'err', text: built.error });
    return;
  }

  const norm = buildNormalization(built);
  if (norm.error) {
    setStatus({ kind: 'err', text: norm.error });
    return;
  }

  const preset = PRESETS[state.aggressiveness];
  setStatus({ kind: 'info', text: 'Solving…' });
  showSolveBusy();                          // spinner + Cancel button (B3)
  // writeUrlState lives in the sibling ui-url-state.js module, which installs
  // AFTER this one — read it from `ui` at call time (runtime), not at install.
  ui.writeUrlState();                       // keep the shareable URL in sync (B1)

  state.altSearchToken++;
  state.altSearchActive = false;
  $('#alt-search-indicator').classList.add('hidden');

  const myToken = ++_solveAndRenderToken;
  // The worker preempts any prior in-flight solve when called again; on the
  // main-thread fallback this is a no-op. Either way we re-check the token
  // after each await so a stale completion doesn't overwrite a newer state.
  const PSW = QD.PrimarySolverWorker;

  const runOne = (opts) => {
    if (PSW && typeof PSW.solve === 'function') return PSW.solve(built, opts);
    // Fallback (unit-test / no-worker environment): sync solve wrapped in a
    // microtask so the path is uniformly async.
    return Promise.resolve().then(() => QD.solveInverseQD(built, opts));
  };

  (async () => {
   // Outer try/finally guarantees the busy indicator (spinner + Cancel) is
   // cleared on every exit path — success, failure, solver error, or
   // supersession — but only by the solve that still owns the token, so a
   // newer in-flight solve keeps its own spinner up (B3).
   try {
    let result;
    try {
      const opts = buildSolverOptions(preset, { findAlternates: false });
      applyNormToOpts(opts, norm);
      // §23: in any PQD mode, let the solver auto-detect the singular ↔
      // non-singular transition (boundary crossing the origin) and re-dispatch
      // to the correct family. Scoped to the four PQD modes; gated by the toggle.
      opts.autoSwitchSingular = state.autoSwitchSingular && /^pqd-/.test(state.mode);
      // 1D: warm-start the drag-end full solve from the last good φ (typically
      // the live-drag result). The solver only uses it when the family matches
      // and it converges; otherwise it falls through to the full pipeline, so
      // this is a pure speedup with no behavior change. φ is plain data
      // (clonePhi shape) → structured-clone-safe across the worker postMessage.
      //
      // CRITICAL: inject the CURRENT normalization gauges (c / α / w₀) into the
      // warm seed via the mode's warmStartUpdate hook FIRST. The solver's
      // warm-start trusts the seed's own gauge, so a stale warmPhi would pin the
      // OLD c / α / w₀ and converge in 0 iterations — making the conformal-radius
      // slider, the α input, and manual w₀ appear inert (state + the solver
      // options changed, but the warm seed silently overrode them).
      //
      // Only pass the seed when it is family/structure-compatible with the target
      // (same compatibility test the live quick-solve uses), so a mode switch or
      // a preset change to a different pole structure starts fresh instead of
      // warm-starting from — and corrupting with the gauge update — an
      // incompatible φ.
      const prevPhi = state.current && state.current.success && state.current.primary
                    ? state.current.primary.phi : null;
      if (prevPhi) {
        const wdesc = modeDescriptor();
        const compatible =
          prevPhi.family === wdesc.familyTag &&
          !!prevPhi.unbounded === !!norm.unbounded &&
          Array.isArray(prevPhi.branches) &&
          prevPhi.branches.length === built.poles.length &&
          prevPhi.branches.every((br, j) =>
            br.A.length === built.poles[j].principal.length);
        if (compatible) {
          const warm = QD.clonePhi(prevPhi);
          if (typeof wdesc.warmStartUpdate === 'function') {
            try { wdesc.warmStartUpdate(warm, norm); } catch (e) { /* keep raw seed */ }
          }
          opts.warmPhi = warm;
        }
      }
      result = await runOne(opts);
      if (myToken !== _solveAndRenderToken) return;   // superseded by a newer call

      // Auto-escalation: if standard pipeline failed, re-run with the
      // exhaustive preset before giving up. Toggleable in the search panel.
      //
      // Auto-escalation is per-family: see MODES[X].autoEscalate. LQDs skip
      // it because non-existence is genuine (Theorem 5.3.2 / 5.6.2 bounds).
      if (modeDescriptor().autoEscalate
          && (!result.success || !result.primary ||
              !(result.primary.univalent && result.primary.identityOK))
          && state.searchOptions.autoEscalate
          && state.aggressiveness !== 'exhaustive') {
        // Item 8: surface the phase transition — the first pass didn't find a
        // valid QD, so we're now widening the search.
        const phaseEl = $('#solve-phase');
        if (phaseEl) phaseEl.textContent = 'escalating to exhaustive search…';
        const exh = buildSolverOptions(PRESETS.exhaustive, { findAlternates: false });
        applyNormToOpts(exh, norm);
        const escalated = await runOne(exh);
        if (myToken !== _solveAndRenderToken) return;
        if (escalated.success) result = escalated;
      }
    } catch (e) {
      if (myToken !== _solveAndRenderToken) return;
      if (e && e.aborted) return;                     // user-initiated cancellation
      setStatus({ kind: 'err', text: 'Solver error: ' + (e && e.message || e) });
      return;
    }
    state.current = result;
    state.current.hData = built;
    state.current.w0Used = norm.w0 || (state.current.primary && state.current.primary.phi && state.current.primary.phi.w0);
    state.current.cUsed  = norm.c;
    state.current.unbounded = !!norm.unbounded;
    state.selectedSolutionIdx = 0;
    publishPrimarySolution();

    if (!result.success) {
      setStatus({
        kind: 'err',
        text: failureGuidance(state.mode) + '\n\n' +
              'Technical detail:\n' +
              '  reason: ' + result.error + '\n' +
              '  attempts: ' + (result.attempts ? result.attempts.length : 0),
      });
      // Item 8: highlight the recovery affordance the guidance points at.
      const th = $('#try-harder-btn');
      if (th) th.classList.add('attention');
      plot.clear();
      $('#alternates-card').classList.add('hidden');
      renderRiemannMap(null);   // hide the φ(z) formula in the Domain-type tile
      // Geom/cusp sections live in the status panel; hide the whole panel on a
      // failed solve (no current solution to summarize).
      updateStatusPanelVisibility();
      // Bounded-PQD failures are often genuine NON-REALIZABILITY (the α-branch
      // folds), not solver weakness. Trace the branch off the critical path and
      // replace the generic message with the realizable-α verdict (failure-only,
      // idle, token-guarded — see scheduleRealizabilityDiagnostic).
      scheduleRealizabilityDiagnostic(built, norm, state.mode, myToken);
      // Try-harder button is always visible — nothing to toggle here.
      return;
    }

    // §23: if the solver auto-switched the PQD regime, reflect the actual
    // family in the UI (radio + cards) WITHOUT re-solving — we already have the
    // correct result. The visible radio flip + the status note in showSolution
    // are the user-facing indicators.
    let altNorm = norm;
    if (result.regimeSwitched && result.primary && result.primary.phi) {
      reflectFamilyMode(result.primary.phi.family);
      altNorm = QD.normFromPhi(result.primary.phi) || norm;
    }

    showSolution(result.primary, built, /*isPrimary=*/true);
    refreshAlternatesPanel();

    startBackgroundAltSearch(built, altNorm);
    runStatusAnalyses();   // geometry / cusps / observables (authoritative, with accuracy)
   } finally {
     if (myToken === _solveAndRenderToken) hideSolveBusy();
   }
  })();
}

// ---------- Solve busy-indicator (spinner + Cancel) ----------------------
// Shown for the duration of a primary solve. The Cancel button aborts the
// warm worker mid-solve and bumps the solve token so any late completion is
// treated as superseded (and therefore never paints).
// Elapsed-time ticker shown next to the "solving…" spinner so a slow solve
// reads as in-progress rather than frozen. Cleared whenever the busy row hides.
let _solveElapsedTimer = null;
function showSolveBusy() {
  const row = $('#solve-busy-row');
  if (row) row.classList.remove('hidden');
  // Item 8: reset the phase label and clear any prior failure cue on the
  // "Try harder" button at the start of a fresh solve.
  const phaseEl = $('#solve-phase');
  if (phaseEl) phaseEl.textContent = 'solving…';
  const th = $('#try-harder-btn');
  if (th) th.classList.remove('attention');
  const elapsedEl = $('#solve-elapsed');
  if (elapsedEl) {
    const t0 = Date.now();
    elapsedEl.textContent = '0.0s';
    if (_solveElapsedTimer) clearInterval(_solveElapsedTimer);
    _solveElapsedTimer = setInterval(() => {
      elapsedEl.textContent = ((Date.now() - t0) / 1000).toFixed(1) + 's';
    }, 100);
  }
}
function hideSolveBusy() {
  const row = $('#solve-busy-row');
  if (row) row.classList.add('hidden');
  if (_solveElapsedTimer) { clearInterval(_solveElapsedTimer); _solveElapsedTimer = null; }
}
function cancelSolve() {
  const PSW = QD.PrimarySolverWorker;
  if (PSW && typeof PSW.cancel === 'function') PSW.cancel();
  // Bump the token so the (rejected) in-flight promise is seen as superseded.
  _solveAndRenderToken++;
  // Also stop any background alternate search tied to this solve.
  state.altSearchToken++;
  state.altSearchActive = false;
  $('#alt-search-indicator').classList.add('hidden');
  hideSolveBusy();
  setStatus({ kind: 'warn', text: 'Solve cancelled.' });
}

// §25: classify the solved Ω by special univalence criteria (convex /
// star-like / spiral-like) AFTER the boundary has painted, off the critical
// path. The checks are cheap (≈360 order-2 Taylor evals) but the idle pass
// keeps the solve snappy; the `myToken === _solveAndRenderToken` guard discards
// results from a superseded solve. Result is stashed on the envelope and
// rendered into the dedicated card.

// Run the three status-panel analyses against the CURRENT solution, off the
// critical path. They share _analysisToken (bumped here) so older in-flight
// idle callbacks are discarded — safe to call repeatedly. With opts.live the
// observables pass is cheaper and skips the heavy accuracy estimate.
function runStatusAnalyses(opts) {
  const cur = state.current;
  if (!cur || !cur.success || !cur.primary || !cur.primary.phi) return;
  const token = ++_analysisToken;
  scheduleGeomClassification(cur.primary, token);
  scheduleCuspClassification(cur.primary, token);
  scheduleObservables(cur.primary, cur.hData, token, opts);
  scheduleSymmetry(cur.primary, token);
}

// #9: detect the domain's symmetry group (for the annotated-phenomena overlay's
// dashed axes + D_n / Z_n badge), off the idle path like the cusp classifier.
// Cached on state.current.symmetry, keyed by phi identity.
function scheduleSymmetry(sol, token) {
  if (!sol || !sol.phi || !QD.detectSymmetry) return;
  const phi = sol.phi;
  const cached = state.current && state.current.symmetry;
  if (cached && cached._phiRef === phi) return;     // already computed for this φ
  const run = () => {
    if (token !== _analysisToken) return;           // superseded by a newer analysis
    let sym;
    try { sym = QD.detectSymmetry(phi); } catch (e) { return; }
    if (token !== _analysisToken) return;
    if (sym) sym._phiRef = phi;
    if (state.current) state.current.symmetry = sym;
    // Repaint so the symmetry axes (drawn from state.current.symmetry) appear.
    if (state.showPhenomena && typeof plot !== 'undefined' && plot) plot.render();
  };
  const idle = window.requestIdleCallback || ((fn) => setTimeout(fn, 0));
  idle(run, { timeout: 250 });
}

// Throttled live refresh during a pole/slider drag (called per quick-solve
// frame). Caps the rate so the cheap analyses track the drag without flooding;
// the drag-end full solve runs the authoritative pass via runStatusAnalyses().
const LIVE_ANALYSIS_MS = 120;
function scheduleLiveAnalysis() {
  const now = Date.now();
  if (now - _liveAnalysisLast < LIVE_ANALYSIS_MS) return;
  _liveAnalysisLast = now;
  runStatusAnalyses({ live: true });
}

function scheduleGeomClassification(sol, token) {
  if (!sol || !sol.phi || !QD.classifyUnivalence) {
    renderGeomProps(null);
    return;
  }
  const run = () => {
    if (token !== _analysisToken) return;           // superseded by a newer analysis
    let geom;
    try {
      geom = QD.classifyUnivalence(sol.phi, { samples: state.samples, univalent: sol.univalent });
    } catch (e) {
      renderGeomProps(null);
      return;
    }
    if (token !== _analysisToken) return;
    if (state.current) { state.current.geomProps = geom; publishPrimarySolution(); }
    renderGeomProps(geom);
  };
  const idle = window.requestIdleCallback || ((fn) => setTimeout(fn, 0));
  idle(run, { timeout: 200 });
}

// §25: render the geometric-properties card from a classifyUnivalence result.
function renderGeomProps(geom) {
  const content = $('#sp-geom-content');
  if (!content) return;
  if (!geom) { content.innerHTML = ''; return; }
  const yn = (v) => v ? '<span class="ok">✓ yes</span>' : '<span class="warn">✗ no</span>';
  const fmt = (x) => (typeof formatExp === 'function') ? formatExp(x) : (+x).toExponential(2);
  // Each criterion shows just "label: ✓/✗" inline; the quantitative detail
  // (the margin min Re(·), or the spiral angle) goes into a hover tooltip
  // (title=) so the panel stays uncluttered.
  const row = (label, ynHtml, tip) =>
    `<span class="geom-row" title="${escapeAttr(tip)}"><span class="key">${label}:</span> ${ynHtml}</span>`;
  const lines = [];

  // Convex (bounded only).
  if (geom.convex && geom.convex.na) {
    lines.push(row('convex', '<span class="key">n/a</span>', 'Convexity is reported only for bounded Ω.'));
  } else if (geom.convex && geom.convex.indeterminate) {
    lines.push(row('convex', '<span class="warn">indeterminate</span>', 'φ′ vanishes on ∂𝔻 (cusp / fold) — convexity indeterminate.'));
  } else if (geom.convex) {
    lines.push(row('convex', yn(geom.convex.is), 'min Re(1 + z·φ″/φ′) = ' + fmt(geom.convex.margin) + '   (convex ⇔ > 0)'));
  }

  // Star-like (w.r.t. center, or ∞ for unbounded).
  const starWrt = geom.bounded ? 'w₀' : '∞';
  lines.push(row('star-like (' + starWrt + ')', yn(geom.starLike.is),
    'min Re(z·φ′/(φ−c)) = ' + fmt(geom.starLike.margin) + '   (star-like ⇔ > 0)'));

  // Spiral-like — the optimal spiral angle λ / arc width goes in the tooltip.
  const sp = geom.spiralLike;
  const spTip = sp.is
    ? 'optimal spiral angle λ ≈ ' + sp.angleDeg.toFixed(1) + '°'
    : 'arg-arc width ' + (sp.arcWidth * 180 / Math.PI).toFixed(0) + '° ≥ 180° (not spiral-like)';
  lines.push(row('spiral-like', yn(sp.is), spTip));

  if (geom.notes && geom.notes.length) {
    lines.push(`<span class="key" style="font-style:italic;">${geom.notes.map(escapeHTML).join('<br>')}</span>`);
  }

  content.innerHTML = lines.join('<br>');
}

// Bounded-PQD realizability diagnostic. When a `pqd-bounded*` solve FAILS, the
// requested-α PQD often simply does not exist: with fixed quadrature data the
// univalent solution branch folds as α grows (the |w|^{2(α−1)} weight shrinks
// the realizable region), so "classically (α=1) solvable" does NOT imply the
// target-α PQD is realizable. QD.diagnosePQDRealizability traces the α-branch
// from α≈1 to the target and locates that fold. It is EXPENSIVE (tens of nested
// Newton solves), so it runs ONLY on failure, in an idle pass, then REPLACES the
// generic failure text with the realizable-α verdict. Token-guarded so a newer
// solve discards a stale result.
function scheduleRealizabilityDiagnostic(hData, norm, mode, token) {
  if (!/^pqd-bounded/.test(mode) || !QD.diagnosePQDRealizability) return;
  const alpha = (norm && norm.alpha) || state.alpha;
  if (!(alpha > 0) || Math.abs(alpha - 1) < 1e-6) return;
  const run = () => {
    if (token !== _solveAndRenderToken) return;       // superseded by a newer solve
    let d;
    try {
      d = QD.diagnosePQDRealizability(hData, {
        alpha,
        w0: norm && norm.w0,
        singular: mode === 'pqd-bounded-singular',
      });
    } catch (e) { return; }                            // keep the generic failure message
    if (token !== _solveAndRenderToken) return;
    const aStr = (+alpha).toPrecision(3).replace(/\.?0+$/, '');
    let verdict;
    if (d.reason === 'fold-below-target') {
      const amax = d.alphaMax.toFixed(2);
      const rel = alpha > 1 ? ('α ≲ ' + amax) : ('α ≳ ' + amax);
      verdict = 'No power-weighted (α = ' + aStr + ') quadrature domain exists for this data. '
        + 'The univalent solution branch folds at α ≈ ' + amax + ' — it is realizable only for '
        + rel + '. Lower α, or move/scale the poles.';
    } else if (d.reason === 'invalid-even-classical') {
      verdict = 'This quadrature data does not define a valid domain even classically (α = 1). '
        + 'Check the residue magnitudes and pole positions.';
    } else if (d.reason === 'non-univalent') {
      verdict = 'A solution exists at α = ' + aStr + ' but its boundary self-intersects '
        + '(not univalent), so it is not a valid domain.';
    } else if (d.reason === 'realizable') {
      verdict = 'A valid α = ' + aStr + ' domain appears to exist but the direct solve missed it — '
        + 'click “Try harder (exhaustive search)”.';
    } else {
      return;
    }
    setStatus({ kind: 'err', text: 'Realizability: ' + verdict });
  };
  const idle = window.requestIdleCallback || ((fn) => setTimeout(fn, 0));
  idle(run, { timeout: 400 });
}

// Boundary cusp detection: classify ∂Ω's cusps + their (p,q) type AFTER the
// boundary has painted, off the critical path (mirrors scheduleGeomClassification).
// The classifier (QD.classifyCusps) finds boundary-adjacent zeros of φ′, reads
// the order m and the (p,q) exponents from φ's exact Taylor coefficients, and
// cross-checks the leading exponent numerically. Result is stashed on the
// envelope (for cross-tab readers) and rendered into the card + plot markers.
function scheduleCuspClassification(sol, token) {
  if (!sol || !sol.phi || !QD.classifyCusps) {
    renderCusps(null);
    return;
  }
  const run = () => {
    if (token !== _analysisToken) return;           // superseded by a newer analysis
    let res;
    try {
      res = QD.classifyCusps(sol.phi, { });
    } catch (e) {
      renderCusps(null);
      return;
    }
    if (token !== _analysisToken) return;
    if (state.current) { state.current.cuspProps = res; publishPrimarySolution(); }
    renderCusps(res);
    // Repaint so the cusp markers (drawn from state.current.cuspProps) appear.
    if (typeof plot !== 'undefined' && plot) plot.render();
  };
  const idle = window.requestIdleCallback || ((fn) => setTimeout(fn, 0));
  idle(run, { timeout: 200 });
}

// Render the boundary-singularities section of the status panel from a
// classifyCusps result.
function renderCusps(res) {
  const content = $('#sp-cusps-content');
  if (!content) return;
  if (!res) { content.innerHTML = ''; return; }
  const cusps = res.cusps || [];
  const lines = [];

  if (cusps.length === 0) {
    lines.push(`<span class="ok">✓ smooth boundary</span> <span class="key">(no cusps)</span>`);
    if (res.notes && res.notes.length) {
      lines.push(`<span class="key" style="font-style:italic;">${res.notes.map(escapeHTML).join('<br>')}</span>`);
    }
  } else {
    const nActual = cusps.filter(c => c.isCusp).length;
    lines.push(`<span class="key">${cusps.length} singular point${cusps.length === 1 ? '' : 's'} on ∂Ω (${nActual} cusp${nActual === 1 ? '' : 's'})</span>`);
    for (const c of cusps) {
      // ● filled = an actual boundary cusp; ○ hollow = incipient (φ′-zero near
      // but not on ∂𝔻 — the type it WOULD have at the bifurcation).
      const glyph = c.isCusp ? '<span class="warn">●</span>' : '<span class="key">○</span>';
      const typeStr = `(${c.type[0]},${c.type[1]})`;
      const where = `θ ≈ ${c.thetaDeg.toFixed(1)}°`;
      const detail = c.isCusp
        ? `<span class="key">m=${c.orderM}</span>`
        : `<span class="key">incipient, d=${Math.abs(c.dist).toFixed(3)}</span>`;
      lines.push(`${glyph} <span class="key">${where} — ${escapeHTML(c.typeLabel || (typeStr + ' cusp'))}</span> · ${detail}`);
    }
  }

  content.innerHTML = lines.join('<br>');
}

// Tier-0 foundational observables: boundary geometry (area / perimeter /
// curvature) + an accuracy estimate, computed AFTER the boundary paints, off the
// critical path (mirrors scheduleCuspClassification). Result is stashed on the
// envelope (state.current.observables) so the curvature overlay can read it, and
// rendered into the "Geometry & accuracy" card. `built` is the solved h-data.
function scheduleObservables(sol, hData, token, opts) {
  opts = opts || {};
  if (!sol || !sol.phi || !QD.boundaryObservables) {
    renderObservables(null);
    return;
  }
  const run = () => {
    if (token !== _analysisToken) return;           // superseded by a newer analysis
    let res;
    try {
      // Live (drag) passes use a lighter sweep and SKIP the heavy accuracy
      // estimate (two identity verifies + a Jacobian): the accuracy line carries
      // its last value forward and refreshes on the drag-end full solve.
      const obs = QD.boundaryObservables(sol.phi, { samples: opts.live ? 512 : 1024 });
      obs._phiRef = sol.phi;              // lets the curvature overlay reuse this compute
      let acc;
      if (opts.live) {
        acc = (state.current && state.current.observables) ? state.current.observables.acc : null;
      } else {
        acc = (hData && QD.estimateAccuracy) ? QD.estimateAccuracy(sol.phi, hData, {}) : null;
      }
      res = { obs, acc };
    } catch (e) {
      renderObservables(null);
      return;
    }
    if (token !== _analysisToken) return;
    if (state.current) { state.current.observables = res; publishPrimarySolution(); }
    renderObservables(res);
    // Repaint so the curvature overlay (drawn from state.current.observables) appears.
    if (typeof plot !== 'undefined' && plot) plot.render();
  };
  const idle = window.requestIdleCallback || ((fn) => setTimeout(fn, 0));
  idle(run, { timeout: 250 });
}

// Render the "Geometry & accuracy" card from a scheduleObservables result.
function renderObservables(res) {
  const content = $('#sp-observables-content');
  if (!content) return;
  if (!res || !res.obs) { content.innerHTML = ''; return; }
  const obs = res.obs, acc = res.acc;
  const num = (x) => {
    if (x == null || !isFinite(x)) return '—';
    const a = Math.abs(x);
    if (x !== 0 && (a >= 1e4 || a < 1e-3)) return (+x).toExponential(3);
    return (+x).toPrecision(4);
  };
  const lines = [];
  lines.push(`<span class="geom-row"><span class="key">area:</span> ${num(obs.area)}</span>`);
  lines.push(`<span class="geom-row"><span class="key">perimeter:</span> ${num(obs.perimeter)}</span>`);
  const cuspDeg = (obs.argMaxCurvatureTheta * 180 / Math.PI).toFixed(1);
  lines.push(`<span class="geom-row" title="max |κ| along ∂Ω at θ ≈ ${cuspDeg}° (κ → ∞ at a cusp)"><span class="key">max curvature:</span> ${num(obs.maxCurvature)}</span>`);
  if (acc && acc.significantDigits != null) {
    const tip = 'identity maxRelDiff  N=' + (acc.relN != null ? (+acc.relN).toExponential(2) : '—')
      + ', 2N=' + (acc.rel2N != null ? (+acc.rel2N).toExponential(2) : '—')
      + (acc.conditionEst != null ? ('   cond ≈ ' + (+acc.conditionEst).toExponential(2)) : '');
    const warn = acc.underResolved ? ' <span class="warn">(under-resolved)</span>' : '';
    lines.push(`<span class="geom-row" title="${escapeAttr(tip)}"><span class="key">accuracy:</span> ≈ ${acc.significantDigits.toFixed(1)} sig. digits${warn}</span>`);
  }
  // Near-cusp honesty (#11): when a φ′ zero approaches |z|=1 the identity check
  // is unreliable (the hole thins) and the geometric criterion governs validity.
  if (acc && acc.nearCusp) {
    const dist = acc.cuspDistance != null ? (+acc.cuspDistance).toExponential(2) : '—';
    const tip = 'A φ′ zero is within ' + dist + ' of |z|=1 — a forming cusp. The '
      + 'quadrature-identity check is unreliable here (the hole thins so interior '
      + 'test points cannot clear ∂Ω); validity is governed by the geometric '
      + 'criterion (univalence + critical modulus).';
    lines.push(`<span class="geom-row" title="${escapeAttr(tip)}"><span class="key">near cusp:</span> <span class="warn">dist ≈ ${dist}</span> — trusting geometry</span>`);
  }
  content.innerHTML = lines.join('<br>');
}

// ---------- Status panel (overlaid bottom-right of the plot) -------------
// The verdict badge + show/hide logic for the transparent #status-panel that
// hosts the geometric properties + boundary singularities. (The Riemann map
// φ(z) lives separately, at the bottom of the Domain-type tile — renderRiemannMap.)

// Map a solved solution to the valid/invalid verdict shown in the panel badge.
function qdValidityBadge(sol) {
  if (!sol)                              return { cls: 'err',  text: '✗ No solution' };
  if (sol.univalent && sol.identityOK)   return { cls: 'ok',   text: '✓ Valid quadrature domain' };
  if (!sol.univalent && !sol.identityOK) return { cls: 'err',  text: '✗ Spurious root (non-univalent + identity fails)' };
  if (!sol.univalent)                    return { cls: 'warn', text: '⚠ Boundary self-intersects (non-univalent)' };
  return { cls: 'warn', text: '⚠ Quadrature identity not satisfied' };
}
function renderValidityBadge(sol) {
  const el = $('#sp-badge');
  if (!el) return;
  const b = qdValidityBadge(sol);
  el.innerHTML = `<span class="${b.cls}">${escapeHTML(b.text)}</span>`;
}

// Show the panel only on the QD tab's inverse view, once there is a current
// solution. Hidden on the Schwarz/param-slice tabs and on the direct view.
function updateStatusPanelVisibility() {
  const panel = $('#status-panel');
  if (!panel) return;
  const tabBtn = document.querySelector('.tab-btn.active');
  const onQdTab = !tabBtn || tabBtn.dataset.tab === 'qd';   // QD is the default tab
  const inverseView = (state.viewMode || 'inverse') === 'inverse';
  const hasSol = !!(state.current && state.current.success);
  panel.classList.toggle('hidden', !(onQdTab && inverseView && hasSol));
}

// §23: silently switch the UI to the mode matching a solved phi's family
// (used after the solver auto-switches the PQD regime). Updates state.mode and
// the mode visuals (which re-sync the compact domain-type control) — but does
// NOT trigger a solve.
const FAMILY_TO_MODE = {
  powerQD:                 'pqd-bounded',
  powerQD_singular:        'pqd-bounded-singular',
  unboundedPQD:            'pqd-unbounded',
  unboundedPQD_singular:   'pqd-unbounded-singular',
};
function reflectFamilyMode(family) {
  const target = FAMILY_TO_MODE[family];
  if (!target || target === state.mode || !MODES[target]) return;
  state.mode = target;
  applyModeVisuals();   // also re-syncs the compact domain-type control
}

// ---------- Display a chosen solution on the plot ------------------------
// sol: the solution to draw (primary or a previewed alternate). isPrimary:
// true only for the primary solve — alternates being previewed pass false so
// they don't reframe the viewport. Auto-fit happens iff (state.autoFit &&
// isPrimary); it is NOT an "autoFit" flag despite some historical call sites.
function showSolution(sol, hData, isPrimary) {
  const boundary = QD.sampleBoundaryAdaptive(sol.phi, state.samples, Math.floor(state.samples * 1.5));
  const boundaryPts = boundary.map(p => p.w);
  const poles = hData.poles.map(p => p.a);

  plot.setData({
    boundaryPts,
    poles,
    w0: sol.phi.unbounded ? null : sol.phi.w0,
    univalent: !!sol.univalent,
    unbounded: !!sol.phi.unbounded,
    hData,
    phi: sol.phi,           // singular-LQD vector field reads q from here
  });

  if (state.autoFit && isPrimary) plot.fit();

  renderRiemannMap(sol.phi);

  // Build status
  const lines = [];
  // §23: if the solver auto-switched the singular ↔ non-singular regime
  // (boundary crossed the origin), lead with a clear indicator. Transient —
  // it clears on the next solve; the mode radio also visibly flips.
  if (state.current && state.current.regimeSwitched) {
    const toSing = state.current.switchedTo === 'singular';
    lines.push(`<span class="warn">⇄ Auto-switched to the ${toSing ? 'singular' : 'non-singular'} regime — the boundary crossed the origin (0 ${toSing ? '∈' : '∉'} Ω)</span>`);
  }
  // The valid/invalid verdict now lives in the status-panel badge (overlaid on
  // the plot, bottom-right); #status keeps the operational detail below.
  renderValidityBadge(sol);
  updateStatusPanelVisibility();
  lines.push(`<span class="key">method:</span> ${escapeHTML(sol.method || '?')}`);
  if (typeof sol.iterations === 'number') {
    lines.push(`<span class="key">Newton iterations:</span> ${sol.iterations}`);
  }
  if (sol.trace) {
    lines.push(`<span class="key">continuation steps:</span> ${sol.trace.length}`);
  }
  lines.push(`<span class="key">Newton residual:</span> ${formatExp(sol.residual)}`);
  lines.push(`<span class="key">degree of φ:</span> ${sol.phi.branches.reduce((a, b) => a + b.A.length, 0)}`);
  if (sol.identity) {
    const v = sol.identity;
    const cls = sol.identityOK ? 'ok' : 'err';
    // Test-function class: per-family verifier sets one of:
    //   v.unbounded     → 1/(w−b)^k for b ∈ K
    //   v.lqdSingular   → monomials w^k vanishing at 0 (k ≥ 1)
    //   default         → monomials w^k including k = 0
    const testClass = describeTestClass(v);
    lines.push(`<span class="key">identity check:</span> <span class="${cls}">max rel diff = ${formatExp(v.maxRelDiff)}</span>` +
               ` <span class="key">(${testClass})</span>`);
  }
  setStatus({ kind: 'raw', html: lines.join('<br>') });

  // Try-harder button is always visible; no per-solution toggle needed.
}

// Build the human-readable test-function-class string from a verifier result.
// Lives here (not on Family) because it's a UI display concern; the verifier
// flags are the source of truth.
function describeTestClass(v) {
  if (v.lqdUnboundedSingular) {
    const nb = v.testPoints ? v.testPoints.length : 0;
    return `w/(w − b)^k for k = 2…${v.maxDeg} at ${nb} test point${nb === 1 ? '' : 's'} in K (vanishing at 0 and ∞)`;
  }
  if (v.lqdUnbounded) {
    return `1/w, 1/w², …, 1/w^${v.maxDeg} (vanishing at ∞; required by L¹(ρ₀))`;
  }
  if (v.unbounded) {
    const nb = v.testPoints ? v.testPoints.length : 1;
    return `1/(w − b)^k for k = 1…${v.maxDeg} at ${nb} test point${nb === 1 ? '' : 's'} in K`;
  }
  if (v.lqdSingular) {
    return `monomials w¹…w^${v.maxDeg} (vanishing at 0; required by L¹(ρ₀))`;
  }
  return `monomials w⁰…w^${v.maxDeg}`;
}


// ---------- Riemann-map formula card ------------------------------------
// Renders (1) symbolic identity, (2) closed-form expression with values
// substituted, (3) parameters table. The LaTeX itself is generated by the
// pure, Node-testable QD.RiemannLatex.build(phi) (riemann-latex.js); this
// function only does the DOM rendering.
// Renders the solved Riemann map into the bottom of the Domain-type tile
// (#dm-riemann). The explicit NUMERICAL closed form is shown inline
// (#dm-riemann-numer); the SYMBOLIC identity is rendered into a sibling div
// (#dm-riemann-sym) that the "?" toggle reveals — it isn't shown by default.
// Hidden entirely until a solve produces a φ.
function renderRiemannMap(phi) {
  const box   = $('#dm-riemann');
  const numer = $('#dm-riemann-numer');
  const sym   = $('#dm-riemann-sym');
  if (!box || !numer) return;
  if (!phi) { box.classList.add('hidden'); return; }

  // Pure LaTeX generation lives in riemann-latex.js (Node-testable).
  const { symbolic, numeric } = QD.RiemannLatex.build(phi);
  numer.innerHTML = '';
  renderKatex(numer, numeric, true);
  if (sym) { sym.innerHTML = ''; renderKatex(sym, symbolic, true); }  // hidden until "?"
  box.classList.remove('hidden');
}

// Render LaTeX `expr` into the given element. Uses KaTeX if available;
// falls back to a plain-text placeholder if the CDN failed to load.
function renderKatex(el, expr, display) {
  if (typeof katex === 'undefined') {
    el.textContent = expr;
    return;
  }
  try {
    katex.render(expr, el, { displayMode: !!display, throwOnError: false });
  } catch (e) {
    el.textContent = expr;
  }
}

// ---------- Alternates panel ---------------------------------------------
function refreshAlternatesPanel() {
  const card = $('#alternates-card');
  const list = $('#alternates-list');
  list.innerHTML = '';

  const all = state.current.success
    ? [state.current.primary, ...(state.current.alternates || [])]
    : [];

  // Show the card whenever we have alternates OR a background search is
  // running, so the "searching…" spinner is visible even before any alt is
  // found. Otherwise hide it.
  if (all.length <= 1 && !state.altSearchActive) {
    card.classList.add('hidden');
    return;
  }
  card.classList.remove('hidden');

  if (all.length <= 1) {
    const note = document.createElement('div');
    note.style.cssText = 'font-size: 11px; color: #777;';
    note.textContent = 'No alternates found yet…';
    list.appendChild(note);
    return;
  }

  all.forEach((sol, i) => {
    const isSel = i === state.selectedSolutionIdx;
    const row = document.createElement('div');
    row.className = 'alt';
    const tag = i === 0 ? 'Primary' : `Alt ${i}`;
    const valid = sol.univalent && sol.identityOK;
    const flag = valid ? '✓' : (sol.univalent && !sol.identityOK ? '?' : '⚠');
    const desc = valid ? 'valid QD'
                       : (!sol.univalent ? 'non-univalent' : 'identity fails');
    row.innerHTML = `
      <span>
        <strong>${tag}</strong>
        <span style="color:#777"> · ${flag} ${desc}</span>
        <span style="color:#777"> · id ${formatExp(sol.identity ? sol.identity.maxRelDiff : null)}</span>
      </span>
      <button class="small ${isSel ? 'primary' : ''}" data-alt-idx="${i}">${isSel ? 'shown' : 'view'}</button>
    `;
    list.appendChild(row);
  });
}

function viewSolutionByIndex(i) {
  if (!state.current || !state.current.success) return;
  const all = [state.current.primary, ...(state.current.alternates || [])];
  if (i < 0 || i >= all.length) return;
  state.selectedSolutionIdx = i;
  showSolution(all[i], state.current.hData, /*isPrimary=*/i === 0);
  refreshAlternatesPanel();
}

// ---------- Background alternate search ---------------------------------
// Runs QD.searchAlternates on the dedicated AUX WORKER (A3) so each pass is
// off the main thread. Previously every chunk ran synchronously via
// setTimeout, janking the 2D plot. Each loop iteration awaits one worker pass,
// applies the acceptance filter, then yields briefly; the `myToken` guard
// stops the loop (and discards a late worker result) once a newer solve or
// search supersedes this one. Falls back to a main-thread microtask when the
// worker is unavailable (file:// origin / no Worker support).
function startBackgroundAltSearch(hData, norm) {
  const preset = PRESETS[state.aggressiveness];
  const so = state.searchOptions;
  const myToken = ++state.altSearchToken;
  state.altSearchActive = true;
  $('#alt-search-indicator').classList.remove('hidden');
  refreshAlternatesPanel();

  const bgChunks   = so.bgChunks   ?? preset.bgAltChunks;
  const keepGoing  = so.keepSearching;
  let chunk = 0;
  // Seed = user override if any, else time-based.
  let seed = so.seed !== null
    ? (so.seed >>> 0)
    : ((Date.now() ^ 0x9E3779B1) >>> 0);

  const PSW = QD.PrimarySolverWorker;
  const runChunk = (known, opts) =>
    (PSW && typeof PSW.searchAlternates === 'function')
      ? PSW.searchAlternates(hData, norm, known, opts)
      : Promise.resolve().then(() => QD.searchAlternates(hData, norm, known, opts));

  const stop = () => {
    if (myToken !== state.altSearchToken) return;   // a newer search owns the UI
    state.altSearchActive = false;
    $('#alt-search-indicator').classList.add('hidden');
    refreshAlternatesPanel();
  };

  (async () => {
    for (;;) {
      if (myToken !== state.altSearchToken) return;             // superseded
      if (!state.current || !state.current.success) { stop(); return; }
      if (!keepGoing && chunk >= bgChunks)          { stop(); return; }
      chunk++;

      let found = [];
      try {
        const known = [state.current.primary, ...(state.current.alternates || [])];
        found = await runChunk(known, buildAltSearchOptions(preset, seed));
      } catch (e) {
        if (e && e.aborted) return;                             // superseded mid-flight
        console.warn('alt search error:', e);
      }
      if (myToken !== state.altSearchToken) return;             // superseded while awaiting
      seed = (seed * 1664525 + 1013904223) >>> 0;

      if (found && found.length > 0) {
        // Acceptance criteria — by default, only valid QDs are shown. Toggle
        // overrides in the panel let the user surface partial / spurious
        // candidates for diagnostic purposes.
        const accept = found.filter(s => {
          if (s.univalent && s.identityOK) return true;
          if (so.showNonUnivalent && !s.univalent) return true;
          if (so.showIdFailing    && s.univalent && !s.identityOK) return true;
          return false;
        });
        if (accept.length > 0) {
          state.current.alternates = (state.current.alternates || []).concat(accept);
          publishPrimarySolution();
          refreshAlternatesPanel();
        }
      }

      await new Promise(r => setTimeout(r, 30));   // gentle yield between passes
    }
  })();
}



    return {
      scheduleSolve, scheduleQuickSolve, solveAndRender, cancelSolve,
      showSolution, refreshAlternatesPanel, viewSolutionByIndex,
      startBackgroundAltSearch, updateStatusPanelVisibility,
    };
  };
})(typeof window !== 'undefined' ? window : globalThis);
