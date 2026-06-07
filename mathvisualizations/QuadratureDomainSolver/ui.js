// =============================================================================
// ui.js -- Frontend for the Quadrature Domain Solver
//
// File layout (search the banner comments / function names below to jump;
// approximate line numbers in brackets drift — names are the stable anchors):
//
//   Preset library            [in ui-presets.js] quadrature-function presets
//   MODE DESCRIPTORS / MODES   [~60] single source of truth for all 10 modes:
//                              family tag, visible cards, preset list, norm
//                              build, vector-field label, auto-escalate flag
//   Aggressiveness presets     [~358] numRestarts / Newton tol / cont. steps
//   State                      [in ui-state.js] the single source of UI truth
//   Search options             [~424] advanced panel: per-phase toggles + overrides
//   Helpers                    [~536] $, debounce, fmtArg, subscripts
//   buildHData / buildW0        [~551] read DOM → solver input
//   Polar slider bookkeeping   [~608]
//   renderPolesList /          [~644] build the pole-control cards
//     renderPolyCoefList              and the poly-coefficient cards
//   quickSolveAndRender        [~798] snappy warm-start path for slider drags
//   solveAndRender             [~874] MAIN solve pipeline (debounced) — see its
//                              own header for the step-by-step contract
//   scheduleGeomClassification / [~992] async geometric-univalence card (§25):
//     renderGeomProps                 convex / star-like / spiral-like status
//   FAMILY_TO_MODE /           [~1055] §23 auto-switch reflector: mirror the
//     reflectFamilyMode               solver-chosen PQD regime into the compact
//                                     domain-type control (via applyModeVisuals)
//   showSolution               [~1071] push a solution to canvas + Riemann display
//   renderRiemannMap           [~1215] DOM render of φ(z) into the Domain-type tile (the
//                              LaTeX itself comes from QD.RiemannLatex.build in
//                              riemann-latex.js — pure + Node-testable)
//   refreshAlternatesPanel     build the alternates UI
//   startBackgroundAltSearch   [~1598] chunked async alternate-hunt after solve
//   class DomainPlot           [~1664] canvas plot: axes, boundary, poles, field
//   DOM wiring                 [~1678] all event listeners
//   applyModeVisuals / setMode [~1946] mode chrome (descriptor-driven) + switch
//   View-mode toggle           [~2095] inverse | direct (HANDOFF #30)
//   mountQolHelp               [~2149] the "?" help popovers on each card
//   Search-options panel wiring[~2436]
//   Direct-view hooks          [~2556] QD.Direct._sendHToInverseTab cross-load
//   Param-slice cross-tab hooks[~2669] QD_UI.snapshotScenario / loadScenarioIntoQdTab
//
// All math is delegated to QD (= window.QD from solver.js); this file only
// translates between DOM state and QD calls.
// =============================================================================
'use strict';

// NB: `complex.js` already declares `Complex` in the shared script lexical
// scope, so we must NOT redeclare it here. Reference it via QD.Complex when
// needed, or via the unaliased global from complex.js.

// ===========================================================================
// Preset library lives in ui-presets.js (P0 split, §1.1).
// ---------------------------------------------------------------------------
// That file attaches each list to a top-level global keeping the original
// const name (QD_PRESETS_BOUNDED, QD_PRESETS_UNBOUNDED, LQD_PRESETS_BOUNDED,
// LQD_PRESETS_BOUNDED_SINGULAR, LQD_PRESETS_UNBOUNDED,
// LQD_PRESETS_UNBOUNDED_SINGULAR) AND mirrors them under window.QD_UI.Presets.*
// for any future namespaced reader. Bare identifiers below resolve to those
// globals — no other ui.js change required by the move.
// ===========================================================================


// ===========================================================================
// Phase-3 UI modularization (item E) — shared injection context.
// ---------------------------------------------------------------------------
// `uiCtx` carries the ui.js closures that extracted modules need (state, the
// descriptor tables, DOM helpers, peer functions). Each extracted module
// exposes a QD_UI.installX(uiCtx) factory; we capture the returns into local
// bindings with their ORIGINAL names so every call site below is unchanged.
// Install points are chosen so each module dependency exists on uiCtx by the
// time its functions actually run. See ui-modes.js / ui-url-state.js for the
// template and ARCHITECTURE.md "factory-injection" for the rationale.
// ===========================================================================
const uiCtx = { state };
// Forward bindings for Phase-3 extracted Inverse-tab modules (item E). Assigned
// by the install calls near the end of this file; referenced by name throughout.
let renderPolesList, renderPolyCoefList;
let modeAllowsPoly, refreshHText, setHTextMsg, parseAndApplyHText;
let scheduleSolve, scheduleQuickSolve, solveAndRender, cancelSolve;
let showSolution, refreshAlternatesPanel, viewSolutionByIndex;
let startBackgroundAltSearch, updateStatusPanelVisibility;

// MODE DESCRIPTORS (R5) + aggressiveness PRESETS live in ui-modes.js, installed
// here (early) so MODES / modeDescriptor / currentPresetList / PRESETS resolve
// for every consumer below. buildNorm reads ui.buildW0 at solve time, so
// uiCtx.buildW0 only needs to be set before the first solve (done after buildW0).
const { MODES, PRESETS, modeDescriptor, currentPresetList } =
  window.QD_UI.installModes(uiCtx);
Object.assign(uiCtx, { MODES, PRESETS, modeDescriptor, currentPresetList });

// ---------- State --------------------------------------------------------
// The `state` object lives in ui-state.js (A1 split). It's still a
// cross-script-realm `const` so the bare identifier resolves here, and
// is also exposed via window.QD_UI.state for explicit-namespace callers.
// See ui-state.js for the field-by-field schema.

// ===========================================================================
// Search options (advanced panel)
// ---------------------------------------------------------------------------
// The panel exposes per-phase toggles and numeric overrides for everything
// in the aggressiveness preset, plus a few solver knobs the preset doesn't
// touch. Blank numeric fields fall back to the preset.
// ===========================================================================

// Read a number from a text/number input. Returns null when blank / NaN
// so the caller knows to fall back to the preset.
function readNumOrNull(sel) {
  const v = $(sel).value.trim();
  if (v === '') return null;
  const n = Number(v);
  return isFinite(n) ? n : null;
}

// Read the entire search-options DOM into state.searchOptions.
function readSearchOptions() {
  const so = state.searchOptions;
  so.phases.direct       = $('#so-phase-direct').checked;
  so.phases.continuation = $('#so-phase-continuation').checked;
  so.phases.multistart   = $('#so-phase-multistart').checked;
  so.phases.diverse      = $('#so-phase-diverse').checked;
  so.phases.deflation    = $('#so-phase-deflation').checked;

  so.numRestarts   = readNumOrNull('#so-num-restarts');
  so.numDiverse    = readNumOrNull('#so-num-diverse');
  so.numDeflation  = readNumOrNull('#so-num-deflation');
  so.bgChunks      = readNumOrNull('#so-bg-chunks');
  so.bgChunkSize   = readNumOrNull('#so-bg-chunk-size');
  so.keepSearching = $('#so-keep-searching').checked;

  so.newtonMaxIter = readNumOrNull('#so-newton-maxiter');
  so.newtonTol     = readNumOrNull('#so-newton-tol');
  so.contTStart    = readNumOrNull('#so-cont-tstart');
  so.contGrow      = readNumOrNull('#so-cont-grow');

  so.deflationAlpha   = readNumOrNull('#so-defl-alpha');
  so.deflationP       = readNumOrNull('#so-defl-p');
  so.deflateFromValid = $('#so-defl-from-valid').checked;

  so.univalenceSamples = readNumOrNull('#so-uni-samples');
  so.identityTol       = readNumOrNull('#so-id-tol');
  so.showNonUnivalent  = $('#so-show-non-univalent').checked;
  so.showIdFailing     = $('#so-show-id-failing').checked;
  so.autoEscalate      = $('#so-auto-escalate').checked;

  so.seed = readNumOrNull('#so-seed');
}

// Clear every override field; checkboxes return to default state.
function resetSearchOptions() {
  ['#so-num-restarts', '#so-num-diverse', '#so-num-deflation',
   '#so-bg-chunks', '#so-bg-chunk-size',
   '#so-newton-maxiter', '#so-newton-tol', '#so-cont-tstart', '#so-cont-grow',
   '#so-defl-alpha', '#so-defl-p',
   '#so-uni-samples', '#so-id-tol', '#so-seed'].forEach(s => { $(s).value = ''; });
  ['#so-phase-direct', '#so-phase-continuation', '#so-phase-multistart',
   '#so-phase-diverse', '#so-phase-deflation', '#so-auto-escalate'
  ].forEach(s => { $(s).checked = true; });
  ['#so-keep-searching', '#so-defl-from-valid',
   '#so-show-non-univalent', '#so-show-id-failing'
  ].forEach(s => { $(s).checked = false; });
  readSearchOptions();
}

// Build the option-bag passed to QD.solveInverseQD. Layers overrides on top
// of an aggressiveness preset.
function buildSolverOptions(preset, { findAlternates = false } = {}) {
  const so = state.searchOptions;
  const opts = {
    numRestarts:      so.numRestarts   ?? preset.numRestarts,
    newton: {
      maxIter:   so.newtonMaxIter ?? preset.newton.maxIter,
      tolerance: so.newtonTol     ?? preset.newton.tolerance,
    },
    continuation: {
      tStart:     so.contTStart ?? preset.continuation.tStart,
      growFactor: so.contGrow   ?? preset.continuation.growFactor,
    },
    univalenceSamples: so.univalenceSamples ?? state.samples,
    identityTol:       so.identityTol       ?? 1e-6,
    findAlternates,
    usePhases:         { ...so.phases },
    deflationAlpha:    so.deflationAlpha ?? 1,
    deflationP:        so.deflationP     ?? 2,
    deflateFromValid:  so.deflateFromValid,
  };
  if (so.numDiverse   !== null) opts.numDiverseSeeds   = so.numDiverse;
  if (so.numDeflation !== null) opts.numDeflationSeeds = so.numDeflation;
  return opts;
}

// Build the option-bag passed to QD.searchAlternates.
function buildAltSearchOptions(preset, seed) {
  const so = state.searchOptions;
  return {
    numRestarts:       so.bgChunkSize ?? preset.bgAltChunkSize,
    seed,
    newton: {
      maxIter:   so.newtonMaxIter ?? preset.newton.maxIter,
      tolerance: so.newtonTol     ?? preset.newton.tolerance,
    },
    univalenceSamples: so.univalenceSamples ?? state.samples,
    identityTol:       so.identityTol       ?? 1e-6,
    deflateFromKnown:  true,
    deflationAlpha:    so.deflationAlpha ?? 1,
    deflationP:        so.deflationP     ?? 2,
  };
}

// ---------- Helpers ------------------------------------------------------
const subs = ['₀','₁','₂','₃','₄','₅','₆','₇','₈','₉'];
const sub  = n => String(n).split('').map(d => subs[+d] || d).join('');

function $(sel, parent = document) { return parent.querySelector(sel); }
function $$(sel, parent = document) { return Array.from(parent.querySelectorAll(sel)); }

function debounce(fn, ms) {
  let t = null;
  return (...args) => {
    if (t) clearTimeout(t);
    t = setTimeout(() => { t = null; fn(...args); }, ms);
  };
}

// ---------- Build hData from state --------------------------------------
function buildHData() {
  const poles = [];
  for (let j = 0; j < state.poles.length; j++) {
    const p = state.poles[j];
    const a = QD.Complex.parse(p.a);
    if (!a) return { error: `Pole ${j+1}: invalid value for a` };
    const principal = [];
    for (let s = 0; s < p.order; s++) {
      const c = QD.Complex.parse(p.residues[s] || '0');
      if (!c) return { error: `Pole ${j+1}: invalid value for C${sub(s+1)}` };
      principal.push(c);
    }
    poles.push({ a, principal });
  }
  // Polynomial part of h. Allowed in any "unbounded-like" mode where the
  // panel is meaningful (classical unbounded + both unbounded-LQD variants).
  // NB: the unbounded-LQD solvers currently DON'T match polynomial-h in
  // their (★) system — see solver-uqd-lqd.js header. The validation below
  // (in solveAndRender via buildNormalization) surfaces a clear error to
  // the user when they attempt to solve with a nonzero LQD polynomial part.
  const polyPart = [];
  if (modeAllowsPoly(state.mode) && state.polyDegree >= 0) {
    for (let l = 0; l <= state.polyDegree; l++) {
      const c = QD.Complex.parse(state.polyCoeffs[l] ?? '0');
      if (!c) return { error: `Poly coef C∞,${l}: invalid value` };
      polyPart.push(c);
    }
  }
  if (poles.length === 0 && polyPart.length === 0) return null;
  return { poles, polyPart };
}

function buildW0(hData) {
  if (state.w0Mode === 'manual') {
    const w0 = QD.Complex.parse(state.w0Manual);
    if (!w0) return { error: 'Invalid value for φ(0)' };
    return { w0 };
  }
  // centroid
  let sumRe = 0, sumIm = 0;
  for (const p of hData.poles) { sumRe += p.a.re; sumIm += p.a.im; }
  return { w0: { re: sumRe / hData.poles.length, im: sumIm / hData.poles.length } };
}
uiCtx.buildW0 = buildW0;  // available before the first solve (ui-modes buildNorm uses it)


// Copy the normalization signal from `norm` into a solver-options object,
// preserving the (unbounded, c) | (w0) | (lqd, w0) distinction. Used at every
// call site that hands off to QD.solveInverseQD.
function applyNormToOpts(opts, norm) {
  modeDescriptor().applyNorm(opts, norm);
  return opts;
}

function buildNormalization(hData) {
  return modeDescriptor().buildNorm(hData, state);
}

// ---------------------------------------------------------------------------
// PrimarySolution publish helper: forward state.current to QD.PrimarySolution
// so cross-tab readers (Schwarz, Sphere, Param-slice) can subscribe instead
// of reaching into state.current directly. ui.js remains the canonical
// writer; this is a published mirror.
// ---------------------------------------------------------------------------
function publishPrimarySolution() {
  if (!QD.PrimarySolution) return;          // tolerate missing module (test contexts)
  QD.PrimarySolution.publish(state.current);
}

// ===========================================================================
// Polar (magnitude / argument) helpers and slider-range bookkeeping
// ===========================================================================

// Slider max for |C| is cached per residue so it persists across re-renders
// and grows automatically when the user types a larger value.
const magSliderMax = {};

function residueKey(poleIdx, s) {
  return `pole-${poleIdx}-residue-${s}`;
}

function magMaxFor(key, mag) {
  // Default starting max is 5; grow when needed; never shrink automatically.
  const current = magSliderMax[key] ?? 5;
  const target = mag > current * 0.95 ? Math.max(current, Math.ceil(mag * 1.5)) : current;
  magSliderMax[key] = target;
  return target;
}

// Format an argument value (radians) as a multiple of π for the slider readout.
function fmtArg(arg) {
  return (arg / Math.PI).toFixed(3) + 'π';
}

// ---------- Pole / poly-coef grid renderers -> ui-pole-grid.js ----------
// renderPolesList / renderPolyCoefList are installed via QD_UI.installPoleGrid
// near the end of this file. The shared slider helpers they use (residueKey /
// magMaxFor / fmtArg / magSliderMax) + escapeHTML / escapeAttr stay below.

function escapeAttr(s) { return String(s).replace(/&/g, '&amp;').replace(/"/g, '&quot;'); }

// ---------- Read controls into state -------------------------------------
// (readPolesFromDOM + markInvalid removed in D8 — both became dead after the
// structured-grid refactor; pole data is now read directly via the grid
// model rather than walking the DOM, and parse-error highlighting is owned
// by #h-text-msg instead of a per-input badge.)

// ---------- Adding / removing poles --------------------------------------
function addPole() {
  state.poles.push({ a: '0', order: 1, residues: ['1'] });
  renderPolesList();
  scheduleSolve();
}
// Add a simple pole (order 1, coefficient 1) at a specific w-plane point.
// Used by the plot's double-click gesture (plot.onAddPole). Mirrors addPole()
// but stamps the clicked position and marks the config as custom (a placed
// pole is a user edit, like a drag).
function addPoleAt(w) {
  state.poles.push({ a: QD.Complex.toString(w, 4), order: 1, residues: ['1'] });
  markAsCustom();
  renderPolesList();
  scheduleSolve();
}
function removePoleAt(idx) {
  state.poles.splice(idx, 1);
  // No defensive "insert a default if empty" — an empty pole list is a valid
  // user intent (e.g. h = polynomial-only in unbounded mode). The solver
  // surfaces "no poles entered" if appropriate. Without this, clicking × on
  // the last remaining pole silently re-inserts a default and looks like
  // "delete did nothing while h-text picked up 1/w".
  renderPolesList();
  scheduleSolve();
}

// ---------- Solve / render / analysis pipeline -> ui-solve.js ------------
// scheduleSolve / scheduleQuickSolve / solveAndRender / cancelSolve /
// showSolution / refreshAlternatesPanel / viewSolutionByIndex /
// startBackgroundAltSearch / updateStatusPanelVisibility (+ the geom/cusp/
// realizability analysis + status-panel badge + Riemann-map rendering) are
// installed via QD_UI.installSolve near the end of this file. The small shared
// helpers escapeHTML / formatExp / setStatus stay here (used by the pipeline
// AND the retained Try-harder / cross-tab handlers).

// escapeHTML: delegates to QD.QoL.escapeHTML (HANDOFF #35 consolidation).
// Falls back to a local impl if qol.js failed to load.
function escapeHTML(s) {
  return (window.QD && window.QD.QoL && window.QD.QoL.escapeHTML)
    ? window.QD.QoL.escapeHTML(s)
    : String(s).replace(/[&<>]/g, c => ({ '&':'&amp;','<':'&lt;','>':'&gt;' }[c]));
}
function formatExp(x) {
  if (x === null || x === undefined || !isFinite(x)) return '—';
  if (x === 0) return '0';
  return x.toExponential(3);
}

function setStatus(s) {
  const el = $('#status');
  if (s.kind === 'raw') { el.innerHTML = s.html; return; }
  const cls = s.kind === 'err' ? 'err' : s.kind === 'warn' ? 'warn' : s.kind === 'ok' ? 'ok' : '';
  el.innerHTML = cls ? `<span class="${cls}">${escapeHTML(s.text)}</span>` : escapeHTML(s.text);
}

// ---------- (Riemann-map + alternates + bg-search) -> ui-solve.js --------

// ===========================================================================
// Canvas plotting — DomainPlot class lives in ui-domain-plot.js (P0 split).
// It receives its ui.js closures (state, modeDescriptor, formatTick, sub) via
// dependency injection. The installer call below resolves the class for use
// in the DOM wiring section further down.
// ===========================================================================


function formatTick(v, step) {
  // Choose precision based on step size
  const digits = Math.max(0, -Math.floor(Math.log10(step)));
  return Number(v.toFixed(digits)).toString();
}

// ---------- Wire everything up ------------------------------------------
// DomainPlot lives in ui-domain-plot.js; install with the ui.js closures
// it needs (state, modeDescriptor, formatTick, sub).
const DomainPlot = window.QD_UI.installDomainPlot({ state, modeDescriptor, formatTick, sub });
const plot = new DomainPlot($('#canvas'), $('#plot-readout'));
// Coalesce a burst of resize events (window drag) into one canvas resize per
// frame (A8). plot.resize() reallocates the canvas backing store + repaints,
// so running it on every raw resize event is wasteful.
let _resizeScheduled = false;
window.addEventListener('resize', () => {
  if (_resizeScheduled) return;
  _resizeScheduled = true;
  requestAnimationFrame(() => { _resizeScheduled = false; plot.resize(); });
});

// HANDOFF #34 (revised): re-render the QD plot whenever the QD tab becomes
// active. The 2D canvas is shared with Schwarz (CPU pyramid + orbit polyline)
// and Param-slice (image data); without this the user sees stale graphics
// from the previous tab until they pan, zoom, or re-solve.
//
// We MUST defer the render to a microtask so it runs AFTER every other
// synchronous tab-changed listener for this dispatch. In particular,
// schwarz-ui.js's exit branch (registered later in script load order)
// clears the 2D canvas via ctx.clearRect — if we render synchronously we
// fire first and Schwarz's clear immediately wipes the freshly-drawn
// pixels (the bug fixed here). Microtasks drain before the browser
// paints, so there is no flicker.
//
// Future tabs adding their own exit-clear would be vulnerable to the
// same listener-order trap; keep this microtask deferral.
document.addEventListener('tab-changed', e => {
  if (!e.detail || e.detail.tab !== 'qd') return;
  queueMicrotask(() => {
    // Stale-tab guard: a rapid qd → schwarz double-click would otherwise
    // briefly paint the QD canvas while the user is already on Schwarz.
    const active = document.querySelector('.tab-btn.active');
    if (!active || active.dataset.tab !== 'qd') return;
    plot.resize();
  });
});

// Click-and-drag of quadrature nodes on the plot. Live updates while
// dragging via quick-solve; on release run the full solver pipeline so
// alternates and background search get a proper pass at the new value.
plot.onPoleDrag = (idx, w) => {
  if (idx < 0 || idx >= state.poles.length) return;
  // Move the marker to the cursor IMMEDIATELY, decoupled from the solve, so the
  // dot never appears stuck while a (possibly slow) live-solve frame is in
  // flight. The boundary catches up when the live result lands.
  plot.setLivePole(idx, w);
  const text = QD.Complex.toString(w, 4);
  state.poles[idx].a = text;
  // Update the matching text input in the side panel (no slider for a_j).
  const aInput = document.querySelector(
    `#poles-list .pole[data-idx="${idx}"] input[data-field="a"]`);
  if (aInput) aInput.value = text;
  markAsCustom();
  scheduleQuickSolve();
};
plot.onPoleDragEnd = () => { scheduleSolve(); };

// Double-click on the plot drops a new simple pole (coefficient 1) at the
// clicked w. Inverse view only — the direct view hides this canvas, and the
// gesture only makes sense for the inverse problem (it edits h's poles).
plot.onAddPole = (w) => {
  if ((state.viewMode || 'inverse') !== 'inverse') return;
  addPoleAt(w);
};


// ---------- Custom h(w) text + modeAllowsPoly -> ui-h-text.js ------------
// modeAllowsPoly / refreshHText / setHTextMsg / parseAndApplyHText are installed
// via QD_UI.installHText near the end of this file (they read renderPolesList /
// renderPolyCoefList + the shared helpers via uiCtx).

// ---------- Preset dropdown ---------------------------------------------
function populatePresetDropdown() {
  const sel = $('#preset-select');
  const list = currentPresetList();
  sel.innerHTML = '<option value="">— custom —</option>' +
    list.map(p => `<option value="${p.id}">${escapeHTML(p.label)}</option>`).join('');
}

function applyPreset(id) {
  const p = currentPresetList().find(x => x.id === id);
  if (!p) return;
  applyConfig(p);
}

// Apply a config object { poles, c, w0, alpha, q, polyCoeffs } to the live state
// for the CURRENT mode (state.mode must already be set). Shared by family-preset
// loads (applyPreset) and thesis-example loads (loadThesisExample, #8).
function applyConfig(p) {
  state.poles = p.poles.map(po => ({
    a: po.a,
    order: po.order,
    residues: po.residues.slice(),
  }));
  if (state.mode === 'bounded') {
    state.w0Mode = 'auto';
    $('input[name="w0mode"][value="auto"]').checked = true;
    $('#w0-manual').disabled = true;
  } else if (state.mode === 'pqd-bounded' || state.mode === 'pqd-bounded-singular') {
    // α may be any real > 0 (α ≠ 1); default 2.
    let presetAlpha = +p.alpha;
    if (!(presetAlpha > 0) || presetAlpha === 1) presetAlpha = 2;
    state.alpha = presetAlpha;
    const inp = $('#alpha-input');
    if (inp) inp.value = String(presetAlpha);
    // Default to 'auto' (Centroid of poles) for ALL bounded PQD families,
    // singular included: buildNorm uses the live centroid (recomputed each solve,
    // so it tracks pole drags and stays interior) and validates w₀ ≠ 0, surfacing
    // a clear error if the centroid degenerates to the origin. The preset's w₀ is
    // pre-filled into the Manual field for users who switch to Manual.
    state.w0Manual = p.w0 || '1';
    $('#w0-manual').value = state.w0Manual;
    state.w0Mode = 'auto';
    $('input[name="w0mode"][value="auto"]').checked = true;
    $('#w0-manual').disabled = true;
  } else if (state.mode === 'lqd-bounded' || state.mode === 'lqd-bounded-singular') {
    // Default to 'auto' (Centroid of poles), like the other bounded families.
    // buildNorm validates w₀ ≠ 0 (LQD needs a non-origin w₀); the preset's w₀ is
    // pre-filled into the Manual field for users who switch to Manual.
    state.w0Mode = 'auto';
    state.w0Manual = p.w0 || '1';
    $('input[name="w0mode"][value="auto"]').checked = true;
    $('#w0-manual').disabled = true;
    $('#w0-manual').value = state.w0Manual;
    // Singular: reset q to '0' by default (the user can dial it via slider).
    if (state.mode === 'lqd-bounded-singular') {
      const presetQ = p.q || '0';
      setQ(presetQ);
    }
  } else {
    if (typeof p.c === 'number') setC(p.c);
    // Unbounded PQD presets carry α (any real > 0, α ≠ 1; default 2).
    if (state.mode === 'pqd-unbounded' || state.mode === 'pqd-unbounded-singular') {
      let presetAlpha = +p.alpha;
      if (!(presetAlpha > 0) || presetAlpha === 1) presetAlpha = 2;
      state.alpha = presetAlpha;
      const inp = $('#alpha-input');
      if (inp) inp.value = String(presetAlpha);
    }
    // Polynomial part: preset may provide polyCoeffs (string[]).
    if (Array.isArray(p.polyCoeffs)) {
      state.polyCoeffs = p.polyCoeffs.slice();
      state.polyDegree = p.polyCoeffs.length - 1;
    } else {
      state.polyDegree = -1;
      state.polyCoeffs = [];
    }
    syncPolyDegreeInput();
  }
  for (const k of Object.keys(magSliderMax)) delete magSliderMax[k];
  renderPolesList();
  renderPolyCoefList();
  scheduleSolve();
}

// Load a thesis example (#8): switch to its family, apply its config, frame the
// view, and turn on the annotated-phenomena overlay. The oracle card + dropdown
// bookkeeping live in ui-thesis.js; this is the ui.js-internal half that needs
// setMode / applyConfig / plot. Exposed via uiCtx.loadThesisExample.
function loadThesisExample(ex) {
  if (!ex || !MODES[ex.mode]) return;
  if (state.mode !== ex.mode) { state.mode = ex.mode; applyModeVisuals(); }
  $('#preset-select').value = '';            // not a family preset
  applyConfig(ex);                           // poles/c/q/alpha/w0/polyCoeffs (+ scheduleSolve)
  if (ex.view && plot && plot.view) {
    plot.view.cx = ex.view.cx; plot.view.cy = ex.view.cy; plot.view.scale = ex.view.scale;
    state.autoFit = false;                   // honor the framed view (don't auto-fit over it)
    plot.render();
  }
  state.showPhenomena = true;                // self-explaining exhibit
  const tog = $('#phenomena-toggle'); if (tog) tog.checked = true;
}

// Programmatic setter for q (complex) that keeps text input, |q| / arg sliders,
// and state in sync. Used by applyPreset and by the q-event handlers.
function setQ(qStr) {
  state.q = qStr;
  const qval = QD.Complex.parse(qStr) || { re: 0, im: 0 };
  const mag = Math.hypot(qval.re, qval.im);
  const arg = Math.atan2(qval.im, qval.re);
  const txt = $('#q-manual');
  if (txt) txt.value = qStr;
  const magS = $('#q-mag-slider');
  const argS = $('#q-arg-slider');
  if (magS) {
    if (mag > +magS.max) magS.max = (mag * 1.5).toFixed(3);
    magS.value = mag;
  }
  if (argS) argS.value = arg;
  const magL = $('#q-mag-val');
  const argL = $('#q-arg-val');
  if (magL) magL.textContent = mag.toFixed(3);
  if (argL) argL.textContent = fmtArg(arg);
}

// Programmatic setter for c that keeps slider, text input, and state in sync.
function setC(c) {
  state.c = c;
  const slider = $('#c-slider');
  const text   = $('#c-manual');
  if (slider) {
    if (c > +slider.max) slider.max = (c * 1.5).toFixed(3);
    slider.value = c;
  }
  if (text) text.value = c;   // the number input doubles as the readout
}

// Update every mode-dependent visual (normalization cards, α-card, hints,
// vector-field label, preset dropdown, polynomial UI) to match the current
// state.mode. Driven entirely by the per-mode descriptor (modeDescriptor() →
// the MODES table), so it scales to all ten modes without per-mode branches
// here — each descriptor declares which cards.{w0,c,poly,q,alpha} it shows.
// Preserves the user's pole data (a_j and C_{j,s}); only swaps chrome.
// Does NOT force the w₀ mode and does NOT trigger a solve — those are setMode's
// job. Factored out of setMode (§23) so the auto-switch reflector can update
// the UI to the solver-chosen family without kicking off a redundant re-solve.
// CANONICAL mode-refresh. INVARIANT: any code path that assigns `state.mode`
// programmatically MUST call this immediately after (or go through setMode, which
// does). It is the single source of truth for mode-dependent card visibility
// (#c-card, the #map-params-card row groups, the poly section), α validation, the
// c input (via setC), the poly list, and the domain-type segmented control (via
// syncDomainModeControl). Bypassing it leaves stale UI until the next refresh
// (the historical _sendHToInverseTab c-card bug). It is idempotent — safe to call
// even when the mode is unchanged.
function applyModeVisuals() {
  const desc = modeDescriptor();
  // Card visibility from descriptor.
  // Map-parameters card merges the φ(0), PQD-α, and singular-LQD-q knobs
  // (index.html #map-params-card). Toggle each row group from its descriptor
  // flag, and the parent card iff ANY group is shown. The conformal radius c
  // stays in its own #c-card (nested under the φ(z) display in the Domain card).
  $('#map-alpha-rows').classList.toggle('hidden', !desc.cards.alpha);
  $('#map-w0-rows').classList.toggle('hidden',    !desc.cards.w0);
  $('#map-q-rows').classList.toggle('hidden',     !desc.cards.q);
  $('#map-params-card').classList.toggle('hidden',
    !(desc.cards.alpha || desc.cards.w0 || desc.cards.q));
  $('#c-card').classList.toggle('hidden',            !desc.cards.c);
  $('#poly-part-section').classList.toggle('hidden', !desc.cards.poly);
  // α rows: visible only in PQD modes. When hidden, force state.alpha back to 1
  // so the next mode-switch doesn't carry a stale PQD config; when shown, ensure
  // state.alpha is a valid PQD value (> 0, ≠ 1; default 2).
  if (desc.cards.alpha) {
    if (!(state.alpha > 0) || state.alpha === 1) state.alpha = 2;
    const inp = $('#alpha-input');
    if (inp) inp.value = String(state.alpha);
  } else if (state.alpha !== 1) {
    state.alpha = 1;
    const inp = $('#alpha-input');
    if (inp) inp.value = '2';            // input default for next visit
  }
  // Hint elements: show only the one this mode names (if any).
  for (const hintId of ['lqd-hint', 'lqd-singular-hint', 'lqd-unbounded-hint', 'lqd-unbounded-singular-hint', 'pqd-hint', 'pqd-singular-hint', 'pqd-unbounded-hint', 'pqd-unbounded-singular-hint']) {
    const el = $('#' + hintId);
    if (el) el.style.display = (desc.hint === hintId) ? '' : 'none';
  }
  const vfExt = $('#vf-external-opt');
  if (vfExt) vfExt.textContent = desc.externalFieldLabel;
  populatePresetDropdown();
  markAsCustom();
  setC(state.c);
  // Sync the polynomial UI with state.polyDegree.
  syncPolyDegreeInput();
  renderPolyCoefList();
  // Reflect the current mode into the compact domain-type control. Runs on
  // EVERY mode change — interactive, §23 auto-switch, preset, URL restore.
  syncDomainModeControl(state.mode);
  // Item 5: a plain-language one-liner of what's being solved.
  const summ = $('#dm-summary');
  if (summ) summ.textContent = modeSummary(state.mode);
}

// Item 5: plain-language description of the active mode, for #dm-summary.
function modeSummary(mode) {
  const d = decomposeMode(mode);
  const weight = d.weight === 'classical' ? 'classical (unweighted)'
    : d.weight === 'pqd' ? 'power-weighted (|w|^(2(α−1)))'
    : 'log-weighted (1/|w|²)';
  const extent = d.domain === 'bounded' ? 'bounded' : 'unbounded (reaches ∞)';
  const sing = d.singular ? ', with the origin inside Ω' : '';
  return `Solving for a ${extent} ${weight} quadrature domain Ω from your h(w)${sing}.`;
}

// ---------- Compact 3-axis domain-type control --------------------------
// The 10 modes factor as {weight} × {bounded|unbounded} × {singular}. Classical
// has no singular variant. composeMode maps the three controls → a MODES key;
// decomposeMode is the inverse; syncDomainModeControl reflects state.mode back
// into the buttons + checkbox (called from applyModeVisuals).
function composeMode(weight, domain, singular) {
  if (weight === 'classical') return domain;            // 'bounded' | 'unbounded'
  return `${weight}-${domain}${singular ? '-singular' : ''}`;
}
function decomposeMode(mode) {
  if (mode === 'bounded' || mode === 'unbounded') {
    return { weight: 'classical', domain: mode, singular: false };
  }
  const m = /^(pqd|lqd)-(bounded|unbounded)(-singular)?$/.exec(mode);
  if (!m) return { weight: 'classical', domain: 'bounded', singular: false };
  return { weight: m[1], domain: m[2], singular: !!m[3] };
}
function syncDomainModeControl(mode) {
  const d = decomposeMode(mode);
  $$('#dm-weight .seg-btn').forEach(b => b.classList.toggle('active', b.dataset.weight === d.weight));
  $$('#dm-domain .seg-btn').forEach(b => b.classList.toggle('active', b.dataset.domain === d.domain));
  const sing = $('#dm-singular');
  if (sing) {
    const classical = d.weight === 'classical';   // Classical has no singular variant
    sing.checked = d.singular;
    sing.disabled = classical;
    sing.parentElement.style.opacity = classical ? '0.45' : '';
  }
}
// Read the three controls and switch to the composed mode.
function applyDomainModeControl() {
  const weight = ($('#dm-weight .seg-btn.active') || {}).dataset
    ? $('#dm-weight .seg-btn.active').dataset.weight : 'classical';
  const domain = ($('#dm-domain .seg-btn.active') || {}).dataset
    ? $('#dm-domain .seg-btn.active').dataset.domain : 'bounded';
  const sing = $('#dm-singular');
  const singular = !!(sing && sing.checked && weight !== 'classical');
  setMode(composeMode(weight, domain, singular));
}

function setMode(newMode) {
  if (!MODES[newMode]) return;
  if (state.mode === newMode) return;
  state.mode = newMode;
  applyModeVisuals();
  // All families default to 'auto' (Centroid of poles) — including LQD and the
  // singular families. buildNorm computes the centroid each solve and validates
  // w₀ ≠ 0 where required (LQD / singular), surfacing a clear error if it
  // degenerates; the user can switch to Manual at any time.
  scheduleSolve();
}

function syncPolyDegreeInput() {
  const inp = $('#poly-degree');
  if (inp) inp.value = state.polyDegree;
}

// Selecting a preset loads it; the user editing anything afterward reverts
// the dropdown to "— custom —".
function markAsCustom() {
  $('#preset-select').value = '';
  // A manual edit (or a family-preset load) means we're no longer showing a
  // thesis example — ui-thesis.js listens for this to clear its oracle card.
  document.dispatchEvent(new CustomEvent('qd-customized'));
}

// Per-pole event delegation. Handles three kinds of `input` events:
//   • text fields for a_j, order, and C_{j,s}
//   • magnitude (|C|) range sliders
//   • argument (arg) range sliders
// Slider changes write back to the C_{j,s} text field and trigger a solve;
// text-field changes update both sliders to match.
$('#poles-list').addEventListener('input', e => {
  const t = e.target;
  const poleDiv = t.closest('.pole');
  if (!poleDiv) return;
  const idx = +poleDiv.dataset.idx;
  const pole = state.poles[idx];

  // Any manual edit reverts the preset dropdown to "— custom —".
  markAsCustom();

  // --- mag / arg sliders ---
  if (t.classList.contains('slider1d-mag') || t.classList.contains('slider1d-arg')) {
    const sIdx = +t.dataset.s;
    const block = t.closest('.residue-block');
    const magSlider = $('.slider1d-mag', block);
    const argSlider = $('.slider1d-arg', block);
    const mag = +magSlider.value;
    const arg = +argSlider.value;
    $('.mag-val', block).textContent = mag.toFixed(3);
    $('.arg-val', block).textContent = fmtArg(arg);
    const c = { re: mag * Math.cos(arg), im: mag * Math.sin(arg) };
    const text = Complex.toString(c, 4);
    pole.residues[sIdx] = text;
    $('.residue', block).value = text;
    scheduleQuickSolve();   // live update during drag
    return;
  }

  // --- text inputs ---
  const field = t.dataset.field;
  if (field === 'a') {
    pole.a = t.value;
  }
  else if (field === 'order') {
    const newOrder = Math.max(1, Math.min(6, +t.value || 1));
    pole.order = newOrder;
    while (pole.residues.length < newOrder) pole.residues.push('0');
    pole.residues.length = newOrder;
    renderPolesList();
  }
  else if (field === 'residue') {
    const sIdx = +t.dataset.s;
    pole.residues[sIdx] = t.value;
    // Sync the two sliders to the parsed value.
    const c = Complex.parse(t.value);
    if (c) {
      const mag = Math.hypot(c.re, c.im);
      const arg = Math.atan2(c.im, c.re);
      const block = t.closest('.residue-block');
      const magSlider = $('.slider1d-mag', block);
      const argSlider = $('.slider1d-arg', block);
      const key = residueKey(idx, sIdx);
      const newMax = magMaxFor(key, mag);
      if (newMax > +magSlider.max) magSlider.max = newMax;
      magSlider.value = mag;
      argSlider.value = arg;
      $('.mag-val', block).textContent = mag.toFixed(3);
      $('.arg-val', block).textContent = fmtArg(arg);
    }
  }
  scheduleSolve();
});
$('#poles-list').addEventListener('click', e => {
  // Use closest() so the handler still fires if the click lands on a child
  // of the remove button (e.g. an inner glyph element).
  const removeBtn = e.target.closest('[data-action="remove"]');
  if (!removeBtn) return;
  const poleDiv = removeBtn.closest('.pole');
  if (!poleDiv) return;
  markAsCustom();
  removePoleAt(+poleDiv.dataset.idx);
});

// When a slider drag ENDS (mouseup), re-run the full solver to refresh
// alternates and kick off the background search again.
$('#poles-list').addEventListener('change', e => {
  if (e.target.classList.contains('slider1d-mag') ||
      e.target.classList.contains('slider1d-arg')) {
    scheduleSolve();
  }
});
$('#add-pole').addEventListener('click', () => { markAsCustom(); addPole(); });

// ---------- View-mode toggle (HANDOFF #30) ------------------------------
// Inverse | direct segmented control at the top of the QD tab. The inverse
// view is the existing QD/LQD UI (wrapped in #qd-inverse-content in
// index.html); the direct view is the former Direct-problem tab UI,
// relocated into #controls-direct (now a sibling of #qd-inverse-content
// inside #controls-qd). Direct UI is lazy-mounted on first switch.
function mountViewToggle() {
  const qdRoot = document.getElementById('controls-qd');
  if (!qdRoot) return;
  const card = document.createElement('section');
  card.id = 'qd-view-toggle';
  card.className = 'card';
  card.innerHTML = `
    <div class="segmented" role="tablist" aria-label="View mode">
      <button class="seg-btn active" data-view="inverse" type="button">inverse</button>
      <button class="seg-btn"        data-view="direct"  type="button">direct</button>
    </div>
  `;
  // Insert as the first child of #controls-qd, BEFORE #qd-inverse-content.
  qdRoot.insertBefore(card, qdRoot.firstChild);
  card.querySelectorAll('.seg-btn').forEach(btn => {
    btn.addEventListener('click', () => setViewMode(btn.dataset.view));
  });
}

function setViewMode(mode) {
  if (mode !== 'inverse' && mode !== 'direct') return;
  if (mode === state.viewMode) return;
  state.viewMode = mode;
  // Toggle segmented-control highlight.
  document.querySelectorAll('#qd-view-toggle .seg-btn').forEach(btn => {
    btn.classList.toggle('active', btn.dataset.view === mode);
  });
  const inv = document.getElementById('qd-inverse-content');
  const dir = document.getElementById('controls-direct');
  if (inv) inv.style.display = (mode === 'inverse') ? '' : 'none';
  if (dir) dir.style.display = (mode === 'direct')  ? '' : 'none';
  if (mode === 'direct') {
    // Lazy-mount Direct UI on first switch.
    if (!state.directMounted && window.QD && QD.Direct && QD.Direct._mountUI) {
      QD.Direct._mountUI();
      state.directMounted = true;
    }
    if (window.QD && QD.Direct && QD.Direct._activate) QD.Direct._activate();
  }
  updateStatusPanelVisibility();   // panel shows only on the inverse view
}

mountViewToggle();

// -----------------------------------------------------------------------------
// QoL: attach "?" help buttons to the inverse-tab cards (HANDOFF #33).
// Static cards in index.html; lazy-mounted cards (Direct, Schwarz, Param-slice)
// wire their own help inside their respective ui modules.
// -----------------------------------------------------------------------------
function mountQolHelp() {
  if (!window.QD || !window.QD.QoL || !window.QD.QoL.attachHelp) return;
  const H = window.QD.QoL.attachHelp;
  // All help-popover prose lives in app/ui-strings.js (QD.Strings.help.*) so it
  // can be edited in one place. Bail quietly if the strings module isn't loaded.
  const help = window.QD.Strings && window.QD.Strings.help;
  if (!help) return;
  const headerOf = (cardSelector) => {
    const card = document.querySelector(cardSelector);
    return card ? card.querySelector('h2') : null;
  };
  // Item 6: an app-level "What is a quadrature domain?" intro, as a "?" next to
  // the title — the missing on-ramp for a newcomer.
  const title = document.querySelector('.app-header-row h1');
  if (title) H(title, help.intro);
  H(headerOf('#domain-mode-card'), help.domainType);
  H(headerOf('#h-card'), help.hCard);
  H(headerOf('#map-params-card'), help.mapParams);
  H(headerOf('#c-card'), help.cCard);
  H(headerOf('#solver-settings-card'), help.solverSettings);
  H(headerOf('#search-options-card'), help.searchOptions);
  H(headerOf('#status-card'), help.status);
  H(document.querySelector('#sp-geom summary'), help.geom);
  H(document.querySelector('#sp-cusps summary'), help.cusps);
  // (The Riemann-map symbolic identity is shown via the "?" toggle next to the
  // numerical φ(z) in the Domain-type tile — no separate help popover.)
  H(headerOf('#alternates-card'), help.alternates);
}
mountQolHelp();

// Copy-link affordance: surface the (already-maintained) URL-hash state as a
// one-click shareable link. Reuses QD.QoL.copyButton (clipboard + toast); the
// hash already encodes mode / h(w) / gauges / active tab via ui-url-state.js.
(function mountCopyLink() {
  const host = $('#copy-link-host');
  if (!host || !(window.QD && QD.QoL && QD.QoL.copyButton)) return;
  const btn = QD.QoL.copyButton(() => location.href,
    { title: 'Copy a shareable link to this configuration' });
  btn.classList.remove('copy-btn');
  btn.classList.add('small');
  btn.textContent = '🔗 Copy link';
  host.appendChild(btn);
})();

// Relocate the advanced "Search options" card to the BOTTOM of the inverse
// sidebar — it's rarely touched, so it shouldn't sit between the everyday
// controls. (Done in JS rather than in markup to keep the source order readable
// and the move trivially reversible.)
{
  const soCard = document.getElementById('search-options-card');
  const host = document.getElementById('qd-inverse-content');
  if (soCard && host) host.appendChild(soCard);
}

// QoL: copy button on the h(w) text input (HANDOFF #33).
(function mountHTextCopyButton() {
  if (!window.QD || !window.QD.QoL || !window.QD.QoL.copyButton) return;
  const parseBtn = document.getElementById('h-parse');
  if (!parseBtn) return;
  const copy = window.QD.QoL.copyButton(() => {
    const inp = document.getElementById('h-text');
    return inp ? inp.value : '';
  }, { title: 'Copy h(w) text' });
  copy.style.marginLeft = '6px';
  parseBtn.parentNode.insertBefore(copy, parseBtn.nextSibling);
})();

// Domain-type toggle
// Compact domain-type control: weight + domain segmented buttons + singular
// checkbox all route through applyDomainModeControl() → setMode().
$$('#dm-weight .seg-btn').forEach(b => b.addEventListener('click', () => {
  $$('#dm-weight .seg-btn').forEach(x => x.classList.toggle('active', x === b));
  // Classical has no singular variant — clear it before composing the mode.
  if (b.dataset.weight === 'classical') { const s = $('#dm-singular'); if (s) s.checked = false; }
  applyDomainModeControl();
}));
$$('#dm-domain .seg-btn').forEach(b => b.addEventListener('click', () => {
  $$('#dm-domain .seg-btn').forEach(x => x.classList.toggle('active', x === b));
  applyDomainModeControl();
}));
{
  const sing = $('#dm-singular');
  if (sing) sing.addEventListener('change', applyDomainModeControl);
}
// Reflect the initial mode into the control NOW so the singular checkbox starts
// disabled under the default classical mode (it was toggleable until the first
// mode change otherwise).
syncDomainModeControl(state.mode);

// Relocate the conformal-radius card (#c-card) to sit directly below the
// Riemann-map formula inside the Domain-type tile. Its mode-driven visibility
// (applyModeVisuals toggles `.hidden`) is id-based, so the move is transparent.
{
  const cCard = document.getElementById('c-card');
  const dmRiemann = document.getElementById('dm-riemann');
  if (cCard && dmRiemann) dmRiemann.insertAdjacentElement('afterend', cCard);
}

// "?" toggle that reveals/hides the symbolic Riemann-map identity inline.
{
  const tog = $('#dm-riemann-sym-toggle');
  const sym = $('#dm-riemann-sym');
  if (tog && sym) tog.addEventListener('click', () => sym.classList.toggle('hidden'));
}

// Preset dropdown
populatePresetDropdown();
// Reflect the default load as a *named* preset rather than "— custom —": the
// initial state (bounded, h = 1/w) is exactly the 'unit-disk' preset, so show
// that label as the starting point. Setting .value to a non-existent option is
// a no-op, so this is safe if the preset list ever changes. (B2: orientation.)
{
  const sel = $('#preset-select');
  if (sel && currentPresetList().some(p => p.id === 'unit-disk')) sel.value = 'unit-disk';
}
$('#preset-select').addEventListener('change', e => {
  if (e.target.value) {
    document.dispatchEvent(new CustomEvent('qd-customized'));   // leaving any thesis example
    applyPreset(e.target.value);
  }
});

// Custom h(w) parse button + Enter-to-parse on the text input.
$('#h-parse').addEventListener('click', () => parseAndApplyHText());
$('#h-text').addEventListener('keydown', e => {
  if (e.key === 'Enter') { e.preventDefault(); parseAndApplyHText(); }
});

// w0 mode (bounded only)
$$('input[name="w0mode"]').forEach(r => r.addEventListener('change', e => {
  state.w0Mode = e.target.value;
  $('#w0-manual').disabled = state.w0Mode !== 'manual';
  scheduleSolve();
}));
$('#w0-manual').addEventListener('input', e => {
  state.w0Manual = e.target.value;
  if (state.w0Mode === 'manual') scheduleSolve();
});

// α input (PQD mode): free-entry number, any real α > 0 (α ≠ 1), plus
// quick-pick buttons (½, 2, 3). Routes Newton to Family.powerQD via the
// pqd-bounded mode's buildNorm/applyNorm.
const alphaInputEl = $('#alpha-input');
if (alphaInputEl) {
  alphaInputEl.addEventListener('input', e => {
    const a = +e.target.value;
    if (a > 0) { state.alpha = a; markAsCustom(); scheduleSolve(); }
  });
}
$$('.alpha-quick').forEach(btn => btn.addEventListener('click', () => {
  const a = +btn.getAttribute('data-alpha');
  if (!(a > 0)) return;
  state.alpha = a;
  if (alphaInputEl) alphaInputEl.value = String(a);
  markAsCustom();
  scheduleSolve();
}));

// Conformal radius c (unbounded only): slider drives quick solve, manual
// text drives full solve.
$('#c-slider').addEventListener('input', e => {
  const c = +e.target.value;
  state.c = c;
  $('#c-manual').value = c;   // the number input IS the readout (no #c-val span)
  markAsCustom();
  scheduleQuickSolve();
});
$('#c-slider').addEventListener('change', () => { scheduleSolve(); });
$('#c-manual').addEventListener('input', e => {
  const c = +e.target.value;
  if (!(c > 0) || !isFinite(c)) return;
  setC(c);
  markAsCustom();
  scheduleSolve();
});

// "Estimate max c" — automatically find c* = the largest conformal radius for
// which a VALID unbounded QD exists (univalent + quadrature identity), via the
// bracket+bisection estimator in solver-cmax.js driven through the solver worker.
// On success: cap the slider range at c* and jump to ≈0.99·c* (the largest clean
// domain). Only meaningful for unbounded families (the #c-card is unbounded-only,
// and the handler re-checks norm.unbounded).
const cEstimateBtn = $('#c-estimate-btn');
if (cEstimateBtn) cEstimateBtn.addEventListener('click', () => {
  const btn      = $('#c-estimate-btn');
  const busy     = $('#c-estimate-busy');
  const busyText = $('#c-estimate-busy-text');
  const resultEl = $('#c-estimate-result');
  const valueEl  = $('#c-estimate-value');
  const noteEl   = $('#c-estimate-note');
  const showResult = (html, note) => {
    if (valueEl) valueEl.innerHTML = html;
    if (noteEl)  noteEl.textContent = note || '';
    if (resultEl) resultEl.classList.remove('hidden');
  };
  if (!QD.estimateMaxConformalRadius) {
    setStatus({ kind: 'err', text: 'Conformal-radius estimator is unavailable.' });
    return;
  }
  btn.disabled = true;
  if (busy) busy.classList.remove('hidden');
  if (busyText) busyText.textContent = 'searching…';
  (async () => {
    try {
      const built = buildHData();
      if (!built || built.error) {
        setStatus({ kind: 'err', text: (built && built.error) || 'No valid input.' });
        return;
      }
      const norm = buildNormalization(built);
      if (norm.error) { setStatus({ kind: 'err', text: norm.error }); return; }
      if (!norm.unbounded) {
        setStatus({ kind: 'warn', text: 'Estimate max c applies only to unbounded families.' });
        return;
      }
      // Reuse the user's solver settings layered on a thorough preset (robust
      // near the existence boundary, so a hard-to-find root isn't mis-read as
      // "no valid QD"). findAlternates off — we only need the primary's validity.
      const preset = PRESETS.thorough || PRESETS.standard;
      const baseOpts = buildSolverOptions(preset, { findAlternates: false });
      applyNormToOpts(baseOpts, norm);

      const PSW = QD.PrimarySolverWorker;
      const solveFn = (PSW && typeof PSW.solve === 'function')
        ? (h, o) => PSW.solve(h, o)
        : (h, o) => QD.solveInverseQD(h, o);

      const res = await QD.estimateMaxConformalRadius(built, baseOpts, solveFn, {
        cStart: state.c,
        progress: (p) => {
          if (busyText && p && isFinite(p.c)) busyText.textContent = 'searching… c≈' + (+p.c).toFixed(3);
        },
      });

      if (!res.found) {
        if (res.reason === 'no-invalid-below-ceiling') {
          showResult('No finite maximum found',
                     'the unbounded QD stays valid up to the search ceiling c ≤ ' + res.ceiling.toFixed(2) + '.');
          setStatus({ kind: 'ok',
            text: 'No critical c found below the ceiling — the unbounded QD remains valid up to c ≤ ' + res.ceiling.toFixed(2) + '.' });
        } else {
          showResult('No valid QD at this scale',
                     'no valid unbounded QD at or below the current c — adjust the quadrature data.');
          setStatus({ kind: 'warn', text: 'No valid unbounded QD found at or below the current scale.' });
        }
        return;
      }

      const cStar = res.cMax;
      // Note the mechanism that set c* (res.mechanism): a CUSP (a φ′ zero reaching
      // |z| = 1, after which the boundary self-overlaps) vs a FOLD / existence limit
      // (the QD branch simply ends, with no cusp). For a cusp, annotate the angle.
      let note = (res.mechanism === 'fold')
        ? 'largest existing domain; beyond c* no valid unbounded QD exists (an existence/fold limit, not a cusp).'
        : 'largest clean domain; beyond c* the boundary cusps and self-overlaps.';
      try {
        if (res.mechanism !== 'fold' && res.phiAtMax && QD.classifyCusps) {
          const cz = QD.classifyCusps(res.phiAtMax);
          if (cz && cz.cusps && cz.cusps.length) {
            note = 'incipient cusp near θ ≈ ' + cz.cusps[0].thetaDeg.toFixed(0) +
                   '°; beyond c* the boundary cusps and self-overlaps.';
          }
        }
      } catch (_) { /* annotation is best-effort */ }

      // c* confidence (#11): how trustworthy the estimate is — a clean cusp with a
      // tight bracket reads high; a soft fold/existence limit reads lower.
      const confPct = (typeof res.confidence === 'number')
        ? Math.round(res.confidence * 100) : null;
      const confTag = (confPct != null)
        ? ' <span class="muted" title="Confidence that c* is correctly located: blends the cusp/fold mechanism cleanliness with the final bracket tightness.">(confidence ' + confPct + '%)</span>'
        : '';
      showResult('Max conformal radius <strong>c* ≈ ' + cStar.toFixed(4) + '</strong>' + confTag, note);

      // Cap the slider range at c* and jump EXACTLY to c* — the extremal domain.
      // c* (= res.cMax) is the largest c the estimator PROVED valid (cLowValid),
      // so a re-solve at exactly c* succeeds; no 0.99 backoff is needed. Set the
      // max FIRST (at full precision, so value == max) so setC's auto-expand
      // doesn't fight it. Re-solving (rather than reusing res.phiAtMax) keeps the
      // solution flowing through the normal pipeline so PrimarySolution subscribers
      // — the oracle and Faber cards — refresh.
      const slider = $('#c-slider');
      const text   = $('#c-manual');
      const capStr = cStar.toFixed(4);
      if (slider) slider.max = String(cStar);
      if (text)   text.max   = String(cStar);
      setC(cStar);
      markAsCustom();
      scheduleSolve();
      setStatus({ kind: 'ok',
        text: 'Estimated max conformal radius c* ≈ ' + capStr +
              '. Slider set to c* — the extremal (cusped/limit) domain; nudge c down slightly for a clean interior.' });
    } catch (e) {
      if (e && e.aborted) return;
      setStatus({ kind: 'err', text: 'Estimate max c error: ' + ((e && e.message) || e) });
    } finally {
      btn.disabled = false;
      if (busy) busy.classList.add('hidden');
    }
  })();
});

// Singular-LQD charge q (complex). Text input drives full solve; |q| / arg
// sliders drive the snappy warm-start path.
const qManualInp = $('#q-manual');
if (qManualInp) qManualInp.addEventListener('input', e => {
  const parsed = QD.Complex.parse(e.target.value);
  if (!parsed) return;
  setQ(e.target.value);
  markAsCustom();
  scheduleSolve();
});
const qMagS = $('#q-mag-slider');
const qArgS = $('#q-arg-slider');
function readQFromSliders() {
  const mag = +qMagS.value, arg = +qArgS.value;
  const c = { re: mag * Math.cos(arg), im: mag * Math.sin(arg) };
  const text = QD.Complex.toString(c, 4);
  state.q = text;
  $('#q-manual').value = text;
  $('#q-mag-val').textContent = mag.toFixed(3);
  $('#q-arg-val').textContent = fmtArg(arg);
}
if (qMagS) {
  qMagS.addEventListener('input', () => { readQFromSliders(); markAsCustom(); scheduleQuickSolve(); });
  qMagS.addEventListener('change', () => scheduleSolve());
}
if (qArgS) {
  qArgS.addEventListener('input', () => { readQFromSliders(); markAsCustom(); scheduleQuickSolve(); });
  qArgS.addEventListener('change', () => scheduleSolve());
}

// Polynomial part of h (unbounded mode)
$('#poly-degree').addEventListener('input', e => {
  const d = parseInt(e.target.value, 10);
  if (isNaN(d) || d < -1 || d > 6) return;
  state.polyDegree = d;
  if (d >= 0) {
    while (state.polyCoeffs.length < d + 1) state.polyCoeffs.push('0');
    state.polyCoeffs.length = d + 1;
  } else {
    // d = -1: clear coeffs but preserve them in case the user wants to come back?
    // For simplicity, just leave state.polyCoeffs as-is (render hides it).
  }
  markAsCustom();
  renderPolyCoefList();
  scheduleSolve();
});

// Per-coef events on the polynomial section.
$('#poly-coefs-list').addEventListener('input', e => {
  const t = e.target;
  const block = t.closest('.residue-block');
  if (!block) return;
  const l = +block.dataset.polyL;

  if (t.classList.contains('slider1d-poly-mag') || t.classList.contains('slider1d-poly-arg')) {
    const magSlider = $('.slider1d-poly-mag', block);
    const argSlider = $('.slider1d-poly-arg', block);
    const mag = +magSlider.value;
    const arg = +argSlider.value;
    $('.poly-mag-val', block).textContent = mag.toFixed(3);
    $('.poly-arg-val', block).textContent = fmtArg(arg);
    const c = { re: mag * Math.cos(arg), im: mag * Math.sin(arg) };
    const text = QD.Complex.toString(c, 4);
    state.polyCoeffs[l] = text;
    $('.poly-coef', block).value = text;
    markAsCustom();
    scheduleQuickSolve();
    return;
  }

  if (t.classList.contains('poly-coef')) {
    state.polyCoeffs[l] = t.value;
    const c = QD.Complex.parse(t.value);
    if (c) {
      const mag = Math.hypot(c.re, c.im);
      const arg = Math.atan2(c.im, c.re);
      const magSlider = $('.slider1d-poly-mag', block);
      const argSlider = $('.slider1d-poly-arg', block);
      const key = `poly-coef-${l}`;
      const newMax = magMaxFor(key, mag);
      if (newMax > +magSlider.max) magSlider.max = newMax;
      magSlider.value = mag;
      argSlider.value = arg;
      $('.poly-mag-val', block).textContent = mag.toFixed(3);
      $('.poly-arg-val', block).textContent = fmtArg(arg);
    }
    markAsCustom();
    scheduleSolve();
  }
});

// Slider drag end → full solve for the polynomial section.
$('#poly-coefs-list').addEventListener('change', e => {
  if (e.target.classList.contains('slider1d-poly-mag') ||
      e.target.classList.contains('slider1d-poly-arg')) {
    scheduleSolve();
  }
});

// solver settings
$('#samples').addEventListener('input', e => {
  state.samples = Math.max(50, Math.min(5000, +e.target.value || 500));
  // Re-render current solution with new sample count.
  if (state.current && state.current.success) {
    const all = [state.current.primary, ...(state.current.alternates || [])];
    showSolution(all[state.selectedSolutionIdx] || all[0], state.current.hData,
                 state.selectedSolutionIdx === 0);
  }
});
$('#aggressiveness').addEventListener('change', e => {
  state.aggressiveness = e.target.value;
  scheduleSolve();
});
$('#auto-fit').addEventListener('change', e => {
  state.autoFit = e.target.checked;
});
$('#auto-switch-singular').addEventListener('change', e => {
  state.autoSwitchSingular = e.target.checked;
  scheduleSolve();   // re-solve so the change takes effect immediately
});
$('#vector-field-mode').addEventListener('change', e => {
  state.vectorFieldMode = e.target.value;
  plot.render();
});
$('#critical-set-toggle').addEventListener('change', e => {
  state.showCriticalSet = e.target.checked;
  plot.render();
});
$('#curvature-toggle').addEventListener('change', e => {
  state.showCurvature = e.target.checked;
  plot.render();
});
$('#phenomena-toggle')?.addEventListener('change', e => {
  state.showPhenomena = e.target.checked;
  plot.render();
});
// Faber-roots overlay toggle. ui-faber.js keeps the card's "Plot roots" checkbox
// in sync with this Layers toggle (and pushes the root payload onto state.faberRoots).
$('#faber-roots-toggle')?.addEventListener('change', e => {
  state.showFaberRoots = e.target.checked;
  const cardBox = $('#faber-show-roots');
  if (cardBox && cardBox.checked !== e.target.checked) {
    cardBox.checked = e.target.checked;
    cardBox.dispatchEvent(new Event('change'));
  }
  plot.render();
});

// ---------- Search-options panel wiring ----------------------------------
// Every field updates state.searchOptions on `input`/`change`. Inputs that
// change solver behavior schedule a fresh solve; display-only toggles
// (showNonUnivalent / showIdFailing) only re-render the alternates panel.
(function wireSearchOptions() {
  const allInputs = [
    '#so-phase-direct', '#so-phase-continuation', '#so-phase-multistart',
    '#so-phase-diverse', '#so-phase-deflation',
    '#so-num-restarts', '#so-num-diverse', '#so-num-deflation',
    '#so-bg-chunks', '#so-bg-chunk-size', '#so-keep-searching',
    '#so-newton-maxiter', '#so-newton-tol', '#so-cont-tstart', '#so-cont-grow',
    '#so-defl-alpha', '#so-defl-p', '#so-defl-from-valid',
    '#so-uni-samples', '#so-id-tol',
    '#so-show-non-univalent', '#so-show-id-failing', '#so-auto-escalate',
    '#so-seed',
  ];
  const displayOnly = new Set(['#so-show-non-univalent', '#so-show-id-failing']);

  allInputs.forEach(sel => {
    const el = $(sel);
    if (!el) return;
    const evt = el.type === 'checkbox' ? 'change' : 'input';
    el.addEventListener(evt, () => {
      readSearchOptions();
      if (displayOnly.has(sel)) {
        refreshAlternatesPanel();
      } else {
        scheduleSolve();
      }
    });
  });

  $('#so-reset').addEventListener('click', () => {
    resetSearchOptions();
    scheduleSolve();
  });

  // Reseed button: just bumps the alt-search token and starts a fresh round
  // using the current primary, without re-solving.
  $('#so-reseed').addEventListener('click', () => {
    if (!state.current || !state.current.success) return;
    // Reconstruct the FULL norm from the solved phi so the alt-search re-selects
    // the same family. The old sparse {unbounded,c}|{w0} dropped alpha/lqd/
    // singular/q and misrouted every non-classical family to boundedQD.
    const norm = QD.normFromPhi(state.current.primary && state.current.primary.phi);
    if (!norm) return;
    startBackgroundAltSearch(state.current.hData, norm);
  });

  // Initial population.
  readSearchOptions();
})();

// Try-harder button: runs the "exhaustive" aggressiveness preset (which has
// a much larger multistart budget, tighter Newton, longer continuation, and
// implicitly more deflation rounds). Useful when the current default-level
// solve has flagged the primary as non-valid (spurious / non-univalent).
// Cancel an in-flight primary solve (B3).
{
  const cancelBtn = $('#solve-cancel-btn');
  if (cancelBtn) cancelBtn.addEventListener('click', cancelSolve);
}

// Status-panel collapse + dock (item 1). Collapse hides the body (badge stays);
// dock moves the whole panel into the sidebar so it never covers the domain.
// Both preferences persist in localStorage and survive reloads.
{
  const panel = $('#status-panel');
  const collapseBtn = $('#sp-collapse');
  const dockBtn = $('#sp-dock');
  const dockHost = $('#status-dock-host');
  const plotArea = $('#plot-area');
  const LS_COLLAPSE = 'qd-status-collapsed';
  const LS_DOCK = 'qd-status-docked';
  let docked = false;
  try {
    state.statusPanelCollapsed = localStorage.getItem(LS_COLLAPSE) === '1';
    docked = localStorage.getItem(LS_DOCK) === '1';
  } catch (e) { /* private mode → defaults */ }

  const applyCollapse = () => {
    if (!panel) return;
    panel.classList.toggle('collapsed', !!state.statusPanelCollapsed);
    if (collapseBtn) collapseBtn.setAttribute('aria-expanded', state.statusPanelCollapsed ? 'false' : 'true');
  };
  const applyDock = () => {
    if (!panel) return;
    panel.classList.toggle('docked', docked);
    if (docked && dockHost && panel.parentNode !== dockHost) dockHost.appendChild(panel);
    else if (!docked && plotArea && panel.parentNode !== plotArea) plotArea.appendChild(panel);
    if (dockBtn) {
      dockBtn.textContent = docked ? '⇥' : '⇤';
      dockBtn.title = docked ? 'Pop the panel back onto the plot' : 'Dock the panel into the sidebar (clear the plot)';
      dockBtn.setAttribute('aria-pressed', docked ? 'true' : 'false');
    }
  };

  if (collapseBtn) collapseBtn.addEventListener('click', () => {
    state.statusPanelCollapsed = !state.statusPanelCollapsed;
    try { localStorage.setItem(LS_COLLAPSE, state.statusPanelCollapsed ? '1' : '0'); } catch (e) { /* ignore */ }
    applyCollapse();
  });
  if (dockBtn) dockBtn.addEventListener('click', () => {
    docked = !docked;
    try { localStorage.setItem(LS_DOCK, docked ? '1' : '0'); } catch (e) { /* ignore */ }
    applyDock();
  });
  applyCollapse();
  applyDock();   // reflect the persisted dock state on load
}

// Item 4: first-run coachmark over the plot — shown once, then remembered.
{
  const coach = $('#plot-coach');
  const dismissBtn = $('#plot-coach-dismiss');
  let seen = true;
  try { seen = localStorage.getItem('qd-coach-seen') === '1'; } catch (e) { /* treat as seen */ }
  const hideCoach = () => {
    if (coach) coach.style.display = 'none';
    try { localStorage.setItem('qd-coach-seen', '1'); } catch (e) { /* ignore */ }
  };
  if (coach && !seen) {
    coach.style.display = 'flex';
    if (dismissBtn) dismissBtn.addEventListener('click', hideCoach);
    setTimeout(hideCoach, 14000);   // never linger
  } else if (coach) {
    coach.style.display = 'none';
  }
}

$('#try-harder-btn').addEventListener('click', () => {
  const btn = $('#try-harder-btn');
  const busy = $('#try-harder-busy');
  btn.disabled = true;
  busy.classList.remove('hidden');
  (async () => {
    try {
      const built = buildHData();
      if (!built || built.error) {
        setStatus({ kind: 'err', text: 'No valid input.' });
        return;
      }
      const norm = buildNormalization(built);
      if (norm.error) {
        setStatus({ kind: 'err', text: norm.error });
        return;
      }
      const preset = PRESETS.exhaustive;
      const opts = buildSolverOptions(preset, { findAlternates: true });
      applyNormToOpts(opts, norm);
      // Route through the worker when available; main-thread fallback runs
      // inline (the spinner above + the async wrap let it still paint).
      const PSW = QD.PrimarySolverWorker;
      const result = (PSW && typeof PSW.solve === 'function')
        ? await PSW.solve(built, opts)
        : QD.solveInverseQD(built, opts);
      state.current = result;
      state.current.hData = built;
      state.current.w0Used = norm.w0 || (state.current.primary && state.current.primary.phi && state.current.primary.phi.w0);
      state.current.cUsed  = norm.c;
      state.current.unbounded = !!norm.unbounded;
      state.selectedSolutionIdx = 0;
      publishPrimarySolution();
      if (!result.success) {
        setStatus({ kind: 'err',
          text: 'Exhaustive search still found no algebraic root.\n  reason: ' + result.error });
        return;
      }
      showSolution(result.primary, built, /*isPrimary=*/false);
      refreshAlternatesPanel();
    } catch (e) {
      if (e && e.aborted) return;
      setStatus({ kind: 'err', text: 'Try-harder error: ' + (e && e.message || e) });
    } finally {
      btn.disabled = false;
      busy.classList.add('hidden');
    }
  })();
});

// plot controls
$('#btn-fit').addEventListener('click', () => plot.fit());
$('#btn-reset').addEventListener('click', () => plot.reset());

// alternates panel: click to view
$('#alternates-list').addEventListener('click', e => {
  const idx = e.target.dataset.altIdx;
  if (idx !== undefined) viewSolutionByIndex(+idx);
});

// ===========================================================================
// Phase-3 UI modularization (item E) — extracted-module installs.
// ---------------------------------------------------------------------------
// uiCtx (declared at the top) carries the ui.js closures the extracted modules
// need. Populate the shared helpers, then install each module and capture its
// exports back into the forward-declared lets with their ORIGINAL names, so
// every call site above is unchanged. Order: pole-grid + h-text (mutually
// referential, resolved via uiCtx at call time) before url-state (whose
// applyUrlState calls parseAndApplyHText). The initial pole/poly render runs
// here, AFTER the installs, so the now-let-bound renderers + h-text mirror are
// available (they were function-declaration-hoisted before this split).
// ===========================================================================
Object.assign(uiCtx, {
  $, $$, sub, debounce, plot,
  escapeHTML, escapeAttr, formatExp, setStatus,
  residueKey, magMaxFor, fmtArg, magSliderMax,
  syncPolyDegreeInput, markAsCustom,
  applyModeVisuals, setC, setQ,
  buildHData, buildW0, buildNormalization, buildSolverOptions,
  buildAltSearchOptions, applyNormToOpts, publishPrimarySolution,
});

// Solve / render / analysis pipeline (installed before pole-grid/h-text so the
// shared scheduleSolve they call via uiCtx is available).
({
  scheduleSolve, scheduleQuickSolve, solveAndRender, cancelSolve,
  showSolution, refreshAlternatesPanel, viewSolutionByIndex,
  startBackgroundAltSearch, updateStatusPanelVisibility,
} = window.QD_UI.installSolve(uiCtx));
Object.assign(uiCtx, {
  scheduleSolve, scheduleQuickSolve, solveAndRender, cancelSolve,
  showSolution, refreshAlternatesPanel, viewSolutionByIndex,
  startBackgroundAltSearch, updateStatusPanelVisibility,
});

({ renderPolesList, renderPolyCoefList } = window.QD_UI.installPoleGrid(uiCtx));
uiCtx.renderPolesList = renderPolesList;
uiCtx.renderPolyCoefList = renderPolyCoefList;

({ modeAllowsPoly, refreshHText, setHTextMsg, parseAndApplyHText } =
  window.QD_UI.installHText(uiCtx));
Object.assign(uiCtx, { modeAllowsPoly, refreshHText, setHTextMsg, parseAndApplyHText });

// URL/hash state (B1) — extracted to ui-url-state.js.
const { writeUrlState, applyUrlState } = window.QD_UI.installUrlState(uiCtx);
uiCtx.writeUrlState = writeUrlState;
uiCtx.applyUrlState = applyUrlState;

// Thesis-example gallery + analytic-oracle card (#8) — ui-thesis.js.
uiCtx.loadThesisExample = loadThesisExample;
if (window.QD_UI && window.QD_UI.installThesis) window.QD_UI.installThesis(uiCtx);

// Faber-polynomials card (UQD) — ui-faber.js. setFaberRoots is the decoupling
// hook: ui-faber.js writes the root payload + renders without touching `plot`.
uiCtx.setFaberRoots = (payload) => {
  state.faberRoots = payload;
  if (payload && !state.showFaberRoots) {
    state.showFaberRoots = true;
    const t = $('#faber-roots-toggle'); if (t) t.checked = true;
  }
  plot.render();
};
if (window.QD_UI && window.QD_UI.installFaber) window.QD_UI.installFaber(uiCtx);

// Initial structured-grid render (relocated from just after the plot setup, so
// the let-bound renderers + modeAllowsPoly exist by the time it runs).
renderPolesList();
renderPolyCoefList();
$('#poly-part-section').classList.toggle('hidden', !modeAllowsPoly(state.mode));
// Item 5: seed the "what you're solving" summary for the initial mode (later
// mode changes refresh it via applyModeVisuals).
{ const summ = $('#dm-summary'); if (summ) summ.textContent = modeSummary(state.mode); }

// Item 7: per-tab subtitle clarifying what each view does and that they all
// operate on the same solved domain (Schwarz / Param-slice show "using h(w)=…").
const TAB_SUBTITLES = {
  qd: 'Inverse problem: Recover the domain Ω and associated conformal map φ from a given quadrature function h(w)',
  schwarz: 'Iterate the Schwarz reflection of your current domain — the fractal tiling set.',
  'param-slice': 'Sweep a parameter and map where valid quadrature domains exist.',
};
function updateTabSubtitle(tab) {
  const el = $('#tab-subtitle');
  if (!el) return;
  let text = TAB_SUBTITLES[tab] || '';
  if (tab && tab !== 'qd') {
    const hInput = $('#h-text');
    const h = (hInput && hInput.value || '').trim();
    if (h) text += '  ·  using h(w) = ' + h;
  }
  el.textContent = text;
}
document.addEventListener('tab-changed', (e) => updateTabSubtitle(e.detail && e.detail.tab));
updateTabSubtitle('qd');

// Keep the URL in sync when the active tab changes (the QD-config writes
// happen via solveAndRender; this covers pure tab switches).
document.addEventListener('tab-changed', () => writeUrlState());
// Show/hide the bottom-right status panel with the QD tab (fires on every tab change).
document.addEventListener('tab-changed', () => updateStatusPanelVisibility());

// Initial solve — restore a shared/bookmarked config from the URL if present,
// otherwise solve the default state (B1). Item 3: a brand-new visitor (no saved
// URL, no prior visit) is greeted with a more compelling default than the bare
// unit disk — the cardioid (one pole, but it shows a cusp), so the first screen
// demonstrates what the app does. Returning visitors keep the plain default.
if (applyUrlState()) {
  scheduleSolve();             // ensure a solve even if no h param was present
} else {
  let firstVisit = false;
  try {
    firstVisit = !localStorage.getItem('qd-seen');
    localStorage.setItem('qd-seen', '1');
  } catch (e) { /* private mode → treat as returning */ }
  if (firstVisit && currentPresetList().some(p => p.id === 'cardioid')) {
    $('#preset-select').value = 'cardioid';
    applyPreset('cardioid');
  } else {
    solveAndRender();
  }
}

// ---------- Hooks for the Direct view (within the QD tab) ----------------------------
// direct-ui.js calls these to (a) push a ∂Ω preview onto the shared canvas,
// (b) send a computed h back to the QD tab's inverse view and switch view modes.
window.QD = window.QD || {};
window.QD.Direct = window.QD.Direct || {};

window.QD.Direct._setPlotBoundary = function (boundaryPts, opts) {
  // Display the user's φ-boundary on the canvas. Accepts an `unbounded`
  // flag in opts so the bounded-vs-unbounded shading convention matches
  // what the inverse solver uses. opts.overlayBoundary, if present, is
  // drawn as a dashed gold curve over the main boundary (used by the
  // round-trip diagnostic to show the inverse-recovered φ).
  opts = opts || {};
  plot.setData({
    boundaryPts,
    poles: [],
    w0: boundaryPts.length ? boundaryPts[0] : { re: 0, im: 0 },
    univalent: !boundarySelfIntersectsSimple(boundaryPts),
    unbounded: !!opts.unbounded,
    overlayBoundary: opts.overlayBoundary || null,
    vfMode: 'off',
    hData: { poles: [] },
    phi: null,
  });
};

window.QD.Direct._setPlotOverlay = function (overlayBoundary) {
  // Append/replace the overlay boundary without disturbing the main one.
  if (!plot.data) return;
  plot.data.overlayBoundary = overlayBoundary || null;
  plot.render();
};

// Cheap O(N²) self-intersection check — sufficient for the preview.
function boundarySelfIntersectsSimple(pts) {
  const N = pts.length;
  if (N < 4) return false;
  for (let i = 0; i < N; i++) {
    const a1 = pts[i], a2 = pts[(i + 1) % N];
    for (let j = i + 2; j < N; j++) {
      if (j === N - 1 && i === 0) continue;
      const b1 = pts[j], b2 = pts[(j + 1) % N];
      if (segmentsIntersect(a1, a2, b1, b2)) return true;
    }
  }
  return false;
}
function segmentsIntersect(p1, p2, p3, p4) {
  function ccw(a, b, c) {
    return (c.im - a.im) * (b.re - a.re) > (b.im - a.im) * (c.re - a.re);
  }
  return ccw(p1, p3, p4) !== ccw(p2, p3, p4) &&
         ccw(p1, p2, p3) !== ccw(p1, p2, p4);
}

window.QD.Direct._sendHToInverseTab = function (hData, opts) {
  // Populate the inverse-tab state from a Direct-computed h (+ the family params
  // carried in opts) and re-render, then switch to the inverse view.
  //
  // opts carries the FULL family identity: { unbounded?, alpha? (⇒ PQD), lqd?,
  // singular?, c?, q? }. We compose the matching inverse mode from it (so an
  // unbounded PQD lands in `pqd-unbounded`, NOT classical `unbounded`), set the
  // family parameters, then go through the canonical refresh.
  //
  // INVARIANT (see applyModeVisuals): any programmatic `state.mode` write MUST be
  // followed by applyModeVisuals() — the single source of truth for which
  // mode-dependent UI (#c-card, the #map-params-card row groups, poly section, …)
  // is visible and for syncing the domain-type control. The earlier bug here set
  // state.mode + syncDomainModeControl only, leaving #c-card hidden until a
  // refresh; applyModeVisuals fixes that and is idempotent.
  opts = opts || {};
  const weight   = opts.alpha ? 'pqd' : (opts.lqd ? 'lqd' : 'classical');
  const domain   = opts.unbounded ? 'unbounded' : 'bounded';
  const singular = !!opts.singular && weight !== 'classical';
  const mode     = composeMode(weight, domain, singular);

  state.poles = hData.poles.map(p => ({
    a: QD.Complex.toString(p.a, 6),
    order: p.principal.length,
    residues: p.principal.map(c => QD.Complex.toString(c, 6)),
  }));

  // Family parameters (set BEFORE applyModeVisuals: it validates α against the
  // mode and repopulates the c input from state.c via setC).
  if (typeof opts.alpha === 'number' && opts.alpha > 0) state.alpha = opts.alpha;
  if (typeof opts.c === 'number' && opts.c > 0) state.c = opts.c;
  if (domain === 'unbounded') {
    const polyPart = hData.polyPart || [];
    state.polyDegree = polyPart.length - 1;
    state.polyCoeffs = polyPart.map(c => QD.Complex.toString(c, 6));
  } else {
    state.polyDegree = -1;
    state.polyCoeffs = [];
  }

  // Apply the mode + refresh ALL mode-dependent UI (cards, α/c population, poly
  // list, domain-type control) in one canonical call.
  state.mode = mode;
  applyModeVisuals();

  // q (LQD-singular origin pole) lives in the now-visible #map-q-rows group.
  if (opts.q) setQ(QD.Complex.toString(opts.q, 6));

  renderPolesList();
  const ps = document.getElementById('preset-select');
  if (ps) ps.value = '';

  // Switch view-mode to inverse and re-solve (HANDOFF #30: the Direct UI is now a
  // view within the QD tab, so this is a view-mode switch rather than a tab switch).
  setViewMode('inverse');
  solveAndRender();
};

// =============================================================================
// Public hooks for the Parameter-slice tab.
//
// The slice tab needs two things from ui.js:
//   • snapshotScenario()              — read out the current { hData, norm, mode }
//   • loadScenarioIntoQdTab(s, mode)  — push a scenario back into state, re-solve,
//                                       and switch to the QD tab.
//
// We expose these on window.QD_UI so param-slice-ui.js can find them
// without having to be loaded after this file.
// =============================================================================
window.QD_UI = window.QD_UI || {};

window.QD_UI.snapshotScenario = function () {
  const built = buildHData();
  if (!built || built.error) return null;
  const norm = buildNormalization(built);
  if (norm.error) return null;
  // Defensive: ensure polyPart is always present (slice UI inspects it).
  if (!built.polyPart) built.polyPart = [];
  return { hData: built, norm, mode: state.mode };
};

window.QD_UI.loadScenarioIntoQdTab = function (scenario, mode) {
  if (!scenario || !scenario.hData) return;
  // Switch mode first (re-runs setMode's side-effects: card visibility,
  // preset dropdown, polynomial panel toggle). setMode is a no-op if the
  // mode hasn't changed.
  if (mode && mode !== state.mode) {
    setMode(mode);   // setMode → applyModeVisuals syncs the compact control
  }
  // Reflect hData into state.poles + state.polyDegree + state.polyCoeffs.
  const hData = scenario.hData;
  state.poles = hData.poles.map(p => ({
    a: QD.Complex.toString(p.a, 6),
    order: p.principal.length,
    residues: p.principal.map(c => QD.Complex.toString(c, 6)),
  }));
  const polyPart = hData.polyPart || [];
  state.polyDegree = polyPart.length - 1;
  state.polyCoeffs = polyPart.map(c => QD.Complex.toString(c, 6));

  // Reflect norm fields (c, q, w0) into state + DOM.
  const norm = scenario.norm || {};
  if (typeof norm.c === 'number' && norm.c > 0) {
    state.c = norm.c;
    const cInput  = $('#c-manual');
    const cSlider = $('#c-slider');
    if (cInput)  cInput.value  = norm.c.toString();
    if (cSlider) cSlider.value = norm.c.toString();
  }
  if (norm.q) {
    state.q = QD.Complex.toString(norm.q, 6);
    const qInput = $('#q-manual');
    if (qInput) qInput.value = state.q;
  }
  if (norm.w0) {
    state.w0Manual = QD.Complex.toString(norm.w0, 6);
    state.w0Mode = 'manual';
    const wManual = $('#w0-manual');
    const wRadio  = document.querySelector('input[name="w0mode"][value="manual"]');
    if (wManual) { wManual.value = state.w0Manual; wManual.disabled = false; }
    if (wRadio)  { wRadio.checked = true; }
  }

  syncPolyDegreeInput();
  renderPolesList();
  renderPolyCoefList();
  $('#poly-part-section').classList.toggle('hidden', !modeAllowsPoly(state.mode));
  markAsCustom();

  // Switch to the QD tab + run the full solver.
  const tabBtn = document.querySelector('.tab-btn[data-tab="qd"]');
  if (tabBtn) tabBtn.click();
  solveAndRender();
};
