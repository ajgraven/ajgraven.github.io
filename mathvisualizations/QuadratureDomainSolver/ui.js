// =============================================================================
// ui.js -- Frontend for the Quadrature Domain Solver
//
// File layout (search for the banner comments to jump):
//
//   Preset library          Quadrature-function presets (bounded / unbounded)
//   Aggressiveness presets  numRestarts / Newton tol / continuation steps
//   State                   The single source of truth for the UI
//   Search options          Advanced panel: per-phase toggles + overrides
//   Helpers                 $, debounce, fmtArg, subscripts
//   buildHData / buildW0    Read DOM → solver input
//   Polar slider bookkeeping
//   renderPolesList,        Build the pole-control cards
//   renderPolyCoefList      Build the poly-coefficient cards
//   solveAndRender,         Main solve pipeline (debounced)
//   quickSolveAndRender     Snappy warm-start path for slider drags
//   showSolution            Push a solution to the canvas + Riemann-map card
//   refreshAlternatesPanel  Build the alternates UI
//   startBackgroundAltSearch  Chunked async alternate-hunt after each solve
//   class DomainPlot        Canvas plot: axes, boundary, poles, vector field
//   DOM wiring              All event listeners
//
// All math is delegated to QD (= window.QD from solver.js); this file only
// translates between DOM state and QD calls.
// =============================================================================
'use strict';

// NB: `complex.js` already declares `Complex` in the shared script lexical
// scope, so we must NOT redeclare it here. Reference it via QD.Complex when
// needed, or via the unaliased global from complex.js.

// ===========================================================================
// Quadrature-function preset library
// ---------------------------------------------------------------------------
// Each preset is a quadrature function h(w) the user can load with one click.
// Selecting a preset replaces the pole data (and, for unbounded presets,
// also c). Any manual edit afterward reverts the dropdown to "— custom —".
//
// Preset shape:
//   {
//     id:         string,    // unique
//     label:      string,    // shown in dropdown
//     poles: [
//       { a: string, order: int, residues: [string,...] },
//       ...
//     ],
//     polyCoeffs: [string,...] | undefined,   // unbounded only
//     c:          number      | undefined,    // unbounded only
//   }
//
// All numeric values are stored as strings so they can pass through
// Complex.parse, which accepts forms like '1', '-0.5', '0.4+0.5i', 'i'.
//
// To ADD a preset: copy any entry, change the values. To REORDER: just
// reorder the array — the dropdown is built in array order.
// ===========================================================================

const QD_PRESETS_BOUNDED = [
  //  id                 label                                                 poles
  { id: 'unit-disk',     label: 'Unit disk:  h = 1/w',
    poles: [ { a: '0', order: 1, residues: ['1'] } ] },

  { id: 'cardioid',      label: 'Cardioid:  h = 1.5/w + 0.5/w²',
    poles: [ { a: '0', order: 2, residues: ['1.5', '0.5'] } ] },

  { id: 'two-point-sym', label: 'Two-point symmetric:  1.5/(w−1) + 1.5/(w+1)',
    poles: [
      { a:  '1', order: 1, residues: ['1.5'] },
      { a: '-1', order: 1, residues: ['1.5'] },
    ] },

  { id: 'triangle', label: 'Equilateral 3-point on unit circle',
    poles: [
      { a:  '1',              order: 1, residues: ['1'] },
      { a: '-0.5+0.8660254i', order: 1, residues: ['1'] },
      { a: '-0.5-0.8660254i', order: 1, residues: ['1'] },
    ] },
];

const QD_PRESETS_UNBOUNDED = [
  //  id                       label                                                  poles  + c (+ polyCoeffs)
  { id: 'unb-1pt-pos',  label: 'One-point positive charge:  h = 1/(w − 2),  c = 2',
    poles: [ { a: '2', order: 1, residues: ['1'] } ],
    c: 2 },

  { id: 'unb-1pt-neg',  label: 'One-point negative charge:  h = −0.5/(w − 2),  c = 0.7',
    poles: [ { a: '2', order: 1, residues: ['-0.5'] } ],
    c: 0.7 },

  { id: 'unb-1pt-imag', label: 'One-point imaginary charge:  h = i/(w − 2),  c = 0.8',
    poles: [ { a: '2', order: 1, residues: ['i'] } ],
    c: 0.8 },

  { id: 'unb-deltoid',  label: 'Deltoid:  h = w²,  c = 0.5',
    poles: [],
    polyCoeffs: ['0', '0', '1'],
    c: 0.5 },

  { id: 'unb-2pt-nonuniq', label: 'Two-point non-uniqueness:  1/(w−1) + 1/(w+1),  c = 0.4',
    poles: [
      { a:  '1', order: 1, residues: ['1'] },
      { a: '-1', order: 1, residues: ['1'] },
    ],
    c: 0.4 },
];

// ===========================================================================
// LQD presets (bounded non-singular log-weighted quadrature domains)
// ---------------------------------------------------------------------------
// Closed-form one-point cases from Theorem 5.3.2: Ω ∈ QD₀(α/(w−w₀)) iff
// 0 < α ≤ π² with Ω = {|ln(w/w₀)|² < α}; double point forms at α = π².
// Each preset includes the required w₀ since the LQD mode is manual-only
// (the centroid default doesn't apply when 0 ∉ Ω̄ is a hard constraint).
// ===========================================================================

const LQD_PRESETS_BOUNDED = [
  { id: 'lqd-1pt-small', label: 'One-pt: α = 0.5 / (w − 1),  w₀ = 1',
    poles: [ { a: '1', order: 1, residues: ['0.5'] } ],
    w0: '1' },

  { id: 'lqd-1pt-medium', label: 'One-pt: α = 2 / (w − 1),  w₀ = 1',
    poles: [ { a: '1', order: 1, residues: ['2'] } ],
    w0: '1' },

  { id: 'lqd-1pt-large', label: 'One-pt: α = 9 / (w − 1),  w₀ = 1  (near critical α = π²)',
    poles: [ { a: '1', order: 1, residues: ['9'] } ],
    w0: '1' },

  { id: 'lqd-1pt-shifted', label: 'Shifted one-pt: α = 0.4 / (w − 2),  w₀ = 2',
    poles: [ { a: '2', order: 1, residues: ['0.4'] } ],
    w0: '2' },

  { id: 'lqd-1pt-complex', label: 'Complex w₀: α = 0.5 / (w − (1+i)),  w₀ = 1+i',
    poles: [ { a: '1+i', order: 1, residues: ['0.5'] } ],
    w0: '1+i' },

  { id: 'lqd-3pt-equi', label: 'Equilateral 3-pt around w₀ = 3 (existence depends on residues)',
    poles: [
      { a:  '3.5',                     order: 1, residues: ['0.2'] },
      { a:  '2.75+0.4330127i',         order: 1, residues: ['0.2'] },
      { a:  '2.75-0.4330127i',         order: 1, residues: ['0.2'] },
    ],
    w0: '3' },
];

// ===========================================================================
// LQD presets — BOUNDED SINGULAR
// ---------------------------------------------------------------------------
// q is dialed separately by the user (sliders in #q-card), defaulting to 0.
// Each preset supplies (poles, w0). At q = 0 the family degenerates to a
// "singular LQD with no log-charge at origin" — 0 ∈ Ω but the residue is
// 0; as |q| grows we move along the Theorem 5.6.2 family.
// ===========================================================================
const LQD_PRESETS_BOUNDED_SINGULAR = [
  { id: 'lqd-s-thm-562',
    label: 'Thm 5.6.2 family: 0.5/(w − 2),  w₀ = 1  (dial q with slider)',
    poles: [ { a: '2', order: 1, residues: ['0.5'] } ],
    w0: '1' },

  { id: 'lqd-s-shifted',
    label: 'Shifted Thm 5.6.2: 0.3/(w − 1.5),  w₀ = 0.6  (dial q)',
    poles: [ { a: '1.5', order: 1, residues: ['0.3'] } ],
    w0: '0.6' },

  { id: 'lqd-s-2pt-sym',
    label: 'Two-pt symmetric: 0.4/(w−1) + 0.4/(w+1),  w₀ = 0.5+0.5i  (dial q)',
    poles: [
      { a:  '1', order: 1, residues: ['0.4'] },
      { a: '-1', order: 1, residues: ['0.4'] },
    ],
    w0: '0.5+0.5i' },
];

// ===========================================================================
// LQD presets — UNBOUNDED non-singular
// ---------------------------------------------------------------------------
// Ω is unbounded with 0 ∉ Ω̄. Conformal radius c = φ'(∞) > 0 (slider).
// For h ≡ 0 (no finite poles), φ(z) = cz and Ω = exterior of disk of
// radius c — the trivial unbounded LQD.
// ===========================================================================
const LQD_PRESETS_UNBOUNDED = [
  { id: 'lqd-u-trivial', label: 'Trivial:  h = 0,  c = 0.5  (Ω = ext. disk radius c)',
    poles: [], c: 0.5 },

  { id: 'lqd-u-1pt', label: 'One-pt:  h = 1/(w − 2),  c = 0.6',
    poles: [ { a: '2', order: 1, residues: ['1'] } ],
    c: 0.6 },

  { id: 'lqd-u-1pt-small-c', label: 'One-pt, small c:  h = 1/(w − 2),  c = 0.3',
    poles: [ { a: '2', order: 1, residues: ['1'] } ],
    c: 0.3 },

  { id: 'lqd-u-2pt-sym', label: 'Two-pt symmetric:  1/(w−2) + 0.6/(w+1.5),  c = 0.4',
    poles: [
      { a:  '2',    order: 1, residues: ['1'] },
      { a: '-1.5',  order: 1, residues: ['0.6'] },
    ],
    c: 0.4 },
];

// ===========================================================================
// LQD presets — UNBOUNDED SINGULAR
// ---------------------------------------------------------------------------
// Ω unbounded with both 0 ∈ Ω and ∞ ∈ Ω.  q (complex, dialable) is the
// residue of h at the origin; q = 0 is allowed (singular LQD with no
// log-charge). Default q = 0; dial via slider in #q-card.
// ===========================================================================
const LQD_PRESETS_UNBOUNDED_SINGULAR = [
  { id: 'lqd-us-1pt', label: 'One-pt:  h = q/w + 1/(w−2),  c = 0.6  (dial q)',
    poles: [ { a: '2', order: 1, residues: ['1'] } ],
    c: 0.6 },

  { id: 'lqd-us-2pt-sym', label: 'Two-pt:  q/w + 1/(w−2) + 0.6/(w+1.5),  c = 0.4',
    poles: [
      { a:  '2',    order: 1, residues: ['1'] },
      { a: '-1.5',  order: 1, residues: ['0.6'] },
    ],
    c: 0.4 },
];

// ===========================================================================
// MODE DESCRIPTORS (R5)
// ---------------------------------------------------------------------------
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
// ===========================================================================
const MODES = {
  'bounded': {
    label: 'Bounded QD',
    familyTag:        undefined,           // legacy: untagged phi (boundedQD)
    cards: { w0: true, c: false, poly: false, q: false },
    hint: null,
    presets:          () => QD_PRESETS_BOUNDED,
    externalFieldLabel: 'External field   w − h̄(w)',
    externalFieldKind:  'qd',              // 'qd' = w − h̄;  'lqd' = ln|w|²/w̄ − h̄
    vectorFieldOriginAbs2Floor: 1e-30,     // origin not in Ω, no special clip
    extraHContrib:    null,                // no extra terms beyond polyPart + finite poles
    autoEscalate:     true,
    requireManualW0:  false,
    buildNorm(hData, state) {
      const w0 = buildW0(hData);
      if (w0.error) return w0;
      return { w0: w0.w0 };
    },
    applyNorm(opts, norm) { opts.w0 = norm.w0; },
    warmStartUpdate(initPhi, norm) { initPhi.w0 = { re: norm.w0.re, im: norm.w0.im }; },
  },
  'unbounded': {
    label: 'Unbounded QD',
    familyTag:        undefined,           // legacy: untagged phi (unboundedQD)
    cards: { w0: false, c: true, poly: true, q: false },
    hint: null,
    presets:          () => QD_PRESETS_UNBOUNDED,
    externalFieldLabel: 'External field   w − h̄(w)',
    externalFieldKind:  'qd',
    vectorFieldOriginAbs2Floor: 1e-30,
    extraHContrib:    null,
    autoEscalate:     true,
    requireManualW0:  false,
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
    externalFieldLabel: 'External field   ln|w|²/w̄ − h̄(w)',
    externalFieldKind:  'lqd',
    vectorFieldOriginAbs2Floor: 1e-30,     // 0 ∉ Ω̄, no special clip
    extraHContrib:    null,
    autoEscalate:     false,                // existence is constrained (Thm 5.3.2)
    requireManualW0:  true,                 // w₀ must be ≠ 0
    buildNorm(hData, state) {
      const w0 = buildW0(hData);
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
    externalFieldLabel: 'External field   ln|w|²/w̄ − h̄(w)',
    externalFieldKind:  'lqd',
    vectorFieldOriginAbs2Floor: 1e-30,    // 0 ∈ K, no special clip
    extraHContrib:    null,
    autoEscalate:     false,
    requireManualW0:  false,
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
    externalFieldLabel: 'External field   ln|w|²/w̄ − h̄(w)',
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
    requireManualW0:  false,
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
    externalFieldLabel: 'External field   ln|w|²/w̄ − h̄(w)',
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
    requireManualW0:  true,
    buildNorm(hData, state) {
      const w0 = buildW0(hData);
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

// ===========================================================================
// Aggressiveness presets
// ---------------------------------------------------------------------------
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
// ===========================================================================

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

// ---------- State --------------------------------------------------------
const state = {
  // h(w) data: array of { a: string, order: int, residues: [string,...] }
  poles: [
    { a: '0', order: 1, residues: ['1'] },
  ],
  mode: 'bounded',           // 'bounded' | 'unbounded'
  c: 0.5,                    // conformal radius (unbounded mode only)
  polyDegree: -1,            // polynomial-part degree m_∞ (-1 = none)
  polyCoeffs: [],            // strings, length = max(polyDegree+1, 0)
  w0Mode: 'auto',            // 'auto' | 'manual' (bounded mode only)
  w0Manual: '0',
  q: '0',                    // residue of h at origin (singular LQD mode only)
  aggressiveness: 'standard',
  samples: 500,
  autoFit: false,        // off by default — slider drags shouldn't reframe the view
  vectorFieldMode: 'off',    // 'off' | 'polya' | 'external'
  showCriticalSet: false,    // overlay w-images of {z : φ'(z) = 0}

  // Solver result
  current: null,             // { primary, alternates, w0, attempts }
  selectedSolutionIdx: 0,    // 0 = primary, 1+ = alternate index

  // Debounce + background search bookkeeping
  solveTimer: null,
  altSearchActive: false,
  altSearchToken: 0,

  // Advanced search-options panel. All numeric fields are blank-meaning-
  // "use preset"; toggles are explicit booleans. Mutated by readSearchOptions.
  searchOptions: {
    phases: {
      direct: true, continuation: true, multistart: true,
      diverse: true, deflation: true,
    },
    numRestarts:       null,   // null → preset
    numDiverse:        null,
    numDeflation:      null,
    bgChunks:          null,
    bgChunkSize:       null,
    keepSearching:     false,

    newtonMaxIter:     null,
    newtonTol:         null,
    contTStart:        null,
    contGrow:          null,

    deflationAlpha:    null,   // null → 1
    deflationP:        null,   // null → 2
    deflateFromValid:  false,

    univalenceSamples: null,   // null → state.samples
    identityTol:       null,   // null → 1e-6
    showNonUnivalent:  false,
    showIdFailing:     false,
    autoEscalate:      true,

    seed:              null,   // null → time-based
  },

  // View mode (HANDOFF #30): the former "Direct problem" tab is folded
  // into this tab as a segmented toggle. 'inverse' shows the existing
  // QD/LQD UI; 'direct' shows the relocated #controls-direct content
  // (mounted lazily by QD.Direct._mountUI on first switch).
  viewMode:       'inverse',     // 'inverse' | 'direct'
  directMounted:  false,
};

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

// ---------- Render the pole controls ------------------------------------
function renderPolesList() {
  const list = $('#poles-list');
  list.innerHTML = '';

  state.poles.forEach((pole, idx) => {
    const div = document.createElement('div');
    div.className = 'pole';
    div.dataset.idx = idx;
    div.innerHTML = `
      <div class="pole-header">
        <span class="pole-num">Pole ${idx + 1}</span>
        <button type="button" class="small danger" data-action="remove" title="Remove this pole">×</button>
      </div>
      <div class="row">
        <label>a${sub(idx+1)} =
          <input type="text" class="cnum" data-field="a" value="${escapeAttr(pole.a)}">
        </label>
      </div>
      <div class="row">
        <label>Order:
          <input type="number" min="1" max="6" value="${pole.order}" data-field="order" style="width: 56px;">
        </label>
      </div>
      <div class="residues"></div>
    `;
    const residuesEl = $('.residues', div);
    for (let s = 0; s < pole.order; s++) {
      const cval = Complex.parse(pole.residues[s] || '0') || { re: 0, im: 0 };
      const mag = Math.hypot(cval.re, cval.im);
      const arg = Math.atan2(cval.im, cval.re);
      const key = residueKey(idx, s);
      const magMax = magMaxFor(key, mag);

      const block = document.createElement('div');
      block.className = 'residue-block';
      block.dataset.s = s;
      block.innerHTML = `
        <div class="residue-row">
          <span class="label-fixed">C${sub(idx+1)}${sub(s+1)}</span>
          =
          <input type="text" class="cnum residue" data-field="residue" data-s="${s}" value="${escapeAttr(pole.residues[s] || '')}">
        </div>
        <div class="slider1d-row">
          <label>|C|</label>
          <input type="range" class="slider1d slider1d-mag" data-s="${s}"
                 min="0" max="${magMax}" step="any" value="${mag}">
          <span class="slider1d-val mag-val">${mag.toFixed(3)}</span>
        </div>
        <div class="slider1d-row">
          <label>arg</label>
          <input type="range" class="slider1d slider1d-arg" data-s="${s}"
                 min="${-Math.PI}" max="${Math.PI}" step="any" value="${arg}">
          <span class="slider1d-val arg-val">${fmtArg(arg)}</span>
        </div>
      `;
      residuesEl.appendChild(block);
    }
    list.appendChild(div);
  });
  if (typeof refreshHText === 'function') refreshHText();
}

// Render the polynomial-part coefficient list. One block per C_{∞,l} for
// l = 0..polyDegree, with magnitude/argument sliders matching the residue
// rows. Visible in any mode where polynomial-h is meaningful (classical
// unbounded + both unbounded-LQD variants — see modeAllowsPoly).
function renderPolyCoefList() {
  const list = $('#poly-coefs-list');
  if (!list) return;
  list.innerHTML = '';
  const deg = state.polyDegree;
  if (!modeAllowsPoly(state.mode) || deg < 0) return;

  // Ensure polyCoeffs has at least deg+1 entries (pad with '0').
  while (state.polyCoeffs.length < deg + 1) state.polyCoeffs.push('0');
  state.polyCoeffs.length = deg + 1;          // truncate any extras

  for (let l = 0; l <= deg; l++) {
    const cval = QD.Complex.parse(state.polyCoeffs[l] || '0') || { re: 0, im: 0 };
    const mag = Math.hypot(cval.re, cval.im);
    const arg = Math.atan2(cval.im, cval.re);
    const key = `poly-coef-${l}`;
    const magMax = magMaxFor(key, mag);
    const block = document.createElement('div');
    block.className = 'residue-block';
    block.dataset.polyL = l;
    block.innerHTML = `
      <div class="residue-row">
        <span class="label-fixed">C<sub>∞,${l}</sub></span>
        =
        <input type="text" class="cnum poly-coef" data-poly-l="${l}" value="${escapeAttr(state.polyCoeffs[l] || '')}">
      </div>
      <div class="slider1d-row">
        <label>|C|</label>
        <input type="range" class="slider1d slider1d-poly-mag" data-poly-l="${l}"
               min="0" max="${magMax}" step="any" value="${mag}">
        <span class="slider1d-val poly-mag-val">${mag.toFixed(3)}</span>
      </div>
      <div class="slider1d-row">
        <label>arg</label>
        <input type="range" class="slider1d slider1d-poly-arg" data-poly-l="${l}"
               min="${-Math.PI}" max="${Math.PI}" step="any" value="${arg}">
        <span class="slider1d-val poly-arg-val">${fmtArg(arg)}</span>
      </div>
    `;
    list.appendChild(block);
  }
  if (typeof refreshHText === 'function') refreshHText();
}

function escapeAttr(s) { return String(s).replace(/&/g, '&amp;').replace(/"/g, '&quot;'); }

// ---------- Read controls into state -------------------------------------
function readPolesFromDOM() {
  state.poles = $$('.pole').map(div => {
    const idx = +div.dataset.idx;
    const a = $('input[data-field="a"]', div).value;
    const order = Math.max(1, Math.min(6, +$('input[data-field="order"]', div).value || 1));
    const residues = $$('input[data-field="residue"]', div).map(inp => inp.value);
    return { a, order, residues };
  });
}

// Highlight invalid inputs that failed to parse.
function markInvalid(errorMsg) {
  $$('.cnum').forEach(inp => inp.classList.remove('invalid'));
  if (!errorMsg) return;
  // We don't currently know which input is bad; mark them all softly. The
  // status panel surfaces the message clearly.
}

// ---------- Adding / removing poles --------------------------------------
function addPole() {
  state.poles.push({ a: '0', order: 1, residues: ['1'] });
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

// ---------- Solving (debounced) ------------------------------------------
const scheduleSolve = debounce(() => { solveAndRender(); }, 250);

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

function quickSolveAndRender() {
  const built = buildHData();
  if (!built || built.error) return;
  const norm = buildNormalization(built);
  if (norm.error) return;

  const preset = PRESETS[state.aggressiveness];
  const unbounded = !!norm.unbounded;
  const desc = modeDescriptor();
  const expectedFamilyTag = desc.familyTag;

  // Select the Family backend so the rest of the quick-solve path is
  // family-agnostic. Warm-start when the previous solution matches family
  // tag, branch structure, and bounded/unbounded mode; otherwise fresh init.
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

  let result;
  try {
    // newtonSolve auto-dispatches Family from initPhi.family — works for QD/LQD.
    result = QD.newtonSolve(initPhi, built, { ...preset.newton, maxIter: 30 });
  } catch (e) { return; }

  if (!result.success) {
    scheduleSolve();   // out of warm-start basin → fall back to debounced full solve
    return;
  }

  const phi = family.canonicalizePhi(result.phi);
  const sol = {
    ...result,
    phi,
    method: 'live',
    univalent: QD.isBoundaryUnivalent(phi, state.samples),
    identity:  family.verifyQuadratureIdentity(phi, built, { numSamples: state.samples }),
  };
  sol.identityOK = sol.identity.maxRelDiff < 1e-6;

  state.current = {
    success: true,
    primary: sol,
    alternates: [],
    hData: built,
    w0Used: norm.w0,
    cUsed:  norm.c,
    unbounded,
    attempts: [],
  };
  state.selectedSolutionIdx = 0;

  showSolution(sol, built, /* autoFit = */ false);
  refreshAlternatesPanel();
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

  state.altSearchToken++;
  state.altSearchActive = false;
  $('#alt-search-indicator').classList.add('hidden');

  setTimeout(() => {
    let result;
    try {
      const opts = buildSolverOptions(preset, { findAlternates: false });
      applyNormToOpts(opts, norm);
      result = QD.solveInverseQD(built, opts);

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
        const exh = buildSolverOptions(PRESETS.exhaustive, { findAlternates: false });
        applyNormToOpts(exh, norm);
        const escalated = QD.solveInverseQD(built, exh);
        if (escalated.success) result = escalated;
      }
    } catch (e) {
      setStatus({ kind: 'err', text: 'Solver error: ' + e.message });
      return;
    }
    state.current = result;
    state.current.hData = built;
    state.current.w0Used = norm.w0;
    state.current.cUsed  = norm.c;
    state.current.unbounded = !!norm.unbounded;
    state.selectedSolutionIdx = 0;

    if (!result.success) {
      setStatus({
        kind: 'err',
        text: 'No solution found.\n' +
              '  reason: ' + result.error + '\n' +
              '  attempts: ' + (result.attempts ? result.attempts.length : 0),
      });
      plot.clear();
      $('#alternates-card').classList.add('hidden');
      $('#riemann-map-card').classList.add('hidden');
      // Try-harder button is always visible — nothing to toggle here.
      return;
    }

    showSolution(result.primary, built, /*autoFit=*/true);
    refreshAlternatesPanel();

    startBackgroundAltSearch(built, norm);
  }, 0);
}

// ---------- Display a chosen solution on the plot ------------------------
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
  const valid = sol.univalent && sol.identityOK;
  if (valid) {
    lines.push(`<span class="ok">✓ Valid quadrature domain</span>`);
  } else if (!sol.univalent && !sol.identityOK) {
    lines.push(`<span class="err">✗ Spurious algebraic root (non-univalent AND identity fails)</span>`);
  } else if (!sol.univalent) {
    lines.push(`<span class="warn">⚠ Algebraic root: boundary self-intersects (non-univalent)</span>`);
  } else {
    lines.push(`<span class="warn">⚠ Algebraic root: quadrature identity not satisfied (likely spurious)</span>`);
  }
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

// ---------- Riemann-map formula card ------------------------------------
// Renders (1) symbolic identity, (2) closed-form expression with values
// substituted, (3) parameters table.
//
// Per-family pieces (R6): RIEMANN_FRAGMENTS keys to family tag (undefined =
// legacy QD/UQD). Each fragment declares:
//   symbolic(sumBody, unb, mInf)        → LaTeX for (1)
//   numericLeader(phi)                   → array of LaTeX rows preceding the Σ
//   numericTrailer(phi)                  → array of LaTeX rows after the Σ
//   extraParameterRows(phi)              → [{name, value}] inserted after w_0
//   stripFirstPlus                       → bool: should leading "+\," be
//                                          dropped from the first sum row?
const RIEMANN_FRAGMENTS = {
  // Bounded classical QD (legacy untagged phi or family === 'boundedQD'):
  // φ(z) = w_0 + Σ … . No leader transforms; sum is appended directly.
  '_boundedQD': {
    symbolic: (sumBody) => String.raw`\varphi(z) \;=\; w_0 \;+\; ${sumBody}`,
    numericLeader: (phi) => [String.raw`\varphi(z) \;\approx\; ${katexCmpx(phi.w0)}`],
    numericTrailer: () => [],
    extraParameterRows: () => [],
    stripFirstPlus: false,
  },
  // Unbounded classical QD:
  // φ(z) = c·z + Σ_l F_l / z^l + Σ … .
  '_unboundedQD': {
    symbolic: (sumBody, unb, mInf) => {
      const polySym = mInf >= 0
        ? String.raw` \;+\; \sum_{l=0}^{m_\infty} \frac{F_{l}}{z^{l}}`
        : '';
      return String.raw`\varphi(z) \;=\; c\, z${polySym} \;+\; ${sumBody}`;
    },
    numericLeader: (phi) => {
      const rows = [String.raw`\varphi(z) \;\approx\; ${Number(phi.c.toFixed(5))}\, z`];
      const polyA = phi.polyA || [];
      for (let l = 0; l < polyA.length; l++) {
        const F = polyA[l];
        if (l === 0) {
          rows.push(String.raw`+\, ${katexCmpxParen(F)}`);
        } else {
          const zl = l === 1 ? 'z' : `z^{${l}}`;
          rows.push(String.raw`+\, \dfrac{${katexCmpxParen(F)}}{${zl}}`);
        }
      }
      return rows;
    },
    numericTrailer: () => [],
    extraParameterRows: () => [],
    stripFirstPlus: false,
  },
  // Bounded non-singular LQD: φ(z) = w_0 · exp(Σ …).
  'boundedLQD': {
    symbolic: (sumBody) =>
      String.raw`\varphi(z) \;=\; w_0 \cdot \exp\!\Bigl(${sumBody}\Bigr)`,
    numericLeader: (phi) => [
      String.raw`\varphi(z) \;\approx\; ${katexCmpx(phi.w0)} \cdot \exp\Bigl(`,
    ],
    numericTrailer: () => [String.raw`\Bigr)`],
    extraParameterRows: () => [],
    stripFirstPlus: true,
  },
  // Unbounded non-singular LQD: φ(z) = c·z · exp(r#(z) − r#(∞)) on 𝔻*.
  // The (− r#(∞)) absorbs the ∞-gauge so the displayed leading coefficient
  // really is c. The numerical expansion shows the actual A_{j,k} the solver
  // returned and the additional − r#(∞) term (computed in closed form).
  'unboundedLQD': {
    symbolic: (sumBody) =>
      String.raw`\varphi(z) \;=\; c\, z \cdot \exp\!\Bigl(${sumBody} \;-\; r_\#(\infty)\Bigr)`,
    numericLeader: (phi) => [
      String.raw`\varphi(z) \;\approx\; ${Number(phi.c.toFixed(5))}\, z \cdot \exp\Bigl(`,
    ],
    numericTrailer: (phi) => [
      String.raw`-\, ${katexCmpxParen(rHashAtInfinityForDisplay(phi))}`,
      String.raw`\Bigr)`,
    ],
    extraParameterRows: (phi) => [
      { name: String.raw`r_\#(\infty)`, value: rHashAtInfinityForDisplay(phi) },
    ],
    stripFirstPlus: true,
  },
  // Unbounded singular LQD: φ(z) = c·|z₀|·z·b_{z₀}(z)·exp(r#(z) − r#(∞)).
  'unboundedLQD_singular': {
    symbolic: (sumBody) =>
      String.raw`\varphi(z) \;=\; c\cdot|z_0|\cdot z\cdot b_{z_0}(z) \cdot \exp\!\Bigl(${sumBody} \;-\; r_\#(\infty)\Bigr),\quad b_{z_0}(z) = -\tfrac{\overline{z_0}}{|z_0|}\cdot\tfrac{z - z_0}{1 - \overline{z_0}\, z}`,
    numericLeader: (phi) => {
      const absZ0 = QD.Complex.abs(phi.z0).toFixed(5);
      const z0Latex = katexCmpxParen(phi.z0);
      return [
        String.raw`\varphi(z) \;\approx\; ${Number(phi.c.toFixed(5))}\cdot ${absZ0}\cdot z\cdot b_{${z0Latex}}(z) \cdot \exp\Bigl(`,
      ];
    },
    numericTrailer: (phi) => [
      String.raw`-\, ${katexCmpxParen(rHashAtInfinityForDisplay(phi))}`,
      String.raw`\Bigr)`,
    ],
    extraParameterRows: (phi) => [
      { name: String.raw`z_0`,        value: phi.z0 },
      { name: String.raw`q`,          value: phi.q },
      { name: String.raw`r_\#(\infty)`, value: rHashAtInfinityForDisplay(phi) },
    ],
    stripFirstPlus: true,
  },
  // Bounded singular LQD: φ(z) = γ · b_{z_0}(z) · exp(Σ …).
  'boundedLQD_singular': {
    symbolic: (sumBody) =>
      String.raw`\varphi(z) \;=\; \gamma \cdot b_{z_0}(z) \cdot \exp\!\Bigl(${sumBody}\Bigr),`
      + String.raw`\quad b_{z_0}(z) = -\tfrac{\overline{z_0}}{|z_0|}\cdot\tfrac{z - z_0}{1 - \overline{z_0}\, z}`,
    numericLeader: (phi) => [
      String.raw`\varphi(z) \;\approx\; ${katexCmpxParen(phi.gamma)} \cdot b_{${katexCmpxParen(phi.z0)}}(z) \cdot \exp\Bigl(`,
    ],
    numericTrailer: () => [String.raw`\Bigr)`],
    extraParameterRows: (phi) => [
      { name: String.raw`z_0`,    value: phi.z0 },
      { name: String.raw`\gamma`, value: phi.gamma },
      { name: String.raw`q`,      value: phi.q },
    ],
    stripFirstPlus: true,
  },
};

function getRiemannFragment(phi) {
  if (phi.family && RIEMANN_FRAGMENTS[phi.family]) return RIEMANN_FRAGMENTS[phi.family];
  // Legacy: family tag absent → QD/UQD by phi.unbounded.
  return phi.unbounded ? RIEMANN_FRAGMENTS['_unboundedQD'] : RIEMANN_FRAGMENTS['_boundedQD'];
}

// r#(∞) — delegate to QD.LqdCommon.rHashAtInfinity (defined in
// solver-lqd-common.js). Kept under this local name so the
// RIEMANN_FRAGMENTS table can reference it without a per-call namespace
// lookup.
function rHashAtInfinityForDisplay(phi) {
  return QD.LqdCommon.rHashAtInfinity(phi);
}

function renderRiemannMap(phi) {
  const card = $('#riemann-map-card');
  const content = $('#riemann-map-content');
  if (!phi) { card.classList.add('hidden'); return; }
  card.classList.remove('hidden');

  const maxOrder = phi.branches.reduce((m, b) => Math.max(m, b.A.length), 1);
  const unb = !!phi.unbounded;
  const polyA = phi.polyA || [];
  const m_inf = polyA.length - 1;
  const frag = getRiemannFragment(phi);

  // --- (1) symbolic identity ---
  const sumBody = maxOrder === 1
    ? String.raw`\sum_{j} \overline{A_{j}}\, \frac{z}{1 - \overline{z_{j}}\, z}`
    : String.raw`\sum_{j,k} \overline{A_{j,k}}\, \frac{z^{k}}{\bigl(1 - \overline{z_{j}}\, z\bigr)^{k}}`;
  const sym = frag.symbolic(sumBody, unb, m_inf);

  // --- (2) closed-form expression with values substituted ---
  // Build the additive sum-of-branches rows. Each row starts with "+\, ".
  const sumTermsRows = [];
  phi.branches.forEach((br) => {
    const zjC = QD.Complex.conj(br.z);
    br.A.forEach((a, k) => {
      const aC = QD.Complex.conj(a);
      const power = k + 1;
      const numerator = power === 1 ? 'z' : `z^{${power}}`;
      const denomCore = String.raw`1 - ${katexCmpxParen(zjC)}\, z`;
      const denom = power === 1 ? `\\bigl(${denomCore}\\bigr)` : `\\bigl(${denomCore}\\bigr)^{${power}}`;
      sumTermsRows.push(String.raw`+\, ${katexCmpxParen(aC)}\, \dfrac{${numerator}}{${denom}}`);
    });
  });

  // Assemble: leader rows + (optionally-de-plussed first sum term + rest) +
  // trailer rows. Families like LQDs that wrap the sum in exp(…) strip the
  // leading "+\," from the first term.
  const numerRows = frag.numericLeader(phi);
  if (sumTermsRows.length === 0 && frag.stripFirstPlus) {
    numerRows.push('0');                 // empty sum → write a literal 0
  } else {
    for (let i = 0; i < sumTermsRows.length; i++) {
      const row = (i === 0 && frag.stripFirstPlus)
        ? sumTermsRows[i].replace(/^\+\\,\s*/, '')
        : sumTermsRows[i];
      numerRows.push(row);
    }
  }
  for (const row of frag.numericTrailer(phi)) numerRows.push(row);
  // First row plain; subsequent rows ampersand-aligned at the leading "+".
  const numerLatex =
    String.raw`\begin{aligned}` +
    numerRows.map((row, i) => i === 0 ? row : `& ${row}`).join(' \\\\[2pt] ') +
    String.raw`\end{aligned}`;

  // --- (3) parameters table ---
  const rows = [];
  if (unb) {
    rows.push({ name: String.raw`c`, value: { re: phi.c, im: 0 } });
  } else {
    rows.push({ name: String.raw`w_0`, value: phi.w0 });
  }
  for (const r of frag.extraParameterRows(phi)) rows.push(r);
  phi.branches.forEach((br, j) => {
    rows.push({ name: String.raw`z_{${j + 1}}`, value: br.z });
    br.A.forEach((a, k) => {
      const sub = br.A.length === 1 ? `${j + 1}` : `${j + 1},${k + 1}`;
      rows.push({ name: String.raw`A_{${sub}}`, value: a });
    });
  });
  // F_l rows (polynomial part of φ at ∞).
  if (unb && m_inf >= 0) {
    for (let l = 0; l <= m_inf; l++) {
      rows.push({ name: String.raw`F_{${l}}`, value: polyA[l] });
    }
  }
  const tableRows = rows.map((r, i) =>
    `<tr><td class="rm-name" data-tex="${escapeAttr(r.name)}"></td>
         <td class="rm-value">${escapeHTML(QD.Complex.toString(r.value, 5))}</td></tr>`
  ).join('');

  content.innerHTML = `
    <div class="rm-sym"></div>
    <details class="rm-section" open>
      <summary>Closed-form expression</summary>
      <div class="rm-numer"></div>
    </details>
    <details class="rm-section">
      <summary>Parameters</summary>
      <table class="rm-params"><tbody>${tableRows}</tbody></table>
    </details>
  `;

  renderKatex($('.rm-sym',   content), sym,        true);
  renderKatex($('.rm-numer', content), numerLatex, true);
  $$('td.rm-name', content).forEach(td => {
    renderKatex(td, td.dataset.tex, false);
  });
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

// Format a complex number as LaTeX-safe text (no math-mode commands; pure
// real/imag/general handled naturally).
function katexCmpx(c) {
  const re = c.re, im = c.im;
  const fmt = (x) => Number(x.toFixed(5)).toString();
  if (Math.abs(im) < 1e-12) return fmt(re);
  if (Math.abs(re) < 1e-12) {
    if (Math.abs(im - 1) < 1e-12) return 'i';
    if (Math.abs(im + 1) < 1e-12) return '-i';
    return fmt(im) + 'i';
  }
  const sign = im >= 0 ? ' + ' : ' - ';
  const iAbs = Math.abs(im);
  const iPart = Math.abs(iAbs - 1) < 1e-12 ? '' : fmt(iAbs);
  return fmt(re) + sign + iPart + 'i';
}

// Wrap in \left(...\right) unless the value is a bare non-negative real.
function katexCmpxParen(c) {
  const s = katexCmpx(c);
  if (/^\d+(\.\d+)?$/.test(s)) return s;
  return String.raw`\left(${s}\right)`;
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

  function tick() {
    if (myToken !== state.altSearchToken) return;       // superseded
    if (!state.current || !state.current.success) {
      state.altSearchActive = false;
      $('#alt-search-indicator').classList.add('hidden');
      refreshAlternatesPanel();
      return;
    }
    if (!keepGoing && chunk >= bgChunks) {
      state.altSearchActive = false;
      $('#alt-search-indicator').classList.add('hidden');
      refreshAlternatesPanel();
      return;
    }
    chunk++;

    let found = [];
    try {
      const known = [state.current.primary, ...(state.current.alternates || [])];
      found = QD.searchAlternates(hData, norm, known,
        buildAltSearchOptions(preset, seed));
    } catch (e) {
      console.warn('alt search error:', e);
    }
    seed = (seed * 1664525 + 1013904223) >>> 0;

    if (found.length > 0) {
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
        refreshAlternatesPanel();
      }
    }

    setTimeout(tick, 30);   // yield
  }

  setTimeout(tick, 80);
}

// ===========================================================================
// Canvas plotting
// ===========================================================================

// Hover hit-test radius for the pole-proximity annotation in the readout
// (HANDOFF #33 / #35). Larger than the click hit-radius (9 px, the default
// `radius` parameter on _hitTestPole) so the cursor doesn't need to be
// pixel-perfect over the pole dot to see its label.
const POLE_HOVER_HIT_RADIUS_PX = 12;

class DomainPlot {
  constructor(canvas, readout) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.readout = readout;
    this.dpr = window.devicePixelRatio || 1;

    this.view = { cx: 0, cy: 0, scale: 100 };  // pixels per unit
    this.data = null;

    // Callbacks for click-drag on quadrature-node dots. Set by ui.js.
    this.onPoleDrag    = null;   // (idx, worldPoint) -> void
    this.onPoleDragEnd = null;   // (idx)            -> void

    this.attachEvents();
    this.resize();
  }

  // Returns the index of the pole dot under (x, y) in CSS pixels, or -1 if
  // none is within the hit-test radius.
  _hitTestPole(x, y, radius = 9) {
    if (!this.data || !this.data.poles) return -1;
    let bestI = -1, bestD2 = radius * radius;
    for (let i = 0; i < this.data.poles.length; i++) {
      const sp = this.toScreen(this.data.poles[i].re, this.data.poles[i].im);
      const dx = sp.x - x, dy = sp.y - y;
      const d2 = dx * dx + dy * dy;
      if (d2 <= bestD2) { bestD2 = d2; bestI = i; }
    }
    return bestI;
  }

  resize() {
    const rect = this.canvas.getBoundingClientRect();
    this.cssW = rect.width;
    this.cssH = rect.height;
    this.canvas.width  = Math.round(this.cssW * this.dpr);
    this.canvas.height = Math.round(this.cssH * this.dpr);
    this.render();
  }

  toScreen(re, im) {
    return {
      x: this.cssW / 2 + (re - this.view.cx) * this.view.scale,
      y: this.cssH / 2 - (im - this.view.cy) * this.view.scale,
    };
  }
  toWorld(x, y) {
    return {
      re: this.view.cx + (x - this.cssW / 2) / this.view.scale,
      im: this.view.cy - (y - this.cssH / 2) / this.view.scale,
    };
  }

  attachEvents() {
    let panning = false, lastX = 0, lastY = 0;
    let draggingPole = -1;          // index, or -1

    // Belt-and-suspenders gate: the QD/LQD inverse tab is the only one that
    // owns the main #canvas as a 2D drawing surface; the Schwarz and Riemann-
    // sphere tabs overlay their own GL canvases.  These handlers must early-
    // return when those other tabs are active, otherwise a drag/wheel in
    // those tabs would trigger this.render() and repaint the QD plot over
    // the other tab's content.  (The sphere tab additionally puts its GL
    // canvas on top with pointer-events:auto so most events never reach
    // this listener — this check is the second line of defense, and the
    // only line of defense for the Schwarz tab whose GL canvas sits below.)
    //
    // Exception: an in-progress pan or pole-drag started on the QD tab is
    // allowed to complete even if the user mid-drag switched tabs (avoids
    // a stuck "panning=true" state).
    function _qdTabActive() {
      const btn = document.querySelector('.tab-btn.active');
      return btn && btn.dataset.tab === 'qd';
    }

    this.canvas.addEventListener('mousedown', e => {
      if (!_qdTabActive()) return;
      const rect = this.canvas.getBoundingClientRect();
      const x = e.clientX - rect.left, y = e.clientY - rect.top;
      // Hit-test a pole first; if hit, take ownership of this drag for the
      // pole rather than starting a pan.
      const hit = this._hitTestPole(x, y);
      if (hit >= 0) {
        draggingPole = hit;
        this.canvas.style.cursor = 'grabbing';
        e.preventDefault();
        return;
      }
      panning = true; lastX = e.clientX; lastY = e.clientY;
      this.canvas.style.cursor = 'grabbing';
    });

    window.addEventListener('mousemove', e => {
      // Allow an already-in-progress drag/pole-drag to finish even after a
      // tab switch.  Brand-new hover/readout updates require the QD tab.
      const inProgressDrag = panning || draggingPole >= 0;
      if (!inProgressDrag && !_qdTabActive()) return;

      const rect = this.canvas.getBoundingClientRect();
      const x = e.clientX - rect.left, y = e.clientY - rect.top;

      if (draggingPole >= 0) {
        const w = this.toWorld(x, y);
        if (this.onPoleDrag) this.onPoleDrag(draggingPole, w);
        return;
      }
      if (panning) {
        const dx = e.clientX - lastX, dy = e.clientY - lastY;
        this.view.cx -= dx / this.view.scale;
        this.view.cy += dy / this.view.scale;
        lastX = e.clientX; lastY = e.clientY;
        this.render();
      }

      // Mouse-coordinate readout + hover cursor (pointer over a pole, grab
      // elsewhere) — only when the cursor is over the canvas.
      if (x >= 0 && x <= this.cssW && y >= 0 && y <= this.cssH) {
        const w = this.toWorld(x, y);
        let text = `w = ${w.re.toFixed(4)} ${w.im >= 0 ? '+' : '-'} ${Math.abs(w.im).toFixed(4)}i`;
        // Append nearby-pole annotation when the cursor is within the
        // hit-test radius (HANDOFF #33). Useful for quickly identifying
        // which pole a residue belongs to as the user scans.
        const hitIdx = this._hitTestPole(x, y, POLE_HOVER_HIT_RADIUS_PX);
        if (hitIdx >= 0) {
          const a = this.data.poles[hitIdx];
          const SUBS = ['₀','₁','₂','₃','₄','₅','₆','₇','₈','₉'];
          const sub = String(hitIdx + 1).split('').map(d => SUBS[+d] || d).join('');
          text += `  ·  near pole a${sub} = ${a.re.toFixed(3)}${a.im >= 0 ? '+' : '-'}${Math.abs(a.im).toFixed(3)}i`;
        }
        this.readout.textContent = text;
        if (!panning && draggingPole < 0) {
          this.canvas.style.cursor = hitIdx >= 0 ? 'pointer' : 'grab';
        }
      }
    });

    window.addEventListener('mouseup', () => {
      if (draggingPole >= 0) {
        const idx = draggingPole;
        draggingPole = -1;
        this.canvas.style.cursor = 'grab';
        if (this.onPoleDragEnd) this.onPoleDragEnd(idx);
        return;
      }
      panning = false;
      this.canvas.style.cursor = 'grab';
    });

    this.canvas.addEventListener('wheel', e => {
      if (!_qdTabActive()) return;
      e.preventDefault();
      const rect = this.canvas.getBoundingClientRect();
      const x = e.clientX - rect.left, y = e.clientY - rect.top;
      const wBefore = this.toWorld(x, y);
      const factor = Math.exp(-e.deltaY * 0.001);
      this.view.scale = Math.max(1e-3, Math.min(1e7, this.view.scale * factor));
      const wAfter = this.toWorld(x, y);
      this.view.cx += wBefore.re - wAfter.re;
      this.view.cy += wBefore.im - wAfter.im;
      this.render();
    }, { passive: false });
  }

  setData(d) {
    this.data = d;
    this.render();
  }
  clear() {
    this.data = null;
    this.render();
  }

  reset() {
    this.view = { cx: 0, cy: 0, scale: 100 };
    this.render();
  }

  fit() {
    if (!this.data || !this.data.boundaryPts || this.data.boundaryPts.length === 0) return;
    let minRe = Infinity, maxRe = -Infinity, minIm = Infinity, maxIm = -Infinity;
    for (const p of this.data.boundaryPts) {
      if (p.re < minRe) minRe = p.re; if (p.re > maxRe) maxRe = p.re;
      if (p.im < minIm) minIm = p.im; if (p.im > maxIm) maxIm = p.im;
    }
    for (const p of this.data.poles) {
      if (p.re < minRe) minRe = p.re; if (p.re > maxRe) maxRe = p.re;
      if (p.im < minIm) minIm = p.im; if (p.im > maxIm) maxIm = p.im;
    }
    const dx = Math.max(1e-6, maxRe - minRe);
    const dy = Math.max(1e-6, maxIm - minIm);
    const pad = 0.15;
    const sx = this.cssW / (dx * (1 + 2*pad));
    const sy = this.cssH / (dy * (1 + 2*pad));
    this.view.scale = Math.min(sx, sy);
    this.view.cx = (minRe + maxRe) / 2;
    this.view.cy = (minIm + maxIm) / 2;
    this.render();
  }

  render() {
    const c = this.ctx;
    c.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    c.clearRect(0, 0, this.cssW, this.cssH);
    c.fillStyle = '#fafafa';
    c.fillRect(0, 0, this.cssW, this.cssH);

    this.drawGrid();
    this.drawAxes();

    // Vector field underlays the domain so the boundary remains crisp.
    if (state.vectorFieldMode !== 'off' && this.data && this.data.hData) {
      this.drawVectorField();
    }

    if (this.data && this.data.boundaryPts && this.data.boundaryPts.length > 0) {
      this.drawBoundary();
    }
    // Optional dashed overlay (e.g. for Direct-tab round-trip diagnostics).
    if (this.data && this.data.overlayBoundary && this.data.overlayBoundary.length > 0) {
      this.drawOverlayBoundary();
    }
    if (this.data && this.data.poles)  this.drawPoles();
    if (this.data && this.data.w0)     this.drawW0();
    // Critical-set image overlay (zeros of φ', mapped to w-plane).  Drawn
    // last so the markers sit on top of everything; lazy-computed and
    // cached on state.current.criticalSet to avoid recomputing per render.
    if (state.showCriticalSet && this.data && this.data.phi) {
      this.drawCriticalSet();
    }
  }

  drawOverlayBoundary() {
    const c = this.ctx;
    const pts = this.data.overlayBoundary;
    c.save();
    c.beginPath();
    const p0 = this.toScreen(pts[0].re, pts[0].im);
    c.moveTo(p0.x, p0.y);
    for (let i = 1; i < pts.length; i++) {
      const p = this.toScreen(pts[i].re, pts[i].im);
      c.lineTo(p.x, p.y);
    }
    c.closePath();
    c.strokeStyle = '#d4a017';                            // gold
    c.lineWidth = 1.8;
    if (c.setLineDash) c.setLineDash([6, 4]);
    c.stroke();
    c.restore();
  }

  // ----- Vector-field overlay: V(w) = conj(h(w)) ---------------------------
  // Samples h on a grid that is anchored to WORLD coordinates (multiples of
  // a "nice" step in the w-plane), so that as the user pans/zooms each arrow
  // stays glued to a specific point in the domain. The step is chosen so
  // that adjacent grid points are roughly 26 px apart on screen at the
  // current zoom, snapped to a 1 / 2 / 5 / 10 sequence times a power of 10.
  //
  // Each arrow's length and opacity scale by tanh(log10(1 + |h|)) so a few
  // large values near poles don't dominate the visual.
  //
  // Screen-direction note: conj(h) = Re(h) - i Im(h) is the vector
  // (Re h, -Im h) in the math plane. Screen y is flipped relative to the
  // imaginary axis, so the screen-space arrow direction is just (Re h, Im h).
  drawVectorField() {
    const c = this.ctx;
    const hData = this.data.hData;
    if (!hData) return;
    // h(w) = poly part + sum of principal parts. Bail only if BOTH are empty;
    // a poly-only unbounded h still has a vector field worth drawing.
    const polyLen = (hData.polyPart && hData.polyPart.length) || 0;
    if (hData.poles.length === 0 && polyLen === 0) return;

    // Pick a nice world-coordinate step targeting ~26 px on screen.
    const targetPx = 26;
    const targetWorld = targetPx / this.view.scale;
    const exp10 = Math.floor(Math.log10(targetWorld));
    const frac = targetWorld / Math.pow(10, exp10);
    let mult;
    if      (frac < 1.5) mult = 1;
    else if (frac < 3)   mult = 2;
    else if (frac < 7)   mult = 5;
    else                 mult = 10;
    const worldStep = mult * Math.pow(10, exp10);
    const stepPx = worldStep * this.view.scale;
    const arrowMaxLen = stepPx * 0.72;

    // Visible world bounds (screen y is flipped vs imaginary axis).
    const wTL = this.toWorld(0, 0);
    const wBR = this.toWorld(this.cssW, this.cssH);
    const minRe = Math.min(wTL.re, wBR.re);
    const maxRe = Math.max(wTL.re, wBR.re);
    const minIm = Math.min(wTL.im, wBR.im);
    const maxIm = Math.max(wTL.im, wBR.im);
    const iMin = Math.floor(minRe / worldStep);
    const iMax = Math.ceil (maxRe / worldStep);
    const jMin = Math.floor(minIm / worldStep);
    const jMax = Math.ceil (maxIm / worldStep);

    c.lineWidth = 1;
    c.lineCap = 'round';

    for (let i = iMin; i <= iMax; i++) {
      for (let j = jMin; j <= jMax; j++) {
        const wRe = i * worldStep;
        const wIm = j * worldStep;
        const screen = this.toScreen(wRe, wIm);
        const sx = screen.x;
        const sy = screen.y;
        // Skip if the grid point is well off-canvas.
        if (sx < -arrowMaxLen || sx > this.cssW + arrowMaxLen ||
            sy < -arrowMaxLen || sy > this.cssH + arrowMaxLen) continue;
        const w = { re: wRe, im: wIm };

        // Skip arrows that fall too close to any pole in screen-space; the
        // field diverges and the visual is meaningless there.
        let skip = false;
        for (const p of hData.poles) {
          const sp = this.toScreen(p.a.re, p.a.im);
          if (Math.hypot(sp.x - sx, sp.y - sy) < 9) { skip = true; break; }
        }
        if (skip) continue;

        // Evaluate h(w) = Σ_{l=0..m_∞} C_{∞,l} w^l + Σ_j Σ_s C_{j,s}/(w − a_j)^s
        let hr = 0, hi = 0;
        let bad = false;
        // Polynomial part of h (unbounded only; otherwise polyPart is empty).
        const polyPart = hData.polyPart || [];
        if (polyPart.length > 0) {
          let wPowRe = 1, wPowIm = 0;                       // w^0
          for (let l = 0; l < polyPart.length; l++) {
            const cR = polyPart[l].re, cI = polyPart[l].im;
            hr += cR * wPowRe - cI * wPowIm;
            hi += cR * wPowIm + cI * wPowRe;
            if (l + 1 < polyPart.length) {
              const nR = wPowRe * w.re - wPowIm * w.im;
              const nI = wPowRe * w.im + wPowIm * w.re;
              wPowRe = nR; wPowIm = nI;
            }
          }
        }
        // Finite-pole part.
        for (const pole of hData.poles) {
          const dx = w.re - pole.a.re, dy = w.im - pole.a.im;
          const d2 = dx*dx + dy*dy;
          if (d2 < 1e-14) { bad = true; break; }
          // dPowRe + dPowIm·i  =  (w - a)^{s+1}  for s = 0, 1, ...
          let dPowRe = dx, dPowIm = dy;
          for (let s = 0; s < pole.principal.length; s++) {
            const cR = pole.principal[s].re, cI = pole.principal[s].im;
            const den = dPowRe*dPowRe + dPowIm*dPowIm;
            if (den < 1e-30) { bad = true; break; }
            // C / (dPowRe + dPowIm·i)  =  C · (dPowRe − dPowIm·i) / den
            hr += (cR * dPowRe + cI * dPowIm) / den;
            hi += (cI * dPowRe - cR * dPowIm) / den;
            // Advance (w-a)^{s+1} -> (w-a)^{s+2}
            const nRe = dPowRe * dx - dPowIm * dy;
            const nIm = dPowRe * dy + dPowIm * dx;
            dPowRe = nRe; dPowIm = nIm;
          }
          if (bad) break;
        }
        if (bad) continue;

        // Per-family extra contributions to h (e.g. q/w pole at origin for
        // singular LQDs). Default is null = no extra terms.
        const desc = modeDescriptor();
        if (desc.extraHContrib) {
          const extra = desc.extraHContrib(w, hData, this.data && this.data.phi, state);
          hr += extra.re;
          hi += extra.im;
        }

        // Pólya field is V = conj(h). External-potential field depends on
        // the QD family:
        //   classical/PQD:   V = w − conj(h)        (∇ of |w|² − 2 Re H)
        //   LQD:             V = ln|w|²/conj(w) − conj(h)
        //                                            (∇ of (1/2)ln²|w|² − 2 Re H)
        // Math-plane vectors (before screen y-flip):
        //   conj(h)      = (Re h, −Im h)
        //   w − conj(h)  = (Re w − Re h, Im w + Im h)
        //   ln|w|²/conj(w) = ln|w|² · (Re w, Im w) / |w|²
        //                  = (Re w · ln|w|²/|w|², Im w · ln|w|²/|w|²)
        // Screen y is flipped vs the imaginary axis, so the screen-direction
        // negates the math y-component.
        let fieldX, fieldY;
        if (state.vectorFieldMode === 'external') {
          if (desc.externalFieldKind === 'lqd') {
            // V = ln|w|² / conj(w) − conj(h). Clip near origin per descriptor
            // (singular LQDs need a larger floor since 0 ∈ Ω).
            const absW2 = w.re * w.re + w.im * w.im;
            if (absW2 < desc.vectorFieldOriginAbs2Floor) continue;
            const logScale = Math.log(absW2) / absW2;
            fieldX = w.re * logScale - hr;
            fieldY = -(w.im * logScale) - hi;
            // ln|w|²/conj(w) = ln|w|² · w/|w|², so its math-Im = (Im w) · ln|w|²/|w|²;
            // screen-y = −(math-Im); and conj(h)'s math-Im = −Im h, screen-y = +Im h.
            // V = (math)  (Re w · L − Re h, Im w · L + Im h)  with L = ln|w|²/|w|²
            // screen      (Re w · L − Re h, −(Im w · L + Im h))
            //           = (Re w · L − Re h, −Im w · L − Im h)
            // (Re-deriving here matches the assignment above.)
          } else {
            fieldX =  w.re - hr;
            fieldY = -w.im - hi;
          }
        } else {
          // Pólya field V = conj(h): screen vector is (Re h, Im h).
          fieldX = hr;
          fieldY = hi;
        }
        const mag = Math.hypot(fieldX, fieldY);
        if (!isFinite(mag) || mag === 0) continue;

        // Length + opacity from tanh(log10(1 + |V|)) — short/dim for small
        // values, saturating for very large ones.
        const sat = Math.tanh(Math.log10(1 + mag));
        const len = arrowMaxLen * (0.25 + 0.75 * sat);
        const alpha = 0.18 + 0.55 * sat;

        const dirX = fieldX / mag;
        const dirY = fieldY / mag;

        // Center the arrow on the grid point.
        const baseX = sx - 0.35 * len * dirX;
        const baseY = sy - 0.35 * len * dirY;
        const tipX  = sx + 0.65 * len * dirX;
        const tipY  = sy + 0.65 * len * dirY;

        c.strokeStyle = `rgba(58, 84, 124, ${alpha.toFixed(3)})`;
        c.beginPath();
        c.moveTo(baseX, baseY);
        c.lineTo(tipX, tipY);
        c.stroke();

        // Arrowhead
        const ahLen = Math.max(3, len * 0.32);
        const ahAng = 0.45;
        const ang = Math.atan2(dirY, dirX);
        c.beginPath();
        c.moveTo(tipX, tipY);
        c.lineTo(tipX - ahLen * Math.cos(ang - ahAng), tipY - ahLen * Math.sin(ang - ahAng));
        c.moveTo(tipX, tipY);
        c.lineTo(tipX - ahLen * Math.cos(ang + ahAng), tipY - ahLen * Math.sin(ang + ahAng));
        c.stroke();
      }
    }
  }

  // Pick a "nice" tick spacing for the current scale
  niceStep() {
    const target = 80 / this.view.scale;            // ~80 px per major grid line
    const exp = Math.floor(Math.log10(target));
    const frac = target / Math.pow(10, exp);
    let step;
    if      (frac < 1.5) step = 1;
    else if (frac < 3)   step = 2;
    else if (frac < 7)   step = 5;
    else                 step = 10;
    return step * Math.pow(10, exp);
  }

  drawGrid() {
    const c = this.ctx;
    const step = this.niceStep();
    const tl = this.toWorld(0, 0);
    const br = this.toWorld(this.cssW, this.cssH);
    const minRe = Math.floor(tl.re / step) * step;
    const maxRe = Math.ceil(br.re / step) * step;
    const minIm = Math.floor(br.im / step) * step;
    const maxIm = Math.ceil(tl.im / step) * step;

    c.strokeStyle = '#e8eaef';
    c.lineWidth = 1;
    c.beginPath();
    for (let r = minRe; r <= maxRe + 1e-9; r += step) {
      const x = this.toScreen(r, 0).x;
      c.moveTo(x, 0); c.lineTo(x, this.cssH);
    }
    for (let i = minIm; i <= maxIm + 1e-9; i += step) {
      const y = this.toScreen(0, i).y;
      c.moveTo(0, y); c.lineTo(this.cssW, y);
    }
    c.stroke();

    // tick labels
    c.fillStyle = '#777';
    c.font = '10px ui-monospace, "SF Mono", Consolas, monospace';
    c.textAlign = 'left';
    c.textBaseline = 'top';
    const y0 = Math.max(2, Math.min(this.cssH - 12, this.toScreen(0, 0).y + 2));
    for (let r = minRe; r <= maxRe + 1e-9; r += step) {
      if (Math.abs(r) < step * 1e-6) continue;
      c.fillText(formatTick(r, step), this.toScreen(r, 0).x + 2, y0);
    }
    c.textAlign = 'left';
    const x0 = Math.max(2, Math.min(this.cssW - 30, this.toScreen(0, 0).x + 2));
    for (let i = minIm; i <= maxIm + 1e-9; i += step) {
      if (Math.abs(i) < step * 1e-6) continue;
      c.fillText(formatTick(i, step) + 'i', x0, this.toScreen(0, i).y + 2);
    }
  }

  drawAxes() {
    const c = this.ctx;
    c.strokeStyle = '#bbb';
    c.lineWidth = 1;
    c.beginPath();
    const yAxisX = this.toScreen(0, 0).x;
    const xAxisY = this.toScreen(0, 0).y;
    c.moveTo(0, xAxisY); c.lineTo(this.cssW, xAxisY);
    c.moveTo(yAxisX, 0); c.lineTo(yAxisX, this.cssH);
    c.stroke();
  }

  drawBoundary() {
    const c = this.ctx;
    const pts = this.data.boundaryPts;
    c.beginPath();
    const p0 = this.toScreen(pts[0].re, pts[0].im);
    c.moveTo(p0.x, p0.y);
    for (let i = 1; i < pts.length; i++) {
      const p = this.toScreen(pts[i].re, pts[i].im);
      c.lineTo(p.x, p.y);
    }
    c.closePath();

    const ok = this.data.univalent;
    if (this.data.unbounded) {
      // Unbounded: shade the bounded complement K (= inside of the boundary
      // curve) in a contrasting muted color and outline ∂Ω.
      c.fillStyle = ok ? 'rgba(180, 195, 220, 0.45)' : 'rgba(220, 180, 180, 0.45)';
      c.fill('evenodd');
      c.strokeStyle = ok ? '#1a3e7a' : '#b53030';
      c.lineWidth = 1.8;
      c.stroke();
    } else {
      // Bounded: shade Ω (= inside of the curve) in the standard tint.
      c.fillStyle   = ok ? 'rgba(86, 119, 168, 0.16)' : 'rgba(181, 48, 48, 0.14)';
      c.fill('evenodd');
      c.strokeStyle = ok ? '#1a3e7a' : '#b53030';
      c.lineWidth = 1.6;
      c.stroke();
    }
  }

  drawPoles() {
    const c = this.ctx;
    c.font = '11px system-ui, sans-serif';
    c.textBaseline = 'middle';
    for (let i = 0; i < this.data.poles.length; i++) {
      const p = this.data.poles[i];
      const s = this.toScreen(p.re, p.im);
      c.beginPath();
      c.arc(s.x, s.y, 5.5, 0, 2*Math.PI);
      c.fillStyle = '#b53030';
      c.fill();
      c.strokeStyle = '#fff'; c.lineWidth = 1.5;
      c.stroke();
      c.fillStyle = '#b53030';
      c.textAlign = 'left';
      c.fillText('a' + sub(i+1), s.x + 7, s.y);
    }
  }

  drawW0() {
    const c = this.ctx;
    const s = this.toScreen(this.data.w0.re, this.data.w0.im);
    c.strokeStyle = '#1a3e7a';
    c.lineWidth = 1.6;
    c.beginPath();
    c.moveTo(s.x - 5, s.y); c.lineTo(s.x + 5, s.y);
    c.moveTo(s.x, s.y - 5); c.lineTo(s.x, s.y + 5);
    c.stroke();
    c.fillStyle = '#1a3e7a';
    c.font = '11px system-ui, sans-serif';
    c.textAlign = 'left';
    c.textBaseline = 'top';
    c.fillText('φ(0)', s.x + 6, s.y + 4);
  }

  // -------------------------------------------------------------------------
  // Critical-set overlay: w-plane images of {z : φ'(z) = 0}.
  //
  // Computed lazily on first request (and on toggle-on after a fresh solve),
  // cached on state.current.criticalSet so subsequent pans/zooms don't pay
  // the Newton cost.  The cache is keyed by reference identity of
  // this.data.phi — if a new solve produces a new phi object, the cache
  // is recomputed automatically.
  //
  // Visual encoding:
  //   severity 'critical' (zero of φ' strictly inside the relevant disk):
  //     red filled disk, 6 px radius.  This is the bad case — φ is
  //     non-univalent.
  //   severity 'near'     (|z| within 0.05 of the unit circle):
  //     orange filled disk, 5 px radius.  Imminent-degeneracy warning.
  //   severity 'safe'     (zero of φ' outside the relevant disk):
  //     small gray hollow circle, 3.5 px radius.  Background info.
  //
  // Each marker carries a 1-letter tag showing severity.
  // -------------------------------------------------------------------------
  drawCriticalSet() {
    const phi = this.data.phi;
    if (!phi || typeof QD === 'undefined' || !QD.findCriticalPoints) return;

    // Look up / refresh cache.
    if (!state.current) return;
    const cached = state.current.criticalSet;
    let cs;
    if (cached && cached._phiRef === phi) {
      cs = cached;
    } else {
      try {
        cs = QD.findCriticalPoints(phi);
      } catch (e) {
        return;       // silent on solver error — overlay is purely diagnostic
      }
      cs._phiRef = phi;
      state.current.criticalSet = cs;
    }
    if (!cs.points || cs.points.length === 0) return;

    const c = this.ctx;
    c.save();
    c.font = '10px ui-monospace, Consolas, monospace';
    c.textBaseline = 'middle';
    c.textAlign    = 'left';

    for (const p of cs.points) {
      const s = this.toScreen(p.w.re, p.w.im);
      // Skip if off-screen by a large margin (saves draw work; markers very
      // far outside the visible region clutter the corner-clip).
      if (s.x < -40 || s.x > this.cssW + 40) continue;
      if (s.y < -40 || s.y > this.cssH + 40) continue;

      let fill, stroke, r, tag;
      switch (p.severity) {
        case 'critical':
          fill   = '#d12d2d';
          stroke = '#ffffff';
          r      = 6.0;
          tag    = '!';
          break;
        case 'near':
          fill   = '#d97706';   // orange
          stroke = '#ffffff';
          r      = 5.0;
          tag    = '~';
          break;
        case 'safe':
        default:
          fill   = null;        // hollow
          stroke = '#888888';
          r      = 3.5;
          tag    = '';
          break;
      }
      c.beginPath();
      c.arc(s.x, s.y, r, 0, 2 * Math.PI);
      if (fill) {
        c.fillStyle = fill;
        c.fill();
      }
      c.strokeStyle = stroke;
      c.lineWidth   = (p.severity === 'safe') ? 1.2 : 1.6;
      c.stroke();

      if (tag) {
        c.fillStyle = (p.severity === 'critical') ? '#d12d2d' : '#a85706';
        c.fillText(tag + ' |z|=' + p.absZ.toFixed(3), s.x + r + 3, s.y);
      }
    }
    c.restore();
  }
}

function formatTick(v, step) {
  // Choose precision based on step size
  const digits = Math.max(0, -Math.floor(Math.log10(step)));
  return Number(v.toFixed(digits)).toString();
}

// ---------- Wire everything up ------------------------------------------
const plot = new DomainPlot($('#canvas'), $('#plot-readout'));
window.addEventListener('resize', () => plot.resize());

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

renderPolesList();
renderPolyCoefList();
$('#poly-part-section').classList.toggle('hidden', !modeAllowsPoly(state.mode));

// Polynomial part of h is meaningful exactly in the three unbounded family
// panels. Keep this predicate centralized so refreshHText / parseAndApplyHText
// agree with what the mode descriptors expose (cards.poly).
function modeAllowsPoly(mode) {
  return mode === 'unbounded' ||
         mode === 'lqd-unbounded' ||
         mode === 'lqd-unbounded-singular';
}

// ---------- Custom h(w) text input --------------------------------------
// The #h-text input is a two-way-coupled mirror of the structured pole grid
// and polynomial-part coefficient list. refreshHText() rebuilds the text
// from current state; parseAndApplyHText() goes the other direction via
// QD.parseH (Phase 1 strict PFD walker → Phase 2 general-rational fallback).
//
// Refresh is called from renderPolesList / renderPolyCoefList / setMode /
// applyPreset so the text mirrors structural state. The per-keystroke
// pole-residue text-field edits don't trigger a refresh (they'd cause
// double-translation churn while the user types); the next solve / preset
// / mode switch syncs the text box.
function refreshHText() {
  const inp = document.getElementById('h-text');
  if (!inp) return;
  try {
    const poles = state.poles.map(po => {
      const a = QD.Complex.parse(po.a) || { re: 0, im: 0 };
      const residues = po.residues.slice(0, po.order).map(r =>
        QD.Complex.parse(r) || { re: 0, im: 0 });
      return { a, order: po.order, residues };
    });
    let polyCoeffs = [];
    if (modeAllowsPoly(state.mode) && state.polyDegree >= 0) {
      polyCoeffs = state.polyCoeffs.slice(0, state.polyDegree + 1).map(s =>
        QD.Complex.parse(s) || { re: 0, im: 0 });
    }
    inp.value = QD.formatH({ poles, polyCoeffs });
    setHTextMsg('');
  } catch (e) {
    // Defensive: never let formatter errors break the panel.
  }
}

function setHTextMsg(msg, kind) {
  const el = document.getElementById('h-text-msg');
  if (!el) return;
  el.textContent = msg || '';
  el.style.color = (kind === 'warn') ? '#9a6a00' : '#b53030';
}

function parseAndApplyHText() {
  const inp = document.getElementById('h-text');
  if (!inp) return;
  const expr = inp.value.trim();
  if (!expr) { setHTextMsg('Enter an expression in w.'); return; }
  let parsed;
  try {
    parsed = QD.parseH(expr, math, { mode: state.mode });
  } catch (e) {
    setHTextMsg(e.message || String(e));
    return;
  }

  // Convert parsed.poles (Complex-typed) back to the state's string form.
  if (parsed.poles.length === 0) {
    // Need at least one row in the grid so the user can extend it.
    state.poles = [{ a: '0', order: 1, residues: ['0'] }];
  } else {
    state.poles = parsed.poles.map(p => ({
      a: QD.Complex.format(p.a),
      order: p.order,
      residues: p.residues.map(c => QD.Complex.format(c)),
    }));
  }

  if (modeAllowsPoly(state.mode)) {
    if (parsed.polyCoeffs.length > 0) {
      state.polyCoeffs = parsed.polyCoeffs.map(c => QD.Complex.format(c));
      state.polyDegree = parsed.polyCoeffs.length - 1;
    } else {
      state.polyDegree = -1;
      state.polyCoeffs = [];
    }
    syncPolyDegreeInput();
  }

  for (const k of Object.keys(magSliderMax)) delete magSliderMax[k];
  renderPolesList();
  renderPolyCoefList();
  markAsCustom();
  if (parsed.warnings && parsed.warnings.length) {
    setHTextMsg('Parsed with warning: ' + parsed.warnings[0], 'warn');
  } else {
    setHTextMsg('');
  }
  scheduleSolve();
}

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
  state.poles = p.poles.map(po => ({
    a: po.a,
    order: po.order,
    residues: po.residues.slice(),
  }));
  if (state.mode === 'bounded') {
    state.w0Mode = 'auto';
    $('input[name="w0mode"][value="auto"]').checked = true;
    $('#w0-manual').disabled = true;
  } else if (state.mode === 'lqd-bounded' || state.mode === 'lqd-bounded-singular') {
    // LQDs are always manual: w₀ must be explicit and ≠ 0.
    state.w0Mode = 'manual';
    state.w0Manual = p.w0 || '1';
    $('input[name="w0mode"][value="manual"]').checked = true;
    $('#w0-manual').disabled = false;
    $('#w0-manual').value = state.w0Manual;
    // Singular: reset q to '0' by default (the user can dial it via slider).
    if (state.mode === 'lqd-bounded-singular') {
      const presetQ = p.q || '0';
      setQ(presetQ);
    }
  } else {
    if (typeof p.c === 'number') setC(p.c);
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
  if (text) text.value = c;
  const lbl = $('#c-val');
  if (lbl) lbl.textContent = c.toFixed(3);
}

// Switch UI between bounded QD, unbounded QD, and bounded LQD modes.
// Preserves pole data (the user's a_j and C_{j,s} stay), just swaps the
// normalization card visible and refreshes the preset list. The three
// modes share most of the UI; the differences are:
//   bounded     → w₀ card visible, polynomial-part hidden, h ∈ Rat₀(Ω)
//   unbounded   → c card visible, polynomial-part available, h ∈ Rat(Ω)
//   lqd-bounded → w₀ card visible (manual only, ≠ 0), polynomial-part hidden,
//                 h ∈ Rat₀(Ω), weight is ρ₀(w) = |w|⁻² instead of 1
function setMode(newMode) {
  if (!MODES[newMode]) return;
  if (state.mode === newMode) return;
  state.mode = newMode;
  const desc = modeDescriptor();
  // Card visibility from descriptor.
  $('#w0-card').classList.toggle('hidden',           !desc.cards.w0);
  $('#c-card').classList.toggle('hidden',            !desc.cards.c);
  $('#poly-part-section').classList.toggle('hidden', !desc.cards.poly);
  $('#q-card').classList.toggle('hidden',            !desc.cards.q);
  // Hint elements: show only the one this mode names (if any).
  for (const hintId of ['lqd-hint', 'lqd-singular-hint']) {
    const el = $('#' + hintId);
    if (el) el.style.display = (desc.hint === hintId) ? '' : 'none';
  }
  // LQD modes require an explicit nonzero w₀ → default w₀ to manual.
  if (desc.requireManualW0 && state.w0Mode === 'auto') {
    state.w0Mode = 'manual';
    const am = $('input[name="w0mode"][value="manual"]');
    if (am) am.checked = true;
    $('#w0-manual').disabled = false;
  }
  const vfExt = $('#vf-external-opt');
  if (vfExt) vfExt.textContent = desc.externalFieldLabel;
  populatePresetDropdown();
  markAsCustom();
  setC(state.c);
  // Sync the polynomial UI with state.polyDegree.
  syncPolyDegreeInput();
  renderPolyCoefList();
  scheduleSolve();
}

function syncPolyDegreeInput() {
  const inp = $('#poly-degree');
  if (inp) inp.value = state.polyDegree;
}

// Selecting a preset loads it; the user editing anything afterward reverts
// the dropdown to "— custom —".
function markAsCustom() { $('#preset-select').value = ''; }

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
  const headerOf = (cardSelector) => {
    const card = document.querySelector(cardSelector);
    return card ? card.querySelector('h2') : null;
  };
  // Domain type
  H(headerOf('#domain-mode-card'),
    `<b>Domain type.</b> Pick which family the inverse solver should target.
     <i>Bounded</i> Ω is the classical case (Riemann map of 𝔻); <i>Unbounded</i>
     Ω uses 𝔻* with a conformal radius c at ∞. <i>LQD</i> variants add a
     log-weighted exponent; <i>singular</i> variants allow a higher-order pole
     of h at the origin (with a free residue q).`);
  // Quadrature function h(w)
  H(headerOf('#h-card'),
    `<b>Quadrature data h(w).</b> Sum of rational and polynomial terms.
     Edit poles + residues structurally below, or paste a math.js expression
     in the textbox at the top. The inverse solver finds Ω whose
     quadrature data matches this h.`);
  H(headerOf('#q-card'),
    `<b>Residue at the origin (q).</b> For <i>singular</i> LQDs, q is a free
     parameter representing the residue of the log-weighted Schwarz function
     at w=0. The solver enforces a closed-form constraint linking q to the
     finite poles and any polynomial part of h.`);
  H(headerOf('#w0-card'),
    `<b>Riemann map center φ(0).</b> The image of the disk center 0 ∈ 𝔻 under
     the Riemann map. Together with c (the conformal radius) this fixes the
     gauge of φ. For bounded families w₀ is a free parameter; for unbounded
     families it is implicit and not editable.`);
  H(headerOf('#c-card'),
    `<b>Conformal radius c = φ'(∞).</b> For unbounded families, scales the
     Riemann map's behaviour at infinity. Together with w₀ this fixes the
     gauge of φ.`);
  H(headerOf('#solver-settings-card'),
    `<b>Solver settings.</b> Choose a preset (Fast / Default / Thorough) or
     fine-tune via the <i>Search options</i> panel. The preset balances
     Newton iterations, identity-check sample count, and how many alternate
     solution branches the solver attempts to find.`);
  H(headerOf('#search-options-card'),
    `<b>Search options.</b> Each phase is a distinct strategy for finding a φ
     consistent with h(w). Direct = single Newton from the initial guess;
     continuation = parameter-homotopy from a related solved scenario;
     multistart = many random seeds; diverse + deflation = explicit
     branch-finding.`);
  H(headerOf('#status-card'),
    `<b>Status.</b> Live readout of the solver: convergence diagnostics,
     identity residual, univalence, and which branches succeeded.`);
  H(headerOf('#riemann-map-card'),
    `<b>Riemann map φ(z).</b> The solved φ as a closed-form rational
     expression (bounded) or rational+polynomial (unbounded). Each pole
     of h corresponds to a pole structure of φ inside 𝔻.`);
  H(headerOf('#alternates-card'),
    `<b>Alternate solutions.</b> When more than one φ satisfies the same h
     (multiple branches), the solver lists them here. Click an alternate to
     promote it to the primary.`);
}
mountQolHelp();

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
$$('input[name="domain-mode"]').forEach(r => r.addEventListener('change', e => {
  setMode(e.target.value);
}));

// Preset dropdown
populatePresetDropdown();
$('#preset-select').addEventListener('change', e => {
  if (e.target.value) applyPreset(e.target.value);
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

// Conformal radius c (unbounded only): slider drives quick solve, manual
// text drives full solve.
$('#c-slider').addEventListener('input', e => {
  const c = +e.target.value;
  state.c = c;
  $('#c-val').textContent = c.toFixed(3);
  $('#c-manual').value = c;
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
$('#vector-field-mode').addEventListener('change', e => {
  state.vectorFieldMode = e.target.value;
  plot.render();
});
$('#critical-set-toggle').addEventListener('change', e => {
  state.showCriticalSet = e.target.checked;
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
    const norm = state.current.unbounded
      ? { unbounded: true, c: state.current.cUsed }
      : { w0: state.current.w0Used };
    startBackgroundAltSearch(state.current.hData, norm);
  });

  // Initial population.
  readSearchOptions();
})();

// Try-harder button: runs the "exhaustive" aggressiveness preset (which has
// a much larger multistart budget, tighter Newton, longer continuation, and
// implicitly more deflation rounds). Useful when the current default-level
// solve has flagged the primary as non-valid (spurious / non-univalent).
$('#try-harder-btn').addEventListener('click', () => {
  const btn = $('#try-harder-btn');
  const busy = $('#try-harder-busy');
  btn.disabled = true;
  busy.classList.remove('hidden');
  // Yield so the spinner can paint, then solve.
  setTimeout(() => {
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
      const result = QD.solveInverseQD(built, opts);
      state.current = result;
      state.current.hData = built;
      state.current.w0Used = norm.w0;
      state.current.cUsed  = norm.c;
      state.current.unbounded = !!norm.unbounded;
      state.selectedSolutionIdx = 0;
      if (!result.success) {
        setStatus({ kind: 'err',
          text: 'Exhaustive search still found no algebraic root.\n  reason: ' + result.error });
        return;
      }
      showSolution(result.primary, built, /*autoFit=*/false);
      refreshAlternatesPanel();
    } finally {
      btn.disabled = false;
      busy.classList.add('hidden');
    }
  }, 30);
});

// plot controls
$('#btn-fit').addEventListener('click', () => plot.fit());
$('#btn-reset').addEventListener('click', () => plot.reset());

// alternates panel: click to view
$('#alternates-list').addEventListener('click', e => {
  const idx = e.target.dataset.altIdx;
  if (idx !== undefined) viewSolutionByIndex(+idx);
});

// Initial solve
solveAndRender();

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
  // Populate the QD/LQD state from hData (and, for unbounded, the conformal
  // radius c) and re-render the inverse-tab sidebar. Then switch tabs.
  opts = opts || {};
  const unbounded = !!opts.unbounded;

  state.poles = hData.poles.map(p => ({
    a: QD.Complex.toString(p.a, 6),
    order: p.principal.length,
    residues: p.principal.map(c => QD.Complex.toString(c, 6)),
  }));

  if (unbounded) {
    // Unbounded mode: set c-slider + polyPart, no manual w0.
    const ubRadio = document.querySelector('input[name="domain-mode"][value="unbounded"]');
    if (ubRadio) ubRadio.checked = true;
    state.mode = 'unbounded';
    if (typeof opts.c === 'number' && opts.c > 0) {
      const cInput = document.getElementById('c-manual');
      const cSlider = document.getElementById('c-slider');
      if (cInput)  cInput.value  = opts.c.toString();
      if (cSlider) cSlider.value = opts.c.toString();
      state.c = opts.c;
    }
    // Populate polyPart fields.
    const polyPart = hData.polyPart || [];
    const polyDeg = polyPart.length - 1;
    const polyDegInput = document.getElementById('poly-degree');
    if (polyDegInput) polyDegInput.value = polyDeg.toString();
    state.polyDegree = polyDeg;
    state.polyCoefs = polyPart.map(c => QD.Complex.toString(c, 6));
    renderPolyCoefList();
    $('#poly-part-section').classList.toggle('hidden', false);
  } else {
    // Bounded mode: set manual w0 to the (single) pole location.
    const boundedRadio = document.querySelector('input[name="domain-mode"][value="bounded"]');
    if (boundedRadio) boundedRadio.checked = true;
    state.mode = 'bounded';
    if (hData.poles.length === 1) {
      document.querySelector('input[name="w0mode"][value="manual"]').checked = true;
      document.getElementById('w0-manual').disabled = false;
      document.getElementById('w0-manual').value = QD.Complex.toString(hData.poles[0].a, 6);
      state.w0Mode = 'manual';
      state.w0Manual = QD.Complex.toString(hData.poles[0].a, 6);
    }
    $('#poly-part-section').classList.toggle('hidden', true);
  }

  renderPolesList();
  document.getElementById('preset-select').value = '';

  // Switch view-mode to inverse and re-solve (HANDOFF #30: the Direct UI is
  // now a view within the QD tab, so this is a view-mode switch rather than
  // a tab switch).
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
    const radio = document.querySelector(`input[name="domain-mode"][value="${mode}"]`);
    if (radio) radio.checked = true;
    setMode(mode);
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
    const cVal    = $('#c-val');
    if (cInput)  cInput.value  = norm.c.toString();
    if (cSlider) cSlider.value = norm.c.toString();
    if (cVal)    cVal.textContent = norm.c.toFixed(3);
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
