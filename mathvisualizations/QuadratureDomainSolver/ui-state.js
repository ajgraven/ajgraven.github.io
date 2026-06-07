// =============================================================================
// ui-state.js -- The single source-of-truth state object for the QD/LQD tab.
//
// Extracted from ui.js by A1 (post-P3 audit). Continues to live at script
// scope (`const state = {...}` in a classic <script>) so ui.js / schwarz-ui.js
// / direct-ui.js etc. can read it via the cross-script lexical environment,
// same as before. Additionally exposed on `window.QD_UI.state` for callers
// that prefer explicit namespacing.
//
// Mutators: ui.js owns writes (see solveAndRender, quickSolveAndRender,
// the various event listeners). Readers across modules subscribe to
// QD.PrimarySolution for the SOLVE result; for non-solve state (mode,
// preset selection, etc.) reads remain via the script-realm `state`
// binding.
//
// Shape:
//   poles[]       — h(w) principal parts (parsed from the structured grid)
//   mode          — one of the MODES keys in ui.js (NB these are MODE keys,
//                   not family tags): 'bounded' | 'unbounded' |
//                   'lqd-bounded' | 'lqd-bounded-singular' |
//                   'lqd-unbounded' | 'lqd-unbounded-singular' |
//                   'pqd-bounded' | 'pqd-bounded-singular' |
//                   'pqd-unbounded' | 'pqd-unbounded-singular'
//   c, polyDegree, polyCoeffs, w0Mode, w0Manual, q, alpha — per-family knobs
//   aggressiveness, samples, autoFit, vectorFieldMode, showCriticalSet
//   current       — last successful solve envelope (mirrors QD.PrimarySolution.get());
//                   see the PrimaryEnvelope typedef in primary-solution.js for the
//                   full field set (primary, alternates, hData, w0Used, cUsed,
//                   geomProps, regimeSwitched, …)
//   selectedSolutionIdx — primary (0) vs alternate (1+) shown on canvas
//   altSearchActive, altSearchToken — bg-search bookkeeping
//   searchOptions — advanced-panel overrides (all numeric fields null → use preset)
//   viewMode, directMounted — HANDOFF #30 inverse|direct toggle bookkeeping
// =============================================================================

'use strict';

// eslint-disable-next-line no-redeclare
const state = {
  // h(w) data: array of { a: string, order: int, residues: [string,...] }
  poles: [
    { a: '0', order: 1, residues: ['1'] },
  ],
  mode: 'bounded',           // a MODES key (see ui.js): 'bounded' | 'unbounded'
                             // | 'lqd-*' | 'pqd-*' (10 in total). NOT a family tag.
  c: 0.5,                    // conformal radius (unbounded mode only)
  polyDegree: -1,            // polynomial-part degree m_∞ (-1 = none)
  polyCoeffs: [],            // strings, length = max(polyDegree+1, 0)
  w0Mode: 'auto',            // 'auto' (centroid of poles) | 'manual'; default for all families
  w0Manual: '0',
  q: '0',                    // residue of h at origin (singular LQD mode only)
  alpha: 1,                  // power for power-weighted QDs (PQDs); α=1 keeps
                             // classical QD behavior, any real α>0, α≠1 routes
                             // to a PQD family (Family.powerQD etc.). Shown in
                             // all four pqd-* modes (the α card), not in others.
  aggressiveness: 'standard',
  samples: 500,
  autoFit: false,            // off by default — slider drags shouldn't reframe the view
  autoSwitchSingular: true,  // §23: auto-switch singular ⇄ non-singular PQD when ∂Ω crosses 0
  vectorFieldMode: 'off',    // 'off' | 'polya' | 'external'
  showCriticalSet: false,    // overlay w-images of {z : φ'(z) = 0}
  showCurvature: false,      // overlay ∂Ω colored by curvature |κ| (cusp heat-strip)
  showPhenomena: false,      // #9 annotated-phenomena overlay (harmonic/curvature peaks + symmetry axes)
  showFaberRoots: false,     // overlay Faber-polynomial roots in the ζ-plane (UQD only)
  faberRoots: null,          // payload from ui-faber.js: { mode:'all'|'single', N, n, sets:[{n,roots,converged}] }

  // Solver result
  current: null,             // last solve envelope (see PrimaryEnvelope typedef
                             // in primary-solution.js for the full field set)
  selectedSolutionIdx: 0,    // 0 = primary, 1+ = alternate index

  // Background search bookkeeping (the solve debounce uses a closure-local
  // timer in ui.js scheduleSolve, not a state field).
  altSearchActive: false,
  altSearchToken: 0,

  // Advanced search-options panel. All numeric fields are blank-meaning-
  // "use preset"; toggles are explicit booleans. Mutated by readSearchOptions.
  searchOptions: {
    phases: {
      direct: true, continuation: true, multistart: true,
      diverse: true, deflation: true,
    },
    numRestarts:       null,
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

  // Transparent solution-summary panel overlaid bottom-right of the plot.
  statusPanelCollapsed: false,   // collapsed → only the verdict badge bar shows
};

// Explicit namespace export for callers that prefer not to rely on the
// cross-script lexical environment.
if (typeof window !== 'undefined') {
  window.QD_UI = window.QD_UI || {};
  window.QD_UI.state = state;
}
