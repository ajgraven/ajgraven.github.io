// @ts-check
// =============================================================================
// qd.mjs -- ES-module façade over the classic-script QD namespace.
//
// The app's solver, math, and UI files are loaded as classic <script> tags
// because (a) the param-slice and primary-solver Workers fetch raw .js
// sources and concatenate them into a Blob at runtime — a pattern that
// pre-dates native module workers, and (b) it keeps the no-build-step
// guarantee. Native ESM is, however, far more pleasant for go-to-definition,
// for downstream tooling (bundlers, type-checkers), and for any external
// consumer that wants `import { Complex, solveInverseQD } from './qd.mjs'`.
//
// This module re-exports the public QD surface as named ESM exports. It is
// the migration bridgehead, not a port: consumers can adopt ESM at their
// own pace by switching `window.QD.solveInverseQD(...)` to
// `import { solveInverseQD } from './qd.mjs'`.
//
// Load order: this module assumes the classic <script> tags have run and
// populated `window.QD`. In `index.html` add `<script type="module"
// src="qd.mjs"></script>` AFTER the existing script block if you want any
// inline module-script consumer to find it.
//
// Migration plan (suggested staging):
//   1. (NOW) This façade — no risk to existing flow.
//   2. (LATER) Convert leaf modules (complex.js, taylor.js, primary-solution.js)
//      to dual exports: keep the existing IIFE branch AND add an `export`
//      block for ESM imports. Both populate the same shapes.
//   3. (LATER) Convert solver-*.js families, paying attention to family-
//      registration order; the Family registry needs an explicit `init()`
//      call once all family files have imported.
//   4. (LATER) Rework the Worker pool to use `new Worker(new URL('./worker.mjs',
//      import.meta.url), { type: 'module' })` and let the worker import the
//      same modules natively — drops the runtime Blob bundling.
//   5. (LATER) UI files (ui.js, schwarz-ui.js, …) — these can stay last since
//      they have the most DOM coupling and the fewest cross-imports.
// =============================================================================

const _global = (typeof window !== 'undefined') ? window
              : (typeof globalThis !== 'undefined') ? globalThis : {};
const _QD = _global.QD;

if (!_QD) {
  throw new Error(
    'qd.mjs: window.QD is undefined. Load the classic-script QD bundle ' +
    '(complex.js, solver.js, family files, …) BEFORE importing this module.'
  );
}

// Math primitives.
export const Complex = _QD.Complex;
export const Taylor  = _QD.Taylor;

// Solver core.
export const evalPhi              = _QD.evalPhi;
export const phiTaylorAt          = _QD.phiTaylorAt;
export const residual             = _QD.residual;
export const residualNorm         = _QD.residualNorm;
export const packPhi              = _QD.packPhi;
export const unpackPhi            = _QD.unpackPhi;
export const newtonSolve          = _QD.newtonSolve;
export const solveInverseQD       = _QD.solveInverseQD;
export const searchAlternates     = _QD.searchAlternates;
export const isBoundaryUnivalent  = _QD.isBoundaryUnivalent;
export const sampleBoundary       = _QD.sampleBoundary;
export const sampleBoundaryAdaptive = _QD.sampleBoundaryAdaptive;
export const binomialCoeff        = _QD.binomialCoeff;
export const selectFamily         = _QD.selectFamily;
export const registerFamily       = _QD.registerFamily;
export const packPhiBySchema      = _QD.packPhiBySchema;
export const unpackPhiBySchema    = _QD.unpackPhiBySchema;
export const applySchemaClamps    = _QD.applySchemaClamps;

// Linear algebra (P1.2 — Householder QR).
export const solveLinearSystem    = _QD.solveLinearSystem;
export const solveLeastSquares    = _QD.solveLeastSquares;
export const houseQR              = _QD.houseQR;
export const numericalJacobian    = _QD.numericalJacobian;

// Inverse-Faber primitives (QD.Faber.*).
export const Faber                = _QD.Faber;

// Direct problem (QD.Direct.*).
export const Direct               = _QD.Direct;

// Schwarz / Sphere / param-slice subsystems (object namespaces).
export const Schwarz              = _QD.Schwarz;
export const Sphere               = _QD.Sphere;

// Custom h(w) text input.
export const parseH               = _QD.parseH;
export const formatH              = _QD.formatH;

// Critical-set kernel.
export const CriticalSet          = _QD.CriticalSet;
export const findCriticalPoints   = _QD.findCriticalPoints;

// Cross-tab envelope (P0.1a) + warm worker (P0.2).
export const PrimarySolution      = _QD.PrimarySolution;
export const PrimarySolverWorker  = _QD.PrimarySolverWorker;

// Family registry — convenience.
export const Family               = _QD.Family;

// The whole QD object — escape hatch for anything not yet enumerated above.
export default _QD;
