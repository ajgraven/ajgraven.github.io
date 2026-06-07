// @ts-nocheck — public-API JSDoc anchors only; full ts-check is a later pass
// =============================================================================
// solver.js -- Inverse-problem solver: shared infrastructure
//
// This file contains everything that is FAMILY-AGNOSTIC: the family registry,
// dispatcher functions for `evalPhi` / `phiTaylorAt` / `residual` / `packPhi`
// / `unpackPhi` / `canonicalizePhi`, the Newton-Raphson driver, deflation, the
// boundary sampler & univalence checker, the top-level solveInverseQD
// orchestrator, the alternates searcher, and the QD namespace export.
//
// Per-family methods (Phi evaluation, residual blocks, identity verifier,
// initial-guess strategies, continuation strategies) live in dedicated files:
//   solver-qd.js                  — Family.boundedQD
//   solver-uqd.js                 — Family.unboundedQD
//   solver-lqd.js                 — Family.boundedLQD            (non-singular)
//   solver-lqd-singular.js        — Family.boundedLQD_singular
//   solver-uqd-lqd.js             — Family.unboundedLQD          (non-singular)
//   solver-uqd-lqd-singular.js    — Family.unboundedLQD_singular
//
// Each family file calls `QD.registerFamily('X')` after populating
// `QD.Family.X`. The dispatcher walks `familyDispatchOrder` (most-specific
// first; classical boundedQD is the catch-all default).
//
// Math background (Graven, Chapters II–V): the inverse problem reduces to a
// real algebraic system in (z_j, A_{j,k}, possibly F_l, z_0, γ, ...). The
// shape varies per family; this file just owns the solver mechanics.
// =============================================================================

// --------- Phi data structure ---------------------------------------------
//
// phi = {
//   family:    string,                            // 'boundedQD' | 'unboundedQD' | …
//   w0:        Complex,                            // φ(0) (bounded)
//   c:         number,                             // φ'(∞) (unbounded)
//   z0, gamma, q: Complex,                         // singular LQDs
//   polyA:     [Complex, ...],                     // F_l (unbounded poly part)
//   branches:  [ { z: Complex, A: [Complex, ...] }, ... ],
//   unbounded: bool,                               // legacy fallback tag
// }
//
// h_data = { poles: [{a, principal:[…]}, …], polyPart: [Complex, …] }
//
// JSDoc anchors for the most-touched public types — IDEs and `// @ts-check`
// in callers will use these to flag silent field-drop bugs (the HANDOFF #26
// /#28 class of failure).
/**
 * @typedef {{ re: number, im: number }} ComplexC
 *
 * @typedef {Object} Pole
 * @property {ComplexC} a
 * @property {ComplexC[]} principal
 *
 * @typedef {Object} HData
 * @property {Pole[]} poles
 * @property {ComplexC[]} [polyPart]
 *
 * @typedef {Object} Branch
 * @property {ComplexC} z
 * @property {ComplexC[]} A
 *
 * @typedef {Object} Phi
 * @property {string} [family]
 * @property {ComplexC} [w0]
 * @property {number} [c]
 * @property {boolean} [unbounded]
 * @property {ComplexC} [z0]
 * @property {ComplexC} [gamma]
 * @property {ComplexC} [q]
 * @property {Branch[]} branches
 * @property {ComplexC[]} [polyA]
 * @property {ComplexC[]} [lqdBeta]   — unbounded-LQD polynomial-h β-correction (HANDOFF #21)
 * @property {ComplexC[]} [lqdGamma]  — higher-order pole-at-origin synth branch (HANDOFF #24)
 *
 * @typedef {Object} SolveResult
 * @property {boolean} success
 * @property {Object} [primary]      — { phi: Phi, univalent, identity, identityOK, … }
 * @property {Object[]} [alternates]
 * @property {string} [error]
 * @property {Object[]} [attempts]
 */
// =============================================================================

// --------- Named numeric constants ----------------------------------------
// Centralized so changes propagate everywhere; literal sites elsewhere are
// being migrated to these names (HANDOFF #27 cleanup).
//
//   ZERO_THRESHOLD   — `Complex.abs2(c) < ZERO_THRESHOLD` ⇒ c is "zero". Used
//                      to detect a = 0 entries (higher-order pole at origin)
//                      and to test residue/pole-location nullity in solvers.
//   DISK_CLAMP_OUT   — schemas with side='out' default cap. Keeps |z| just
//                      OUTSIDE the unit disk for unbounded families.
//   DISK_CLAMP_IN    — schemas with side='in' default cap. Keeps |z| just
//                      INSIDE the unit disk for bounded families.
//   Z0_MAX_RADIUS    — upper bound on |z₀| for unbounded singular LQDs to
//                      keep Newton out of the z₀ → ∞ degenerate limit.
//   DEFAULT_FD_EPS   — forward-difference step for the numerical Jacobian.
const ZERO_THRESHOLD = 1e-20;
const DISK_CLAMP_OUT = 1.0001;
const DISK_CLAMP_IN  = 0.9999;
const Z0_MAX_RADIUS  = 1000;
const DEFAULT_FD_EPS = 1e-7;

// --------- Generic clone: must propagate every family-specific field ------
function clonePhi(phi) {
  return {
    family: phi.family,
    unbounded: !!phi.unbounded,
    c: phi.c,
    w0: phi.w0 ? Complex.clone(phi.w0) : { re: 0, im: 0 },
    polyA: phi.polyA ? phi.polyA.map(Complex.clone) : [],
    // lqdBeta: polynomial-h extension for unbounded LQDs (Σ β_l / z^l inside
    // the exp). Cloned through so warm-start across c-continuation steps
    // preserves the β's.
    lqdBeta: phi.lqdBeta ? phi.lqdBeta.map(Complex.clone) : [],
    // lqdGamma: higher-order-pole-at-origin extension for unbounded singular
    // LQDs (Σ c_l / z^l inside the exp, distinct from lqdBeta).
    lqdGamma: phi.lqdGamma ? phi.lqdGamma.map(Complex.clone) : [],
    z0:    phi.z0    ? Complex.clone(phi.z0)    : undefined,
    gamma: phi.gamma ? Complex.clone(phi.gamma) : undefined,
    q:     phi.q     ? Complex.clone(phi.q)     : undefined,
    // alpha: PQD power (Family.powerQD / powerQD_singular). Plain number;
    // must survive cloning or warm-start/capture silently drops it → NaN powers.
    alpha: phi.alpha,
    branches: phi.branches.map(br => ({
      z: Complex.clone(br.z),
      A: br.A.map(Complex.clone)
    }))
  };
}

// --------- Family dispatch helpers ----------------------------------------
// Each generic primitive (evalPhi, phiTaylorAt, residual, packPhi, unpackPhi,
// canonicalizePhi) is a thin dispatcher to Family[phi.family].method.
// Legacy fallback: if phi has no .family tag, use phi.unbounded to pick
// 'unboundedQD' vs 'boundedQD'. This keeps callers that hand-construct phi
// objects (older tests, the unbounded-init path) working.
//
// B10: WeakMap-cache the resolved Family object per phi reference. In the
// param-slice hot path _resolveFamily fires thousands of times per pixel
// (once per evalPhi / phiTaylorAt call inside each Newton iteration).
// Caching avoids the string-OR-ternary + property lookup once per phi.
const _familyCache = new WeakMap();
function _resolveFamily(phi) {
  let fam = _familyCache.get(phi);
  if (fam !== undefined) return fam;
  const name = phi.family || (phi.unbounded ? 'unboundedQD' : 'boundedQD');
  fam = Family[name];
  if (!fam) throw new Error("solver.js: family not registered: " + name);
  _familyCache.set(phi, fam);
  return fam;
}

function evalPhi(z, phi) {
  return _resolveFamily(phi).evalPhi(z, phi);
}

function phiTaylorAt(z0, phi, L) {
  return _resolveFamily(phi).phiTaylorAt(z0, phi, L);
}

function residual(phi, hData, options) {
  return _resolveFamily(phi).residual(phi, hData, options);
}

function packPhi(phi)                  { return _resolveFamily(phi).packPhi(phi); }
function unpackPhi(v, template)        { return _resolveFamily(template).unpackPhi(v, template); }
function canonicalizePhi(phi)          { return _resolveFamily(phi).canonicalizePhi(phi); }

// --------- Residual / Jacobian / linear solve -----------------------------
function residualNorm(F) {
  let s = 0;
  for (let i = 0; i < F.length; i++) s += F[i] * F[i];
  return Math.sqrt(s);
}

// ---------------------------------------------------------------------------
// Real Householder QR.  Used as the primitive for both solveLinearSystem
// (square) and solveLeastSquares (overdetermined).  Householder QR is
// backward-stable; for tall systems it solves the least-squares problem
// directly without forming the (much more ill-conditioned) normal equations
// A^T A x = A^T b.
//
// Returns { R, betas, diag, rank, condEst }:
//   * R       — upper-triangular (n × n).
//   * betas   — Householder reflector scalars (length min(m,n)).
//   * diag    — pre-reflector diagonal of R (length min(m,n)); when any
//               entry is below `singularTol`, rank is reduced by one.
//   * rank    — numerical rank (count of |diag[k]| ≥ singularTol).
//   * condEst — a cheap estimate of the condition number:
//                 max|diag[k]| / min|diag[k]| over the retained pivots.
//
// A and the optional rhs vector b are NOT mutated; the routine works on a
// scratch copy.  Callers wanting the Q · b transformation only (without
// triangular-solve) get it via `applyQt(b)` on the returned handle.
// ---------------------------------------------------------------------------
const QR_SINGULAR_TOL = 1e-13;

/**
 * Householder QR of a real m × n matrix (m ≥ n).
 * @param {number[][]} A
 * @returns {{ R: number[][], betas: number[], diag: number[],
 *             rank: number, condEst: number,
 *             applyQt: (b: number[]) => number[],
 *             backSolve: (rhs: number[]) => number[],
 *             m: number, n: number }}
 */
function houseQR(A) {
  const m = A.length;
  if (m === 0) throw new Error("houseQR: empty matrix");
  const n = A[0].length;
  if (m < n) throw new Error("houseQR: under-determined (m < n) not supported");
  // Working copy. R is stored in-place in the upper triangle of `R`.
  const R = A.map(row => row.slice());
  const betas = new Array(n).fill(0);
  const diag = new Array(n).fill(0);
  // Workspace for column k.
  for (let k = 0; k < n; k++) {
    // Compute the 2-norm of column k from row k downward.
    let normSq = 0;
    for (let i = k; i < m; i++) normSq += R[i][k] * R[i][k];
    const norm = Math.sqrt(normSq);
    if (norm === 0) { diag[k] = 0; betas[k] = 0; continue; }
    // Householder vector v: v = x - sign(x_k) * ||x|| * e_k. Sign chosen so
    // that x_k + sign(x_k) * ||x|| has the same sign as x_k (avoids cancellation).
    const sign = R[k][k] >= 0 ? 1 : -1;
    const alpha = -sign * norm;
    diag[k] = alpha;
    // v_k = x_k - alpha (in-place), so v_k = x_k + sign*norm > 0 when sign chosen as above.
    R[k][k] = R[k][k] - alpha;
    // beta = 2 / (v^T v). Recompute v^T v from updated v.
    let vNormSq = R[k][k] * R[k][k];
    for (let i = k + 1; i < m; i++) vNormSq += R[i][k] * R[i][k];
    const beta = (vNormSq > 0) ? 2 / vNormSq : 0;
    betas[k] = beta;
    // Apply reflector to remaining columns: A[:,j] -= beta * (v^T A[:,j]) * v.
    for (let j = k + 1; j < n; j++) {
      let dot = 0;
      for (let i = k; i < m; i++) dot += R[i][k] * R[i][j];
      const f = beta * dot;
      for (let i = k; i < m; i++) R[i][j] -= f * R[i][k];
    }
    // R[k][k] keeps v_k (the full Householder vector entry, sign and all)
    // so applyQt() can read it directly. diag[k] = alpha is the post-reflector
    // diagonal entry used by backSolve. Strictly-below-diagonal slots keep
    // v's tail. Reconstructing v_k from beta + tail would lose the sign on
    // x_k < 0; storing v_k explicitly avoids that.
  }
  // Rank + condition-number estimate from the diagonal of R.
  let absMax = 0, absMin = Infinity, rank = 0;
  for (let k = 0; k < n; k++) {
    const a = Math.abs(diag[k]);
    if (a >= QR_SINGULAR_TOL) {
      rank++;
      if (a > absMax) absMax = a;
      if (a < absMin) absMin = a;
    }
  }
  const condEst = (rank === n && absMin > 0) ? absMax / absMin : Infinity;

  // Apply Q^T to a vector b. v_k is stored directly in R[k][k], so no
  // sign-reconstruction needed.
  function applyQt(bIn) {
    const b = bIn.slice();
    for (let k = 0; k < n; k++) {
      const beta = betas[k];
      if (beta === 0) continue;
      // dot = v^T b restricted to rows k..m-1 (v_k stored at R[k][k]).
      let dot = 0;
      for (let i = k; i < m; i++) dot += R[i][k] * b[i];
      const f = beta * dot;
      for (let i = k; i < m; i++) b[i] -= f * R[i][k];
    }
    return b;
  }

  // Back-substitute on the upper-triangular R for a length-n right-hand side
  // (the leading n entries of Q^T b for least-squares; or the supplied b for
  // square systems). diag[i] holds the post-reflector diagonal entry of R;
  // the strict-upper-triangle is read from R[i][j] for j > i. Throws if R
  // is singular within QR_SINGULAR_TOL.
  function backSolve(rhsLeadingN) {
    const x = new Array(n);
    for (let i = n - 1; i >= 0; i--) {
      let s = rhsLeadingN[i];
      for (let j = i + 1; j < n; j++) s -= R[i][j] * x[j];
      if (Math.abs(diag[i]) < QR_SINGULAR_TOL) throw new Error("singular");
      x[i] = s / diag[i];
    }
    return x;
  }

  return { R, betas, diag, rank, condEst, applyQt, backSolve, m, n };
}

/**
 * Square real system A x = b. Routed through Householder QR for numerical
 * stability and a meaningful "singular" gate.
 * @param {number[][]} A   n × n real matrix
 * @param {number[]}   b   length-n right-hand side
 * @returns {number[]}     length-n solution; throws "singular" if R is rank-deficient
 */
// Above this estimated condition number, run iterative refinement after the QR
// solve (C5). Chosen well below the ~1e8 scale where Float64 QR starts losing
// digits, so refinement kicks in before accuracy degrades but not on
// well-conditioned systems where it would just waste budget.
const ILL_COND_REFINE_THRESHOLD = 1e6;

// Near-cusp accuracy (#11): when the residual Jacobian is ill-conditioned (the
// signature of a forming cusp / fold) the single forward-difference Jacobian and
// single refinement step cap how far Newton can drive the residual. We raise both
// ceilings — but ONLY when conditioning warrants it, so well-conditioned solves
// keep their current cost and byte-for-byte results.
//   • MAX_REFINE_STEPS — bound on residual-correction steps (was effectively 1).
//   • CENTRAL_FD_EPS    — step for central differences (≈ macheps^{1/3}); their
//     O(eps²) truncation gives ~1e-10 Jacobian entries vs ~1e-7 for forward.
//   • CENTRAL_DIFF_COND_TRIGGER — once a step's QR condEst exceeds this, switch
//     the Jacobian to central differences for the remaining iterations.
const MAX_REFINE_STEPS          = 3;
const CENTRAL_FD_EPS            = 1e-5;
const CENTRAL_DIFF_COND_TRIGGER = 1e5;

function solveLinearSystem(A, b, maxRefine = MAX_REFINE_STEPS) {
  const n = A.length;
  if (n === 0 || A[0].length !== n) throw new Error("solveLinearSystem: not square");
  const qr = houseQR(A);
  const Qtb = qr.applyQt(b);
  let x = qr.backSolve(Qtb.slice(0, n));
  // C5: iterative refinement when the system is ill-conditioned. Cheap (one
  // matvec + one Qt-apply + one backSolve per step) — buys back accuracy on
  // cond(A) > 1e6 systems. Skipped entirely when condEst is well-behaved.
  if (qr.condEst > ILL_COND_REFINE_THRESHOLD && isFinite(qr.condEst)) {
    x = _qrIterativeRefine(A, b, x, qr, maxRefine);
  }
  return x;
}

// Bounded iterative refinement: repeat r = b - A·x; dx = QR-solve(r); x ← x + dx
// until the correction stops shrinking or `maxSteps` is reached. One step
// recovers ~½ digit on a cond≈1e6 system; the extra steps matter near a cusp
// where the residual Jacobian climbs toward 1e8+. Operates on a copy of x.
function _qrIterativeRefine(A, b, x, qr, maxSteps) {
  const m = A.length;
  const n = A[0].length;
  const steps = (maxSteps > 0) ? maxSteps : 1;
  let out = x.slice();
  let prevCorr = Infinity;
  for (let s = 0; s < steps; s++) {
    const r = new Array(m);
    for (let i = 0; i < m; i++) {
      let acc = b[i];
      for (let j = 0; j < n; j++) acc -= A[i][j] * out[j];
      r[i] = acc;
    }
    const Qtr = qr.applyQt(r);
    let dx;
    try { dx = qr.backSolve(Qtr.slice(0, n)); }
    catch (_) { break; }   // singular — refinement is a no-op
    let corr = 0;
    for (let k = 0; k < n; k++) { out[k] += dx[k]; corr += dx[k] * dx[k]; }
    corr = Math.sqrt(corr);
    if (!(corr < prevCorr)) break;   // corrections no longer shrinking → stop
    prevCorr = corr;
  }
  return out;
}

// Least-squares solve that also reports the QR condition estimate, so callers
// (newtonSolve) can adapt the Jacobian scheme to ill-conditioning without a
// second factorization. Mirrors solveLeastSquares' refinement gate.
// Returns { x, condEst }; throws "singular" exactly as backSolve does.
function leastSquaresWithCond(A, b, maxRefine = MAX_REFINE_STEPS) {
  const n = A[0].length;
  const qr = houseQR(A);
  const Qtb = qr.applyQt(b);
  let x = qr.backSolve(Qtb.slice(0, n));
  if (qr.condEst > ILL_COND_REFINE_THRESHOLD && isFinite(qr.condEst)) {
    x = _qrIterativeRefine(A, b, x, qr, maxRefine);
  }
  return { x, condEst: qr.condEst };
}

/**
 * Solve min ‖A x − b‖₂ via Householder QR — directly, without forming the
 * normal equations A^T A x = A^T b (which would square the condition number
 * and was the major weakness of the previous implementation).
 * @param {number[][]} A   m × n real matrix with m ≥ n
 * @param {number[]}   b   length-m
 * @returns {number[]}     least-squares solution (length n)
 */
function solveLeastSquares(A, b, maxRefine = MAX_REFINE_STEPS) {
  const m = A.length;
  if (m === 0) throw new Error("solveLeastSquares: empty matrix");
  const n = A[0].length;
  if (m < n) throw new Error("solveLeastSquares: under-determined (m < n) not supported");
  // Square case: same Q^T-then-backsolve path as solveLinearSystem.
  // Tall case: same routine, but the leading n entries of Q^T b are the
  // ones we back-substitute against; the trailing m − n entries are the
  // (orthogonal) residual.
  const qr = houseQR(A);
  const Qtb = qr.applyQt(b);
  let x = qr.backSolve(Qtb.slice(0, n));
  // C5: iterative refinement on ill-conditioned LS too. For LS, the
  // refinement step is exact only if rank(A) == n (the tall block of QR
  // already separates fit from residual); we skip it when condEst flags
  // near-rank-deficiency.
  if (qr.condEst > ILL_COND_REFINE_THRESHOLD && isFinite(qr.condEst)) {
    x = _qrIterativeRefine(A, b, x, qr, maxRefine);
  }
  return x;
}

// Numerical Jacobian via finite differences.
//
// Forward differences (default) need a baseline F(v). Most callers (newtonSolve)
// have just evaluated F at the current iterate to test convergence; passing it
// in as `F0Opt` saves one residual evaluation per Newton step — roughly a
// 1 / (n + 1) speedup, where n = packPhi length (typically 10–20 for QDs).
//
// Central differences (mode === 'central', #11) cost 2n evaluations (F0 is
// unused) but their truncation error is O(eps²) instead of O(eps), giving
// ~1e-10 Jacobian entries vs ~1e-7. newtonSolve switches to this mode
// automatically once a step's condEst flags ill-conditioning (a forming cusp),
// where the more accurate Jacobian lets Newton drive the residual far lower.
function numericalJacobian(v, evalF, eps = DEFAULT_FD_EPS, F0Opt, mode = 'forward') {
  const n = v.length;
  if (mode === 'central') {
    let m = -1;
    let J = null;
    for (let j = 0; j < n; j++) {
      const vPlus = v.slice();  vPlus[j] += eps;
      const vMinus = v.slice(); vMinus[j] -= eps;
      const Fp = evalF(vPlus);
      const Fm = evalF(vMinus);
      if (J === null) { m = Fp.length; J = Array.from({ length: m }, () => new Array(n)); }
      const inv2e = 1 / (2 * eps);
      for (let i = 0; i < m; i++) J[i][j] = (Fp[i] - Fm[i]) * inv2e;
    }
    return J || [];
  }
  const F0 = F0Opt || evalF(v);
  const m = F0.length;
  const J = Array.from({ length: m }, () => new Array(n));
  for (let j = 0; j < n; j++) {
    const vPlus = v.slice();
    vPlus[j] += eps;
    const Fp = evalF(vPlus);
    for (let i = 0; i < m; i++) J[i][j] = (Fp[i] - F0[i]) / eps;
  }
  return J;
}

// --------- Newton-Raphson with Armijo backtracking + deflation ------------
function newtonSolve(phi_init, hData, options = {}) {
  // Family is resolved from the initial phi's tag (or the legacy unbounded
  // flag). Every per-step pack/unpack/residual is routed through the same
  // family for consistency.
  const fam = _resolveFamily(phi_init);

  const {
    maxIter = 80,
    tolerance = 1e-10,
    finiteDiffEps = DEFAULT_FD_EPS,
    armijoFactor = 1e-4,
    backtrackMax = 25,
    minStep = 1e-12,
    enforceInDisk = fam.enforceInDisk,
    enforceOutDisk = fam.enforceOutDisk,
    jacobianFn = numericalJacobian,
    deflationRoots = [],
    deflationAlpha = 1,
    deflationP     = 2,
    // Near-cusp accuracy (#11). jacobianMode:
    //   'auto'    (default) — forward differences until a step's QR condEst
    //              exceeds CENTRAL_DIFF_COND_TRIGGER (a forming cusp / fold),
    //              then central differences for the rest of the solve.
    //   'central' — central differences from the first step (the c*-polish path).
    //   'forward' — never upgrade (legacy behaviour, for A/B and perf-critical use).
    // centralDiffEps / maxRefine tune the high-accuracy path; defaults leave
    // well-conditioned solves unchanged.
    jacobianMode  = 'auto',
    centralDiffEps = CENTRAL_FD_EPS,
    maxRefine      = MAX_REFINE_STEPS,
  } = options;

  const template = phi_init;
  let v = fam.packPhi(phi_init);
  // C1: singular-Jacobian recovery telemetry. Each entry records
  // { iter, condEst, noiseScale, Fnorm } so the UI / alternates panel can
  // surface "Newton perturbed by 3.2e-8 at iter 4" rather than a silent jolt.
  const recoveryEvents = [];

  // Near-cusp accuracy (#11): adaptive Jacobian scheme. `centralMode` flips on
  // once conditioning (or a stalling line search) signals a forming cusp/fold,
  // then stays on for the rest of the solve. `canUpgrade` is false when the
  // caller pinned the mode explicitly ('central' or 'forward').
  let centralMode = (jacobianMode === 'central');
  const canUpgrade = (jacobianMode === 'auto');

  const evalFRaw = (vec) => fam.residual(fam.unpackPhi(vec, template), hData);
  const deflationEta = (vec) => {
    if (deflationRoots.length === 0) return 1;
    let eta = 1;
    for (const r of deflationRoots) {
      let d2 = 0;
      for (let k = 0; k < vec.length; k++) { const d = vec[k] - r[k]; d2 += d*d; }
      const d = Math.sqrt(d2);
      if (d < 1e-12) return Infinity;
      eta *= (1 + deflationAlpha / Math.pow(d, deflationP));
    }
    return eta;
  };
  const evalF = (vec) => {
    const F = evalFRaw(vec);
    const eta = deflationEta(vec);
    if (eta === 1) return F;
    if (!isFinite(eta)) return F.map(() => 1e30);
    return F.map(f => f * eta);
  };

  let F, Fnorm;
  try {
    F = evalF(v);
    Fnorm = residualNorm(F);
  } catch (e) {
    return { success: false, error: "Initial residual failed: " + e.message };
  }

  for (let iter = 0; iter < maxIter; iter++) {
    if (Fnorm < tolerance) {
      return { success: true, phi: fam.unpackPhi(v, template), iterations: iter, residual: Fnorm,
               recoveryEvents: recoveryEvents.length ? recoveryEvents : undefined };
    }
    let J;
    // Pass the just-computed F to the Jacobian builder so it doesn't re-evaluate
    // (forward mode only; central differences ignore it). centralMode raises the
    // Jacobian accuracy once conditioning warrants it.
    const jacEps = centralMode ? centralDiffEps : finiteDiffEps;
    const jacModeArg = centralMode ? 'central' : 'forward';
    try { J = jacobianFn(v, evalF, jacEps, F, jacModeArg); }
    catch (e) {
      return { success: false, error: "Jacobian failed: " + e.message, phi: fam.unpackPhi(v, template), iterations: iter,
               recoveryEvents: recoveryEvents.length ? recoveryEvents : undefined };
    }

    let delta;
    try {
      // condEst comes back from the same factorization used to solve, so the
      // central-difference upgrade costs no extra QR.
      const sol = leastSquaresWithCond(J, F.map(x => -x), maxRefine);
      delta = sol.x;
      if (canUpgrade && !centralMode && isFinite(sol.condEst) &&
          sol.condEst > CENTRAL_DIFF_COND_TRIGGER) {
        centralMode = true;   // forming cusp/fold — sharpen the Jacobian from here
      }
    } catch (e) {
      // C1: scale the singular-recovery perturbation by the current
      // residual magnitude and a coarse condition estimate (from the
      // failing QR), so we don't apply a 1e-9 jolt to an already-tiny
      // residual or a microscopic one to an obviously-divergent state.
      // Floor and ceiling keep the perturbation bounded to a useful range.
      let condEst = 1e8;
      try {
        const qrProbe = houseQR(J);
        if (isFinite(qrProbe.condEst)) condEst = qrProbe.condEst;
      } catch (_) { /* keep default condEst */ }
      const noiseScale = Math.max(1e-12, Math.min(1e-4, Fnorm * 1e-3 * Math.sqrt(condEst / 1e8)));
      let nudgedV = v.map(x => x + (Math.random() - 0.5) * noiseScale);
      try {
        const Fnudged = evalF(nudgedV);
        const Jnudged = jacobianFn(nudgedV, evalF, finiteDiffEps);
        delta = solveLinearSystem(Jnudged, Fnudged.map(x => -x));
        v = nudgedV; F = Fnudged; Fnorm = residualNorm(F); J = Jnudged;
        // Stash recovery telemetry on the result so the UI can surface it.
        recoveryEvents.push({ iter, condEst, noiseScale, Fnorm });
      } catch (e2) {
        return { success: false, error: "Singular Jacobian (recovery failed)",
                 phi: fam.unpackPhi(v, template), iterations: iter, residual: Fnorm,
                 recoveryEvents };
      }
    }

    let alpha = 1.0;
    let accepted = false;
    let v_new, F_new, Fnorm_new;

    for (let bt = 0; bt < backtrackMax; bt++) {
      v_new = v.map((x, i) => x + alpha * delta[i]);

      // Disk-clamp invariant (C4): Newton iterates for branches' z_j must
      // stay on the correct side of |z|=1 — strictly inside for bounded
      // families (φ : 𝔻 → Ω) and strictly outside for unbounded ones
      // (φ : 𝔻* → Ω). Crossing the unit circle puts the iterate in a
      // disjoint branch of the parametric φ and Newton can never recover.
      // This block enforces |z| ≤ DISK_CLAMP_IN or |z| ≥ DISK_CLAMP_OUT
      // before the residual is re-evaluated; it's load-bearing. The
      // node-test suite has a regression test that fails if this clamp
      // is removed.
      if (enforceInDisk) {
        for (let j = 0; j < template.branches.length; j++) {
          const re = v_new[2 * j], im = v_new[2 * j + 1];
          const r = Math.hypot(re, im);
          if (r > DISK_CLAMP_IN) {
            const scl = DISK_CLAMP_IN / r;
            v_new[2 * j] = re * scl; v_new[2 * j + 1] = im * scl;
          }
        }
      } else if (enforceOutDisk) {
        for (let j = 0; j < template.branches.length; j++) {
          const re = v_new[2 * j], im = v_new[2 * j + 1];
          const r = Math.hypot(re, im);
          if (r < DISK_CLAMP_OUT) {
            const scl = r > 1e-12 ? DISK_CLAMP_OUT / r : DISK_CLAMP_OUT;
            v_new[2 * j] = re * scl; v_new[2 * j + 1] = im * scl;
          }
        }
      }

      try { F_new = evalF(v_new); Fnorm_new = residualNorm(F_new); }
      catch (e) { alpha *= 0.5; continue; }

      if (Fnorm_new <= (1 - armijoFactor * alpha) * Fnorm) { accepted = true; break; }
      alpha *= 0.5;
      if (alpha < minStep) break;
    }

    if (!accepted) {
      return { success: false, error: "Line search failed at iter " + iter,
               phi: fam.unpackPhi(v, template), iterations: iter, residual: Fnorm };
    }

    // Stall signal (#11): a heavily-backtracked step (alpha ≤ 1/4) is the other
    // near-cusp signature — the forward-difference Jacobian is too coarse to take
    // a full Newton step. Sharpen the Jacobian for the remaining iterations.
    if (canUpgrade && !centralMode && alpha <= 0.25) centralMode = true;

    v = v_new; F = F_new; Fnorm = Fnorm_new;
  }

  return { success: false, error: "Max iterations exceeded",
           phi: fam.unpackPhi(v, template), iterations: maxIter, residual: Fnorm };
}

// --------- Homotopy helper: scale poles toward w_0 ------------------------
function scaleHDataPoles(hData, t, w0) {
  return {
    poles: hData.poles.map(p => ({
      a: { re: w0.re + t * (p.a.re - w0.re), im: w0.im + t * (p.a.im - w0.im) },
      principal: p.principal.map(Complex.clone),
    })),
  };
}

// --------- Homotopy helper: scale residue strength --------------------------
// Returns hData with every pole's principal coefficients scaled by `s`
// (pole LOCATIONS and any polyPart preserved). The s=0 limit collapses the
// quadrature data to the trivial domain (a near-disk for PQDs); s=1 is the
// full problem. Used as the residue-strength continuation parameter — a
// genuine deformation even for a single pole at w_0 (unlike scaleHDataPoles,
// which is vacuous when w_0 coincides with the only pole).
function scaleHDataResidues(hData, s) {
  return {
    poles: hData.poles.map(p => ({
      a: Complex.clone(p.a),
      principal: p.principal.map(c => Complex.scale(c, s)),
    })),
    polyPart: hData.polyPart ? hData.polyPart.map(c => Complex.scale(c, s)) : undefined,
  };
}

// --------- Boundary sampling + univalence ---------------------------------
function isBoundaryUnivalent(phi, samples = 500) {
  // Families with a continuous-arg sweep sampler (currently the powerQD
  // families, whose φ = (R#)^{1/α} needs αth-root branch tracking) get the
  // SAME boundary far more cheaply via the sweep than via N independent
  // `evalPhi` calls: `evalPhi_PQD` re-runs a per-point arg-continuation walk
  // (argContAt, ~24 R# evals each) at every sample, whereas the sweep unwraps
  // arg INCREMENTALLY between adjacent θ (~1 R# eval each) — a ~20× cut that
  // dominates per-pixel cost in the param-slice tab. maxExtra=0 ⇒ a uniform
  // N-point ring (no curvature refinement), which is what the self-intersection
  // test wants; sampleBoundaryAdaptive is WeakMap-cached by phi, so this also
  // reuses the boundary the plot already computed. Classical families have no
  // override and keep the cheap per-point `sampleBoundary` (unchanged).
  const fam = _resolveFamily(phi);
  let pts;
  if (fam && fam.sampleBoundary) {
    pts = sampleBoundaryAdaptive(phi, samples, 0).map(s => s.w);
  } else {
    pts = sampleBoundary(phi, samples);
  }
  return !boundarySelfIntersects(pts);
}

function sampleBoundary(phi, N) {
  const pts = new Array(N);
  for (let i = 0; i < N; i++) {
    const theta = (2 * Math.PI * i) / N;
    pts[i] = evalPhi({ re: Math.cos(theta), im: Math.sin(theta) }, phi);
  }
  return pts;
}

// Proper-crossing test for two segments p1p2, p3p4 (used by the boundary
// self-intersection / univalence check).
const SEG_PARALLEL_EPS = 1e-15;   // |cross(d1,d2)| below this ⇒ treat as parallel
const SEG_ENDPOINT_EPS = 1e-9;    // exclude near-endpoint touches: require the
                                  // crossing parameters t,u ∈ (ε, 1−ε), so shared
                                  // endpoints of adjacent edges don't count as a
                                  // self-intersection.
function segmentsCross(p1, p2, p3, p4) {
  const d1x = p2.re - p1.re, d1y = p2.im - p1.im;
  const d2x = p4.re - p3.re, d2y = p4.im - p3.im;
  const denom = d1x * d2y - d1y * d2x;
  if (Math.abs(denom) < SEG_PARALLEL_EPS) return false;
  const tx = p3.re - p1.re, ty = p3.im - p1.im;
  const t = (tx * d2y - ty * d2x) / denom;
  const u = (tx * d1y - ty * d1x) / denom;
  return t > SEG_ENDPOINT_EPS && t < 1 - SEG_ENDPOINT_EPS
      && u > SEG_ENDPOINT_EPS && u < 1 - SEG_ENDPOINT_EPS;
}

// Reference O(N²) self-intersection: tests every non-adjacent edge pair of the
// closed polyline pts[0]→…→pts[N-1]→pts[0]. Kept as the small-N fast path, the
// fallback for degenerate geometry, and the correctness oracle in node-test.
function boundarySelfIntersectsBruteForce(pts) {
  const N = pts.length;
  for (let i = 0; i < N; i++) {
    const i2 = (i + 1) % N;
    for (let j = i + 2; j < N; j++) {
      const j2 = (j + 1) % N;
      if (j2 === i) continue;
      if (segmentsCross(pts[i], pts[i2], pts[j], pts[j2])) return true;
    }
  }
  return false;
}

// O(N log N)-typical self-intersection via uniform spatial-grid bucketing.
// Each edge is hashed into the grid cells its bounding box overlaps; only edges
// that share a cell are tested with the exact `segmentsCross` predicate, so the
// verdict is IDENTICAL to the brute-force version — just far fewer pair tests
// (~O(N) for a boundary whose edges are of comparable length, which the uniform
// θ-sweep produces). The all-pairs version is O(N²) and was the dominant cost
// of `isBoundaryUnivalent` at high sample counts (e.g. the param-slice tab).
function boundarySelfIntersects(pts) {
  const N = pts.length;
  if (N < 4) return false;
  // Small N: the grid's setup overhead isn't worth it; brute force is instant.
  if (N <= 32) return boundarySelfIntersectsBruteForce(pts);

  // Bounding box + mean edge length (for the cell size).
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  let perim = 0;
  for (let i = 0; i < N; i++) {
    const p = pts[i];
    if (!Number.isFinite(p.re) || !Number.isFinite(p.im)) {
      return boundarySelfIntersectsBruteForce(pts);   // NaN/∞ ⇒ exact path
    }
    if (p.re < minx) minx = p.re; if (p.re > maxx) maxx = p.re;
    if (p.im < miny) miny = p.im; if (p.im > maxy) maxy = p.im;
    const q = pts[(i + 1) % N];
    perim += Math.hypot(q.re - p.re, q.im - p.im);
  }
  const w = maxx - minx, h = maxy - miny;
  const meanEdge = perim / N;
  // Degenerate (all points coincident / collinear-zero-area or zero edges):
  // fall back to the exact path rather than divide by zero.
  if (!(meanEdge > 0) || (!(w > 0) && !(h > 0))) return boundarySelfIntersectsBruteForce(pts);

  // Cell ≈ 2× mean edge ⇒ a handful of edges per cell. Clamp grid dimensions so
  // a pathological aspect ratio (one huge edge) can't allocate a giant grid.
  const cell = meanEdge * 2;
  const cols = Math.max(1, Math.min(2048, Math.floor(w / cell) + 1));
  const rows = Math.max(1, Math.min(2048, Math.floor(h / cell) + 1));
  const cw = w > 0 ? w / cols : 1;
  const ch = h > 0 ? h / rows : 1;

  // Bucket each edge by every cell its bbox overlaps; test against edges already
  // in that cell. `adjacent` skips the two edges sharing a vertex with edge i
  // (i±1, and the 0↔N-1 wrap) — matching the brute-force loop exactly. A pair
  // sharing >1 cell may be tested more than once; that's harmless (we return on
  // the first true crossing) and bounded by edges-per-cell.
  const buckets = new Map();
  for (let i = 0; i < N; i++) {
    const a = pts[i], b = pts[(i + 1) % N];
    const cx0 = Math.min(a.re, b.re), cx1 = Math.max(a.re, b.re);
    const cy0 = Math.min(a.im, b.im), cy1 = Math.max(a.im, b.im);
    const gx0 = Math.max(0, Math.floor((cx0 - minx) / cw));
    const gx1 = Math.min(cols - 1, Math.floor((cx1 - minx) / cw));
    const gy0 = Math.max(0, Math.floor((cy0 - miny) / ch));
    const gy1 = Math.min(rows - 1, Math.floor((cy1 - miny) / ch));
    for (let gy = gy0; gy <= gy1; gy++) {
      for (let gx = gx0; gx <= gx1; gx++) {
        const key = gy * cols + gx;
        let list = buckets.get(key);
        if (list) {
          for (let t = 0; t < list.length; t++) {
            const j = list[t];
            const adjacent = (i === j + 1) || (j === i + 1) || (i === 0 && j === N - 1) || (j === 0 && i === N - 1);
            if (adjacent) continue;
            if (segmentsCross(pts[i], pts[(i + 1) % N], pts[j], pts[(j + 1) % N])) return true;
          }
          list.push(i);
        } else {
          buckets.set(key, [i]);
        }
      }
    }
  }
  return false;
}

// Adaptive boundary refinement. Splits the longest edge until either the
// budget (maxExtra) is exhausted or the max-edge / mean-edge ratio drops
// below 3. Cache splitJ BEFORE the index-update loop — see comment.
//
// B4: WeakMap-cache results keyed by (phi, baseSamples, maxExtra). Pan
// and zoom on the QD plot don't change phi, just the viewport — but
// sampleBoundaryAdaptive used to re-evaluate ~500 φ-calls each render.
// With the cache, only the first call per phi pays.
const _boundaryCache = new WeakMap();    // phi → Map<budgetKey, samples>
function sampleBoundaryAdaptive(phi, baseSamples = 500, maxExtra = 1500) {
  const key = baseSamples + 'x' + maxExtra;
  let inner = _boundaryCache.get(phi);
  if (inner) {
    const hit = inner.get(key);
    if (hit) return hit;
  }
  // Family-specific override (Q1.3): families whose Riemann map needs
  // branch-tracking on ∂𝔻 (currently powerQD, where φ = (R#)^{1/α} and
  // a per-point principal αth root has an artificial atan2 cut) provide
  // their own sampler. Others fall through to the generic per-sample
  // path which calls `evalPhi` independently at each θ.
  const fam = _resolveFamily(phi);
  const out = fam.sampleBoundary
    ? fam.sampleBoundary(phi, baseSamples, maxExtra)
    : _sampleBoundaryAdaptiveImpl(phi, baseSamples, maxExtra);
  if (!inner) { inner = new Map(); _boundaryCache.set(phi, inner); }
  inner.set(key, out);
  return out;
}

// Curvature-aware boundary refinement (§22). Given a coarse θ-sorted ring
// `initialPts = [{theta, w, contArg?}]` (INCLUDING a θ=2π closure as the last
// entry, equal to the first point) and a family midpoint evaluator
// `evalMid(thetaMid, leftPt) → {theta, w, contArg?}`, recursively bisect each
// edge while the true curve's midpoint bows farther than `relTol·diam` from the
// chord midpoint (the sagitta). Concentrates points where the curve bends
// (e.g. a PQD boundary swinging in toward the origin) instead of equalizing
// arc length the way the old `maxLen < 3·meanLen` rule did. `leftPt` is passed
// to `evalMid` so families that need continuous-arg branch tracking (the αth
// root in PQDs) can unwrap from the left neighbour. Returns the refined ring
// WITHOUT the closure duplicate; θ strictly increasing, no duplicate points.
function refineBoundaryByDeviation(initialPts, evalMid, opts = {}) {
  const n0 = initialPts.length;
  if (n0 < 2) return initialPts.slice();

  // Domain scale = coarse bounding-box diagonal (tolerances are relative to it).
  let minRe = Infinity, maxRe = -Infinity, minIm = Infinity, maxIm = -Infinity;
  for (const p of initialPts) {
    if (p.w.re < minRe) minRe = p.w.re;
    if (p.w.re > maxRe) maxRe = p.w.re;
    if (p.w.im < minIm) minIm = p.w.im;
    if (p.w.im > maxIm) maxIm = p.w.im;
  }
  let diam = Math.hypot(maxRe - minRe, maxIm - minIm);
  if (!(diam > 0)) diam = 1;

  const relTol       = opts.relTol       != null ? opts.relTol       : 1.5e-3;
  const minChordFrac = opts.minChordFrac != null ? opts.minChordFrac : 2e-4;
  const maxDepth     = opts.maxDepth     != null ? opts.maxDepth     : 12;
  const maxPoints    = opts.maxPoints    != null ? opts.maxPoints    : (n0 + 1500);
  const tol      = relTol * diam;          // max allowed sagitta
  const minChord = minChordFrac * diam;    // floor: stop subdividing tiny chords

  const out = [];
  const budget = { n: n0 };                // running total (incl. the closure)

  // In-order (left subtree, node, right subtree) recursion ⇒ θ strictly
  // increasing. Inserts the evaluated midpoint only when the chord is a poor
  // approximation; the caller pushes the right endpoint.
  function bisect(a, b, depth) {
    if (budget.n >= maxPoints || depth >= maxDepth) return;
    const chord = Complex.abs(Complex.sub(b.w, a.w));
    if (chord <= minChord) return;
    const thMid = 0.5 * (a.theta + b.theta);
    const m = evalMid(thMid, a);
    const devRe = m.w.re - 0.5 * (a.w.re + b.w.re);
    const devIm = m.w.im - 0.5 * (a.w.im + b.w.im);
    if (Math.hypot(devRe, devIm) <= tol) return;   // chord already flat enough
    budget.n++;
    bisect(a, m, depth + 1);
    out.push(m);
    bisect(m, b, depth + 1);
  }

  out.push(initialPts[0]);
  for (let i = 0; i < n0 - 1; i++) {
    bisect(initialPts[i], initialPts[i + 1], 0);
    out.push(initialPts[i + 1]);
  }
  return out.slice(0, out.length - 1);   // drop the closure duplicate
}

function _sampleBoundaryAdaptiveImpl(phi, baseSamples = 500, maxExtra = 1500) {
  const N0 = baseSamples;
  const pts = [];
  for (let i = 0; i < N0; i++) {
    const theta = (2 * Math.PI * i) / N0;
    pts.push({ theta, w: evalPhi({ re: Math.cos(theta), im: Math.sin(theta) }, phi) });
  }
  pts.push({ theta: 2 * Math.PI, w: { ...pts[0].w } });
  const evalMid = (thMid) => ({
    theta: thMid,
    w: evalPhi({ re: Math.cos(thMid), im: Math.sin(thMid) }, phi),
  });
  return refineBoundaryByDeviation(pts, evalMid, { maxPoints: N0 + maxExtra });
}

// --------- Binomial helper (used by QD/UQD/LQD identity verifiers) -------
function binomialCoeff(n, k) {
  if (k < 0 || k > n) return 0;
  let r = 1;
  for (let i = 0; i < k; i++) r = r * (n - i) / (i + 1);
  return r;
}

// --------- Diverse initial guess (shared QD/UQD seed generator) -----------
// Goes beyond the "disk-guess + Gaussian perturbation" approach by sampling
// |z_j|, |A_{j,k}|, |F_l| log-uniformly with uniform arg, plus a Joukowski-
// flavored bias for unbounded polynomial cases.
function diverseInitialGuess(hData, norm, rng, r = 0) {
  const unbounded = !!(norm && norm.unbounded);
  const c = unbounded ? norm.c : 1;
  const logUni = (lo, hi) => Math.exp(Math.log(lo) + rng() * (Math.log(hi) - Math.log(lo)));
  const uni    = (lo, hi) => lo + rng() * (hi - lo);

  const zMin = unbounded ? 1.05 : 0.05;
  const zMax = unbounded ? 30   : 0.95;

  const branches = hData.poles.map((pole, jIdx) => {
    const z_mag = logUni(zMin, zMax);
    const a_arg = Complex.arg(pole.a);
    const z_arg = a_arg + uni(-Math.PI, Math.PI) * (0.5 + 0.5 * Math.tanh(r / 5));
    const z = { re: z_mag * Math.cos(z_arg), im: z_mag * Math.sin(z_arg) };
    const A = pole.principal.map(() => {
      const m = logUni(0.01, 10.0);
      const a = uni(0, 2 * Math.PI);
      return { re: m * Math.cos(a), im: m * Math.sin(a) };
    });
    return { z, A };
  });

  const polyA = [];
  if (unbounded && hData.polyPart && hData.polyPart.length > 0) {
    const mInf = hData.polyPart.length - 1;
    for (let l = 0; l <= mInf; l++) {
      let m, a;
      if (l === 1 && mInf >= 1) {
        const sign = rng() < 0.5 ? 1 : -1;
        m = c * (0.2 + 0.6 * rng());
        a = sign > 0 ? uni(-0.4, 0.4) : uni(Math.PI - 0.4, Math.PI + 0.4);
      } else {
        m = logUni(0.01, 5.0);
        a = uni(0, 2 * Math.PI);
      }
      polyA.push({ re: m * Math.cos(a), im: m * Math.sin(a) });
    }
  }
  if (unbounded) {
    return { unbounded: true, c, w0: undefined, polyA, branches };
  } else {
    return { unbounded: false, c: undefined, w0: Complex.clone(norm.w0), polyA: [], branches };
  }
}

// --------- Schema runtime --------------------------------------------------
// Each family may declare a `schema` array describing the layout of its
// packed real vector. The runtime synthesizes packPhi / unpackPhi / clamp
// from the declaration, so adding a new family's unknowns becomes a schema
// edit rather than six lines of bespoke pack/unpack.
//
// Entry kinds (all currently in use by the LQD-singular / unbounded-LQD
// families):
//
//   { kind: 'complex',     name, clamp?: {side, cap, minR, maxR} }
//       — a top-level complex stored at phi[name]
//   { kind: 'complexList', name }
//       — a list of complex stored at phi[name]; length taken from
//         template[name].length
//   { kind: 'branchesZ',   clamp? }
//       — per-branch .z complex; creates the branch objects
//   { kind: 'branchesA' }
//       — per-branch .A list; length per branch taken from template
//
// Clamp options (only `complex` and `branchesZ` entries support `clamp`):
//   side='in'  : keep |c| ≤ cap; optional `minR` pushes outward to ≥ minR.
//   side='out' : keep |c| ≥ cap; optional `maxR` pulls inward to ≤ maxR.
//
// Convention: branchesZ MUST appear before branchesA in the schema (the
// branch objects are created when branchesZ is unpacked; A entries just
// fill them in).
//
// The schema is OPTIONAL — families that don't declare one keep their
// hand-written packPhi / unpackPhi. The schema-driven clamp is invoked
// at the end of unpackPhiBySchema, in addition to Newton's per-step
// enforceInDisk / enforceOutDisk clamp on branches.
// =============================================================================
function _packEntry(entry, phi, v) {
  if (entry.kind === 'complex') {
    const c = phi[entry.name] || { re: 0, im: 0 };
    v.push(c.re, c.im);
  } else if (entry.kind === 'complexList') {
    for (const c of (phi[entry.name] || [])) v.push(c.re, c.im);
  } else if (entry.kind === 'branchesZ') {
    for (const br of phi.branches) v.push(br.z.re, br.z.im);
  } else if (entry.kind === 'branchesA') {
    for (const br of phi.branches) for (const a of br.A) v.push(a.re, a.im);
  } else {
    throw new Error("packPhiBySchema: unknown kind: " + entry.kind);
  }
}

function packPhiBySchema(phi, schema) {
  const v = [];
  for (const entry of schema) _packEntry(entry, phi, v);
  return v;
}

function unpackPhiBySchema(v, template, schema, postProcess) {
  const phi = {
    family: template.family,
    unbounded: !!template.unbounded,
  };
  if (template.w0 !== undefined) phi.w0 = Complex.clone(template.w0);
  if (template.c  !== undefined) phi.c  = template.c;
  if (template.q  !== undefined) phi.q  = Complex.clone(template.q);
  if (template.polyA && !schema.some(e => e.name === 'polyA')) {
    phi.polyA = template.polyA.map(Complex.clone);
  } else {
    phi.polyA = [];
  }
  phi.branches = [];

  let idx = 0;
  for (const entry of schema) {
    if (entry.kind === 'complex') {
      phi[entry.name] = { re: v[idx], im: v[idx + 1] };
      idx += 2;
    } else if (entry.kind === 'complexList') {
      const len = (template[entry.name] || []).length;
      phi[entry.name] = [];
      for (let k = 0; k < len; k++) {
        phi[entry.name].push({ re: v[idx], im: v[idx + 1] });
        idx += 2;
      }
    } else if (entry.kind === 'branchesZ') {
      for (let j = 0; j < template.branches.length; j++) {
        phi.branches.push({ z: { re: v[idx], im: v[idx + 1] }, A: [] });
        idx += 2;
      }
    } else if (entry.kind === 'branchesA') {
      for (let j = 0; j < template.branches.length; j++) {
        const len = template.branches[j].A.length;
        for (let k = 0; k < len; k++) {
          phi.branches[j].A.push({ re: v[idx], im: v[idx + 1] });
          idx += 2;
        }
      }
    }
  }
  applySchemaClamps(phi, schema);
  if (postProcess) postProcess(phi);
  return phi;
}

function applySchemaClamps(phi, schema) {
  // Clamp a complex into the side-of-unit-disk region declared by `cl`:
  //   side='in':   |c| ≤ cap     (𝔻);   optional minR pushes outward to ≥ minR
  //   side='out':  |c| ≥ cap     (𝔻*);  optional maxR pushes inward to ≤ maxR
  // The minR/maxR options handle the "degenerate-limit" cases (e.g. z_0 → 0
  // for bounded singular LQD, z_0 → ∞ for unbounded singular LQD) by
  // keeping Newton inside the well-defined regime.
  const clampComplex = (c, cl) => {
    if (!cl) return;
    const r = Math.hypot(c.re, c.im);
    if (cl.side === 'in') {
      const cap = cl.cap ?? DISK_CLAMP_IN;
      if (r > cap) { const s = cap / r; c.re *= s; c.im *= s; }
      if (cl.minR !== undefined && r < cl.minR) {
        const s = cl.minR / Math.max(r, 1e-15);
        c.re *= s; c.im *= s;
      }
    } else if (cl.side === 'out') {
      const cap = cl.cap ?? DISK_CLAMP_OUT;
      if (r < cap) {
        const s = r > 1e-12 ? cap / r : cap;
        c.re *= s; c.im *= s;
      }
      if (cl.maxR !== undefined && r > cl.maxR) {
        const s = cl.maxR / r;
        c.re *= s; c.im *= s;
      }
    }
  };
  for (const entry of schema) {
    if (entry.kind === 'complex' && entry.clamp) {
      clampComplex(phi[entry.name], entry.clamp);
    } else if (entry.kind === 'branchesZ' && entry.clamp) {
      for (const br of phi.branches) clampComplex(br.z, entry.clamp);
    }
  }
}

// --------- Family registry ------------------------------------------------
// Each family file populates Family.X = { ... } and calls registerFamily('X')
// to be inserted at the head of the dispatch order. selectFamily walks
// the order most-specific-first; boundedQD is the catch-all default.
const Family = {};
const familyDispatchOrder = [];

function selectFamily(opts) {
  for (const name of familyDispatchOrder) {
    const f = Family[name];
    if (f && f.matches && f.matches(opts)) return f;
  }
  return Family.boundedQD;       // belt-and-suspenders fallback
}

// Reconstruct the dispatch-complete `norm` from a solved phi — the inverse of
// normalizeOpts for routing/seeding purposes. Used by the reseed / background
// alt-search path, which must re-select the SAME family and re-seed it.
// Carries every field selectFamily inspects ({alpha, unbounded, lqd, singular})
// plus the values the seed functions read ({w0, c, q}); lqd/singular are
// derived from the family tag, the rest read straight off phi. Returns null for
// a missing phi.
function normFromPhi(phi) {
  if (!phi) return null;
  const norm = {};
  if (phi.w0)            norm.w0 = Complex.clone(phi.w0);
  if (phi.c != null)     norm.c = phi.c;
  if (phi.q)             norm.q = Complex.clone(phi.q);
  if (phi.alpha != null) norm.alpha = phi.alpha;
  if (phi.unbounded)     norm.unbounded = true;
  const fam = phi.family || '';
  if (fam.indexOf('LQD') >= 0)      norm.lqd = true;
  if (fam.indexOf('singular') >= 0) norm.singular = true;
  return norm;
}

// §23: is the origin inside the domain Ω of a solved phi? Used to detect the
// singular ↔ non-singular PQD transition (the boundary ∂Ω passing through 0).
// The sampled boundary is a bounded closed loop ∂Ω = ∂K for both bounded and
// unbounded families; we even-odd ray-cast (0,0) against it (count edges that
// cross the positive x-axis). For unbounded families Ω is the unbounded
// EXTERIOR of that loop, so 0 ∈ Ω ⟺ 0 lies OUTSIDE the loop — hence the
// `phi.unbounded` inversion. (Self-contained ray-cast: QD.Schwarz.pointInPolygon
// is not in the worker bundle.)
function originInsideOmega(phi, N = 256) {
  if (!phi) return false;
  const pts = sampleBoundaryAdaptive(phi, N, Math.floor(N * 1.5)).map(p => p.w);
  const m = pts.length;
  if (m < 3) return false;
  let inside = false;
  for (let i = 0, j = m - 1; i < m; j = i++) {
    const yi = pts[i].im, yj = pts[j].im;
    if ((yi > 0) !== (yj > 0)) {
      const xCross = pts[i].re + (0 - yi) / (yj - yi) * (pts[j].re - pts[i].re);
      if (xCross > 0) inside = !inside;
    }
  }
  return phi.unbounded ? !inside : inside;
}

// §c* — Robust test-point selection for the unbounded-family quadrature-identity
// verifiers. Each verifier checks the identity at test functions f(w)=1/(w−b)^k
// for points b in the bounded complement K (the "hole" Ω encircles). A b chosen
// too close to ∂Ω makes the LHS contour integral ∮ f·… dw under-resolved, and a
// b too close to a pole a_j of h makes the RHS residue terms C/(a_j−b)^{k+s−1}
// blow up; either way that one point's relative error spuriously dominates
// maxRelDiff and a GENUINE quadrature domain reads as identity-failing. That false
// negative is the root cause of the c* under-estimate (HANDOFF): the
// `centroid + 0.18·maxDev` placement the non-singular verifiers used was
// geometry-blind, so as the gauge c grows and the hole distorts a test point
// drifts onto a pole / the boundary.
//
// We pick b's that are provably inside K (even-odd ray-cast against the sampled
// boundary polygon) and rank them by CLEARANCE = min(distance to ∂Ω, distance to
// the nearest pole of h), keeping the best-cleared `numTestPoints`. This is the
// generalisation of the ranking the SINGULAR unbounded verifiers already used
// (verifyQuadratureIdentity_UQDLS / _UPQDS), now shared so all families agree.
//
//   polygonPts : Array<{re,im}>  — sampled ∂Ω (= ∂K), one closed loop
//   poles      : Array<{a:{re,im}, ...}> — h's finite poles (avoided by b)
//   opts.numTestPoints (3), opts.radii ([0.1..0.6]), opts.nAngles (12),
//   opts.avoidOriginEps (0; set >0 when 0 ∈ Ω so the origin isn't picked)
// Returns up to numTestPoints {re,im}; [] when K is too thin to clear any point
// (the caller then reports the identity as indeterminate rather than a false OK).
function chooseHoleTestPoints(polygonPts, poles, opts) {
  opts = opts || {};
  const want   = opts.numTestPoints || 3;
  const radii  = opts.radii || [0.1, 0.2, 0.3, 0.45, 0.6];
  const nAng   = opts.nAngles || 12;
  const origEps = opts.avoidOriginEps || 0;
  const m = polygonPts.length;
  if (m < 3) return [];

  let cx = 0, cy = 0;
  for (const w of polygonPts) { cx += w.re; cy += w.im; }
  cx /= m; cy /= m;
  let maxDev = 0;
  for (const w of polygonPts) {
    const d = Math.hypot(w.re - cx, w.im - cy);
    if (d > maxDev) maxDev = d;
  }
  if (!(maxDev > 0)) return [];

  // Even-odd ray-cast: is (x,y) inside the boundary polygon (= inside K)?
  const inside = (x, y) => {
    let cr = 0;
    for (let i = 0, j = m - 1; i < m; j = i++) {
      const yi = polygonPts[i].im, yj = polygonPts[j].im;
      if ((yi > y) !== (yj > y)) {
        const t = (y - yi) / (yj - yi);
        if (polygonPts[i].re + t * (polygonPts[j].re - polygonPts[i].re) > x) cr++;
      }
    }
    return (cr % 2) === 1;
  };
  const distBoundary = (x, y) => {
    let mn = Infinity;
    for (const w of polygonPts) {
      const d = Math.hypot(w.re - x, w.im - y);
      if (d < mn) mn = d;
    }
    return mn;
  };
  const distPole = (x, y) => {
    let mn = Infinity;
    for (const p of (poles || [])) {
      if (!p || !p.a) continue;
      const d = Math.hypot(p.a.re - x, p.a.im - y);
      if (d < mn) mn = d;
    }
    return mn;   // Infinity when h has no finite poles (poly-only h)
  };

  const cand = [{ re: cx, im: cy }];
  for (const frac of radii) {
    for (let i = 0; i < nAng; i++) {
      const a = 2 * Math.PI * i / nAng, r = frac * maxDev;
      cand.push({ re: cx + r * Math.cos(a), im: cy + r * Math.sin(a) });
    }
  }
  const ranked = [];
  for (const b of cand) {
    if (origEps > 0 && Math.hypot(b.re, b.im) < origEps) continue;  // origin ∈ Ω
    if (!inside(b.re, b.im)) continue;
    const clr = Math.min(distBoundary(b.re, b.im), distPole(b.re, b.im));
    if (clr > 0) ranked.push({ b, clr });
  }
  ranked.sort((p, q) => q.clr - p.clr);
  return ranked.slice(0, want).map(r => r.b);
}

function registerFamily(name) {
  if (familyDispatchOrder.indexOf(name) === -1) {
    familyDispatchOrder.unshift(name);
  }
}

// --------- Top-level solver -----------------------------------------------
// Two-phase:
//   PHASE A (primary): direct → continuation → multistart → diverse →
//     deflation. The first valid QD wins; otherwise we return the best
//     candidate so the user can see what we found.
//   PHASE B (alternates): additional restarts; keep all structurally
//     distinct valid QDs.
// §23 — Automatic singular ↔ non-singular PQD regime switching.
//
// `solveInverseQD` is the public entry (called by the UI, the warm worker, and
// the main-thread fallback). When `options.autoSwitchSingular` is set on a
// bounded OR unbounded PQD request, it delegates to `solvePQDWithAutoSwitch`,
// which solves the requested regime and — only if that regime is wrong (invalid
// result, or the origin's Ω-membership disagrees with the requested singular
// flag) — solves the conjugate regime once and returns it tagged. The decision
// is purely geometric (0 ∈ Ω ⟺ singular), so it settles and cannot ping-pong.
// Every other caller (including the recursive bootstrapW0_PQD classical solve)
// hits `_solveOnce` directly, so there is no recursion or behavior change.
function isPQDOpts(options) {
  const a = options && options.alpha;
  return Number.isFinite(a) && a > 0 && a !== 1 && !options.lqd;
}
function _solveResultValid(res) {
  return !!(res && res.success && res.primary &&
            res.primary.univalent && res.primary.identityOK);
}
function solvePQDWithAutoSwitch(hData, options) {
  const resF = _solveOnce(hData, options);
  const wantSingular = !!options.singular;
  // Requested regime is self-consistent → done in one solve (the common case).
  if (_solveResultValid(resF) &&
      originInsideOmega(resF.primary.phi) === wantSingular) {
    return resF;
  }
  // Otherwise try the conjugate regime exactly once.
  const g = Object.assign({}, options);
  delete g.autoSwitchSingular;
  g.singular = !wantSingular;
  if (!g.singular) delete g.w0;   // non-singular bootstraps its own interior w₀
  const resG = _solveOnce(hData, g);
  if (_solveResultValid(resG) &&
      originInsideOmega(resG.primary.phi) === g.singular) {
    resG.regimeSwitched = true;
    resG.switchedTo = g.singular ? 'singular' : 'nonsingular';
    return resG;
  }
  // Neither regime is cleanly valid (a measure-zero tangency) → preserve the
  // requested regime's result so the user sees the original failure.
  return resF;
}

function solveInverseQD(hData, options = {}) {
  if (options.autoSwitchSingular && isPQDOpts(options)) {
    return solvePQDWithAutoSwitch(hData, options);
  }
  return _solveOnce(hData, options);
}

function _solveOnce(hData, options = {}) {
  const family = selectFamily(options);
  let norm;
  try { norm = family.normalizeOpts(options, hData); }
  catch (e) { return { success: false, error: "solveInverseQD: " + e.message, attempts: [] }; }
  const w0 = norm.w0 || null;
  const c  = norm.c  || null;

  const numRestarts       = options.numRestarts ?? 8;
  const numDiverseSeeds   = options.numDiverseSeeds   ?? Math.max(numRestarts, 12);
  const numDeflationSeeds = options.numDeflationSeeds ?? Math.max(numRestarts,  8);
  const univalenceSamples = options.univalenceSamples ?? 500;
  const findAlternates    = options.findAlternates !== false;
  const newtonOpts        = options.newton ?? {};
  const contOpts          = options.continuation ?? {};
  const identityTol       = options.identityTol ?? 1e-6;
  const identityCheck     = options.identityCheck !== false;

  const usePhases   = options.usePhases ?? {};
  const useDirect       = usePhases.direct       !== false;
  const useContinuation = usePhases.continuation !== false && options.useContinuation !== false;
  const useMultistart   = usePhases.multistart   !== false;
  const useDiverse      = usePhases.diverse      !== false;
  const useDeflation    = usePhases.deflation    !== false;

  const deflationAlpha   = options.deflationAlpha ?? 1;
  const deflationP       = options.deflationP     ?? 2;
  const deflateFromValid = !!options.deflateFromValid;

  // Identity sampling is DECOUPLED from univalenceSamples: the quadrature-identity
  // contour integral needs far more nodes than the polygon self-intersection test
  // (sharply-peaked 1/(w−b)^k integrands, esp. for high-order / origin poles), so a
  // genuine QD is not mis-flagged identity-failing as the gauge c grows (the c*
  // under-estimate root cause). Callers that need rigor pass identitySamples (the c*
  // estimator does); otherwise it tracks univalenceSamples and the UNBOUNDED-family
  // verifiers apply their own ≥1500 floor internally — so bounded families and the
  // param-slice fast preset keep their cheap, low-N identity check unchanged.
  const identitySamples = options.identitySamples ?? univalenceSamples;
  const freshInit = () => family.initialGuess(hData, norm);
  const attachIdentity = (sol) => {
    if (!identityCheck) return sol;
    // adaptiveSamples (#11) flows through so callers can opt out of the near-cusp
    // node escalation (the c* estimator does — it gates on geometry near the
    // cusp). Undefined ⇒ default-on in the verifier.
    sol.identity = family.verifyQuadratureIdentity(sol.phi, hData,
      { numSamples: identitySamples, adaptiveSamples: options.adaptiveSamples });
    sol.identityOK = sol.identity.maxRelDiff < identityTol;
    return sol;
  };
  const isValidQD = (sol) => sol.univalent && (identityCheck ? sol.identityOK : true);

  const attempts = [];
  const candidates = [];
  const evalCandidate = (sol, method) => {
    sol.method = method;
    sol.phi = family.canonicalizePhi(sol.phi);
    sol.univalent = isBoundaryUnivalent(sol.phi, univalenceSamples);
    attachIdentity(sol);
    candidates.push(sol);
    return sol;
  };

  let primary = null;
  // Warm-start seed (drag-end fast path, 1D). When the caller supplies the
  // previous solution's phi via options.warmPhi (e.g. the last live-drag
  // result), try it as the very first Newton seed: a valid result lets us skip
  // the entire multistart pipeline. Gated to a structurally-matching family
  // (the deterministic freshInit().family acts as the family-tag probe) so the
  // seed and the candidate verifier agree — under autoSwitch this naturally
  // confines the seed to the regime it came from. Any miss falls through to the
  // normal phases below, so behavior is byte-identical to the no-seed path when
  // warmPhi is absent or unusable.
  const warmPhi = options.warmPhi || null;
  if (warmPhi) {
    let seed = null;
    try {
      const probe = freshInit();
      if (probe && warmPhi.family === probe.family) seed = clonePhi(warmPhi);
    } catch (_) { seed = null; }
    if (seed) {
      const warm = newtonSolve(seed, hData, newtonOpts);
      attempts.push({ method: "warm-start", success: warm.success, residual: warm.residual });
      if (warm.success) {
        evalCandidate(warm, "warm-start");
        if (isValidQD(candidates[candidates.length - 1])) primary = candidates[candidates.length - 1];
      }
    }
  }
  if (!primary && useDirect) {
    const direct = newtonSolve(freshInit(), hData, newtonOpts);
    attempts.push({ method: "direct", success: direct.success, residual: direct.residual });
    if (direct.success) {
      evalCandidate(direct, "direct");
      if (isValidQD(candidates[candidates.length - 1])) primary = candidates[candidates.length - 1];
    }
  }
  if (!primary && useContinuation) {
    const cont = family.continuationSolve(hData, norm, { ...contOpts, newton: newtonOpts });
    attempts.push({ method: "continuation", success: cont.success, residual: cont.residual, trace: cont.trace });
    if (cont.success) evalCandidate(cont, "continuation");
    if (candidates.length > 0) {
      const last = candidates[candidates.length - 1];
      if (isValidQD(last)) primary = last;
    }
  }
  if (!primary && useMultistart) {
    const rng = mulberry32(0xC0FFEE);
    for (let r = 0; r < numRestarts; r++) {
      const init = family.perturbedInitialGuess(hData, norm, rng, r);
      const res = newtonSolve(init, hData, newtonOpts);
      attempts.push({ method: "primary-restart-" + r, success: res.success, residual: res.residual });
      if (res.success) {
        evalCandidate(res, "primary-restart-" + r);
        if (isValidQD(candidates[candidates.length - 1])) { primary = candidates[candidates.length - 1]; break; }
      }
    }
  }
  if (!primary && useDiverse) {
    const rng = mulberry32(0xD1F1ED);
    for (let r = 0; r < numDiverseSeeds; r++) {
      const init = family.diverseInitialGuess(hData, norm, rng, r);
      const res = newtonSolve(init, hData, newtonOpts);
      attempts.push({ method: "diverse-" + r, success: res.success, residual: res.residual });
      if (res.success) {
        evalCandidate(res, "diverse-" + r);
        if (isValidQD(candidates[candidates.length - 1])) { primary = candidates[candidates.length - 1]; break; }
      }
    }
  }
  if (!primary && useDeflation && candidates.length > 0) {
    const roots = candidates
      .filter(s => (deflateFromValid || !isValidQD(s)) && s.phi)
      .map(s => family.packPhi(s.phi));
    if (roots.length > 0) {
      const rng = mulberry32(0xDEF1A7E);
      for (let r = 0; r < numDeflationSeeds; r++) {
        const init = r < numDeflationSeeds / 2
          ? family.diverseInitialGuess(hData, norm, rng, r)
          : family.perturbedInitialGuess(hData, norm, rng, r);
        const res = newtonSolve(init, hData, {
          ...newtonOpts, deflationRoots: roots, deflationAlpha, deflationP,
        });
        attempts.push({ method: "deflated-" + r, success: res.success, residual: res.residual });
        if (res.success) {
          evalCandidate(res, "deflated-" + r);
          if (isValidQD(candidates[candidates.length - 1])) { primary = candidates[candidates.length - 1]; break; }
        }
      }
    }
  }
  if (!primary && candidates.length > 0) {
    candidates.sort((a, b) => {
      const va = isValidQD(a), vb = isValidQD(b);
      if (va !== vb) return va ? -1 : 1;
      if (a.univalent !== b.univalent) return a.univalent ? -1 : 1;
      const ai = a.identity ? a.identity.maxRelDiff : Infinity;
      const bi = b.identity ? b.identity.maxRelDiff : Infinity;
      if (ai !== bi) return ai - bi;
      return (a.residual ?? Infinity) - (b.residual ?? Infinity);
    });
    primary = candidates[0];
  }
  if (!primary) {
    return { success: false, error: "No algebraic root found by direct, continuation, or multistart",
             attempts, w0, c, unbounded: !!norm.unbounded };
  }

  const solutions = [primary];
  if (findAlternates) {
    const rng = mulberry32(0xBEEFC0DE);
    for (let r = 0; r < numRestarts; r++) {
      const init = family.perturbedInitialGuess(hData, norm, rng, r);
      const res = newtonSolve(init, hData, newtonOpts);
      attempts.push({ method: "alt-" + r, success: res.success, residual: res.residual });
      if (res.success) {
        res.phi = family.canonicalizePhi(res.phi);
        res.univalent = isBoundaryUnivalent(res.phi, univalenceSamples);
        attachIdentity(res);
        const isNew = solutions.every(s => !phisEquivalent(s.phi, res.phi));
        if (isNew) solutions.push({ ...res, method: "restart" });
      }
    }
  }
  solutions.sort((a, b) => {
    const va = isValidQD(a), vb = isValidQD(b);
    if (va !== vb) return va ? -1 : 1;
    if (a.univalent !== b.univalent) return a.univalent ? -1 : 1;
    const ai = a.identity ? a.identity.maxRelDiff : Infinity;
    const bi = b.identity ? b.identity.maxRelDiff : Infinity;
    if (ai !== bi) return ai - bi;
    return (a.residual ?? Infinity) - (b.residual ?? Infinity);
  });

  return {
    success: true,
    primary: solutions[0],
    alternates: solutions.slice(1),
    attempts,
    w0, c, unbounded: !!norm.unbounded,
  };
}

// --------- Background alternates search -----------------------------------
function searchAlternates(hData, norm, knownSolutions, options = {}) {
  const {
    numRestarts       = 16,
    seed              = 0xBEEF0001,
    newton            = {},
    univalenceSamples = 500,
    identityTol       = 1e-6,
    diverseFraction   = 0.5,
    deflateFromKnown  = true,
    deflationAlpha    = 1,
    deflationP        = 2,
  } = options;

  const family = selectFamily(norm);
  const rng = mulberry32(seed);
  const found = [];
  const deflationRoots = deflateFromKnown
    ? knownSolutions.filter(s => s && s.phi).map(s => family.packPhi(s.phi))
    : [];

  for (let r = 0; r < numRestarts; r++) {
    const useDiverse = rng() < diverseFraction;
    const init = useDiverse
      ? family.diverseInitialGuess(hData, norm, rng, r)
      : family.perturbedInitialGuess(hData, norm, rng, r % 5);
    const res = newtonSolve(init, hData, {
      ...newton, deflationRoots, deflationAlpha, deflationP,
    });
    if (res.success) {
      res.phi = family.canonicalizePhi(res.phi);
      res.univalent = isBoundaryUnivalent(res.phi, univalenceSamples);
      res.identity = family.verifyQuadratureIdentity(res.phi, hData, { numSamples: univalenceSamples });
      res.identityOK = res.identity.maxRelDiff < identityTol;
      const all = knownSolutions.concat(found);
      const isNew = all.every(s => !phisEquivalent(s.phi, res.phi));
      if (isNew) found.push({ ...res, method: "background-restart" });
    }
  }
  return found;
}

// Quick deterministic RNG.
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Compare two solved phis up to branch reordering.
function phisEquivalent(a, b, tol = 1e-4) {
  if (a.branches.length !== b.branches.length) return false;
  const used = new Set();
  for (const ba of a.branches) {
    let bestI = -1, bestD = Infinity;
    for (let i = 0; i < b.branches.length; i++) {
      if (used.has(i)) continue;
      const bb = b.branches[i];
      if (ba.A.length !== bb.A.length) continue;
      let d = Complex.abs(Complex.sub(ba.z, bb.z));
      for (let k = 0; k < ba.A.length; k++) {
        d += Complex.abs(Complex.sub(ba.A[k], bb.A[k]));
      }
      if (d < bestD) { bestD = d; bestI = i; }
    }
    if (bestI < 0 || bestD > tol) return false;
    used.add(bestI);
  }
  return Complex.abs(Complex.sub(a.w0 || {re:0,im:0}, b.w0 || {re:0,im:0})) < tol;
}

// Realizability diagnostic for bounded PQDs (α-homotopy fold tracer). Thin
// pass-through to QD.PqdCommon.diagnosePQDRealizability, which is attached by
// solver-pqd-common.js (loaded AFTER this file) — so resolve it at call time.
// DIAGNOSTIC ONLY: runs many nested Newton solves; call only on a FAILED
// bounded-PQD solve to explain non-realizability (see the helper's header).
function diagnosePQDRealizability(hData, opts) {
  const pc = _exports.PqdCommon;
  if (!pc || !pc.diagnosePQDRealizability) {
    throw new Error("diagnosePQDRealizability: solver-pqd-common.js not loaded");
  }
  return pc.diagnosePQDRealizability(hData, opts);
}

// Live (drag) solve core (Tier-2 of the pole-drag perf work). This is the
// cheap warm-start path the inverse tab runs on every drag frame, factored out
// of ui-solve.js's quickSolveAndRender so it can run INSIDE the primary-solver
// worker (off the main thread) as well as on the main thread when no Worker is
// available. Given a caller-chosen initial phi — a warm-started clone of the
// previous solution, or a fresh family.initialGuess — it runs ONE bounded
// Newton solve (maxIter passed via opts.newton) plus reduced-sample univalence
// + identity checks. Returns a plain, structured-clone-safe object the UI turns
// into a 'live' solution. Deliberately NO multistart / continuation / deflation
// / alternates — that's solveInverseQD's job on drag-end.
//
//   opts.newton           -> newtonSolve options (e.g. { ...preset.newton, maxIter: 30 })
//   opts.numSamples       -> sample budget for univalence + identity (default 96)
//   opts.wantOriginInside -> also report originInsideOmega(phi) (PQD regime check)
function liveSolveStep(hData, initPhi, opts = {}) {
  // A usable seed is any phi with at least one branch; the family is resolved
  // from initPhi.family (with the boundedQD/unboundedQD legacy fallback in
  // _resolveFamily), so we deliberately do NOT require an explicit .family tag.
  if (!initPhi || !Array.isArray(initPhi.branches) || initPhi.branches.length === 0) {
    return { success: false };
  }
  const newtonOpts = opts.newton || {};
  const numSamples = opts.numSamples || 96;
  const family = _resolveFamily(initPhi);
  let res;
  try {
    res = newtonSolve(initPhi, hData, newtonOpts);
  } catch (e) {
    return { success: false, error: String((e && e.message) || e) };
  }
  if (!res.success) return { success: false, residual: res.residual };
  const phi = family.canonicalizePhi(res.phi);
  const univalent = isBoundaryUnivalent(phi, numSamples);
  const identity = family.verifyQuadratureIdentity(phi, hData, { numSamples });
  const identityOK = identity.maxRelDiff < 1e-6;
  const originInside = opts.wantOriginInside ? originInsideOmega(phi) : undefined;
  return {
    success: true,
    phi,
    univalent,
    identity,
    identityOK,
    originInside,
    residual: res.residual,
    iterations: res.iterations,
  };
}

// --------- Exports --------------------------------------------------------
// Family files may attach additional helpers to QD after this object is set
// (e.g. QD.diskInitialGuess from solver-qd.js, QD.unboundedInitialGuess from
// solver-uqd.js). Search the family files for `QD.foo = ...` if you need the
// full list at any point — it isn't centralized here on purpose, since the
// optional re-exports are part of each family's surface.
const _exports = {
  Complex, Taylor,
  // Dispatchers + shared
  evalPhi, phiTaylorAt, residual, residualNorm,
  packPhi, unpackPhi, canonicalizePhi,
  clonePhi, phisEquivalent,
  newtonSolve, scaleHDataPoles, scaleHDataResidues,
  solveLinearSystem, solveLeastSquares, houseQR, numericalJacobian,
  isBoundaryUnivalent, sampleBoundary, sampleBoundaryAdaptive, refineBoundaryByDeviation,
  boundarySelfIntersects, boundarySelfIntersectsBruteForce, segmentsCross,
  chooseHoleTestPoints,
  binomialCoeff, diverseInitialGuess,
  solveInverseQD, searchAlternates, mulberry32,
  // Live (drag) solve core — one warm Newton + reduced-sample checks. Runs in
  // the primary-solver worker and the main-thread fallback (see ui-solve.js).
  liveSolveStep,
  // §23 — auto-switch helpers (origin-in-Ω detection for the singular ↔
  // non-singular PQD transition; the switch itself is internal to solveInverseQD).
  originInsideOmega,
  // Bounded-PQD realizability diagnostic (α-homotopy fold tracer; failure path).
  diagnosePQDRealizability,
  // Family registry (populated by each solver-{qd,uqd,lqd,...}.js).
  Family, selectFamily, registerFamily, normFromPhi,
  // Schema runtime — opt-in pack/unpack/clamp from declarative schema.
  packPhiBySchema, unpackPhiBySchema, applySchemaClamps,
  // Named numeric constants.
  ZERO_THRESHOLD, DISK_CLAMP_OUT, DISK_CLAMP_IN, Z0_MAX_RADIUS, DEFAULT_FD_EPS,
};
if (typeof window !== 'undefined') window.QD = _exports;
if (typeof module !== 'undefined' && module.exports) module.exports = _exports;
