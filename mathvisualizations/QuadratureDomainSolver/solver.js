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
//   solver-qd.js              — Family.boundedQD
//   solver-uqd.js             — Family.unboundedQD
//   solver-lqd.js             — Family.boundedLQD (non-singular)
//   solver-lqd-singular.js    — Family.boundedLQD_singular
//   (future) solver-uqd-lqd.js / solver-uqd-lqd-singular.js
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
// =============================================================================

// --------- Generic clone: must propagate every family-specific field ------
function clonePhi(phi) {
  return {
    family: phi.family,
    unbounded: !!phi.unbounded,
    c: phi.c,
    w0: phi.w0 ? Complex.clone(phi.w0) : { re: 0, im: 0 },
    polyA: phi.polyA ? phi.polyA.map(Complex.clone) : [],
    z0:    phi.z0    ? Complex.clone(phi.z0)    : undefined,
    gamma: phi.gamma ? Complex.clone(phi.gamma) : undefined,
    q:     phi.q     ? Complex.clone(phi.q)     : undefined,
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
function _resolveFamily(phi) {
  const name = phi.family || (phi.unbounded ? 'unboundedQD' : 'boundedQD');
  const fam = Family[name];
  if (!fam) throw new Error("solver.js: family not registered: " + name);
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

function solveLinearSystem(A, b) {
  // Standard Gauss-Jordan with partial pivoting. Square systems only.
  const n = A.length;
  if (n === 0 || A[0].length !== n) throw new Error("solveLinearSystem: not square");
  const M = A.map((row, i) => row.concat([b[i]]));
  for (let k = 0; k < n; k++) {
    let pivot = k;
    let pivotMax = Math.abs(M[k][k]);
    for (let i = k + 1; i < n; i++) {
      const v = Math.abs(M[i][k]);
      if (v > pivotMax) { pivotMax = v; pivot = i; }
    }
    if (pivotMax < 1e-15) throw new Error("singular");
    if (pivot !== k) { const tmp = M[pivot]; M[pivot] = M[k]; M[k] = tmp; }
    const div = M[k][k];
    for (let j = k; j <= n; j++) M[k][j] /= div;
    for (let i = 0; i < n; i++) {
      if (i === k) continue;
      const f = M[i][k];
      if (f === 0) continue;
      for (let j = k; j <= n; j++) M[i][j] -= f * M[k][j];
    }
  }
  return M.map(row => row[n]);
}

function solveLeastSquares(A, b) {
  // Solve (A^T A) x = A^T b via Gauss-Jordan. Works for tall (over-determined)
  // systems including the +1 row most of our family residuals carry.
  const m = A.length;
  const n = A[0].length;
  const AT = Array.from({ length: n }, (_, j) => Array.from({ length: m }, (_, i) => A[i][j]));
  const AtA = Array.from({ length: n }, () => new Array(n).fill(0));
  const Atb = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let s = 0;
      for (let k = 0; k < m; k++) s += AT[i][k] * A[k][j];
      AtA[i][j] = s;
    }
    let s = 0;
    for (let k = 0; k < m; k++) s += AT[i][k] * b[k];
    Atb[i] = s;
  }
  return solveLinearSystem(AtA, Atb);
}

// Numerical Jacobian via forward differences.
function numericalJacobian(v, evalF, eps = 1e-7) {
  const F0 = evalF(v);
  const m = F0.length;
  const n = v.length;
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
    finiteDiffEps = 1e-7,
    armijoFactor = 1e-4,
    backtrackMax = 25,
    minStep = 1e-12,
    enforceInDisk = fam.enforceInDisk,
    enforceOutDisk = fam.enforceOutDisk,
    jacobianFn = numericalJacobian,
    deflationRoots = [],
    deflationAlpha = 1,
    deflationP     = 2,
  } = options;

  const template = phi_init;
  let v = fam.packPhi(phi_init);

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
      return { success: true, phi: fam.unpackPhi(v, template), iterations: iter, residual: Fnorm };
    }
    let J;
    try { J = jacobianFn(v, evalF, finiteDiffEps); }
    catch (e) {
      return { success: false, error: "Jacobian failed: " + e.message, phi: fam.unpackPhi(v, template), iterations: iter };
    }

    let delta;
    try {
      delta = solveLeastSquares(J, F.map(x => -x));
    } catch (e) {
      const noiseScale = 1e-9;
      let nudgedV = v.map(x => x + (Math.random() - 0.5) * noiseScale);
      try {
        const Fnudged = evalF(nudgedV);
        const Jnudged = jacobianFn(nudgedV, evalF, finiteDiffEps);
        delta = solveLinearSystem(Jnudged, Fnudged.map(x => -x));
        v = nudgedV; F = Fnudged; Fnorm = residualNorm(F); J = Jnudged;
      } catch (e2) {
        return { success: false, error: "Singular Jacobian (recovery failed)",
                 phi: fam.unpackPhi(v, template), iterations: iter, residual: Fnorm };
      }
    }

    let alpha = 1.0;
    let accepted = false;
    let v_new, F_new, Fnorm_new;

    for (let bt = 0; bt < backtrackMax; bt++) {
      v_new = v.map((x, i) => x + alpha * delta[i]);

      if (enforceInDisk) {
        for (let j = 0; j < template.branches.length; j++) {
          const re = v_new[2 * j], im = v_new[2 * j + 1];
          const r = Math.hypot(re, im);
          const cap = 0.9999;
          if (r > cap) {
            const scl = cap / r;
            v_new[2 * j] = re * scl; v_new[2 * j + 1] = im * scl;
          }
        }
      } else if (enforceOutDisk) {
        for (let j = 0; j < template.branches.length; j++) {
          const re = v_new[2 * j], im = v_new[2 * j + 1];
          const r = Math.hypot(re, im);
          const cap = 1.0001;
          if (r < cap) {
            const scl = r > 1e-12 ? cap / r : cap;
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

// --------- Boundary sampling + univalence ---------------------------------
function isBoundaryUnivalent(phi, samples = 500) {
  return !boundarySelfIntersects(sampleBoundary(phi, samples));
}

function sampleBoundary(phi, N) {
  const pts = new Array(N);
  for (let i = 0; i < N; i++) {
    const theta = (2 * Math.PI * i) / N;
    pts[i] = evalPhi({ re: Math.cos(theta), im: Math.sin(theta) }, phi);
  }
  return pts;
}

function segmentsCross(p1, p2, p3, p4) {
  const d1x = p2.re - p1.re, d1y = p2.im - p1.im;
  const d2x = p4.re - p3.re, d2y = p4.im - p3.im;
  const denom = d1x * d2y - d1y * d2x;
  if (Math.abs(denom) < 1e-15) return false;
  const tx = p3.re - p1.re, ty = p3.im - p1.im;
  const t = (tx * d2y - ty * d2x) / denom;
  const u = (tx * d1y - ty * d1x) / denom;
  return t > 1e-9 && t < 1 - 1e-9 && u > 1e-9 && u < 1 - 1e-9;
}

function boundarySelfIntersects(pts) {
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

// Adaptive boundary refinement. Splits the longest edge until either the
// budget (maxExtra) is exhausted or the max-edge / mean-edge ratio drops
// below 3. Cache splitJ BEFORE the index-update loop — see comment.
function sampleBoundaryAdaptive(phi, baseSamples = 500, maxExtra = 1500) {
  const N0 = baseSamples;
  const pts = [];
  for (let i = 0; i < N0; i++) {
    const theta = (2 * Math.PI * i) / N0;
    pts.push({ theta, w: evalPhi({ re: Math.cos(theta), im: Math.sin(theta) }, phi) });
  }
  pts.push({ theta: 2 * Math.PI, w: { ...pts[0].w } });

  const edges = [];
  for (let i = 0; i < pts.length - 1; i++) {
    edges.push({ i, j: i + 1, len: Complex.abs(Complex.sub(pts[i + 1].w, pts[i].w)) });
  }
  let extra = 0;
  while (extra < maxExtra) {
    let maxI = 0, maxLen = edges[0].len;
    for (let i = 1; i < edges.length; i++) {
      if (edges[i].len > maxLen) { maxLen = edges[i].len; maxI = i; }
    }
    const meanLen = edges.reduce((s, e) => s + e.len, 0) / edges.length;
    if (maxLen < 3 * meanLen || maxLen < 1e-3) break;

    const e = edges[maxI];
    const thMid = (pts[e.i].theta + pts[e.j].theta) / 2;
    const z = { re: Math.cos(thMid), im: Math.sin(thMid) };
    const newPt = { theta: thMid, w: evalPhi(z, phi) };

    // Cache splitJ BEFORE the loop: e.j aliases edges[maxI].j and gets
    // incremented inside the loop when ed === e, which would corrupt the
    // comparison for subsequent edges.
    const splitJ = e.j;
    pts.splice(splitJ, 0, newPt);
    for (const ed of edges) {
      if (ed.i >= splitJ) ed.i++;
      if (ed.j >= splitJ) ed.j++;
    }
    const leftLen  = Complex.abs(Complex.sub(newPt.w, pts[e.i].w));
    const rightLen = Complex.abs(Complex.sub(pts[e.j].w, newPt.w));
    edges.splice(maxI, 1,
      { i: e.i,    j: splitJ, len: leftLen  },
      { i: splitJ, j: e.j,    len: rightLen });
    extra++;
  }
  return pts.slice(0, pts.length - 1);
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

// --------- Schema runtime (R3) --------------------------------------------
// Each family may declare a `schema` describing the layout of its packed
// real vector. The runtime synthesizes packPhi / unpackPhi / clamp from
// the declaration so adding a new family's unknowns is a schema edit, not
// six lines of pack/unpack and a separate clamp call inside Newton.
//
// Schema entry kinds:
//   { kind: 'complex',  name: 'z0',    clamp?: { side:'in'|'out', cap, minR } }
//   { kind: 'complex',  name: 'gamma' }
//   { kind: 'branches', name: 'branches',
//     fields: [ { kind: 'complex', name: 'z', clamp: {...} },
//               { kind: 'complexList', name: 'A', len: 'm_j' } ] }
//   { kind: 'complexList', name: 'polyA', len: 'm_inf' }     // optional, present iff template.polyA?.length
//
// The schema is OPTIONAL — families that don't declare one keep their
// hand-written packPhi/unpackPhi. (Currently only future unbounded LQDs
// would benefit most; existing families continue to work as-is.) The
// schema-driven clamp is wired in via the family's enforceInDisk flag
// plus an optional `extraClamps` callback for non-z_j coordinates.
// =============================================================================

// Schema entry kinds:
//   { kind: 'complex',     name, clamp?: {side, cap, minR} }      — top-level complex
//   { kind: 'complexList', name }                                  — list of complex (length from template[name])
//   { kind: 'branchesZ',   clamp? }                                — per-branch .z complex
//   { kind: 'branchesA' }                                          — per-branch .A list (length from template.branches[j].A)
//
// Convention: branchesZ MUST appear before branchesA. (z's create the
// branch objects; A's just fill them in.)
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
      const cap = cl.cap ?? 0.9999;
      if (r > cap) { const s = cap / r; c.re *= s; c.im *= s; }
      if (cl.minR !== undefined && r < cl.minR) {
        const s = cl.minR / Math.max(r, 1e-15);
        c.re *= s; c.im *= s;
      }
    } else if (cl.side === 'out') {
      const cap = cl.cap ?? 1.0001;
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
function solveInverseQD(hData, options = {}) {
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

  const freshInit = () => family.initialGuess(hData, norm);
  const attachIdentity = (sol) => {
    if (!identityCheck) return sol;
    sol.identity = family.verifyQuadratureIdentity(sol.phi, hData, { numSamples: univalenceSamples });
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
  if (useDirect) {
    const direct = newtonSolve(freshInit(), hData, newtonOpts);
    attempts.push({ method: "direct", success: direct.success, residual: direct.residual });
    if (direct.success) evalCandidate(direct, "direct");
    if (candidates.length > 0 && isValidQD(candidates[0])) primary = candidates[0];
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

// --------- Exports --------------------------------------------------------
const _exports = {
  Complex, Taylor,
  // Dispatchers + shared
  evalPhi, phiTaylorAt, residual, residualNorm,
  packPhi, unpackPhi, canonicalizePhi,
  clonePhi, phisEquivalent,
  newtonSolve, scaleHDataPoles,
  solveLinearSystem, solveLeastSquares, numericalJacobian,
  isBoundaryUnivalent, sampleBoundary, sampleBoundaryAdaptive,
  binomialCoeff, diverseInitialGuess,
  solveInverseQD, searchAlternates, mulberry32,
  // Family registry (populated by solver-{qd,uqd,lqd,lqd-singular}.js).
  Family, selectFamily, registerFamily,
  // Schema runtime (R3) — opt-in pack/unpack/clamp from declarative schema.
  packPhiBySchema, unpackPhiBySchema, applySchemaClamps,
  // Convenience re-exports populated by family files (kept for back-compat
  // with existing callers and tests). Each family adds its bits when it
  // loads — see solver-qd.js / solver-uqd.js.
  // QD.diskInitialGuess, QD.perturbedInitialGuess, QD.computeTargetA,
  // QD.continuationSolve, QD.verifyQuadratureIdentity         (← solver-qd.js)
  // QD.unboundedInitialGuess, QD.perturbedUnboundedInitialGuess,
  // QD.continuationInC, QD.phiLaurentAtInfinity, QD.computeTargetF,
  // QD.verifyQuadratureIdentityUnbounded                      (← solver-uqd.js)
};
if (typeof window !== 'undefined') window.QD = _exports;
if (typeof module !== 'undefined' && module.exports) module.exports = _exports;
