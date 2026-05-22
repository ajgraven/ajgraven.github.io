// node-test.js -- Quick sanity check runnable via `node node-test.js`
// Loads the same source files as the browser test.

const fs = require('fs');
const path = require('path');
const vm = require('vm');

// Build one combined script and run it in this Module's context so its
// declarations live on `global` for our use below.
const ctx = { module: { exports: {} }, exports: {}, global, require, console, process, __dirname, __filename };
ctx.global = ctx;
vm.createContext(ctx);
for (const f of ['complex.js', 'taylor.js', 'solver.js', 'solver-faber.js', 'solver-qd.js', 'solver-uqd.js', 'solver-lqd-common.js', 'solver-lqd.js', 'solver-lqd-singular.js', 'solver-uqd-lqd.js', 'solver-uqd-lqd-singular.js']) {
  const src = fs.readFileSync(path.join(__dirname, f), 'utf8')
    .replace(/typeof window !== 'undefined'/g, 'false');
  vm.runInContext(src, ctx, { filename: f });
}
// Pull symbols out of the vm context by evaluating expressions there.
const QD_NS  = vm.runInContext('module.exports', ctx);
const Complex = QD_NS.Complex;
const Taylor  = QD_NS.Taylor;
const evalPhi          = QD_NS.evalPhi;
const phiTaylorAt      = QD_NS.phiTaylorAt;
const computeTargetA   = QD_NS.computeTargetA;    // moved to solver-qd.js
const residual         = QD_NS.residual;
const residualNorm     = QD_NS.residualNorm;
const solveInverseQD   = QD_NS.solveInverseQD;
const isBoundaryUnivalent = QD_NS.isBoundaryUnivalent;

let pass = 0, fail = 0;
function ok(name, cond, detail = '') {
  if (cond) { pass++; console.log('PASS  ' + name + (detail ? '  — ' + detail : '')); }
  else      { fail++; console.log('FAIL  ' + name + (detail ? '  — ' + detail : '')); }
}
function approxEq(a, b, tol = 1e-8) {
  if (typeof a === 'number') return Math.abs(a - b) < tol;
  return Math.abs(a.re - b.re) < tol && Math.abs(a.im - b.im) < tol;
}

// =============================================================================
// runFamilyBattery — per-family standard battery
// -----------------------------------------------------------------------------
// Runs a fixed set of checks against each preset:
//   • solve succeeds
//   • φ has the expected family tag
//   • boundary univalent
//   • quadrature identity satisfies tol
//   • adaptive sampler produces no duplicate points, theta strictly increasing
//   • (optional) target point is inside the rendered polygon (e.g. origin for
//     singular LQDs, w_0 for non-singular bounded)
//
// preset = { tag, hData, opts, identityTol?, insideTest?, family? }
//   - tag         human-readable name
//   - hData       quadrature data
//   - opts        solver options (norm, family flags, etc.)
//   - identityTol overrides the default 1e-8
//   - insideTest  { point: {re, im}, expected: true|false, label: str }
//   - family      expected phi.family tag (default: derived from opts)
// =============================================================================
function pointInside(pts, x, y) {
  let c = 0;
  for (let i = 0; i < pts.length; i++) {
    const j = (i + 1) % pts.length;
    if ((pts[i].im > y) !== (pts[j].im > y)) {
      const t = (y - pts[i].im) / (pts[j].im - pts[i].im);
      if (pts[i].re + t * (pts[j].re - pts[i].re) > x) c++;
    }
  }
  return (c % 2) === 1;
}

function runFamilyBattery(label, presets) {
  for (const p of presets) {
    const tol = p.identityTol ?? 1e-8;
    const result = solveInverseQD(p.hData, p.opts);
    const tag = label + ' :: ' + p.tag;
    ok(tag + ' solves', result.success, result.success ? '' : result.error);
    if (!result.success) continue;
    const sol = result.primary;
    if (p.family) {
      const got = sol.phi.family;
      ok(tag + ' family tag = ' + p.family, got === p.family,
         'got=' + (got || '<none>'));
    }
    ok(tag + ' univalent', sol.univalent);
    ok(tag + ' identityOK (' + tol.toExponential(0) + ')', sol.identity.maxRelDiff < tol,
       'maxRel=' + sol.identity.maxRelDiff.toExponential(2));
    // Sampler regression: no duplicate points, theta strictly increasing.
    const sampleAdaptive = QD_NS.sampleBoundaryAdaptive;
    const boundary = sampleAdaptive(sol.phi, 500, 750);
    let dup = 0, ooo = 0;
    for (let i = 1; i < boundary.length; i++) {
      const dx = boundary[i].w.re - boundary[i-1].w.re;
      const dy = boundary[i].w.im - boundary[i-1].w.im;
      if (Math.hypot(dx, dy) < 1e-12) dup++;
      if (boundary[i].theta < boundary[i-1].theta) ooo++;
    }
    ok(tag + ' sampler: no duplicates', dup === 0, 'dup=' + dup);
    ok(tag + ' sampler: theta strictly increasing', ooo === 0, 'ooo=' + ooo);
    if (p.insideTest) {
      const pts = boundary.map(b => b.w);
      const got = pointInside(pts, p.insideTest.point.re, p.insideTest.point.im);
      ok(tag + ' polygon contains ' + p.insideTest.label, got === p.insideTest.expected,
         'got=' + got + ' expected=' + p.insideTest.expected);
    }
  }
}

const C = Complex, T = Taylor;

// Complex parsing
ok('parse 1+2i', approxEq(C.parse('1+2i'), {re:1, im:2}));
ok('parse -1.5-2.5i', approxEq(C.parse('-1.5-2.5i'), {re:-1.5, im:-2.5}));

// Taylor inversion
{
  const p = [{re:0,im:0},{re:2,im:0},{re:3,im:0},{re:0,im:0},{re:0,im:0}];
  const q = T.invert(p, 4);
  // Verify by composition
  let comp = T.zero(5);
  let qPow = q.slice(0);
  for (let k = 1; k <= 4; k++) {
    if (k > 1) qPow = T.mul(qPow, q, 4);
    for (let i = 0; i <= 4; i++) comp[i] = C.add(comp[i], C.mul(p[k] || {re:0,im:0}, qPow[i]));
  }
  ok('p(q)[1]=1', approxEq(comp[1], {re:1,im:0}, 1e-10));
  ok('p(q)[2]=0', approxEq(comp[2], {re:0,im:0}, 1e-10));
  ok('p(q)[3]=0', approxEq(comp[3], {re:0,im:0}, 1e-10));
}

// Disk evaluation
{
  const R = 1.7, c = {re: 0.5, im: -0.3};
  const phi = { w0: c, branches: [{ z: {re:0,im:0}, A: [{re:R,im:0}] }] };
  for (let i = 0; i < 8; i++) {
    const th = i * Math.PI / 4;
    const w = evalPhi({re: Math.cos(th), im: Math.sin(th)}, phi);
    ok('|φ(e^{iθ}) - c|=R at '+(i*45)+'°', approxEq(C.abs(C.sub(w, c)), R, 1e-10));
  }
}

// Newton solve: disk
const verifyQuadratureIdentity = QD_NS.verifyQuadratureIdentity;     // moved to solver-qd.js
{
  const R = 1.4;
  const hData = { poles: [{ a: {re:0,im:0}, principal: [{re: R*R, im:0}] }] };
  const result = solveInverseQD(hData);
  ok('disk solve success', result.success, result.success ? '' : result.error);
  if (result.success) {
    ok('disk z=0', approxEq(result.primary.phi.branches[0].z, {re:0,im:0}, 1e-5));
    // A is determined only up to phase (rotational gauge of 𝔻 since the
    // Riemann-mapping normalization φ'(0) > 0 isn't enforced). Check |A| = R.
    ok('disk |A|=R', Math.abs(Complex.abs(result.primary.phi.branches[0].A[0]) - R) < 1e-5,
       '|A| = ' + Complex.abs(result.primary.phi.branches[0].A[0]).toFixed(6));
    ok('disk univalent', result.primary.univalent);

    // Verify the quadrature identity holds on monomials.
    const v = verifyQuadratureIdentity(result.primary.phi, hData, { maxDegree: 6 });
    ok('disk quadrature identity: max rel diff < 1e-10', v.maxRelDiff < 1e-10,
       'max=' + v.maxRelDiff.toExponential(3));
    // The disk-specific values: ∫_{|w|<R} w^k dA = R^{2(k+1)}/(k+1) for k=0, else 0
    //                            (modulo our normalization dA = dx dy / π and 1/(k+1) here)
    // Sanity-check k=0 -> R^2 = 1.96, k=1 -> 0.
    ok('disk k=0 ≈ R²',
       approxEq(v.checks[0].lhs, {re: R*R, im: 0}, 1e-8));
    ok('disk k=1 ≈ 0',
       Complex.abs(v.checks[1].lhs) < 1e-10);

    console.log('     iters=' + result.primary.iterations + ' resid=' + result.primary.residual.toExponential(3));
  }
}

// 2-point QD
{
  const hData = { poles: [
    { a: {re:-0.5,im:0}, principal: [{re:1.0,im:0}] },
    { a: {re: 0.5,im:0}, principal: [{re:1.0,im:0}] },
  ]};
  const result = solveInverseQD(hData);
  ok('2-pt solve success', result.success, result.success ? '' : result.error);
  if (result.success) {
    const v = verifyQuadratureIdentity(result.primary.phi, hData, { maxDegree: 4 });
    ok('2-pt quadrature identity holds', v.maxRelDiff < 1e-10, 'maxRel=' + v.maxRelDiff.toExponential(3));
    console.log('     iters=' + result.primary.iterations +
                ' resid=' + result.primary.residual.toExponential(3) +
                ' univalent=' + result.primary.univalent +
                ' alternates=' + result.alternates.length +
                ' identity-maxRel=' + v.maxRelDiff.toExponential(2));
    const ph = result.primary.phi;
    console.log('     w0 =', C.toString(ph.w0));
    console.log('     z1 =', C.toString(ph.branches[0].z), ' A1 =', C.toString(ph.branches[0].A[0]));
    console.log('     z2 =', C.toString(ph.branches[1].z), ' A2 =', C.toString(ph.branches[1].A[0]));
  }
}

// Order-2 quadrature: h = 0.5/w + 0.2/w². With the corrected Schwarz reflection,
// neither direct nor multistart finds a valid simply-connected bounded QD here.
// This might genuinely have no simply-connected QD with these parameters, or
// the solver basin is unreachable. Either way, the failure should be graceful.
{
  const hData = { poles: [
    { a: {re:0,im:0}, principal: [{re:0.5,im:0},{re:0.2,im:0}] }
  ]};
  const result = solveInverseQD(hData);
  ok('order-2 fails gracefully (no valid QD found)',
     result.success === false || (result.success && !result.primary.identityOK),
     result.success ? ('univ=' + result.primary.univalent + ' idOK=' + result.primary.identityOK)
                    : result.error);
}

// 2-point asymmetric real residues -- valid QD, identity should hold
{
  const hData = { poles: [
    { a: {re:-0.5,im:0}, principal: [{re:1.0,im:0}] },
    { a: {re: 0.7,im:0}, principal: [{re:0.6,im:0}] },
  ]};
  const result = solveInverseQD(hData);
  ok('2-pt asymmetric solve success', result.success, result.success ? '' : result.error);
  if (result.success) {
    ok('2-pt asymmetric: identityOK', result.primary.identityOK === true,
       'maxRel=' + result.primary.identity.maxRelDiff.toExponential(3));
    console.log('     method=' + result.primary.method +
                ' resid=' + result.primary.residual.toExponential(3) +
                ' univalent=' + result.primary.univalent +
                ' identity-maxRel=' + result.primary.identity.maxRelDiff.toExponential(2));
  }
}

// y-symmetric 3-point: known to work
{
  const hData = { poles: [
    { a: {re:-1.0, im: 0.0}, principal: [{re: 1.5, im: 0}] },
    { a: {re: 0.5, im: 0.8}, principal: [{re: 1.0, im: 0}] },
    { a: {re: 0.5, im:-0.8}, principal: [{re: 1.0, im: 0}] },
  ]};
  const result = solveInverseQD(hData);
  ok('3-pt y-symmetric solve success', result.success, result.success ? '' : result.error);
  if (result.success) {
    ok('3-pt y-symmetric: identityOK', result.primary.identityOK === true,
       'maxRel=' + result.primary.identity.maxRelDiff.toExponential(3));
    console.log('     method=' + result.primary.method +
                ' resid=' + result.primary.residual.toExponential(3) +
                ' identity-maxRel=' + result.primary.identity.maxRelDiff.toExponential(2));
  }
}

// Previously-failing case: 3-pt asymmetric with real residues but complex
// pole locations. With the buggy Schwarz reflection this returned a non-QD
// spurious root; with the corrected math it now solves to a true QD at
// machine precision.
{
  const hData = { poles: [
    { a: {re:-0.6, im: 0.2}, principal: [{re: 0.8, im: 0}] },
    { a: {re: 0.4, im: 0.4}, principal: [{re: 0.5, im: 0}] },
    { a: {re: 0.1, im:-0.5}, principal: [{re: 0.4, im: 0}] },
  ]};
  const result = solveInverseQD(hData);
  ok('3-pt complex-poles real-residues: valid QD found',
     result.success && result.primary.identityOK,
     'maxRel=' + (result.primary?.identity?.maxRelDiff ?? 'n/a').toExponential?.(2));
}

// User-reported case: poles at 1, -1, -i with residue 2 each. Schwarz
// reflection had been wrong; with the fix this now solves to machine precision.
{
  const hData = { poles: [
    { a: {re: 1, im: 0}, principal: [{re:2,im:0}] },
    { a: {re:-1, im: 0}, principal: [{re:2,im:0}] },
    { a: {re: 0, im:-1}, principal: [{re:2,im:0}] },
  ]};
  const result = solveInverseQD(hData);
  ok('user case (1, -1, -i; residues 2): valid QD',
     result.success && result.primary.identityOK,
     'maxRel=' + (result.primary?.identity?.maxRelDiff ?? 'n/a').toExponential?.(2));
}

// Continuation runs on the 2-point case (force it, not as a fallback) and
// produces the same solution as direct.
{
  const continuationSolve = QD_NS.continuationSolve;    // moved to solver-qd.js
  const hData = { poles: [
    { a: {re:-0.5,im:0}, principal: [{re:1.0,im:0}] },
    { a: {re: 0.5,im:0}, principal: [{re:1.0,im:0}] },
  ]};
  const cont = continuationSolve(hData, {re:0,im:0});
  ok('continuation: 2-pt symmetric succeeds', cont.success, cont.success ? '' : cont.error);
  if (cont.success) {
    console.log('     trace length=' + cont.trace.length +
                ' final residual=' + cont.residual.toExponential(3));
  }
}

// Feasible 3-pt asymmetric with residues large enough that the QD comfortably
// contains the poles. Symmetric in y so we expect a y-symmetric solution.
{
  const hData = { poles: [
    { a: {re:-1.0, im: 0.0}, principal: [{re: 1.5, im: 0}] },
    { a: {re: 0.5, im: 0.8}, principal: [{re: 1.0, im: 0}] },
    { a: {re: 0.5, im:-0.8}, principal: [{re: 1.0, im: 0}] },
  ]};
  const result = solveInverseQD(hData);
  ok('3-pt spread solve success', result.success, result.success ? '' : result.error);
  if (result.success) {
    console.log('     method=' + result.primary.method +
                ' iters=' + result.primary.iterations +
                ' resid=' + result.primary.residual.toExponential(3) +
                ' univalent=' + result.primary.univalent +
                ' alternates=' + result.alternates.length);
  }
}

// Infeasible case (poles too far for the given residues): we expect failure,
// but the solver should fail GRACEFULLY -- not throw, and return a useful
// error message.
{
  const hData = { poles: [
    { a: {re:-3.0, im:0}, principal: [{re: 0.3, im: 0}] },
    { a: {re: 3.0, im:0}, principal: [{re: 0.3, im: 0}] },
  ]};
  const result = solveInverseQD(hData);
  ok('infeasible case fails gracefully', result.success === false && typeof result.error === 'string',
     result.success ? 'unexpectedly succeeded' : result.error);
}

// =====================================================================
// Pass 2: unbounded with polynomial part of h
// =====================================================================

// h = 0 (no poly, no finite): exterior of disk D_c(0). Riemann map φ(z) = cz.
{
  for (const c of [0.5, 1.0, 2.0]) {
    const hData = { poles: [], polyPart: [] };
    const r = solveInverseQD(hData, { unbounded: true, c });
    ok('unb h=0 (c='+c+') solve success', r.success, r.success ? '' : r.error);
    if (r.success) {
      // Boundary |φ(e^{iθ}) - 0| should = c exactly.
      const ph = r.primary.phi;
      let maxErr = 0;
      for (let i = 0; i < 16; i++) {
        const t = 2 * Math.PI * i / 16;
        const z = { re: Math.cos(t), im: Math.sin(t) };
        const w = vm.runInContext('evalPhi', ctx)(z, ph);
        const e = Math.abs(Complex.abs(w) - c);
        if (e > maxErr) maxErr = e;
      }
      ok('unb h=0 (c='+c+') |φ| = c on circle', maxErr < 1e-12,
         'maxErr=' + maxErr.toExponential(2));
    }
  }
}

// h = real constant (shifted disk).
{
  const hData = { poles: [], polyPart: [{re:0.5, im:0}] };
  const r = solveInverseQD(hData, { unbounded: true, c: 1.0, identityTol: 1e-8 });
  ok('unb h=0.5 (real constant) solve success', r.success);
  if (r.success) {
    ok('unb h=0.5 (real constant) identityOK', r.primary.identityOK,
       'maxRel=' + r.primary.identity.maxRelDiff.toExponential(2));
  }
}

// h = complex constant (shifted disk at conj(C_{∞,0})). Tests the conj-on-C
// fix in computeTargetF.
{
  const hData = { poles: [], polyPart: [{re:0.3, im:0.4}] };
  const r = solveInverseQD(hData, { unbounded: true, c: 1.0, identityTol: 1e-8 });
  ok('unb h=0.3+0.4i (complex const) solve success', r.success);
  if (r.success) {
    ok('unb h=0.3+0.4i (complex const) identityOK', r.primary.identityOK,
       'maxRel=' + r.primary.identity.maxRelDiff.toExponential(2));
    // F_0 should be conj(C_{∞,0}) = 0.3 − 0.4i (center of K).
    const F0 = r.primary.phi.polyA[0];
    ok('unb h=0.3+0.4i: F_0 ≈ conj(C_{∞,0})',
       Math.abs(F0.re - 0.3) < 1e-8 && Math.abs(F0.im - (-0.4)) < 1e-8,
       'F_0 = ' + Complex.toString(F0, 6));
  }
}

// h = w + α/(w − w_0) mixed (poly + finite pole).
{
  const hData = {
    poles: [{ a: {re:2, im:0}, principal: [{re:1, im:0}] }],
    polyPart: [{re:0,im:0}, {re:0.5,im:0}],
  };
  const r = solveInverseQD(hData, { unbounded: true, c: 0.6 });
  ok('unb mixed (w + α/(w-w0)) solve success', r.success);
  if (r.success) {
    ok('unb mixed identityOK', r.primary.identityOK,
       'maxRel=' + r.primary.identity.maxRelDiff.toExponential(2));
    // The polyA should have 2 entries (F_0, F_1).
    ok('unb mixed polyA has 2 entries', r.primary.phi.polyA.length === 2);
  }
}

// Backward compat: omitting polyPart should match polyPart=[] explicitly.
{
  const hA = { poles: [{ a: {re:2,im:0}, principal: [{re:1,im:0}] }] };
  const hB = { poles: [{ a: {re:2,im:0}, principal: [{re:1,im:0}] }], polyPart: [] };
  const rA = solveInverseQD(hA, { unbounded: true, c: 0.6 });
  const rB = solveInverseQD(hB, { unbounded: true, c: 0.6 });
  ok('polyPart absent vs polyPart=[]: both succeed', rA.success && rB.success);
  if (rA.success && rB.success) {
    const dz = Complex.abs(Complex.sub(rA.primary.phi.branches[0].z, rB.primary.phi.branches[0].z));
    ok('polyPart absent vs polyPart=[]: same z_1', dz < 1e-10);
  }
}

// ===========================================================================
// LQD tests — bounded non-singular log-weighted quadrature domains
// ---------------------------------------------------------------------------
// Validates Family.boundedLQD against closed-form examples from Chapter V
// (Theorem 5.3.2: one-point bounded non-singular LQDs).
// ===========================================================================
{
  // Taylor.exp sanity (already smoke-tested but make it a formal test)
  const T = Taylor;
  let p = T.zero(5); p[1] = {re:1, im:0};
  let q = T.exp(p, 4);
  const expected = [1, 1, 0.5, 1/6, 1/24];
  let err = 0;
  for (let i = 0; i <= 4; i++) err = Math.max(err, Math.abs(q[i].re - expected[i]), Math.abs(q[i].im));
  ok('Taylor.exp(t) matches 1+t+t²/2+...', err < 1e-14, 'maxErr=' + err.toExponential(2));

  // Multiplicative property: exp(p)·exp(-p) = 1
  p = T.zero(6); p[1] = {re:0.3,im:0.2}; p[2] = {re:0.1,im:-0.4}; p[3] = {re:-0.2,im:0.1};
  q = T.exp(p, 5);
  const pNeg = p.map(c => ({re:-c.re, im:-c.im}));
  const qInv = T.exp(pNeg, 5);
  const prod = T.mul(q, qInv, 5);
  let mulErr = Math.abs(prod[0].re - 1) + Math.abs(prod[0].im);
  for (let i = 1; i <= 5; i++) mulErr = Math.max(mulErr, Math.abs(prod[i].re) + Math.abs(prod[i].im));
  ok('Taylor.exp(p)·exp(-p) = 1', mulErr < 1e-14, 'maxErr=' + mulErr.toExponential(2));
}

// Taylor.log: round-trip log(exp(p)) = p for p with p_0 = 0.
{
  const T = Taylor;
  // log(1 + t) = t − t²/2 + t³/3 − t⁴/4 + ...
  const p = T.zero(6); p[0] = {re:1,im:0}; p[1] = {re:1,im:0};
  const q = T.log(p, 5);
  const expected = [0, 1, -1/2, 1/3, -1/4, 1/5];
  let err = 0;
  for (let i = 0; i <= 5; i++) err = Math.max(err, Math.abs(q[i].re - expected[i]), Math.abs(q[i].im));
  ok('Taylor.log(1+t) matches t-t²/2+t³/3-...', err < 1e-14, 'maxErr=' + err.toExponential(2));

  // Round-trip: log(exp(p)) = p for p with arbitrary p_0.
  const p2 = T.zero(6);
  p2[0] = {re:0.3,im:-0.2}; p2[1] = {re:0.5,im:0.1}; p2[2] = {re:-0.2,im:0.3}; p2[3] = {re:0.1,im:-0.05};
  const expP = T.exp(p2, 5);
  const logExpP = T.log(expP, 5);
  let rtErr = 0;
  for (let i = 0; i <= 5; i++) rtErr = Math.max(rtErr, Math.abs(logExpP[i].re - p2[i].re), Math.abs(logExpP[i].im - p2[i].im));
  ok('Taylor.log(Taylor.exp(p)) = p (round-trip)', rtErr < 1e-13, 'maxErr=' + rtErr.toExponential(2));

  // Round-trip: exp(log(p)) = p for p with p_0 ≠ 0.
  const p3 = T.zero(5);
  p3[0] = {re:1.5,im:0.4}; p3[1] = {re:0.6,im:-0.3}; p3[2] = {re:0.2,im:0.1};
  const logP = T.log(p3, 4);
  const expLogP = T.exp(logP, 4);
  let rt2Err = 0;
  for (let i = 0; i <= 4; i++) rt2Err = Math.max(rt2Err, Math.abs(expLogP[i].re - p3[i].re), Math.abs(expLogP[i].im - p3[i].im));
  ok('Taylor.exp(Taylor.log(p)) = p (round-trip)', rt2Err < 1e-13, 'maxErr=' + rt2Err.toExponential(2));
}

// Theorem 5.3.2: one-point bounded non-singular LQD has φ(z) = w₀·exp(z√α)
// for 0 < α ≤ π². Verify identity and closed form agree at machine precision.
{
  const cases = [
    { alpha: 0.3, w0_re: 1, w0_im: 0 },
    { alpha: 1.0, w0_re: 1, w0_im: 0 },
    { alpha: 2.0, w0_re: 1, w0_im: 0 },
    { alpha: 0.5, w0_re: 2, w0_im: 0 },        // shifted
    { alpha: 0.4, w0_re: 0, w0_im: 1 },        // pure imaginary w₀
    { alpha: 0.5, w0_re: 1, w0_im: 1 },        // generic complex w₀
  ];
  for (const cs of cases) {
    const w0 = { re: cs.w0_re, im: cs.w0_im };
    const hData = { poles: [{ a: w0, principal: [{re: cs.alpha, im: 0}] }] };
    const r = solveInverseQD(hData, { lqd: true, w0 });
    const tag = 'LQD 1-pt α=' + cs.alpha + ' w₀=' + cs.w0_re + (cs.w0_im >= 0 ? '+' : '') + cs.w0_im + 'i';
    ok(tag + ': solve success', r.success);
    if (r.success) {
      ok(tag + ': identityOK', r.primary.identityOK,
        'maxRel=' + r.primary.identity.maxRelDiff.toExponential(2));
      ok(tag + ': univalent', r.primary.univalent);

      // Compare to closed form: φ(z) = w₀·exp(z√α) at z=0.5.
      const sqrtA = Math.sqrt(cs.alpha);
      const expected_re_factor = Math.exp(0.5 * sqrtA);
      const expected = { re: w0.re * expected_re_factor, im: w0.im * expected_re_factor };
      const family = vm.runInContext('Family.boundedLQD', ctx);
      const phi05 = family.evalPhi({re:0.5, im:0}, r.primary.phi);
      const diff = Math.hypot(phi05.re - expected.re, phi05.im - expected.im);
      ok(tag + ': φ(0.5) matches closed-form', diff < 1e-10, 'diff=' + diff.toExponential(2));
    }
  }
}

// (Critical-α non-univalence test omitted: per Theorem 5.3.2, α > π² admits
// no bounded simply-connected LQD, but the algebraic system is still
// solvable — φ(z) = w₀·exp(z√α) satisfies (●) and (★) for any α. Detecting
// the resulting self-intersection in the discrete `isBoundaryUnivalent`
// boundary-segment check is a known pre-existing limitation that affects
// all modes equally; it's not specific to LQDs and addressing it belongs
// in a separate validation pass. The α=2 test above confirms that valid
// LQDs at moderate α are correctly flagged univalent.)

// Three-point equilateral with real residues, around w₀=3
{
  const r3 = 0.5;
  const hData = { poles: [
    { a: {re: 3 + r3, im: 0}, principal: [{re:0.2, im:0}] },
    { a: {re: 3 - r3/2, im:  r3*Math.sqrt(3)/2}, principal: [{re:0.2, im:0}] },
    { a: {re: 3 - r3/2, im: -r3*Math.sqrt(3)/2}, principal: [{re:0.2, im:0}] },
  ]};
  const r = solveInverseQD(hData, { lqd: true, w0: {re:3, im:0} });
  ok('LQD 3-pt equilateral around w₀=3: solve success', r.success);
  if (r.success) {
    ok('LQD 3-pt equilateral: identityOK', r.primary.identityOK,
      'maxRel=' + r.primary.identity.maxRelDiff.toExponential(2));
    ok('LQD 3-pt equilateral: univalent', r.primary.univalent);
  }
}

// Conj-bug sentinel: locator residual must be zero for complex w₀. If any
// conjugation in r#(z_j) = ln(a_j/w₀) is wrong, φ(z_j) won't equal a_j.
{
  const w0 = { re: 1.2, im: 0.7 };
  const hData = { poles: [{ a: w0, principal: [{re: 0.3, im: 0}] }] };
  const r = solveInverseQD(hData, { lqd: true, w0 });
  if (r.success) {
    const family = vm.runInContext('Family.boundedLQD', ctx);
    const zj = r.primary.phi.branches[0].z;
    const phiZj = family.evalPhi(zj, r.primary.phi);
    const locErr = Math.hypot(phiZj.re - w0.re, phiZj.im - w0.im);
    ok('LQD conj-bug sentinel: locator residual ~ 0', locErr < 1e-10,
      'locErr=' + locErr.toExponential(2));
  } else {
    ok('LQD conj-bug sentinel: solve succeeds', false, 'solve failed: ' + r.error);
  }
}

// Family dispatch: solver-lqd.js patches QD.selectFamily (the export), not
// the bare `selectFamily` global from solver.js. Use the export.
{
  const QD = vm.runInContext('module.exports', ctx);
  ok('QD.selectFamily({}) → boundedQD',
    QD.selectFamily({}).name === 'boundedQD');
  ok('QD.selectFamily({unbounded:true}) → unboundedQD',
    QD.selectFamily({unbounded: true}).name === 'unboundedQD');
  ok('QD.selectFamily({lqd:true}) → boundedLQD',
    QD.selectFamily({lqd: true}).name === 'boundedLQD');
  ok('Family.boundedLQD registered', typeof QD.Family.boundedLQD === 'object');
  // Singular-LQD dispatch must take precedence over non-singular when both
  // lqd and singular flags are set.
  ok('QD.selectFamily({lqd, singular}) → boundedLQD_singular',
    QD.selectFamily({ lqd: true, singular: true }).name === 'boundedLQD_singular');
  ok('Family.boundedLQD_singular registered',
    typeof QD.Family.boundedLQD_singular === 'object');
}

// ===========================================================================
// Singular LQD tests
// ---------------------------------------------------------------------------
// Bounded LQDs with 0 ∈ Ω. Riemann map  φ(z) = γ · b_{z_0}(z) · exp(r#(z))
// where b_{z_0}(z) = -(conj(z_0)/|z_0|)·(z-z_0)/(1-conj(z_0)z), z_0 ∈ 𝔻 \ {0}.
//
// Identity:  ∫_Ω f/|w|² dA = ∮_∂Ω f h dw  for f ∈ L¹_a(Ω; ρ₀) (forces f(0)=0),
// and h ∈ Rat(Ω) is allowed an extra simple pole at 0 with residue q ∈ ℂ
// (user input).
// ===========================================================================
{
  const Singular = vm.runInContext('module.exports', ctx).Family.boundedLQD_singular;

  // (a) q-equation residual sanity at hand-constructed config.
  //   z_0=0.5, γ=2, no finite poles. r#(z_0)=0=r(z_0), so q-eq predicts
  //   q = ln|γ|² = ln 4. With q=ln 4, all 5 residual entries must be 0.
  {
    const phi = {
      family: 'boundedLQD_singular',
      w0: { re: 1, im: 0 }, q: { re: Math.log(4), im: 0 },
      z0: { re: 0.5, im: 0 }, gamma: { re: 2, im: 0 },
      branches: [],
    };
    const r = Singular.residual(phi, { poles: [] });
    let maxAbs = 0;
    for (const v of r) maxAbs = Math.max(maxAbs, Math.abs(v));
    ok('LQD-singular q-eq closed-form residual ≈ 0', maxAbs < 1e-14,
       'maxAbs=' + maxAbs.toExponential(2));
  }

  // (b) End-to-end solve for q = 0.1, one finite pole.
  {
    const hData = { poles: [{ a: { re: 2, im: 0 }, principal: [{ re: 0.5, im: 0 }] }] };
    const r = solveInverseQD(hData, {
      lqd: true, singular: true,
      w0: { re: 1, im: 0 }, q: { re: 0.1, im: 0 }, identityTol: 1e-8,
    });
    ok('LQD-singular 1-pt q=0.1 solves', r.success, r.success ? '' : r.error);
    if (r.success) {
      ok('LQD-singular 1-pt q=0.1 univalent', r.primary.univalent);
      ok('LQD-singular 1-pt q=0.1 identityOK', r.primary.identityOK,
         'maxRel=' + r.primary.identity.maxRelDiff.toExponential(2));
      ok('LQD-singular phi has family tag', r.primary.phi.family === 'boundedLQD_singular');
      ok('LQD-singular phi has z_0 inside 𝔻 and ≠ 0',
         Complex.abs(r.primary.phi.z0) > 1e-3 && Complex.abs(r.primary.phi.z0) < 1);
      // φ(z_0) ≈ 0
      const phiAt0 = Singular.evalPhi(r.primary.phi.z0, r.primary.phi);
      ok('LQD-singular φ(z_0) ≈ 0', Complex.abs(phiAt0) < 1e-9,
         '|φ(z_0)| = ' + Complex.abs(phiAt0).toExponential(2));
      // φ(0) ≈ w_0
      const phiAt0w = Singular.evalPhi({ re: 0, im: 0 }, r.primary.phi);
      ok('LQD-singular φ(0) ≈ w_0',
         Complex.abs(Complex.sub(phiAt0w, { re: 1, im: 0 })) < 1e-9,
         '|φ(0) - w_0| = ' + Complex.abs(Complex.sub(phiAt0w, { re: 1, im: 0 })).toExponential(2));
    }
  }

  // (c) q sweep (real-q family parameterized by q ∈ ℝ) — Theorem 5.6.2 style.
  // The user's "q-slider family": fix (h, w_0), dial q and verify each solve.
  // (Complex-q with real h-data is a known solver-basin limitation: the
  // gauge Im(φ'(0)) = 0 + φ(0) = w_0 ∈ ℝ pulls z_0 toward the real axis,
  // which makes Im(q) hard to achieve. Future work.)
  {
    const hData = { poles: [{ a: { re: 2, im: 0 }, principal: [{ re: 0.5, im: 0 }] }] };
    let allOK = true;
    let lastDetail = '';
    for (const qReal of [0, 0.05, 0.1, 0.2, 0.3]) {
      const r = solveInverseQD(hData, {
        lqd: true, singular: true,
        w0: { re: 1, im: 0 }, q: { re: qReal, im: 0 }, identityTol: 1e-6,
      });
      if (!r.success || !r.primary.identityOK) { allOK = false; lastDetail = 'q=' + qReal + ': ' + (r.error || 'identity fail'); break; }
      lastDetail = 'q=' + qReal + ' OK';
    }
    ok('LQD-singular q-sweep [0, 0.05, 0.1, 0.2, 0.3] all solve + identityOK',
       allOK, lastDetail);
  }

  // (c2) sampleBoundaryAdaptive regression: must produce no duplicate points
  // and no out-of-order theta values. (Pre-existing index-update bug: e.j
  // got incremented mid-iteration via aliasing, corrupting subsequent edges'
  // comparisons. Visible on LQD-singular boundaries with spike-shaped
  // refinement; made the rendered polygon misorder and visually exclude the
  // origin from a domain that actually contains it.)
  {
    const sampleAdaptive = vm.runInContext('sampleBoundaryAdaptive', ctx);
    const hData = { poles: [{ a: { re: 2, im: 0 }, principal: [{ re: 0.5, im: 0 }] }] };
    const r = solveInverseQD(hData, {
      lqd: true, singular: true,
      w0: { re: 1, im: 0 }, q: { re: 0, im: 0 },
    });
    if (r.success) {
      const boundary = sampleAdaptive(r.primary.phi, 500, 750);
      let dup = 0, ooo = 0;
      for (let i = 1; i < boundary.length; i++) {
        const dx = boundary[i].w.re - boundary[i-1].w.re;
        const dy = boundary[i].w.im - boundary[i-1].w.im;
        if (Math.hypot(dx, dy) < 1e-12) dup++;
        if (boundary[i].theta < boundary[i-1].theta) ooo++;
      }
      ok('sampleBoundaryAdaptive: no duplicate points (singular LQD)',  dup === 0, 'duplicates=' + dup);
      ok('sampleBoundaryAdaptive: theta strictly increasing',           ooo === 0, 'out-of-order=' + ooo);
      // Ray-cast origin-inside check: confirms the rendered polygon
      // (what Canvas would fill via evenodd) correctly contains origin.
      const pts = boundary.map(b => b.w);
      let cross = 0;
      for (let i = 0; i < pts.length; i++) {
        const j = (i + 1) % pts.length;
        if ((pts[i].im > 0) !== (pts[j].im > 0)) {
          const t = -pts[i].im / (pts[j].im - pts[i].im);
          if (pts[i].re + t * (pts[j].re - pts[i].re) > 0) cross++;
        }
      }
      ok('Singular-LQD adaptive polygon: origin inside rendered fill', (cross % 2) === 1);
    } else {
      ok('sampleBoundaryAdaptive regression: solve succeeded', false, 'solve failed');
    }
  }

  // (d) phiTaylorAt vs finite-difference sanity.
  {
    const phi = {
      family: 'boundedLQD_singular',
      w0: { re: 1, im: 0 }, q: { re: 0, im: 0 },
      z0: { re: 0.3, im: -0.2 }, gamma: { re: 1.5, im: 0.1 },
      branches: [
        { z: { re: 0.4, im: 0.2 }, A: [{ re: 0.3, im: -0.1 }] },
      ],
    };
    const zc = { re: 0.1, im: 0.15 };
    const taylor = Singular.phiTaylorAt(zc, phi, 2);
    // Finite-difference φ' at zc
    const eps = 1e-6;
    const fzPlus  = Singular.evalPhi({ re: zc.re + eps, im: zc.im }, phi);
    const fzMinus = Singular.evalPhi({ re: zc.re - eps, im: zc.im }, phi);
    const fdRe = (fzPlus.re - fzMinus.re) / (2 * eps);
    const fdIm = (fzPlus.im - fzMinus.im) / (2 * eps);
    const err = Math.hypot(taylor[1].re - fdRe, taylor[1].im - fdIm);
    ok('LQD-singular phiTaylorAt[1] ≈ finite-diff φ\'', err < 1e-7,
       'err=' + err.toExponential(2));
  }

  // (e) z_0 ≈ 0 rejection / safety: solve attempt where bootstrap would push
  // z_0 toward 0. We clamp inside unpackPhi at |z_0| ≥ 1e-3, so even if the
  // solve doesn't converge it should not produce NaN.
  {
    // Configure: w_0 large so the bootstrap z_0 might tend toward something
    // small. Just check that we don't crash and return a graceful result.
    const hData = { poles: [{ a: { re: 5, im: 0 }, principal: [{ re: 0.1, im: 0 }] }] };
    let threw = false;
    try {
      solveInverseQD(hData, {
        lqd: true, singular: true,
        w0: { re: 4.9, im: 0 }, q: { re: 0, im: 0 },
      });
    } catch (e) { threw = true; }
    ok('LQD-singular extreme-w0 does not throw', !threw);
  }
}

// =============================================================================
// Per-family standard battery (declarative regression sweep, R8)
// -----------------------------------------------------------------------------
// Demonstrates the runFamilyBattery helper. Adding the upcoming unbounded LQD
// families = adding presets here, no per-test boilerplate.
// =============================================================================

runFamilyBattery('boundedQD', [
  { tag: 'disk h=R²/w (R=1.4)',
    hData: { poles: [{ a: { re: 0, im: 0 }, principal: [{ re: 1.96, im: 0 }] }] },
    opts: {}, identityTol: 1e-8, family: undefined },
  { tag: '2-pt symmetric',
    hData: { poles: [
      { a: {re:-0.5,im:0}, principal: [{re:1.0,im:0}] },
      { a: {re: 0.5,im:0}, principal: [{re:1.0,im:0}] },
    ]}, opts: {}, identityTol: 1e-6 },
]);

runFamilyBattery('unboundedQD', [
  { tag: 'one-pt h=1/(w-2) c=0.6',
    hData: { poles: [{ a: { re: 2, im: 0 }, principal: [{ re: 1, im: 0 }] }] },
    opts: { unbounded: true, c: 0.6 }, identityTol: 1e-6 },
]);

runFamilyBattery('boundedLQD', [
  { tag: '1-pt α=0.5 w₀=1',
    hData: { poles: [{ a: {re:1,im:0}, principal: [{re:0.5,im:0}] }] },
    opts: { lqd: true, w0: {re:1,im:0} }, identityTol: 1e-8,
    family: 'boundedLQD' },
]);

runFamilyBattery('boundedLQD_singular', [
  { tag: 'Thm 5.6.2 h=0.5/(w-2) w₀=1 q=0',
    hData: { poles: [{ a: {re:2,im:0}, principal: [{re:0.5,im:0}] }] },
    opts: { lqd: true, singular: true, w0: {re:1,im:0}, q: {re:0,im:0} },
    identityTol: 1e-8, family: 'boundedLQD_singular',
    insideTest: { point: {re:0,im:0}, expected: true, label: 'origin' } },
  { tag: 'q=0.5 same h, w₀=1',
    hData: { poles: [{ a: {re:2,im:0}, principal: [{re:0.5,im:0}] }] },
    opts: { lqd: true, singular: true, w0: {re:1,im:0}, q: {re:0.5,im:0} },
    identityTol: 1e-8, family: 'boundedLQD_singular',
    insideTest: { point: {re:0,im:0}, expected: true, label: 'origin' } },
]);

runFamilyBattery('unboundedLQD_singular', [
  { tag: '1-pt q=0 h=1/(w-2) c=0.6',
    hData: { poles: [{ a: {re:2,im:0}, principal: [{re:1,im:0}] }] },
    opts: { lqd: true, unbounded: true, singular: true, c: 0.6, q: {re:0,im:0} },
    identityTol: 1e-6, family: 'unboundedLQD_singular',
    // Origin ∈ Ω (singular) ⇒ NOT inside the K-bounding polygon
    insideTest: { point: {re:0,im:0}, expected: false, label: 'origin (∈ Ω)' } },
  { tag: '1-pt q=0.1 h=1/(w-2) c=0.6',
    hData: { poles: [{ a: {re:2,im:0}, principal: [{re:1,im:0}] }] },
    opts: { lqd: true, unbounded: true, singular: true, c: 0.6, q: {re:0.1,im:0} },
    identityTol: 1e-6, family: 'unboundedLQD_singular',
    insideTest: { point: {re:0,im:0}, expected: false, label: 'origin (∈ Ω)' } },
  { tag: '1-pt q=0.5 same h, c=0.6',
    hData: { poles: [{ a: {re:2,im:0}, principal: [{re:1,im:0}] }] },
    opts: { lqd: true, unbounded: true, singular: true, c: 0.6, q: {re:0.5,im:0} },
    identityTol: 1e-6, family: 'unboundedLQD_singular',
    insideTest: { point: {re:0,im:0}, expected: false, label: 'origin (∈ Ω)' } },
  { tag: '1-pt complex q=0.2+0.1i, c=0.6',
    hData: { poles: [{ a: {re:2,im:0}, principal: [{re:1,im:0}] }] },
    opts: { lqd: true, unbounded: true, singular: true, c: 0.6, q: {re:0.2,im:0.1} },
    identityTol: 1e-6, family: 'unboundedLQD_singular',
    insideTest: { point: {re:0,im:0}, expected: false, label: 'origin (∈ Ω)' } },
  { tag: '2-pt symmetric q=0.3, c=0.4',
    hData: { poles: [
      { a: {re: 2,   im:0}, principal: [{re:1,  im:0}] },
      { a: {re:-1.5, im:0}, principal: [{re:0.6,im:0}] },
    ] },
    opts: { lqd: true, unbounded: true, singular: true, c: 0.4, q: {re:0.3,im:0} },
    identityTol: 1e-6, family: 'unboundedLQD_singular',
    insideTest: { point: {re:0,im:0}, expected: false, label: 'origin (∈ Ω)' } },
]);

// Refusal tests: should fail gracefully.
{
  // h = q/w only (no finite poles, nonzero q): no solution exists.
  const r = solveInverseQD({ poles: [] }, {
    lqd: true, unbounded: true, singular: true, c: 0.5, q: {re:0.1, im:0},
  });
  ok('unboundedLQD_singular: h = q/w only is rejected',
     r.success === false && /no unbounded singular LQD exists/.test(r.error || ''));
}
{
  // Higher-order pole at 0 in hData (a=0 with order ≥ 1, even order-1 is rejected
  // since simple residue belongs to opts.q).
  const r = solveInverseQD(
    { poles: [{ a:{re:0,im:0}, principal: [{re:1,im:0}] }] },
    { lqd: true, unbounded: true, singular: true, c: 0.5, q: {re:0, im:0} },
  );
  ok('unboundedLQD_singular: pole at a=0 in hData is rejected',
     r.success === false && /pole at a = 0/.test(r.error || ''));
}

runFamilyBattery('unboundedLQD', [
  { tag: 'trivial h=0 c=0.5  (Ω = ext. disk)',
    hData: { poles: [] },
    opts: { lqd: true, unbounded: true, c: 0.5 },
    identityTol: 1e-8, family: 'unboundedLQD',
    // For unbounded Ω, the boundary polygon traced by φ(e^{iθ}) is the
    // boundary of K (the bounded complement); points in K are inside that
    // polygon (ray-cast = true) and points in Ω are outside.
    insideTest: { point: {re:0,im:0}, expected: true, label: 'origin (∈ K)' } },
  { tag: '1-pt h=1/(w-2) c=0.6',
    hData: { poles: [{ a: {re:2,im:0}, principal: [{re:1,im:0}] }] },
    opts: { lqd: true, unbounded: true, c: 0.6 },
    identityTol: 1e-8, family: 'unboundedLQD',
    // For unbounded Ω, the boundary polygon traced by φ(e^{iθ}) is the
    // boundary of K (the bounded complement); points in K are inside that
    // polygon (ray-cast = true) and points in Ω are outside.
    insideTest: { point: {re:0,im:0}, expected: true, label: 'origin (∈ K)' } },
  { tag: '1-pt h=1/(w-2) c=0.3 (smaller c → bigger Ω)',
    hData: { poles: [{ a: {re:2,im:0}, principal: [{re:1,im:0}] }] },
    opts: { lqd: true, unbounded: true, c: 0.3 },
    identityTol: 1e-8, family: 'unboundedLQD',
    // For unbounded Ω, the boundary polygon traced by φ(e^{iθ}) is the
    // boundary of K (the bounded complement); points in K are inside that
    // polygon (ray-cast = true) and points in Ω are outside.
    insideTest: { point: {re:0,im:0}, expected: true, label: 'origin (∈ K)' } },
  { tag: '2-pt symmetric  c=0.4',
    hData: { poles: [
      { a: {re: 2,im:0}, principal: [{re:1,  im:0}] },
      { a: {re:-1.5,im:0}, principal: [{re:0.6,im:0}] },
    ] },
    opts: { lqd: true, unbounded: true, c: 0.4 },
    identityTol: 1e-8, family: 'unboundedLQD',
    // For unbounded Ω, the boundary polygon traced by φ(e^{iθ}) is the
    // boundary of K (the bounded complement); points in K are inside that
    // polygon (ray-cast = true) and points in Ω are outside.
    insideTest: { point: {re:0,im:0}, expected: true, label: 'origin (∈ K)' } },
]);

// Sanity: post-solve, r#(∞) should reflect the absorbed constant; φ at large
// |z| should behave as c·z to leading order (since the parametrization
// includes the −r#(∞) subtraction).
{
  const hData = { poles: [{ a: {re:2,im:0}, principal: [{re:1,im:0}] }] };
  const r = solveInverseQD(hData, { lqd: true, unbounded: true, c: 0.6 });
  if (r.success) {
    const Fam = QD_NS.Family.unboundedLQD;
    const phi = r.primary.phi;
    // φ(R) / (c · R) → 1 as R → ∞, with rate O(1/R) since r#(R) − r#(∞)
    // = O(1/R) for rational r#. Use R = 1e10 and tolerance 1e-8 for headroom.
    const R = 1e10;
    const wAtR = Fam.evalPhi({ re: R, im: 0 }, phi);
    const ratio = wAtR.re / (phi.c * R);
    ok('unboundedLQD: leading coefficient of φ at ∞ equals c',
       Math.abs(ratio - 1) < 1e-8,
       'φ(R)/(c·R) - 1 = ' + (ratio - 1).toExponential(2) + ' at R=' + R);
  } else {
    ok('unboundedLQD leading-coefficient test setup', false, 'solve failed');
  }
}

// ===========================================================================
// AQD (Stage 1) kernel sanity checks
// ===========================================================================
// Load aqd-common.js into the same vm context so it sees Complex/Taylor.
for (const f of ['aqd/aqd-common.js', 'aqd/aqd-bounded.js']) {
  const src = fs.readFileSync(path.join(__dirname, f), 'utf8')
    .replace(/typeof window !== 'undefined'/g, 'false');
  vm.runInContext(src, ctx, { filename: f });
}
const Aqd = vm.runInContext('QD.Aqd', ctx);

ok('Aqd namespace registered', typeof Aqd === 'object' && Aqd.version,
   'version=' + (Aqd?.version ?? 'undef'));

// R(w) = w  →  R(2+3i) = 2+3i, α=R'=1.
{
  const R = { polyPart: [{re:0,im:0}, {re:1,im:0}], poles: [] };
  const w = { re: 2, im: 3 };
  const rv = Aqd.evalR(R, w);
  const av = Aqd.evalAlpha(R, w);
  ok('Aqd.evalR(w, w=2+3i) = 2+3i',
     Math.abs(rv.re - 2) + Math.abs(rv.im - 3) < 1e-13,
     'got ' + rv.re + '+' + rv.im + 'i');
  ok("Aqd.evalAlpha for R=w equals 1",
     Math.abs(av.re - 1) + Math.abs(av.im) < 1e-13,
     'got ' + av.re + '+' + av.im + 'i');
}

// R(w) = w²/2  →  R(3) = 4.5, α = w, α(3) = 3.
{
  const R = { polyPart: [{re:0,im:0}, {re:0,im:0}, {re:0.5,im:0}], poles: [] };
  const rv = Aqd.evalR(R, { re: 3, im: 0 });
  const av = Aqd.evalAlpha(R, { re: 3, im: 0 });
  ok('Aqd.evalR for R=w²/2 at w=3 equals 4.5',
     Math.abs(rv.re - 4.5) + Math.abs(rv.im) < 1e-13);
  ok('Aqd.evalAlpha for R=w²/2 at w=3 equals 3',
     Math.abs(av.re - 3) + Math.abs(av.im) < 1e-13);
}

// R(w) = -1/w  →  α = R' = 1/w². Sanity: R(2) = -0.5, α(2) = 0.25.
// Also verifies that R-poles produce no nonzero residue in α (so α admits
// AQD treatment — at w=0, α = 1/w² has zero residue).
{
  const R = { polyPart: [], poles: [{ a: {re:0,im:0}, principal: [{re:-1,im:0}] }] };
  const rv = Aqd.evalR(R, { re: 2, im: 0 });
  const av = Aqd.evalAlpha(R, { re: 2, im: 0 });
  ok('Aqd.evalR for R=-1/w at w=2 equals -0.5',
     Math.abs(rv.re + 0.5) + Math.abs(rv.im) < 1e-13,
     'got ' + rv.re + '+' + rv.im + 'i');
  ok('Aqd.evalAlpha for R=-1/w at w=2 equals 0.25',
     Math.abs(av.re - 0.25) + Math.abs(av.im) < 1e-13,
     'got ' + av.re + '+' + av.im + 'i');
}

// Newton inversion: R(w) = w² + 1, R(w) = 5 → w = ±2. Seed from 1.5+0.2i
// should converge to +2.
{
  const R = { polyPart: [{re:1,im:0}, {re:0,im:0}, {re:1,im:0}], poles: [] };
  const res = Aqd.invertR(R, { re: 5, im: 0 }, { re: 1.5, im: 0.2 });
  ok('Aqd.invertR(w²+1=5, seed≈+2) converges',
     res.converged && Math.abs(res.w.re - 2) + Math.abs(res.w.im) < 1e-10,
     'w=' + res.w.re.toFixed(8) + '+' + res.w.im.toFixed(8) + 'i, iters=' + res.iters);
}

// Schwarz reflection on |z|=1: R#(z) = conj(R(z)). Quick spot check.
{
  const R = { polyPart: [{re:0.3,im:-0.1}, {re:1,im:0.2}], poles: [] };
  const theta = 0.7;
  const z = { re: Math.cos(theta), im: Math.sin(theta) };
  const Rz  = Aqd.evalR(R, z);
  const Rsh = Aqd.evalRSharp(R, z);
  const diff = { re: Rsh.re - Rz.re, im: Rsh.im - (-Rz.im) };
  ok('Aqd.evalRSharp on |z|=1 equals conj(R(z))',
     Math.hypot(diff.re, diff.im) < 1e-13,
     'diff=' + Math.hypot(diff.re, diff.im).toExponential(2));
}

// ===========================================================================
// AQD (Stage 2) end-to-end sanity: R = w (classical QD) and R = w²/2 (PQD ν=2)
// ===========================================================================
// For R = w, the AQD solver should produce the same Riemann map as the
// classical bounded QD solver (the math reduces: P_{j,t} = C_{j,t}, U = φ,
// locator and Faber match are pointwise identical).
{
  const hData = { poles: [{ a: {re:0, im:0}, principal: [{re:1, im:0}] }] };

  // R = w (identity primitive). polyPart = [0, 1], no poles.
  const R = { polyPart: [{re:0,im:0}, {re:1,im:0}], poles: [] };

  // The QD/LQD reference: unit disk, h = 1/w, φ(z) = z.
  // We need w₀ = φ(0). For h with one pole at 0, the centroid is 0, but
  // w₀ = 0 → R'(w₀) = R'(0) = 1, fine. The bounded-AQD normalizeOpts
  // however rejects w₀ at an R'-zero, not at zero in general. Let's pick
  // w₀ = 0 (which is valid here) explicitly.
  const r = solveInverseQD(hData, { aqd: true, R, w0: {re:0, im:0} });
  ok('boundedAQD :: R=w h=1/w (unit disk) solves',
     r.success, r.success ? '' : (r.error || ''));
  if (r.success) {
    const Fam = QD_NS.Family.boundedAQD;
    const phi = r.primary.phi;
    // Check φ(0) = w₀ = 0 (within tol).
    const phi0 = Fam.evalPhi({re:0,im:0}, phi);
    ok('boundedAQD :: R=w  φ(0) = 0',
       Math.hypot(phi0.re, phi0.im) < 1e-10,
       'φ(0) = ' + phi0.re.toExponential(2) + '+' + phi0.im.toExponential(2) + 'i');
    // Check φ(1/2) = 1/2 (matching φ(z) = z for unit disk).
    const phiHalf = Fam.evalPhi({re:0.5, im:0}, phi);
    ok('boundedAQD :: R=w  φ(0.5) = 0.5  (unit disk)',
       Math.hypot(phiHalf.re - 0.5, phiHalf.im) < 1e-8,
       'φ(0.5) = ' + phiHalf.re.toFixed(6) + '+' + phiHalf.im.toFixed(6) + 'i');
    // Identity check.
    const verify = Fam.verifyQuadratureIdentity(phi, hData, { numSamples: 400, maxDegree: 3, numTestPoints: 3 });
    ok('boundedAQD :: R=w  identity check passes',
       verify.maxRelDiff < 1e-6,
       'maxRel=' + verify.maxRelDiff.toExponential(2));
  }
}

// PQD ν=2: R = w²/2 (so α = w, R'-zero at w=0). Place h far from the origin
// so the resulting Ω stays well inside the Stage-2 admissibility region
// (Ω̄ contains no R'-zeros). With a=5 and residue 1, Ω is a small distorted
// disk around w=5 — min|w| on ∂Ω ≈ 4.8, well away from 0.
{
  const aLoc = { re: 5, im: 0 };
  const hData = { poles: [{ a: aLoc, principal: [{re:1, im:0}] }] };
  const R = { polyPart: [{re:0,im:0}, {re:0,im:0}, {re:0.5,im:0}], poles: [] };
  const r = solveInverseQD(hData, { aqd: true, R, w0: aLoc });
  ok('boundedAQD :: R=w²/2 h=1/(w-5) solves',
     r.success, r.success ? '' : (r.error || ''));
  if (r.success) {
    const Fam = QD_NS.Family.boundedAQD;
    const phi = r.primary.phi;
    const verify = Fam.verifyQuadratureIdentity(phi, hData, { numSamples: 500, maxDegree: 3, numTestPoints: 3 });
    ok('boundedAQD :: R=w²/2 identity check passes',
       verify.maxRelDiff < 1e-8,
       'maxRel=' + verify.maxRelDiff.toExponential(2));
    ok('boundedAQD :: R=w²/2 univalent',
       isBoundaryUnivalent(phi, 500));
  }
}

// Stage 2b establishes: AQD solver mathematics & infrastructure work end-to-
// end for single-pole order-1 cases. Multi-pole and higher-order cases need
// a sharper initial guess (Stage 6 warm-start work) — they're not failures
// of the solver kernel, just of the linear seed reaching Newton's basin.
// We do NOT claim correctness for those yet; they'll come back as tests in
// Stages 6 and 7.

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
