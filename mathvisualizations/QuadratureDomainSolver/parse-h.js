// =============================================================================
// parse-h.js -- Custom-text-input parser for quadrature functions h(w).
//
// Exports:
//   QD.parseH(expr, math, opts) → { poles, polyCoeffs }
//   QD.formatH({ poles, polyCoeffs }) → string
//
// Lets the user type h(w) as a single expression like
//     1.5/w + 0.5/w^2
//     1/(w-2) + (1+i)/(w-2)^2 + 0.3*w^2
//     1/(w^2 - 1)
// and converts it to the structured form the inverse solver consumes:
//
//   poles:      [ { a: Complex, order: int, residues: [Complex_1, ..., Complex_order] }, ... ]
//   polyCoeffs: [ C_∞,0, C_∞,1, ..., C_∞,m ]    // empty when no polynomial part
//
// Two passes:
//
//   Phase 1 (strict per-summand walker, EXACT)
//     Splits the expression on top-level + / − into "atoms". Each atom must be
//     one of:
//        • a polynomial monomial   C · w^k                         (k ≥ 0)
//        • a pole atom             C / (w − a)^k    or  C / w^k    (k ≥ 1)
//     For pole atoms we recover `a` symbolically from the monic-shifted
//     denominator (a = −Q_{k−1} / (k · Q_k)) and verify the whole denominator
//     matches (w − a)^k within tol. So `1/w^k` gives a = 0 exactly, `1/(w−2)^3`
//     gives a = 2 exactly, etc.
//
//   Phase 2 (general-rational fallback)
//     Triggered when an atom isn't classifiable (e.g. `1/(w^2 − 1)`). The
//     expression is cross-multiplied into a single P(w)/Q(w), long-divided
//     for the polynomial part, Q is factored via Durand–Kerner, the roots
//     are clustered with a relative-tolerance heuristic, and the principal
//     part at each pole is recovered by shift-and-series-divide.
//
// In bounded mode (opts.mode === 'bounded' or 'lqd-*' families), a nonzero
// polynomial part throws — the caller should switch to unbounded.
//
// All numeric internals use the plain {re, im} Complex format from complex.js.
// =============================================================================

(function (global) {
  // Namespace plumbing mirrors the other solver modules: browser stashes on
  // window.QD; node-test's vm context exposes the QD namespace as
  // module.exports (set by solver.js). We attach parseH / formatH onto QD.
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' && module.exports ? module.exports : (global.QD || (global.QD = {})));

  // Complex resolution: the script tag for complex.js puts Complex on window
  // (browser) or sets module.exports = Complex which solver.js then re-exposes
  // as QD.Complex. Either way, we want the plain Complex object.
  const C =
    (typeof Complex !== 'undefined' ? Complex : null)
    || (QD && QD.Complex)
    || (typeof window !== 'undefined' && window.Complex)
    || (global && global.Complex);
  if (!C) throw new Error("parse-h.js: Complex namespace not found");

  // ===========================================================================
  // Polynomial-in-w helpers (ascending-power Complex[]).
  // ===========================================================================
  function polyZero()           { return [{ re: 0, im: 0 }]; }
  function polyOne()            { return [{ re: 1, im: 0 }]; }
  function polyConst(c)         { return [{ re: c.re, im: c.im }]; }
  function polyVar()            { return [{ re: 0, im: 0 }, { re: 1, im: 0 }]; }      // w
  function polyTrim(p) {
    while (p.length > 1 && Math.hypot(p[p.length - 1].re, p[p.length - 1].im) < 1e-14) p.pop();
    return p;
  }
  function polyAdd(a, b) {
    const n = Math.max(a.length, b.length);
    const out = new Array(n);
    for (let i = 0; i < n; i++) {
      const av = i < a.length ? a[i] : { re: 0, im: 0 };
      const bv = i < b.length ? b[i] : { re: 0, im: 0 };
      out[i] = C.add(av, bv);
    }
    return polyTrim(out);
  }
  function polyNeg(a) { return a.map(c => C.neg(c)); }
  function polyMul(a, b) {
    if (a.length === 0 || b.length === 0) return polyZero();
    const out = new Array(a.length + b.length - 1);
    for (let i = 0; i < out.length; i++) out[i] = { re: 0, im: 0 };
    for (let i = 0; i < a.length; i++) {
      for (let j = 0; j < b.length; j++) {
        out[i + j] = C.add(out[i + j], C.mul(a[i], b[j]));
      }
    }
    return polyTrim(out);
  }
  function polyScale(a, s) { return a.map(c => C.mul(c, s)); }
  function polyPow(a, n) {
    let out = polyOne();
    for (let i = 0; i < n; i++) out = polyMul(out, a);
    return out;
  }
  // Long-divide P / Q → { quotient, remainder } so that P = quotient*Q + remainder.
  function polyDivMod(P, Q) {
    P = P.slice(); Q = polyTrim(Q.slice());
    const dq = Q.length - 1;
    const lc = Q[dq];
    if (Math.hypot(lc.re, lc.im) < 1e-300) throw new Error("polyDivMod: divisor is zero");
    const quotient = [];
    while (P.length - 1 >= dq && P.length > 1) {
      const dp = P.length - 1;
      const k = dp - dq;
      const factor = C.div(P[dp], lc);
      quotient[k] = factor;
      for (let i = 0; i <= dq; i++) {
        P[k + i] = C.sub(P[k + i], C.mul(factor, Q[i]));
      }
      P.pop();
    }
    for (let i = 0; i < quotient.length; i++) if (!quotient[i]) quotient[i] = { re: 0, im: 0 };
    return { quotient: quotient.length ? polyTrim(quotient) : polyZero(), remainder: polyTrim(P) };
  }
  // Substitute w = a + u; return polynomial P(a+u) in u.
  function polyShift(P, a) {
    // Build (a + u)^k incrementally.
    let acc = polyZero();
    let aPlusU = [{ re: a.re, im: a.im }, { re: 1, im: 0 }];     // (a + u)
    let power = polyOne();                                       // (a+u)^0
    for (let k = 0; k < P.length; k++) {
      acc = polyAdd(acc, polyScale(power, P[k]));
      if (k + 1 < P.length) power = polyMul(power, aPlusU);
    }
    return acc;
  }
  // Truncated series division: num / den as Taylor series in u up to degree maxDeg.
  // Requires den[0] != 0. Returns coefficients [t_0, t_1, ..., t_maxDeg].
  function seriesDivide(num, den, maxDeg) {
    if (Math.hypot(den[0].re, den[0].im) < 1e-300) {
      throw new Error("seriesDivide: denominator has zero constant term");
    }
    const d0Inv = C.inv(den[0]);
    const t = new Array(maxDeg + 1);
    for (let k = 0; k <= maxDeg; k++) {
      let s = (k < num.length) ? { re: num[k].re, im: num[k].im } : { re: 0, im: 0 };
      const jMax = Math.min(k, den.length - 1);
      for (let j = 1; j <= jMax; j++) {
        s = C.sub(s, C.mul(den[j], t[k - j]));
      }
      t[k] = C.mul(s, d0Inv);
    }
    return t;
  }

  // ===========================================================================
  // Math.js → Complex constant evaluation (for sub-expressions free of w).
  // ===========================================================================
  function mjxToComplex(val) {
    if (val == null) throw new Error("null math.js value");
    if (typeof val === 'number') return { re: val, im: 0 };
    if (typeof val === 'object') {
      if ('re' in val && 'im' in val) return { re: +val.re, im: +val.im };
      if (typeof val.toJSON === 'function') {
        const j = val.toJSON();
        if (j && 're' in j && 'im' in j) return { re: +j.re, im: +j.im };
      }
    }
    throw new Error("can't convert math.js value to Complex: " + String(val));
  }

  // ===========================================================================
  // General-rational walker (port of direct/parseRationalInZ to variable `w`).
  // Returns { num, den } as polynomials in w.
  // ===========================================================================
  function accRat(node, math) {
    if (!node) throw new Error("null AST node");
    if (node.isConstantNode) return { num: [mjxToComplex(node.value)], den: polyOne() };
    if (node.isSymbolNode) {
      if (node.name === 'w') return { num: polyVar(),               den: polyOne() };
      if (node.name === 'i') return { num: [{ re: 0, im: 1 }],      den: polyOne() };
      throw new Error("unknown symbol '" + node.name + "' (only w and i are allowed)");
    }
    if (node.isParenthesisNode) return accRat(node.content, math);
    if (node.isOperatorNode) {
      const op = node.op;
      if (op === '+') {
        let acc = accRat(node.args[0], math);
        for (let i = 1; i < node.args.length; i++) acc = addR(acc, accRat(node.args[i], math));
        return acc;
      }
      if (op === '-') {
        if (node.args.length === 1) {
          const r = accRat(node.args[0], math);
          return { num: polyNeg(r.num), den: r.den };
        }
        let acc = accRat(node.args[0], math);
        for (let i = 1; i < node.args.length; i++) {
          const r = accRat(node.args[i], math);
          acc = addR(acc, { num: polyNeg(r.num), den: r.den });
        }
        return acc;
      }
      if (op === '*') {
        let acc = accRat(node.args[0], math);
        for (let i = 1; i < node.args.length; i++) acc = mulR(acc, accRat(node.args[i], math));
        return acc;
      }
      if (op === '/') {
        const lhs = accRat(node.args[0], math);
        const rhs = accRat(node.args[1], math);
        if (rhs.num.length === 0 || (rhs.num.length === 1 && Math.hypot(rhs.num[0].re, rhs.num[0].im) < 1e-300)) {
          throw new Error("division by zero");
        }
        return mulR(lhs, { num: rhs.den, den: rhs.num });
      }
      if (op === '^') {
        const en = node.args[1];
        if (!en.isConstantNode) throw new Error("exponent in '^' must be a literal integer");
        const k = (typeof en.value === 'number') ? en.value : Number(en.value);
        if (!Number.isInteger(k)) throw new Error("exponent must be an integer (got " + en.value + ")");
        const base = accRat(node.args[0], math);
        if (k === 0) return { num: polyOne(), den: polyOne() };
        if (k > 0)  return powR(base, k);
        const pos = powR(base, -k);
        if (pos.num.length === 0 || (pos.num.length === 1 && Math.hypot(pos.num[0].re, pos.num[0].im) < 1e-300)) {
          throw new Error("division by zero (negative exponent of zero subexpression)");
        }
        return { num: pos.den, den: pos.num };
      }
      throw new Error("unsupported operator: " + op);
    }
    if (node.isFunctionNode) {
      // Only allow function calls with constant args (no `w`).
      const argRats = node.args.map(a => accRat(a, math));
      const allConst = argRats.every(r => r.num.length <= 1 && r.den.length <= 1);
      if (!allConst) throw new Error("function " + node.name + " requires constant arguments");
      let val;
      try { val = node.evaluate(); }
      catch (e) { throw new Error("could not evaluate " + node.name + "(...): " + (e.message || e)); }
      return { num: [mjxToComplex(val)], den: polyOne() };
    }
    throw new Error("unsupported AST node: " + (node.type || "(unknown)"));
  }
  function addR(a, b) {
    return { num: polyAdd(polyMul(a.num, b.den), polyMul(b.num, a.den)),
             den: polyMul(a.den, b.den) };
  }
  function mulR(a, b) {
    return { num: polyMul(a.num, b.num), den: polyMul(a.den, b.den) };
  }
  function powR(r, k) {
    let n = r.num, d = r.den;
    let an = polyOne(), ad = polyOne();
    for (let i = 0; i < k; i++) { an = polyMul(an, n); ad = polyMul(ad, d); }
    return { num: an, den: ad };
  }

  // ===========================================================================
  // Phase 1: per-summand strict PFD walker.
  // ===========================================================================
  // flattenSum: returns [{ sign: +1/-1, node }, ...]. Splits at top-level + and -.
  function flattenSum(node, signIn) {
    signIn = signIn || 1;
    if (node.isParenthesisNode) return flattenSum(node.content, signIn);
    if (node.isOperatorNode && (node.op === '+' || node.op === '-')) {
      if (node.op === '-' && node.args.length === 1) {
        return flattenSum(node.args[0], -signIn);
      }
      const out = flattenSum(node.args[0], signIn);
      for (let i = 1; i < node.args.length; i++) {
        const sub = flattenSum(node.args[i], node.op === '+' ? signIn : -signIn);
        for (const s of sub) out.push(s);
      }
      return out;
    }
    return [{ sign: signIn, node }];
  }

  // Try to classify a single summand strictly. Returns one of:
  //   { kind: 'poly', degree: k, coeff: Complex }       polynomial monomial C·w^k
  //   { kind: 'pole', a: Complex, order: k, residue: Complex }   C / (w - a)^k
  // Throws on un-classifiable atoms.
  function classifySummand(node, math) {
    const rat = accRat(node, math);
    const num = polyTrim(rat.num);
    const den = polyTrim(rat.den);

    // Polynomial atom: den constant.
    if (den.length === 1) {
      const dInv = C.inv(den[0]);
      const p = num.map(c => C.mul(c, dInv));
      // Must have a single nonzero coefficient.
      let idx = -1;
      for (let k = 0; k < p.length; k++) {
        if (Math.hypot(p[k].re, p[k].im) > 1e-14) {
          if (idx !== -1) throw new Error("not a pure monomial");
          idx = k;
        }
      }
      if (idx === -1) return { kind: 'poly', degree: 0, coeff: { re: 0, im: 0 } };
      return { kind: 'poly', degree: idx, coeff: p[idx] };
    }

    // Pole atom: den has degree ≥ 1; need num constant (degree 0).
    if (num.length !== 1) throw new Error("numerator over (w-a)^k must be constant for strict PFD");
    const k = den.length - 1;
    const lc = den[k];
    // Monic-normalize.
    const lcInv = C.inv(lc);
    const dMon = den.map(c => C.mul(c, lcInv));
    // Candidate root: a = -dMon[k-1] / k.
    const a = C.scale(C.neg(dMon[k - 1]), 1 / k);
    // Verify dMon == (w - a)^k.
    const linear = [C.neg(a), { re: 1, im: 0 }];
    const expanded = polyPow(linear, k);
    if (expanded.length !== dMon.length) throw new Error("denominator is not (w-a)^k");
    for (let i = 0; i < expanded.length; i++) {
      const dr = expanded[i].re - dMon[i].re;
      const di = expanded[i].im - dMon[i].im;
      const scale = 1 + Math.hypot(dMon[i].re, dMon[i].im);
      if (Math.hypot(dr, di) > 1e-10 * scale) throw new Error("denominator is not (w-a)^k");
    }
    const residue = C.mul(num[0], lcInv);
    return { kind: 'pole', a, order: k, residue };
  }

  // ===========================================================================
  // Phase 2: general-rational decomposition via Durand–Kerner.
  // ===========================================================================
  // We need polynomialRoots + groupRootsByMultiplicity. Pull from QD.Direct
  // when available (browser/node) — both expose them.
  function getDirect() {
    if (QD && QD.Direct) return QD.Direct;
    if (typeof window !== 'undefined' && window.QD && window.QD.Direct) return window.QD.Direct;
    if (global && global.QD && global.QD.Direct) return global.QD.Direct;
    throw new Error("parse-h: QD.Direct (polynomialRoots) not available — load direct/direct-common.js first");
  }

  function phase2Decompose(P, Q, opts) {
    opts = opts || {};
    const clusterTol = opts.clusterTol || 1e-6;
    const warnings = [];

    // Long-divide for polynomial part.
    let polyPart = [];
    let remainder = P;
    if (Q.length === 1) {
      // den constant; whole thing is polynomial
      const dInv = C.inv(Q[0]);
      const all = P.map(c => C.mul(c, dInv));
      // Trim leading zeros at HIGH end already done; the polyPart length is full.
      polyPart = all.length === 1 && Math.hypot(all[0].re, all[0].im) < 1e-14 ? [] : all;
      return { poles: [], polyCoeffs: polyPart, warnings };
    }
    if (P.length >= Q.length) {
      const dm = polyDivMod(P, Q);
      polyPart = dm.quotient;
      remainder = dm.remainder;
      // If quotient is just [0], drop it.
      if (polyPart.length === 1 && Math.hypot(polyPart[0].re, polyPart[0].im) < 1e-14) polyPart = [];
    }

    // If remainder is zero, no poles.
    if (remainder.length === 0 || (remainder.length === 1 && Math.hypot(remainder[0].re, remainder[0].im) < 1e-14)) {
      return { poles: [], polyCoeffs: polyPart, warnings };
    }

    // Factor Q via Durand–Kerner.
    const Direct = getDirect();
    const rawRoots = Direct.polynomialRoots(Q);
    const groups = Direct.groupRootsByMultiplicity(rawRoots, clusterTol);

    // Cluster-spread warning: per group, recompute the spread of raw roots
    // assigned to it and warn if any spread is uncomfortably close to clusterTol.
    // (We don't get the assignment back from groupRootsByMultiplicity, so re-cluster
    // here for diagnostic purposes.)
    for (const g of groups) {
      if (g.multiplicity < 2) continue;
      let maxSpread = 0;
      for (const r of rawRoots) {
        const d = Math.hypot(r.re - g.root.re, r.im - g.root.im);
        if (d < clusterTol) maxSpread = Math.max(maxSpread, d);
      }
      if (maxSpread > 0.1 * clusterTol) {
        warnings.push("pole at " + C.format(g.root) + " has cluster spread " +
                      maxSpread.toExponential(2) + " (close to tol " +
                      clusterTol.toExponential(0) + ") — multiplicity may be off by one");
      }
    }

    // For each pole (a, m): shift Q by a, divide out u^m, series-divide P_shifted / denQuot
    // up to degree m-1 to read off C_s = T_{m-s}.
    const poles = [];
    for (const g of groups) {
      const a = g.root, m = g.multiplicity;
      const Pshift = polyShift(P, a);
      const Qshift = polyShift(Q, a);
      // Divide Qshift by u^m → coefficients [Q_m, Q_{m+1}, ...].
      // (Lower coefficients should be ≈ 0 by construction of m.)
      if (Qshift.length <= m) {
        warnings.push("internal: Qshift too short for pole at " + C.format(a));
        continue;
      }
      const denQuot = Qshift.slice(m);
      const T = seriesDivide(Pshift, denQuot, m - 1);
      const residues = new Array(m);
      for (let s = 1; s <= m; s++) residues[s - 1] = T[m - s];
      poles.push({ a, order: m, residues });
    }

    return { poles, polyCoeffs: polyPart, warnings };
  }

  // ===========================================================================
  // Public API: parseH and formatH.
  // ===========================================================================
  function parseH(expr, math, opts) {
    opts = opts || {};
    const mode = opts.mode || 'unbounded';      // 'bounded' / 'unbounded' / lqd-*
    // Polynomial part of h is allowed exactly in the three unbounded family
    // panels (classical unbounded + unbounded LQD ± singular). Bounded modes
    // — including bounded LQDs — must have polyCoeffs = [].
    const allowPoly = (mode === 'unbounded' ||
                       mode === 'lqd-unbounded' ||
                       mode === 'lqd-unbounded-singular');

    if (!math || !math.parse) throw new Error("parseH: math.js required");
    if (typeof expr !== 'string' || !expr.trim()) throw new Error("empty expression");

    let root;
    try { root = math.parse(expr); }
    catch (e) { throw new Error("parse error: " + (e.message || e)); }

    // ---------- Phase 1: strict per-summand walk ----------
    const summands = flattenSum(root, 1);
    let strictOK = true;
    const polyAccum = new Map();           // degree → Complex
    const poleAccum = new Map();           // "a.re|a.im|order" → { a, order, residue }

    try {
      for (const s of summands) {
        const cls = classifySummand(s.node, math);
        if (cls.kind === 'poly') {
          const k = cls.degree;
          const c = polyAccum.get(k) || { re: 0, im: 0 };
          const add = s.sign > 0 ? cls.coeff : C.neg(cls.coeff);
          polyAccum.set(k, C.add(c, add));
        } else {
          const key = poleKey(cls.a, cls.order);
          const cur = poleAccum.get(key);
          const add = s.sign > 0 ? cls.residue : C.neg(cls.residue);
          if (cur) cur.residue = C.add(cur.residue, add);
          else      poleAccum.set(key, { a: cls.a, order: cls.order, residue: add });
        }
      }
    } catch (e) {
      strictOK = false;
    }

    let result;
    if (strictOK) {
      // Assemble polyCoeffs[] (dense; length = max-degree + 1).
      let maxDeg = -1;
      for (const k of polyAccum.keys()) if (k > maxDeg) maxDeg = k;
      const polyCoeffs = [];
      for (let k = 0; k <= maxDeg; k++) {
        polyCoeffs.push(polyAccum.get(k) || { re: 0, im: 0 });
      }
      // Drop ALL-zero polynomial part.
      let allZero = true;
      for (const c of polyCoeffs) if (Math.hypot(c.re, c.im) > 1e-14) { allZero = false; break; }
      const finalPoly = allZero ? [] : polyCoeffs;

      // Group pole atoms by location a → assemble principal-part residue arrays.
      // Multiple orders at the same a contribute different s indices.
      const byLoc = new Map();
      for (const v of poleAccum.values()) {
        const lk = locKey(v.a);
        let entry = byLoc.get(lk);
        if (!entry) { entry = { a: v.a, residuesByOrder: new Map() }; byLoc.set(lk, entry); }
        const cur = entry.residuesByOrder.get(v.order) || { re: 0, im: 0 };
        entry.residuesByOrder.set(v.order, C.add(cur, v.residue));
      }
      const poles = [];
      for (const entry of byLoc.values()) {
        let order = 0;
        for (const k of entry.residuesByOrder.keys()) if (k > order) order = k;
        const residues = [];
        let anyNonzero = false;
        for (let s = 1; s <= order; s++) {
          const c = entry.residuesByOrder.get(s) || { re: 0, im: 0 };
          residues.push(c);
          if (Math.hypot(c.re, c.im) > 1e-14) anyNonzero = true;
        }
        if (anyNonzero) poles.push({ a: entry.a, order, residues });
      }
      result = { poles, polyCoeffs: finalPoly, warnings: [] };
    } else {
      // ---------- Phase 2: general-rational fallback ----------
      const rat = accRat(root, math);
      const num = polyTrim(rat.num);
      const den = polyTrim(rat.den);
      if (num.length === 0 || (num.length === 1 && Math.hypot(num[0].re, num[0].im) < 1e-14)) {
        return { poles: [], polyCoeffs: [], warnings: [] };
      }
      result = phase2Decompose(num, den, opts);
    }

    // Mode enforcement: bounded modes reject any polynomial part.
    if (!allowPoly && result.polyCoeffs.length > 0) {
      // Tolerate exactly-zero polyCoeffs even if length > 0 (defensive).
      let anyNonzero = false;
      for (const c of result.polyCoeffs) if (Math.hypot(c.re, c.im) > 1e-14) { anyNonzero = true; break; }
      if (anyNonzero) {
        throw new Error("polynomial part of h is only valid in unbounded mode " +
                        "(switch the mode or remove the w^k terms)");
      }
      result.polyCoeffs = [];
    }

    return result;
  }

  function poleKey(a, order) {
    return roundForKey(a.re) + '|' + roundForKey(a.im) + '|' + order;
  }
  function locKey(a) {
    return roundForKey(a.re) + '|' + roundForKey(a.im);
  }
  function roundForKey(x) {
    if (x === 0) return '0';
    // Identify two `a` values as the same if they agree to 12 significant digits.
    return x.toPrecision(12);
  }

  // ===========================================================================
  // formatH: structured h → math.js source string.
  // ===========================================================================
  function formatH(h) {
    const terms = [];
    const poles = h.poles || [];
    const polyCoeffs = h.polyCoeffs || [];

    // Polynomial part: C_∞,l · w^l, ascending l.
    for (let l = 0; l < polyCoeffs.length; l++) {
      const c = polyCoeffs[l];
      if (!c || Math.hypot(c.re, c.im) < 1e-14) continue;
      terms.push(formatMonomial(c, l));
    }
    // Pole parts: each pole a, each residue order s.
    for (const p of poles) {
      const a = p.a;
      for (let s = 1; s <= p.residues.length; s++) {
        const C0 = p.residues[s - 1];
        if (!C0 || Math.hypot(C0.re, C0.im) < 1e-14) continue;
        terms.push(formatPoleTerm(C0, a, s));
      }
    }
    if (terms.length === 0) return '0';
    return terms.join(' + ').replace(/\+ -/g, '- ');
  }

  function needsParensForCoeff(s) {
    // Wrap complex literals with an interior sign in parens, but NOT plain
    // numbers (including those with a leading minus).
    return /[+\-]/.test(s.slice(1));
  }

  function formatMonomial(c, k) {
    if (k === 0) return C.format(c);
    const cs = C.format(c);
    const wPart = (k === 1) ? 'w' : 'w^' + k;
    if (cs === '1')  return wPart;
    if (cs === '-1') return '-' + wPart;
    return (needsParensForCoeff(cs) ? '(' + cs + ')' : cs) + '*' + wPart;
  }

  function formatPoleTerm(c, a, s) {
    // Denominator: w, (w-a), w^s, (w-a)^s, depending on a and s.
    let denStr;
    if (Math.hypot(a.re, a.im) < 1e-14) {
      denStr = (s === 1) ? 'w' : 'w^' + s;
    } else {
      // (w - a) — write as `(w - <a>)` with a's sign absorbed when negative.
      const aStr = C.format(a);
      const inner = needsParensForCoeff(aStr) ? 'w - (' + aStr + ')'
                                              : (aStr.startsWith('-') ? 'w + ' + aStr.slice(1)
                                                                      : 'w - ' + aStr);
      denStr = (s === 1) ? '(' + inner + ')' : '(' + inner + ')^' + s;
    }
    const cs = C.format(c);
    if (cs === '1')  return '1/'  + denStr;
    if (cs === '-1') return '-1/' + denStr;
    return (needsParensForCoeff(cs) ? '(' + cs + ')' : cs) + '/' + denStr;
  }

  // ===========================================================================
  // Export.
  // ===========================================================================
  // Attach onto QD without clobbering existing keys (Complex, Taylor, Direct, …).
  QD.parseH  = parseH;
  QD.formatH = formatH;
})(typeof globalThis !== 'undefined' ? globalThis : this);
