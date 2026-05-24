// =============================================================================
// direct-common.js -- Direct problem: given a Riemann map φ, compute the
// quadrature function h such that φ(𝔻) ∈ QD(h).
//
// Supports BOUNDED classical QD with:
//   • polynomial φ           (boundedQD)
//   • rational  φ = P(z)/Q(z) with Q non-vanishing on 𝔻̄  (boundedQDRational)
//   • numerical fallback for arbitrary analytic-in-𝔻̄ φ   (numericalBoundedQD)
//
// And UNBOUNDED classical QD with Laurent-at-∞ φ:
//   φ(z) = c·z + F_0 + F_1/z + … (handled by unboundedQD).
//
// All variants produce hData of the same shape used by the inverse-problem
// solver:  hData = { poles: [{a, principal: [...]}, ...], polyPart: [...] }.
//
// DERIVATION (bounded polynomial case; see Graven thesis §3.2 for the full
// argument and §4.3 / Ch.6 for the rational / AQD generalizations):
//
//   • On ∂Ω, the Schwarz function σ(w) = w̄. By Green's theorem the bounded
//     classical QD identity  ∫_Ω f dA = (1/(2i)) ∮ f·h dw  is satisfied with
//     h equal to the principal-part representative of σ.
//   • On |z| = 1: σ∘φ(z) = conj(φ(z)) = φ#(z) = Σ_{l≥1} conj(c_l)·z^{−l}.
//   • Analytic continuation: σ(w) = φ#(φ⁻¹(w)) extends into Ω with a single
//     pole at w₀ = φ(0) (the φ-preimage of z = 0, where φ# is singular).
//   • Local Laurent at w = w₀ + ζ:
//
//       σ(w₀ + ζ) = Σ_{l≥1} conj(c_l) · ψ̃(ζ)^{−l},
//
//     where ψ̃(ζ) is the formal Taylor inverse of phiTilde = [0, c_1, …, c_n].
//   • Factor ψ̃(ζ) = ψ̃[1]·ζ·u(ζ) with u(0) = 1, giving
//
//       C_k = Σ_{l ≥ k} conj(c_l) · c_1^l · [ζ^{l−k}] u(ζ)^{−l},   k = 1..n.
//
//   The rational variant uses the same forward formula at every pole of
//   R̃(z) inside 𝔻; the per-pole computation is factored into
//   forwardLocalPrincipal.
//
// All Taylor primitives used (invert, reciprocal, mul, pow, compose) live
// in taylor.js. The polynomial root-finder (Durand–Kerner) and the
// reverse-conjugate / Taylor-around-z₀ polynomial helpers live in this file.
// =============================================================================
'use strict';

(function (global) {

  // Pick the QD namespace consistently with the other solver files. In the
  // browser everything lives on window.QD. In node-test.js, the solver files
  // attach to `module.exports` and we follow suit so we can see QD.Faber etc.
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : (global.QD || (global.QD = {})));
  const Direct = QD.Direct || (QD.Direct = {});
  Direct.version = '0.1.0-mvp';

  const C = (typeof Complex !== 'undefined') ? Complex
          : (typeof global.Complex !== 'undefined' ? global.Complex : null);
  const T = (typeof Taylor !== 'undefined') ? Taylor
          : (typeof global.Taylor !== 'undefined' ? global.Taylor : null);
  if (!C || !T) {
    throw new Error("direct-common.js: complex.js and taylor.js must be loaded first");
  }

  // ===========================================================================
  // boundedQD: polynomial φ → hData (single pole of order n at w₀).
  // ---------------------------------------------------------------------------
  //   coeffs:  [c_0, c_1, c_2, ..., c_n]   Complex array (length n+1, n ≥ 1)
  // Returns:
  //   {
  //     hData:    { poles: [ { a: w_0, principal: [C_1, ..., C_n] } ] },
  //     w0:       c_0,
  //     degree:   n,
  //     warnings: [ ... ]   (e.g., univalence checks)
  //   }
  // Throws if degree < 1 or c_1 = 0.
  // ===========================================================================
  function boundedQD(coeffs) {
    if (!coeffs || coeffs.length < 2) {
      throw new Error("Direct.boundedQD: need at least 2 coefficients (degree ≥ 1)");
    }
    const n = coeffs.length - 1;
    const c1 = coeffs[1];
    if (C.abs2(c1) < 1e-30) {
      throw new Error("Direct.boundedQD: c_1 ≈ 0; φ is not locally univalent at z = 0");
    }
    const w0 = C.clone(coeffs[0]);

    // phiTilde = [0, c_1, ..., c_n], for Taylor.invert.
    const phiTilde = T.zero(n + 1);
    for (let i = 1; i <= n; i++) phiTilde[i] = C.clone(coeffs[i]);

    // ψ̃(ζ) = Taylor.invert(phiTilde, n)
    const psi = T.invert(phiTilde, n);

    // u(ζ) = ψ̃(ζ) / (ψ̃[1] · ζ)
    // Coefficients: u[i] = ψ̃[i+1] / ψ̃[1]   for i = 0..n-1
    // (u[0] = 1 by construction.)
    const psi1Inv = C.inv(psi[1]);
    const u = T.zero(n);                                  // length n: indices 0..n-1
    for (let i = 0; i < n; i++) {
      u[i] = C.mul(psi[i + 1], psi1Inv);
    }

    // For each l = 1..n we need u(ζ)^{−l} up to order ζ^{l−1}.
    // Computing once for the highest l = n is enough; we extract the lower-l
    // versions by repeated multiplication of u^{−1}.
    const uInv = T.reciprocal(u, n - 1);                  // 1/u up to ζ^{n-1}

    // Build u^{-l} for l = 1..n as a list of Taylor series, each truncated to
    // length max(n-1, l-1). Length n-1 suffices because we only need
    // [ζ^{l-k}] u^{-l} for k = 1..l, i.e., max index l-1 ≤ n-1.
    const uPowNeg = [null];                               // uPowNeg[l] = u^{-l}
    uPowNeg[1] = T.truncate(uInv, n - 1);
    for (let l = 2; l <= n; l++) {
      uPowNeg[l] = T.mul(uPowNeg[l - 1], uInv, n - 1);
    }

    // c_1^l prefactors
    const c1Pow = [{ re: 1, im: 0 }];                     // c1Pow[l] = c_1^l
    for (let l = 1; l <= n; l++) c1Pow.push(C.mul(c1Pow[l - 1], c1));

    // Assemble: C_k = Σ_{l ≥ k} conj(c_l) · c_1^l · [ζ^{l−k}] u^{−l}
    const principal = new Array(n);
    for (let k = 1; k <= n; k++) {
      let acc = { re: 0, im: 0 };
      for (let l = k; l <= n; l++) {
        const idx = l - k;
        if (idx >= uPowNeg[l].length) continue;
        const term = C.mul(C.conj(coeffs[l]), c1Pow[l]);
        acc = C.add(acc, C.mul(term, uPowNeg[l][idx]));
      }
      principal[k - 1] = acc;
    }

    // Trim trailing-near-zero principal entries — preserve the leading C_n
    // (we promised a pole of order n). The trim is a courtesy for display.
    let mEff = n;
    while (mEff > 1 && C.abs(principal[mEff - 1]) < 1e-14 * C.abs(principal[0])) {
      mEff--;
    }
    const trimmedPrincipal = principal.slice(0, mEff);

    const warnings = [];
    // Sanity: does φ have a critical point inside 𝔻? Roughly check φ'(0) = c_1 and
    // a quick |φ'| sweep at a few z's. Full univalence check is deferred.
    // (Boundary-univalence check happens at visualization time.)
    if (n >= 2) {
      // φ'(z) = Σ_{l ≥ 1} l·c_l z^{l-1}. Sample at |z| = 0.99 in a few directions.
      let minAbsDeriv = Infinity;
      for (let k = 0; k < 8; k++) {
        const theta = 2 * Math.PI * k / 8;
        const z = { re: 0.99 * Math.cos(theta), im: 0.99 * Math.sin(theta) };
        // Horner-ish evaluation of φ'(z)
        let v = C.scale(coeffs[n], n);
        for (let l = n - 1; l >= 1; l--) {
          v = C.add(C.mul(v, z), C.scale(coeffs[l], l));
        }
        const av = C.abs(v);
        if (av < minAbsDeriv) minAbsDeriv = av;
      }
      if (minAbsDeriv < 1e-3) {
        warnings.push("φ'(z) approaches 0 inside 𝔻 (min ≈ " + minAbsDeriv.toExponential(2) + "); univalence likely fails");
      }
    }

    return {
      hData: {
        poles: [{ a: w0, principal: trimmedPrincipal }],
      },
      w0,
      degree: n,
      warnings,
    };
  }

  // ===========================================================================
  // Boundary samples of φ on z = e^{iθ}, for live ∂Ω preview.
  // ===========================================================================
  function sampleBoundaryPolynomial(coeffs, N) {
    const pts = new Array(N);
    for (let n = 0; n < N; n++) {
      const theta = 2 * Math.PI * n / N;
      const z = { re: Math.cos(theta), im: Math.sin(theta) };
      // Horner: φ(z) = c_0 + z·(c_1 + z·(c_2 + ...))
      let v = C.clone(coeffs[coeffs.length - 1]);
      for (let l = coeffs.length - 2; l >= 0; l--) {
        v = C.add(C.mul(v, z), coeffs[l]);
      }
      pts[n] = v;
    }
    return pts;
  }

  // ===========================================================================
  // parsePolynomialInZ: thin wrapper around parseRationalInZ that enforces
  // the polynomial-and-locally-univalent contract expected by boundedQD.
  // ---------------------------------------------------------------------------
  //   parsePolynomialInZ(astOrString, mathLib) → Complex[]
  //
  // Accepts the same expression grammar as parseRationalInZ, but rejects
  // expressions whose result has a non-trivial denominator (division by a
  // subexpression containing z), constant expressions (no z dependence),
  // degrees above the polynomial cap (12), or c₁ = 0.
  //
  // For full rational expressions, callers should use parseRationalInZ
  // directly and dispatch on Array.isArray.
  // ===========================================================================
  function parsePolynomialInZ(astOrString, math) {
    if (!math || !math.parse) throw new Error("parsePolynomialInZ: math.js required");
    const r = parseRationalInZ(astOrString, math, { maxDegree: 12 });
    if (!Array.isArray(r)) {
      throw new Error("division by non-constant subexpression (use parseRationalInZ for rational φ)");
    }
    if (r.length < 2) {
      // Trimmed to a constant. Distinguish the all-zero case (c₁ = 0 after
      // expansion, e.g. "0*z") from the genuinely-constant case ("i", "3").
      if (r.length === 1 && Math.hypot(r[0].re, r[0].im) < 1e-14) {
        throw new Error("c₁ = 0; φ not locally univalent at 0 (expression reduces to zero)");
      }
      throw new Error("expression has no z-dependence");
    }
    if (Math.hypot(r[1].re, r[1].im) < 1e-14) {
      throw new Error("c₁ = 0; φ not locally univalent at 0");
    }
    return r;
  }

  function mjxToComplexImpl(v) {
    if (typeof v === 'number') return { re: v, im: 0 };
    if (v && typeof v.re === 'number' && typeof v.im === 'number') {
      return { re: v.re, im: v.im };
    }
    if (v && typeof v.toJSON === 'function') {
      const j = v.toJSON();
      if (j && j.mathjs === 'Complex') return { re: j.re, im: j.im };
    }
    if (v && v.constructor && v.constructor.name === 'Complex') {
      return { re: v.re, im: v.im };
    }
    if (v && v.constructor && v.constructor.name === 'BigNumber') {
      return { re: Number(v), im: 0 };
    }
    throw new Error('cannot convert math.js value to Complex: ' + String(v));
  }

  // ===========================================================================
  // parseRationalInZ: extend the parser to handle ARBITRARY rational
  // expressions in z. Returns either a polynomial array (if the parsed
  // expression has trivial denominator [1]) or a RationalForm {num, den}.
  // ---------------------------------------------------------------------------
  // Grammar walker operating on RationalForms at every AST node, so
  // expressions like  z/(z-2) + 1/(z-3)  reduce cleanly to a single P/Q
  // via polynomial cross-multiplication.
  //
  //   parseRationalInZ(astOrString, math) → Complex[]  |  {num: Complex[], den: Complex[]}
  //
  // Throws on:
  //   • Symbol other than z or i.
  //   • Negative integer exponent (would be a rational; the parser handles
  //     positive integer exponents only — for rational use literal division).
  //   • Function call with non-constant arguments.
  //   • Degree (numerator OR denominator) exceeds maxDegree (default 32).
  //
  // Note: c_1-validity (≠0) is the CALLER's responsibility; this parser does
  // not enforce it (parsePolynomialInZ does for the polynomial path).
  // ===========================================================================
  function parseRationalInZ(astOrString, math, options) {
    options = options || {};
    const maxDegree = options.maxDegree || 32;

    if (!math || !math.parse) throw new Error("parseRationalInZ: math.js required");
    let node;
    if (typeof astOrString === 'string') {
      const expr = astOrString.trim();
      if (!expr) throw new Error("empty expression");
      try { node = math.parse(expr); }
      catch (e) { throw new Error('parse error: ' + (e.message || e)); }
    } else {
      node = astOrString;
    }

    const rat = accumulateRationalImpl(node, math);
    trimPolyInPlace(rat.num);
    trimPolyInPlace(rat.den);
    if (rat.num.length - 1 > maxDegree) {
      throw new Error("numerator degree " + (rat.num.length - 1) + " exceeds cap (" + maxDegree + ")");
    }
    if (rat.den.length - 1 > maxDegree) {
      throw new Error("denominator degree " + (rat.den.length - 1) + " exceeds cap (" + maxDegree + ")");
    }
    if (rat.num.length === 0) {
      throw new Error("expression evaluates to 0");
    }

    // If denominator is the constant 1 (after dividing through), simplify.
    if (rat.den.length === 1 && Math.hypot(rat.den[0].re - 1, rat.den[0].im) < 1e-13) {
      return rat.num;     // pure polynomial result
    }
    // If denominator is a nonzero constant, push it into numerator and return polynomial.
    if (rat.den.length === 1) {
      const dInv = C.inv(rat.den[0]);
      const out = rat.num.map(c => C.mul(c, dInv));
      return out;
    }
    // Otherwise: genuine rational. Normalize so denominator's leading coeff is 1.
    const denLead = rat.den[rat.den.length - 1];
    const dInv = C.inv(denLead);
    rat.num = rat.num.map(c => C.mul(c, dInv));
    rat.den = rat.den.map(c => C.mul(c, dInv));
    return { num: rat.num, den: rat.den };
  }

  // Recursive walker → RationalForm {num, den}.
  function accumulateRationalImpl(node, math) {
    if (!node) throw new Error('null AST node');

    if (node.isConstantNode) {
      return { num: [mjxToComplexImpl(node.value)], den: [{ re: 1, im: 0 }] };
    }
    if (node.isSymbolNode) {
      if (node.name === 'z') {
        return { num: [{re:0,im:0}, {re:1,im:0}], den: [{re:1,im:0}] };
      }
      if (node.name === 'i') {
        return { num: [{re:0,im:1}], den: [{re:1,im:0}] };
      }
      throw new Error("unknown symbol '" + node.name + "' (only z and i are allowed)");
    }
    if (node.isParenthesisNode) {
      return accumulateRationalImpl(node.content, math);
    }
    if (node.isOperatorNode) {
      const op = node.op;
      if (op === '+') {
        let acc = accumulateRationalImpl(node.args[0], math);
        for (let i = 1; i < node.args.length; i++) {
          const rhs = accumulateRationalImpl(node.args[i], math);
          acc = addRat(acc, rhs);
        }
        return acc;
      }
      if (op === '-') {
        if (node.args.length === 1) {
          const r = accumulateRationalImpl(node.args[0], math);
          return { num: r.num.map(c => C.neg(c)), den: r.den };
        }
        let acc = accumulateRationalImpl(node.args[0], math);
        for (let i = 1; i < node.args.length; i++) {
          const rhs = accumulateRationalImpl(node.args[i], math);
          acc = addRat(acc, { num: rhs.num.map(c => C.neg(c)), den: rhs.den });
        }
        return acc;
      }
      if (op === '*') {
        let acc = accumulateRationalImpl(node.args[0], math);
        for (let i = 1; i < node.args.length; i++) {
          const rhs = accumulateRationalImpl(node.args[i], math);
          acc = mulRat(acc, rhs);
        }
        return acc;
      }
      if (op === '/') {
        const lhs = accumulateRationalImpl(node.args[0], math);
        const rhs = accumulateRationalImpl(node.args[1], math);
        if (rhs.num.length === 0 || (rhs.num.length === 1 && C.abs(rhs.num[0]) < 1e-300)) {
          throw new Error("division by zero");
        }
        return mulRat(lhs, { num: rhs.den, den: rhs.num });    // a/b · b/a
      }
      if (op === '^') {
        const expNode = node.args[1];
        if (!expNode.isConstantNode) throw new Error("exponent in '^' must be a literal integer");
        const k = (typeof expNode.value === 'number') ? expNode.value : Number(expNode.value);
        if (!Number.isInteger(k)) {
          throw new Error("exponent must be an integer (got " + expNode.value + ")");
        }
        const base = accumulateRationalImpl(node.args[0], math);
        if (k === 0) return { num: [{re:1,im:0}], den: [{re:1,im:0}] };
        if (k > 0) return powRat(base, k);
        // k < 0: 1/base^|k|
        const positive = powRat(base, -k);
        if (positive.num.length === 0 || (positive.num.length === 1 && C.abs(positive.num[0]) < 1e-300)) {
          throw new Error("division by zero (negative exponent of zero subexpression)");
        }
        return { num: positive.den, den: positive.num };
      }
      throw new Error("unsupported operator: " + op);
    }
    if (node.isFunctionNode) {
      // Only allow function calls with fully constant args (no z).
      const argRats = node.args.map(a => accumulateRationalImpl(a, math));
      const allConst = argRats.every(r =>
        r.num.length <= 1 && r.den.length <= 1);
      if (!allConst) {
        throw new Error("function " + node.name + " requires constant arguments");
      }
      let val;
      try { val = node.evaluate(); }
      catch (e) { throw new Error("could not evaluate " + node.name + "(...): " + (e.message || e)); }
      return { num: [mjxToComplexImpl(val)], den: [{re:1,im:0}] };
    }
    throw new Error('unsupported AST node: ' + (node.type || '(unknown)'));
  }

  // Polynomial operations (ascending-power Complex[]):
  function addPolys(a, b) {
    const n = Math.max(a.length, b.length);
    const out = new Array(n);
    for (let i = 0; i < n; i++) {
      const av = (i < a.length) ? a[i] : { re: 0, im: 0 };
      const bv = (i < b.length) ? b[i] : { re: 0, im: 0 };
      out[i] = C.add(av, bv);
    }
    return out;
  }
  function mulPolys(a, b) {
    const n = a.length, m = b.length;
    const out = new Array(n + m - 1);
    for (let i = 0; i < out.length; i++) out[i] = { re: 0, im: 0 };
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < m; j++) {
        out[i + j] = C.add(out[i + j], C.mul(a[i], b[j]));
      }
    }
    return out;
  }
  function trimPolyInPlace(a) {
    while (a.length > 1 && C.abs(a[a.length - 1]) < 1e-15) a.pop();
  }

  // Rational operations: keep them in {num, den} form. No common-factor
  // reduction yet (the kernel handles roots regardless, so reduction is
  // a polish issue).
  function addRat(a, b) {
    return {
      num: addPolys(mulPolys(a.num, b.den), mulPolys(b.num, a.den)),
      den: mulPolys(a.den, b.den),
    };
  }
  function mulRat(a, b) {
    return {
      num: mulPolys(a.num, b.num),
      den: mulPolys(a.den, b.den),
    };
  }
  function powRat(r, k) {
    if (k === 0) return { num: [{re:1,im:0}], den: [{re:1,im:0}] };
    let accNum = r.num, accDen = r.den;
    for (let i = 1; i < k; i++) {
      accNum = mulPolys(accNum, r.num);
      accDen = mulPolys(accDen, r.den);
    }
    return { num: accNum, den: accDen };
  }

  // ===========================================================================
  // Inverse direction: format a Complex[] back to canonical math.js source.
  // ===========================================================================
  function polynomialToString(coeffs) {
    if (!coeffs || coeffs.length === 0) return '0';
    const terms = [];
    for (let k = 0; k < coeffs.length; k++) {
      const c = coeffs[k];
      if (!c || (Math.abs(c.re) < 1e-15 && Math.abs(c.im) < 1e-15)) continue;
      terms.push(formatTerm(c, k));
    }
    if (terms.length === 0) return '0';
    return terms.join(' + ').replace(/\+ -/g, '- ');
  }

  function formatTerm(c, k) {
    const cs = formatComplex(c);
    if (k === 0) return cs;
    const zPart = (k === 1) ? 'z' : 'z^' + k;
    if (cs === '1') return zPart;
    if (cs === '-1') return '-' + zPart;
    // Wrap complex literals in parens when followed by z
    const needsParens = /[+\-]/.test(cs.slice(1));    // sign in interior → complex
    return (needsParens ? '(' + cs + ')' : cs) + '*' + zPart;
  }

  // Thin wrapper around Complex.format for use by polynomialToString /
  // formatTerm. Preserves the integer-snap + ±i short-form behavior.
  function formatComplex(c) {
    return C.format(c);
  }

  // ===========================================================================
  // unboundedQD: Laurent-at-infinity φ → hData.
  // ---------------------------------------------------------------------------
  // Input:
  //   c    real positive  (conformal radius;  φ'(∞) = c)
  //   F    Complex[]      [F_0, F_1, ..., F_{m-1}]
  //                       so that φ(z) = c·z + F_0 + F_1/z + ... + F_{m-1}/z^{m-1}
  //                       (m = F.length;  m = 0 ⇒ φ = c·z, exterior of disk)
  //
  // Output:
  //   {
  //     hData:   { poles: [...], polyPart: [C_∞,0, ..., C_∞,m-1] },
  //     finitePoleHandled: boolean,
  //     warnings: [...]
  //   }
  //
  // POLYNOMIAL PART (always computed, for any m):
  //   Back-substitute the dual of inverseFaberAtInfinity. For each l from
  //   m−1 down to 0:
  //     conj(F_l) − Σ_{l' > l} conj(C_∞,l') · [u^{l'−l}] g(u)^{l'}    [in conj]
  //   then  C_∞,l = (conj(...)) / c^l.
  //   (Diagonal entry [u^0] g^l = c^l makes this triangular.)
  //
  // FINITE POLES (handled only in simple cases):
  //   • m = 0  →  single pole at w = 0, residue c²        (exterior of disk |w|=c)
  //   • m = 1  →  single pole at w = F_0, residue c²       (exterior of disk centered at F_0)
  //   • m ≥ 2 with F_1 = F_2 = ... = 0  →  same as m = 1
  //   • Otherwise: σ(w) has branch-cut structure in K and Ω is generically
  //     NOT a classical QD; we leave finitePoles = [] and emit a warning.
  //     Such cases would require an unbounded-rational φ ansatz (not yet
  //     implemented; the bounded rational kernel boundedQDRational handles
  //     the bounded analog).
  // ===========================================================================
  function unboundedQD(c, F) {
    if (typeof c !== 'number' || c <= 0 || !isFinite(c)) {
      throw new Error("Direct.unboundedQD: c must be a positive real number");
    }
    F = F || [];
    const m = F.length;

    // ---- Polynomial part: back-substitute the triangular system. ----
    const polyPart = new Array(m);
    if (m > 0) {
      // g(u) = c + F_0·u + F_1·u² + ... + F_{m-1}·u^m   (length m+1)
      const g = T.zero(m + 1);
      g[0] = { re: c, im: 0 };
      for (let i = 1; i <= m; i++) g[i] = C.clone(F[i - 1]);

      // Precompute g^l for l = 0..m-1, each truncated to length m.
      const gPow = [T.zero(m)];
      gPow[0][0] = { re: 1, im: 0 };
      for (let l = 1; l < m; l++) {
        gPow[l] = T.mul(gPow[l - 1], g, m - 1);
      }

      for (let l = m - 1; l >= 0; l--) {
        let known = { re: 0, im: 0 };
        for (let lp = l + 1; lp < m; lp++) {
          const idx = lp - l;
          if (idx >= gPow[lp].length) continue;
          const M = gPow[lp][idx];
          // Term = conj(C_∞,lp) · conj(M)
          const term = C.mul(C.conj(polyPart[lp]), C.conj(M));
          known = C.add(known, term);
        }
        // conj(F_l) − known = conj(C_∞,l) · c^l, so C_∞,l = conj(diff) / c^l.
        const diff = C.sub(F[l], known);
        polyPart[l] = C.scale(C.conj(diff), 1 / Math.pow(c, l));
      }
    }

    // ---- Finite poles ----
    const warnings = [];
    const finitePoles = [];
    let finitePoleHandled = true;
    if (m === 0) {
      finitePoles.push({ a: { re: 0, im: 0 }, principal: [{ re: c * c, im: 0 }] });
    } else {
      // Simple case: F_l = 0 for all l ≥ 1.
      let allZero = true;
      for (let l = 1; l < m; l++) {
        if (C.abs(F[l]) > 1e-14) { allZero = false; break; }
      }
      if (allZero) {
        finitePoles.push({ a: C.clone(F[0]), principal: [{ re: c * c, im: 0 }] });
      } else {
        finitePoleHandled = false;
        warnings.push("F_l ≠ 0 for some l ≥ 1: h's finite poles are not computed " +
                      "(σ is generically non-rational; this Ω is unlikely to be a " +
                      "classical QD without an unbounded-rational φ ansatz, which is " +
                      "not yet implemented).");
      }
    }

    return {
      hData: { poles: finitePoles, polyPart },
      c,
      F,
      finitePoleHandled,
      warnings,
    };
  }

  // Boundary samples of an unbounded-Laurent φ on z = e^{iθ}.
  function sampleBoundaryLaurent(c, F, N) {
    const pts = new Array(N);
    const m = F.length;
    for (let n = 0; n < N; n++) {
      const theta = 2 * Math.PI * n / N;
      const z = { re: Math.cos(theta), im: Math.sin(theta) };
      // φ(z) = c·z + Σ_l F_l / z^l. On |z|=1, 1/z = conj(z).
      let w = C.scale(z, c);
      if (m > 0) {
        let zInvPow = { re: 1, im: 0 };                    // z^{-l}, start at l=0
        for (let l = 0; l < m; l++) {
          w = C.add(w, C.mul(F[l], zInvPow));
          // Next power: zInvPow ← zInvPow / z = zInvPow * conj(z) on |z|=1.
          // (For numerical sampling on |z|=1 only.)
          zInvPow = { re: zInvPow.re * z.re + zInvPow.im * z.im,
                      im: zInvPow.im * z.re - zInvPow.re * z.im };
        }
      }
      pts[n] = w;
    }
    return pts;
  }

  // ===========================================================================
  // numericalBoundedQD: free-form φ → hData via DFT + polynomial truncation.
  // ---------------------------------------------------------------------------
  // For φ(z) analytic in 𝔻̄ (not necessarily polynomial), the Fourier expansion
  // on |z|=1 reads φ(e^{iθ}) = Σ_{k≥0} c_k e^{ikθ}  (the c_k are the Taylor
  // coefficients of φ at z=0; the negative-frequency coefficients are zero by
  // analyticity).
  //
  // Algorithm:
  //   1. Sample φ at N points on |z|=1.
  //   2. Discrete Fourier transform (naive O(N·K) — N=256 is plenty for
  //      smooth φ; no FFT dependency).
  //   3. Extract c_k for k = 0..maxOrder. Check |c_{-k}| ≈ 0 as an
  //      analyticity diagnostic (large values ⇒ φ is NOT analytic in 𝔻̄, so
  //      the inferred polynomial is meaningless).
  //   4. Truncate: drop trailing coefficients below `tol`.
  //   5. Call the symbolic boundedQD with the truncated coefficient list.
  //
  // For a polynomial φ of degree ≤ maxOrder, the result is EXACT (DFT
  // recovers the exact c_k for a band-limited signal). For non-polynomial
  // analytic φ (e.g. exp(z), 1/(z+2), log(1+z)), the result is the QD
  // associated to the polynomial truncation, with a diagnostic indicating
  // how far we truncated.
  //
  // Caveat: this is for the BOUNDED classical case only. For unbounded /
  // LQD / AQD numerical fallback see future stages.
  // ===========================================================================
  function numericalBoundedQD(phiFn, options) {
    options = options || {};
    const N        = options.numSamples || 256;
    const maxOrder = options.maxOrder   || 12;
    const tol      = options.tol        || 1e-8;

    // 1. Sample φ on |z| = 1.
    const samples = new Array(N);
    for (let n = 0; n < N; n++) {
      const theta = 2 * Math.PI * n / N;
      const z = { re: Math.cos(theta), im: Math.sin(theta) };
      samples[n] = phiFn(z);
      if (!isFinite(samples[n].re) || !isFinite(samples[n].im)) {
        throw new Error("φ returned a non-finite value at z = e^{iθ}, θ ≈ " + theta.toFixed(3));
      }
    }

    // 2. DFT helper: φ_k = (1/N) Σ_n φ(z_n) e^{-ikθ_n}
    function dft(k) {
      let re = 0, im = 0;
      for (let n = 0; n < N; n++) {
        const theta = 2 * Math.PI * n / N;
        const cosT = Math.cos(-k * theta);
        const sinT = Math.sin(-k * theta);
        re += samples[n].re * cosT - samples[n].im * sinT;
        im += samples[n].re * sinT + samples[n].im * cosT;
      }
      return { re: re / N, im: im / N };
    }

    // 3. Taylor coefficients c_k for k = 0..maxOrder.
    const c = new Array(maxOrder + 1);
    for (let k = 0; k <= maxOrder; k++) c[k] = dft(k);

    // Analyticity diagnostic: |c_{-k}| should be ≈ 0 for k > 0.
    let analyticityScore = 0;
    for (let k = 1; k <= Math.min(10, Math.floor(N / 2)); k++) {
      const v = dft(-k);
      const mag = Math.hypot(v.re, v.im);
      if (mag > analyticityScore) analyticityScore = mag;
    }

    // 4. Truncate: drop trailing near-zero coefficients (keep at least degree 1).
    let truncOrder = maxOrder;
    const scale = Math.max(...c.map(x => Math.hypot(x.re, x.im)));
    const cutoff = Math.max(tol, scale * tol);
    while (truncOrder > 1 && Math.hypot(c[truncOrder].re, c[truncOrder].im) < cutoff) {
      truncOrder--;
    }
    const cTrunc = c.slice(0, truncOrder + 1);

    // Validate c_1 ≠ 0. For non-analytic φ this commonly fails (e.g.,
    // φ = conj(z) has all positive-frequency coefficients zero). Return a
    // soft diagnostic rather than throwing.
    if (Math.hypot(cTrunc[1].re, cTrunc[1].im) < 1e-12) {
      return {
        hData: { poles: [] },
        w0: cTrunc[0],
        taylorCoeffs: c,
        truncationOrder: 0,
        analyticityScore,
        polynomialSuffices: false,
        warnings: [
          "Inferred c_1 ≈ 0; φ is not locally univalent at z=0. " +
          (analyticityScore > tol
            ? "(Negative-frequency Fourier mass = " + analyticityScore.toExponential(2) +
              " — φ is not analytic in 𝔻̄, so DFT recovery is meaningless.)"
            : "(Smooth interior critical point at 0; h cannot be defined.)"),
        ],
      };
    }

    // 5. Call the symbolic boundedQD.
    const result = boundedQD(cTrunc);

    // Aggregate diagnostics + warnings.
    const warnings = result.warnings.slice();
    if (analyticityScore > tol) {
      warnings.push("Non-zero negative-frequency Fourier coefficients (max " +
        analyticityScore.toExponential(2) +
        "); φ does not appear to be analytic in 𝔻̄. The computed h is meaningless for non-analytic φ.");
    }
    if (truncOrder === maxOrder) {
      const tailMag = Math.hypot(c[maxOrder].re, c[maxOrder].im);
      if (tailMag > tol) {
        warnings.push("Polynomial truncation at degree " + maxOrder +
          " is approximate (|c_" + maxOrder + "| ≈ " + tailMag.toExponential(2) +
          " > tol). Increase maxOrder if higher precision is needed.");
      }
    }

    return {
      hData:             result.hData,
      w0:                cTrunc[0],
      taylorCoeffs:      c,
      truncationOrder:   truncOrder,
      analyticityScore,
      polynomialSuffices: analyticityScore < tol,
      warnings,
    };
  }

  // ===========================================================================
  // evalH: evaluate hData at a complex point w.
  //   h(w) = Σ_l polyPart[l] · w^l   +   Σ_j Σ_s C_{j,s} / (w − a_j)^s
  // Returns Complex {re, im}.
  // ===========================================================================
  function evalH(hData, w) {
    let vre = 0, vim = 0;
    const poly = hData.polyPart || [];
    if (poly.length > 0) {
      // Horner: acc = poly[L]; for l = L-1..0: acc = acc·w + poly[l]
      let ar = poly[poly.length - 1].re, ai = poly[poly.length - 1].im;
      for (let l = poly.length - 2; l >= 0; l--) {
        const nr = ar * w.re - ai * w.im + poly[l].re;
        const ni = ar * w.im + ai * w.re + poly[l].im;
        ar = nr; ai = ni;
      }
      vre += ar; vim += ai;
    }
    for (const pole of (hData.poles || [])) {
      const dr = w.re - pole.a.re, di = w.im - pole.a.im;
      const d2 = dr * dr + di * di;
      if (d2 < 1e-30) continue;
      let invR = dr / d2, invI = -di / d2;        // 1/(w − a)
      const stepR = invR, stepI = invI;
      for (let s = 0; s < pole.principal.length; s++) {
        const Cs = pole.principal[s];
        vre += Cs.re * invR - Cs.im * invI;
        vim += Cs.re * invI + Cs.im * invR;
        if (s + 1 < pole.principal.length) {
          const nr = invR * stepR - invI * stepI;
          const ni = invR * stepI + invI * stepR;
          invR = nr; invI = ni;
        }
      }
    }
    return { re: vre, im: vim };
  }

  // ===========================================================================
  // verifyBoundaryIdentity: check that  h(φ(z)) − conj(φ(z))  is analytic in
  // 𝔻 (composed with φ), via Fourier negative-frequency mass on |z|=1.
  // ---------------------------------------------------------------------------
  // For any classical QD, the Schwarz function σ(w) = w̄ on ∂Ω extends
  // meromorphically into Ω with poles at the quadrature nodes. h is the
  // SUM OF PRINCIPAL PARTS of σ at those nodes (modulo any analytic-in-Ω
  // part that's absorbed into the polyPart for unbounded shapes, or dropped
  // for bounded shapes). So h ≠ conj(w) pointwise in general; rather,
  //
  //     R(w) := σ(w) − h(w)        is analytic in Ω,
  //
  // and on ∂Ω we have R = conj(w) − h(w). Pulling back via φ:  R∘φ is
  // analytic in 𝔻. As a Fourier series on |z|=1, an analytic-in-𝔻 function
  // has ONLY non-negative-frequency terms.
  //
  // So the diagnostic is:
  //   Take Δ(θ) := h(φ(e^{iθ})) − conj(φ(e^{iθ})).
  //   Compute its discrete Fourier coefficients ĉ_k.
  //   Report   negMass = √(Σ_{k<0} |ĉ_k|²)   — should be ≈ 0 for any valid QD.
  //
  // For BOUNDED mode  with c_0 ≠ 0, Δ has a nonzero zero-mode (analytic
  //   constant) but negMass ≈ 0.
  // For UNBOUNDED mode where h includes the polyPart-constant, Δ ≡ 0 and
  //   both negMass and zeroMass are ≈ 0.
  // For non-QD shapes (e.g., unbounded φ with F_l ≠ 0 for l ≥ 1, or non-
  //   analytic numerical φ), negMass is significantly non-zero.
  //
  // Returns:
  //   { negMass, zeroMass, posMass, N, samples: Δ-values }
  // ===========================================================================
  function verifyBoundaryIdentity(hData, boundaryPts, options) {
    options = options || {};
    const N = boundaryPts.length;
    const K = Math.min(Math.floor(N / 2) - 1, options.maxFreq || 24);

    // Δ(θ_n) = h(φ(e^{iθ_n})) − conj(φ(e^{iθ_n})).
    const delta = new Array(N);
    let scale = 0;
    for (let n = 0; n < N; n++) {
      const w = boundaryPts[n];
      if (!isFinite(w.re) || !isFinite(w.im)) {
        delta[n] = { re: 0, im: 0 };
        continue;
      }
      const hv = evalH(hData, w);
      delta[n] = { re: hv.re - w.re, im: hv.im + w.im };       // h − conj(w)
      const s = Math.max(Math.hypot(hv.re, hv.im), Math.hypot(w.re, w.im));
      if (s > scale) scale = s;
    }

    // Naive DFT for k = -K..K.
    let negMass = 0, posMass = 0, zeroMass = 0;
    for (let k = -K; k <= K; k++) {
      let re = 0, im = 0;
      for (let n = 0; n < N; n++) {
        const theta = 2 * Math.PI * n / N;
        const cosT = Math.cos(-k * theta);
        const sinT = Math.sin(-k * theta);
        re += delta[n].re * cosT - delta[n].im * sinT;
        im += delta[n].re * sinT + delta[n].im * cosT;
      }
      const mag2 = (re * re + im * im) / (N * N);
      if      (k < 0) negMass  += mag2;
      else if (k > 0) posMass  += mag2;
      else            zeroMass = mag2;
    }
    return {
      negMass:  Math.sqrt(negMass),
      posMass:  Math.sqrt(posMass),
      zeroMass: Math.sqrt(zeroMass),
      scale,
      N,
      maxFreq: K,
    };
  }

  // ===========================================================================
  // Polynomial root finder (Durand–Kerner / Weierstrass simultaneous iteration).
  // ---------------------------------------------------------------------------
  // For a monic-normalized polynomial p(z) = z^n + a_{n-1} z^{n-1} + ... + a_0
  // we iterate
  //     r_i ← r_i − p(r_i) / ∏_{j ≠ i} (r_i − r_j)
  // until all updates fall below tol or iterCap is reached. Initial guesses
  // are the standard "spread points" 0.4·(0.9 + 0.9i)^k on a circle around
  // the polynomial's centroid — these break the symmetry that traps the
  // method on multiple roots.
  //
  // Handles complex coefficients, complex roots, multi-roots (with small loss
  // of precision per multiplicity), and degenerate cases. Sufficient for
  // degrees up to ~30 in this app (well above what the user is likely to need).
  //
  //   polynomialRoots(coeffs) → Complex[]
  //     coeffs in ascending-power order: [a_0, a_1, ..., a_n], a_n ≠ 0.
  //   Returns n roots (possibly repeated within tol).
  // ===========================================================================
  function polynomialRoots(coeffsIn, options) {
    options = options || {};
    const iterCap = options.iterCap || 200;
    const tol     = options.tol || 1e-13;

    const coeffs = coeffsIn.slice();
    // Strip trailing zeros.
    while (coeffs.length > 1 && C.abs(coeffs[coeffs.length - 1]) < 1e-300) {
      coeffs.pop();
    }
    const n = coeffs.length - 1;
    if (n <= 0) return [];

    // Normalize to monic.
    const an = coeffs[n];
    if (C.abs(an) < 1e-300) throw new Error("polynomialRoots: leading coefficient is zero");
    const anInv = C.inv(an);
    const a = coeffs.map(c => C.mul(c, anInv));   // a[n] = 1

    // Degree 1: trivial.
    if (n === 1) return [C.neg(a[0])];

    // Degree 2: closed-form (more accurate than Durand-Kerner here).
    if (n === 2) {
      // z² + b·z + c = 0  ⇒  z = (-b ± √(b² − 4c)) / 2
      const b = a[1], c = a[0];
      const disc = C.sub(C.mul(b, b), C.scale(c, 4));
      const sq = csqrt(disc);
      const z1 = C.scale(C.add(C.neg(b), sq), 0.5);
      const z2 = C.scale(C.sub(C.neg(b), sq), 0.5);
      return [z1, z2];
    }

    // Initial guesses: spread on a circle of radius R around polynomial centroid.
    // R = 1 + max_k |a_k| (Cauchy's bound on root magnitude). Centroid = −a_{n−1}/n.
    let R = 1;
    for (let k = 0; k < n; k++) R = Math.max(R, 1 + C.abs(a[k]));
    const cent = C.scale(a[n - 1], -1 / n);
    const roots = new Array(n);
    for (let i = 0; i < n; i++) {
      const ang = 2 * Math.PI * (i + 0.25) / n;     // off-axis to avoid symmetry traps
      roots[i] = {
        re: cent.re + 0.4 * R * Math.cos(ang),
        im: cent.im + 0.4 * R * Math.sin(ang),
      };
    }

    // Durand-Kerner iteration.
    for (let it = 0; it < iterCap; it++) {
      let maxDelta = 0;
      const next = new Array(n);
      for (let i = 0; i < n; i++) {
        const pi = evalPolyAscending(a, roots[i]);
        // ∏_{j ≠ i} (r_i − r_j)
        let denom = { re: 1, im: 0 };
        for (let j = 0; j < n; j++) {
          if (j === i) continue;
          denom = C.mul(denom, C.sub(roots[i], roots[j]));
        }
        if (C.abs(denom) < 1e-300) {
          // Coincident estimates: nudge slightly.
          next[i] = { re: roots[i].re + 1e-7, im: roots[i].im + 1e-7 };
          maxDelta = Math.max(maxDelta, 1e-7);
          continue;
        }
        const delta = C.div(pi, denom);
        next[i] = C.sub(roots[i], delta);
        const dm = C.abs(delta);
        if (dm > maxDelta) maxDelta = dm;
      }
      for (let i = 0; i < n; i++) roots[i] = next[i];
      if (maxDelta < tol) break;
    }
    return roots;
  }

  // Square root for complex numbers (principal branch).
  function csqrt(z) {
    const r = Math.hypot(z.re, z.im);
    if (r < 1e-300) return { re: 0, im: 0 };
    const u = Math.sqrt((r + z.re) / 2);
    const v = Math.sqrt((r - z.re) / 2) * Math.sign(z.im || 1);
    return { re: u, im: v };
  }

  // Horner-style eval of polynomial in ASCENDING-power form (a[0] + a[1]·z + ...).
  function evalPolyAscending(a, z) {
    let v = C.clone(a[a.length - 1]);
    for (let k = a.length - 2; k >= 0; k--) {
      v = C.add(C.mul(v, z), a[k]);
    }
    return v;
  }

  // Group roots within tolerance and return [{root, multiplicity}, ...]
  function groupRootsByMultiplicity(roots, tol) {
    tol = tol || 1e-6;
    const groups = [];
    for (const r of roots) {
      let found = false;
      for (const g of groups) {
        if (Math.hypot(r.re - g.root.re, r.im - g.root.im) < tol) {
          // Merge by averaging — accumulates multiplicity.
          const m = g.multiplicity;
          g.root = { re: (g.root.re * m + r.re) / (m + 1),
                     im: (g.root.im * m + r.im) / (m + 1) };
          g.multiplicity++;
          found = true;
          break;
        }
      }
      if (!found) groups.push({ root: { re: r.re, im: r.im }, multiplicity: 1 });
    }
    return groups;
  }

  Direct.polynomialRoots          = polynomialRoots;
  Direct.evalPolyAscending        = evalPolyAscending;
  Direct.groupRootsByMultiplicity = groupRootsByMultiplicity;

  // ===========================================================================
  // boundedQDRational: bounded classical QD direct problem for rational φ.
  // ---------------------------------------------------------------------------
  // INPUT
  //   P, Q :  Complex[] in ascending-power order (P[0] + P[1]·z + ... + P[p]·z^p)
  //   φ(z) = P(z) / Q(z), assumed analytic on 𝔻̄ (so Q has no zeros in 𝔻̄).
  //
  // OUTPUT
  //   {
  //     hData: { poles: [{a: w_j, principal: [C_{j,1}, ..., C_{j,k_j}]}, ...] },
  //     poleData: [{z: z_j, w: w_j, multiplicity}, ...],
  //     warnings: [...],
  //   }
  //
  // MATH (derivation in the user-facing plan accepted before implementation)
  //
  //   φ#(z) := conj(φ(1/conj(z)))  =  z^{q−p} · P̃(z) / Q̃(z)
  //
  // where p = deg P, q = deg Q, and X̃(z) = Σ_k conj(X_{deg X − k}) z^k is the
  // reverse-conjugate of polynomial X.
  //
  // On |z| = 1, σ ∘ φ = φ#. Analytic continuation into 𝔻 has poles at:
  //   • z = 0 with multiplicity (p − q), if p > q. (Maps to w_0 = φ(0).)
  //   • z = 1/conj(r_i) for each root r_i of Q, with multiplicity = mult(r_i).
  //     (Maps to w_j = φ(z_j) ∈ Ω.)
  //
  // For each such pole z_j with multiplicity k_j we:
  //   1. Extract the local Laurent of R̃ at z_j (principal-part coefficients d).
  //   2. Compute the local Taylor of φ at z_j (phiTilde, constant term zero).
  //   3. Apply QD.Faber.inverseFaberAtPole(d, phiTilde) to convert d → A,
  //      where the A's are the principal-part coefficients of σ at w_j in
  //      powers of (w − w_j).
  //
  // The result is hData with one principal-part entry per pole of R̃.
  // ===========================================================================
  function boundedQDRational(P, Q, options) {
    options = options || {};
    const validateTol = options.validateTol || 1e-6;

    // Trim trailing zeros (sanity).
    P = trimTrailingZeros(P);
    Q = trimTrailingZeros(Q);
    if (!Q || Q.length === 0) throw new Error("Direct.boundedQDRational: Q is the zero polynomial");
    if (!P || P.length === 0) throw new Error("Direct.boundedQDRational: P is the zero polynomial");

    const p = P.length - 1;
    const q = Q.length - 1;

    // ---- VALIDATE: Q must have no zeros in 𝔻̄. ----
    let qRoots = [];
    if (q >= 1) {
      qRoots = polynomialRoots(Q);
      for (const r of qRoots) {
        const mag = Math.hypot(r.re, r.im);
        if (mag <= 1 + validateTol) {
          throw new Error("Direct.boundedQDRational: Q has a root at z = " +
            r.re.toFixed(6) + (r.im >= 0 ? '+' : '') + r.im.toFixed(6) + 'i ' +
            '(|z| ≈ ' + mag.toFixed(4) + ' ≤ 1); φ is not analytic on the closed unit disk.');
        }
      }
    }
    // (q = 0 ⇒ Q is a nonzero constant; trivially no zeros.)

    // ---- BUILD R̃(z) = N(z) / D(z) ----
    //   q >= p :  N = z^{q−p} · P̃,            D = Q̃
    //   q <  p :  N = P̃,                       D = z^{p−q} · Q̃
    const Ptil = reverseConjugate(P);
    const Qtil = reverseConjugate(Q);
    let N, D;
    if (q >= p) {
      N = shiftPolynomialUp(Ptil, q - p);
      D = Qtil.slice();
    } else {
      N = Ptil.slice();
      D = shiftPolynomialUp(Qtil, p - q);
    }

    // ---- FIND POLES of R̃ inside 𝔻 (= roots of D). ----
    // D's zeros are: possibly z=0 (if p > q, with multiplicity p−q) and the
    // inverted roots of Q (= roots of Qtil = 1/conj(r_i)).
    const polesOfR = [];
    if (p > q) {
      polesOfR.push({ z: { re: 0, im: 0 }, multiplicity: p - q });
    }
    if (q >= 1) {
      // Inverted Q-roots: z_j = 1/conj(r_j).
      const inverted = qRoots.map(r => {
        const m2 = r.re * r.re + r.im * r.im;
        return { re: r.re / m2, im: r.im / m2 };
      });
      const groups = groupRootsByMultiplicity(inverted, 1e-7);
      for (const g of groups) polesOfR.push({ z: g.root, multiplicity: g.multiplicity });
    }

    if (polesOfR.length === 0) {
      // q = p = 0, i.e., φ is a constant. Degenerate.
      return {
        hData: { poles: [] },
        poleData: [],
        warnings: ['φ is a constant; h is identically zero.'],
      };
    }

    // ---- PER-POLE EXTRACTION ----
    const hPoles = [];
    const poleData = [];
    const warnings = [];

    for (const pole of polesOfR) {
      const zj = pole.z;
      const kj = pole.multiplicity;

      // 1. Local Taylor of N at z_j up to order kj − 1 (need k_j terms for F).
      //    We also need Taylor of D at z_j up to order 2·k_j or so (need
      //    coefficients beyond t^{k_j} for the reciprocal to give the
      //    truncated F up to order k_j − 1).
      //    Concretely: F(t) = N_T(t) / [t^{−k_j} · D_T(t)], so we need D_T's
      //    coefficients from index k_j up to index 2·k_j − 1.
      const Ntayl = polyTaylorAt(N, zj, kj - 1);             // length k_j
      const Dtayl = polyTaylorAt(D, zj, 2 * kj - 1);         // length 2·k_j

      // 2. D̃(t) = D(z_j + t) / t^{k_j} = [D_T[k_j], D_T[k_j+1], ..., D_T[2k_j-1]].
      //    Must have D̃(0) = D_T[k_j] ≠ 0 (else multiplicity was wrong).
      const Dtilde = new Array(kj);
      for (let i = 0; i < kj; i++) Dtilde[i] = Dtayl[kj + i] || { re: 0, im: 0 };
      if (C.abs(Dtilde[0]) < 1e-12) {
        warnings.push("Numerical issue: D's Taylor coefficient at t^" + kj +
                      " is near zero at pole z=" + complexFmt(zj) +
                      "; root multiplicity may be wrong.");
        continue;
      }

      // 3. F(t) = N_T(t) / D̃(t), Taylor up to order k_j − 1.
      const DtildeInv = T.reciprocal(Dtilde, kj - 1);
      const F = T.mul(Ntayl, DtildeInv, kj - 1);              // length k_j

      // 4. Principal part of R̃ at z_j: d_m = F[k_j − m] for m = 1..k_j.
      const d = new Array(kj);
      for (let m = 1; m <= kj; m++) d[m - 1] = F[kj - m];

      // 5. Image w_j = φ(z_j) = P(z_j) / Q(z_j).
      const Pzj = evalPolyAscending(P, zj);
      const Qzj = evalPolyAscending(Q, zj);
      if (C.abs(Qzj) < 1e-14) {
        warnings.push("Numerical issue: Q is near zero at the pole z=" + complexFmt(zj) +
                      "; image is ill-defined.");
        continue;
      }
      const wj = C.div(Pzj, Qzj);

      // 6. Local Taylor of φ at z_j up to order k_j (need φ', φ''/2!, ..., φ^(k_j)/k_j!).
      //    φ = P/Q. Taylor of φ at z_j = (Taylor of P at z_j) · (Taylor of 1/Q at z_j).
      const Pt = polyTaylorAt(P, zj, kj);                    // length k_j+1
      const Qt = polyTaylorAt(Q, zj, kj);                    // length k_j+1
      const QtInv = T.reciprocal(Qt, kj);                    // 1/Q Taylor
      const phiT = T.mul(Pt, QtInv, kj);                     // φ Taylor, length k_j+1

      // 7. phiTilde for the forward principal-part computation:
      //      phiTilde[0] = 0  (drop constant; locator absorbed by w_j)
      //      phiTilde[i] = i-th Taylor coefficient of φ at z_j for i ≥ 1.
      const phiTilde = T.zero(kj + 1);
      for (let i = 1; i <= kj; i++) phiTilde[i] = C.clone(phiT[i]);

      // 8. Forward principal-part conversion (the same primitive boundedQD
      //    uses for the polynomial φ case, applied here per R̃-pole).
      const A = forwardLocalPrincipal(d, phiTilde);

      hPoles.push({ a: wj, principal: A });
      poleData.push({ z: zj, w: wj, multiplicity: kj });
    }

    return {
      hData: { poles: hPoles },
      poleData,
      warnings,
    };
  }

  // ---------------------------------------------------------------------------
  // Polynomial helpers (in ascending-power Complex[] form).
  // ---------------------------------------------------------------------------

  // Trim trailing zero coefficients (keep at least one entry).
  function trimTrailingZeros(p) {
    if (!p) return p;
    let n = p.length;
    while (n > 1 && C.abs(p[n - 1]) < 1e-300) n--;
    return p.slice(0, n);
  }

  // X̃(z) = Σ_k conj(X[deg X − k]) z^k. Reverse the coefficient list and
  // conjugate each entry.
  function reverseConjugate(p) {
    const n = p.length;
    const out = new Array(n);
    for (let k = 0; k < n; k++) out[k] = C.conj(p[n - 1 - k]);
    return out;
  }

  // z^m · p(z): prepend m zeros to the coefficient list.
  function shiftPolynomialUp(p, m) {
    if (m <= 0) return p.slice();
    const out = new Array(m + p.length);
    for (let i = 0; i < m; i++) out[i] = { re: 0, im: 0 };
    for (let i = 0; i < p.length; i++) out[m + i] = C.clone(p[i]);
    return out;
  }

  // Taylor expansion of polynomial p at z = z0, up to order L.
  // Returns [p(z_0), p'(z_0)/1!, ..., p^{(L)}(z_0)/L!], length L+1.
  //
  // Repeated synthetic-division: dividing by (z − z_0) repeatedly produces
  // successive remainders that are exactly the Taylor coefficients at z_0.
  function polyTaylorAt(p, z0, L) {
    const n = p.length - 1;
    const out = new Array(L + 1);
    let work = p.slice();                              // mutable copy
    for (let k = 0; k <= L; k++) {
      if (work.length === 0) { out[k] = { re: 0, im: 0 }; continue; }
      // Synthetic division of `work` by (z − z_0): result is quotient (length-1)
      // and remainder = work[0]'s replacement (which equals work_evaluated_at_z0).
      // Standard inner loop with z0 as the test point.
      const m = work.length;
      const q = new Array(m - 1);                      // quotient
      let rem = work[m - 1];                           // working accumulator
      for (let i = m - 2; i >= 0; i--) {
        q[i] = C.clone(rem);
        rem = C.add(work[i], C.mul(rem, z0));
      }
      out[k] = rem;                                    // = p(z_0) on first pass, etc.
      work = q;
    }
    return out;
  }

  function complexFmt(c) {
    return c.re.toFixed(4) + (c.im >= 0 ? '+' : '') + c.im.toFixed(4) + 'i';
  }

  // ===========================================================================
  // forwardLocalPrincipal: principal-part of σ at a local pole.
  // ---------------------------------------------------------------------------
  // Given:
  //   d         : Complex[]    residues d_1, ..., d_m  of  R̃(z_j + t) in t.
  //                            (R̃(z_j + t) = d_m/t^m + ... + d_1/t + regular)
  //   phiTilde  : Complex[]    length m+1, phiTilde[0] = 0,
  //                            phiTilde[i] = i-th Taylor coefficient of φ at z_j.
  //
  // Returns:
  //   C : Complex[m]   such that  σ(w_j + ζ) − regular = Σ_{k=1..m} C_k · ζ^{-k}
  //                    where w_j = φ(z_j).
  //
  // Formula:
  //   ψ̃(ζ) = Taylor-inverse of phiTilde (so φ̃(ψ̃(ζ)) = ζ; ψ̃[0]=0, ψ̃[1]=1/c_1).
  //   u(ζ) := ψ̃(ζ) / (ψ̃[1] · ζ),  u(0) = 1.
  //   C_k = Σ_{l ≥ k}  d_l · (1/ψ̃[1])^l · [ζ^{l−k}] u(ζ)^{−l}
  //       = Σ_{l ≥ k}  d_l · c_1^l       · [ζ^{l−k}] u(ζ)^{−l},
  // where c_1 = phiTilde[1] = φ'(z_j) (since ψ̃[1] = 1/c_1).
  //
  // This is the SAME primitive used by boundedQD (polynomial φ case); the
  // only difference is that the polynomial case had d_l = conj(c_l) (residues
  // of φ# at z=0).
  // ===========================================================================
  function forwardLocalPrincipal(d, phiTilde) {
    const m = d.length;
    if (m === 0) return [];
    if (phiTilde.length < m + 1) {
      throw new Error("forwardLocalPrincipal: phiTilde must have length ≥ m+1 (got " +
                      phiTilde.length + ", need " + (m + 1) + ")");
    }
    if (C.abs(phiTilde[1]) < 1e-14) {
      throw new Error("forwardLocalPrincipal: phiTilde[1] = 0; φ has a critical point at the pole");
    }

    // ψ̃ = Taylor inverse of phiTilde, length m+1.
    const psi = T.invert(phiTilde, m);
    const psi1Inv = C.inv(psi[1]);                 // = c_1 = phiTilde[1]
    const c1 = psi1Inv;

    // u(ζ) = ψ̃(ζ) / (ψ̃[1] · ζ), as a Taylor of length m:  u[i] = ψ̃[i+1] · c_1.
    const u = T.zero(m);
    for (let i = 0; i < m; i++) {
      u[i] = C.mul(psi[i + 1], psi1Inv);
    }

    // u^{-l} for l = 1..m, each truncated to ζ-degree m−1.
    const uInv = T.reciprocal(u, m - 1);
    const uPowNeg = [null];
    uPowNeg[1] = T.truncate(uInv, m - 1);
    for (let l = 2; l <= m; l++) {
      uPowNeg[l] = T.mul(uPowNeg[l - 1], uInv, m - 1);
    }

    // c_1^l for l = 0..m.
    const c1Pow = [{ re: 1, im: 0 }];
    for (let l = 1; l <= m; l++) c1Pow.push(C.mul(c1Pow[l - 1], c1));

    // C_k = Σ_{l ≥ k}  d_l · c_1^l · [ζ^{l−k}] u^{−l}.
    const out = new Array(m);
    for (let k = 1; k <= m; k++) {
      let acc = { re: 0, im: 0 };
      for (let l = k; l <= m; l++) {
        const idx = l - k;
        if (idx >= uPowNeg[l].length) continue;
        const term = C.mul(d[l - 1], c1Pow[l]);
        acc = C.add(acc, C.mul(term, uPowNeg[l][idx]));
      }
      out[k - 1] = acc;
    }
    return out;
  }

  Direct.boundedQDRational  = boundedQDRational;
  Direct.forwardLocalPrincipal = forwardLocalPrincipal;
  Direct.reverseConjugate   = reverseConjugate;
  Direct.shiftPolynomialUp  = shiftPolynomialUp;
  Direct.polyTaylorAt       = polyTaylorAt;
  Direct.trimTrailingZeros  = trimTrailingZeros;

  Direct.boundedQD                = boundedQD;
  Direct.unboundedQD              = unboundedQD;
  Direct.numericalBoundedQD       = numericalBoundedQD;
  Direct.evalH                    = evalH;
  Direct.verifyBoundaryIdentity   = verifyBoundaryIdentity;
  Direct.sampleBoundaryPolynomial = sampleBoundaryPolynomial;
  Direct.sampleBoundaryLaurent    = sampleBoundaryLaurent;
  Direct.parsePolynomialInZ       = parsePolynomialInZ;
  Direct.parseRationalInZ         = parseRationalInZ;
  Direct.polynomialToString       = polynomialToString;

  if (typeof module !== 'undefined' && module.exports) {
    module.exports.Direct = Direct;
  }

}(typeof window !== 'undefined' ? window
   : typeof global   !== 'undefined' ? global
   : globalThis));
