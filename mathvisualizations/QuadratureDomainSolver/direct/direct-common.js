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
// WEIGHTED families take the rational KERNEL (R#/r#) rather than φ directly:
//   • bounded power/log     boundedPowerQD / boundedLogQD                (0∉Ω)
//   • bounded singular       boundedPowerQDSingular / boundedLogQDSingular (0∈Ω)
//     — see the "SINGULAR weighted forward kernels" section header.
//   • unbounded power/log + singular  unboundedPowerQD[Singular] /
//     unboundedLogQD[Singular] (∞∈Ω) — see the "UNBOUNDED weighted forward
//     kernels (Theorem 4.3.7)" section header for the full math + conventions.
//
// All variants produce hData of the same shape used by the inverse-problem
// solver:  hData = { poles: [{a, principal: [...]}, ...], polyPart: [...] }
// (singular variants additionally return an origin term separately: r₀ for
// bounded-PQD, q for bounded-/unbounded-LQD — unbounded-PQD-singular has none).
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

  // ===========================================================================
  // WEIGHTED-FAMILY FORWARD KERNELS (φ → h) — bounded non-singular LQD + PQD.
  // ---------------------------------------------------------------------------
  // Classically rational φ ⟺ QD, so the Direct tab takes a rational φ. For the
  // weighted families the RATIONAL class is the KERNEL, not φ:
  //   • PQD:  φ = (R#(z))^{1/α},  R# rational  ⟺  power-weighted QD
  //   • LQD:  φ = w₀·exp(r#(z)),  r# rational  ⟺  log-weighted QD
  // So these kernels take the rational kernel K (= R# / r#), reuse the INVERSE
  // solver's parametrization, and read h off by INVERTING the same (★) chain the
  // inverse solver encodes — guaranteed-correct, no new weighted-Faber math.
  //
  // Per pole the (★) chain C_j → A_j is LINEAR and upper-triangular:
  //   PQD:  A_{j,k} = α·inverseFaberAtPole( modifiedResidues_PQD(C_j,α), φ̃_j )[k]
  //   LQD:  A_{j,k} =   inverseFaberAtPole( LqdModifiedResidues(C_j,a_j), φ̃_j )[k]
  // We build that small matrix by PROBING the existing (tested) forward functions
  // on unit residue vectors, then back-substitute for C_j — so no inverse formula
  // is hand-derived. A_j (the kernel's Möbius-branch coeffs) is recovered from the
  // kernel's principal part at ζ_p = 1/conj(z_j) by another small triangular solve.
  // ===========================================================================

  // Upper-triangular complex back-substitution: solve M·x = b where M[i][j] = 0
  // for j < i (rows/cols 0-indexed, length n). Returns x (Complex[n]).
  function solveUpperTriComplex(M, b, n) {
    const x = new Array(n);
    for (let i = n - 1; i >= 0; i--) {
      let acc = C.clone(b[i]);
      for (let j = i + 1; j < n; j++) acc = C.sub(acc, C.mul(M[i][j], x[j]));
      if (C.abs(M[i][i]) < 1e-300) throw new Error("solveUpperTriComplex: singular diagonal at " + i);
      x[i] = C.div(acc, M[i][i]);
    }
    return x;
  }

  // Principal-part residues of a rational N/D at a root z0 of D of multiplicity
  // mult: returns [d_1, …, d_mult] with d_s = coeff of (z − z0)^{-s}. Mirrors the
  // per-pole extraction inside boundedQDRational (steps 1–4).
  function localPrincipalResidues(N, D, z0, mult) {
    const Ntayl = polyTaylorAt(N, z0, mult - 1);            // length mult
    const Dtayl = polyTaylorAt(D, z0, 2 * mult - 1);        // length 2·mult
    const Dtilde = new Array(mult);
    for (let i = 0; i < mult; i++) Dtilde[i] = Dtayl[mult + i] || { re: 0, im: 0 };
    if (C.abs(Dtilde[0]) < 1e-12) throw new Error("localPrincipalResidues: wrong multiplicity at z=" + complexFmt(z0));
    const DtildeInv = T.reciprocal(Dtilde, mult - 1);
    const F = T.mul(Ntayl, DtildeInv, mult - 1);            // length mult
    const d = new Array(mult);
    for (let s = 1; s <= mult; s++) d[s - 1] = F[mult - s];
    return d;
  }

  // Möbius-branch coefficients A_{j,1..m} of a kernel from its principal-part
  // residues dR at ζ_p = 1/conj(z_j). The branch basis term
  //   basis_k(z) = z^k / (1 − conj(z_j) z)^k
  // has, at ζ_p, principal part  Σ_s B[s][k] (z−ζ_p)^{-s}  with (k ≥ s)
  //   B[s][k] = (−1)^k · binom(k, k−s) · ζ_p^{k+s}
  // (since 1 − conj(z_j) z = −conj(z_j)(z − ζ_p) and conj(z_j) = 1/ζ_p). Then
  //   dR_s = Σ_{k≥s} conj(A_{j,k}) · B[s][k],  an upper-triangular solve for conj(A).
  function branchCoeffsFromResidues(dR, zeta, m) {
    // Integer binomials binom(k, k−s) for 1 ≤ s ≤ k ≤ m.
    const nCr = (nn, rr) => { let r = 1; for (let i = 0; i < rr; i++) r = r * (nn - i) / (i + 1); return Math.round(r); };
    // ζ^p for p = 0..2m.
    const zpow = [{ re: 1, im: 0 }];
    for (let p = 1; p <= 2 * m; p++) zpow.push(C.mul(zpow[p - 1], zeta));
    const B = [];
    for (let s = 1; s <= m; s++) {
      const row = new Array(m).fill(null).map(() => ({ re: 0, im: 0 }));
      for (let k = s; k <= m; k++) {
        const sign = (k % 2 === 0) ? 1 : -1;
        row[k - 1] = C.scale(zpow[k + s], sign * nCr(k, k - s));
      }
      B.push(row);
    }
    const conjA = solveUpperTriComplex(B, dR, m);            // conj(A_{j,k})
    return conjA.map(x => C.conj(x));                        // A_{j,k}
  }

  // Real-scale a Taylor series (multiply every coefficient by a real factor).
  function scaleTaylorReal(tay, s) { return tay.map(c => ({ re: c.re * s, im: c.im * s })); }

  // boundedPowerQD: rational R#(z) + α → hData (bounded power-weighted QD, 0∉Ω).
  //   Rhash:  { num: Complex[], den: Complex[] }  (ascending-power), or Complex[]
  //           (polynomial ⇒ den = [1], degenerate: no finite poles).
  //   alpha:  real > 0, ≠ 1.
  // Returns { hData:{poles:[{a,principal}]}, poleData, w0, alpha, warnings }.
  function boundedPowerQD(Rhash, alpha, options) {
    options = options || {};
    const validateTol = options.validateTol || 1e-6;
    if (!(alpha > 0) || Math.abs(alpha - 1) < 1e-12) {
      throw new Error("Direct.boundedPowerQD: need α > 0, α ≠ 1 (α = 1 is the classical QD)");
    }
    let num, den;
    if (Array.isArray(Rhash)) { num = trimTrailingZeros(Rhash.slice()); den = [{ re: 1, im: 0 }]; }
    else { num = trimTrailingZeros(Rhash.num.slice()); den = trimTrailingZeros(Rhash.den.slice()); }

    // R# must be analytic AND non-vanishing on 𝔻̄: den roots and num roots all
    // strictly outside the closed disk (the αth root is single-valued only then).
    const evalRat = (z) => C.div(evalPolyAscending(num, z), evalPolyAscending(den, z));
    if (den.length > 1) for (const r of polynomialRoots(den)) {
      if (Math.hypot(r.re, r.im) <= 1 + validateTol)
        throw new Error("Direct.boundedPowerQD: R# has a pole at |z| ≤ 1; not analytic on 𝔻̄.");
    }
    if (num.length > 1) for (const r of polynomialRoots(num)) {
      if (Math.hypot(r.re, r.im) <= 1 + validateTol)
        throw new Error("Direct.boundedPowerQD: R# has a zero at |z| ≤ 1; (R#)^{1/α} is not single-valued.");
    }

    // w₀ = φ(0) = (R#(0))^{1/α} on the principal branch (the gauge anchor).
    const R0 = evalRat({ re: 0, im: 0 });
    const w0 = C.cpow(R0, 1 / alpha);
    const anchorArg0 = Math.atan2(R0.im, R0.re);            // = α·arg(w0)

    // Poles ζ_p of R# (roots of den, |ζ_p| > 1) ⇒ node preimages z_j = 1/conj(ζ_p).
    const hPoles = [], poleData = [], warnings = [];
    if (den.length <= 1) {
      return { hData: { poles: [] }, poleData: [], w0, alpha,
               warnings: ['R# is a polynomial (no finite poles) — degenerate; h is constant.'] };
    }
    const groups = groupRootsByMultiplicity(polynomialRoots(den), 1e-7);
    const evalRHashRaw = (z /*, phi */) => evalRat(z);

    for (const g of groups) {
      const zeta = g.root, m = g.multiplicity;
      const m2 = zeta.re * zeta.re + zeta.im * zeta.im;
      const zj = { re: zeta.re / m2, im: -zeta.im / m2 };    // 1/conj(ζ_p)

      // A_j: kernel Möbius-branch coeffs from R#'s principal part at ζ_p.
      const dR = localPrincipalResidues(num, den, zeta, m);
      const A = branchCoeffsFromResidues(dR, zeta, m);

      // a_j = φ(z_j) on the anchored αth-root branch (continuous from φ(0)=w₀).
      const argZj = QD.PqdCommon.argContAt(null, zj, evalRHashRaw, anchorArg0, { re: 0, im: 0 });
      const Rzj = evalRat(zj);
      const mag = Math.pow(C.abs2(Rzj), 0.5 / alpha);
      const aj = { re: mag * Math.cos(argZj / alpha), im: mag * Math.sin(argZj / alpha) };

      // φ̃_j = Taylor of (φ(z_j+t) − a_j), anchored: φ = exp((1/α)·logR#) with the
      // log-constant overridden to the anchored value (so the WHOLE series sits on
      // the correct sheet — same trick as phiTaylorAt_PQD).
      const NT = polyTaylorAt(num, zj, m), DT = polyTaylorAt(den, zj, m);
      const RT = T.mul(NT, T.reciprocal(DT, m), m);          // Taylor of R# at z_j
      const L = T.log(RT, m);
      L[0] = { re: L[0].re, im: argZj };                     // anchored log constant
      const phiFull = T.exp(scaleTaylorReal(L, 1 / alpha), m); // φ Taylor, phiFull[0]=a_j
      const phiTilde = T.zero(m + 1);
      for (let i = 1; i <= m; i++) phiTilde[i] = C.clone(phiFull[i]);

      // C_j: invert the (★) chain. Build the upper-triangular map C → A by probing
      // the tested forward functions on unit residue vectors, then back-substitute.
      const Mmat = [];
      for (let k = 0; k < m; k++) Mmat.push(new Array(m));
      for (let c = 0; c < m; c++) {
        const e = new Array(m).fill(null).map((_, i) => ({ re: i === c ? 1 : 0, im: 0 }));
        const D = QD.modifiedResidues_PQD({ a: aj, principal: e }, alpha);
        const Aprobe = QD.Faber.inverseFaberAtPole(D, phiTilde).map(x => C.scale(x, alpha));
        for (let k = 0; k < m; k++) Mmat[k][c] = Aprobe[k];
      }
      const Cj = solveUpperTriComplex(Mmat, A, m);

      hPoles.push({ a: aj, principal: Cj });
      poleData.push({ z: zj, w: aj, multiplicity: m });
    }

    return { hData: { poles: hPoles }, poleData, w0, alpha, warnings };
  }

  // boundedLogQD: rational r#(z) + w₀ → hData (bounded log-weighted QD, 0∉Ω̄).
  //   φ = w₀·exp(r#(z)). exp is entire, so there is no branch/anchoring issue.
  function boundedLogQD(rhash, w0, options) {
    options = options || {};
    const validateTol = options.validateTol || 1e-6;
    let num, den;
    if (Array.isArray(rhash)) { num = trimTrailingZeros(rhash.slice()); den = [{ re: 1, im: 0 }]; }
    else { num = trimTrailingZeros(rhash.num.slice()); den = trimTrailingZeros(rhash.den.slice()); }
    if (!w0 || C.abs(w0) < 1e-14) throw new Error("Direct.boundedLogQD: w₀ (= φ(0)) must be nonzero");

    const evalRat = (z) => C.div(evalPolyAscending(num, z), evalPolyAscending(den, z));
    if (den.length > 1) for (const r of polynomialRoots(den)) {
      if (Math.hypot(r.re, r.im) <= 1 + validateTol)
        throw new Error("Direct.boundedLogQD: r# has a pole at |z| ≤ 1; not analytic on 𝔻̄.");
    }

    const hPoles = [], poleData = [], warnings = [];
    if (den.length <= 1) {
      return { hData: { poles: [] }, poleData: [], w0, warnings: ['r# is a polynomial (no finite poles) — degenerate.'] };
    }
    const cexp = (z) => { const e = Math.exp(z.re); return { re: e * Math.cos(z.im), im: e * Math.sin(z.im) }; };
    const groups = groupRootsByMultiplicity(polynomialRoots(den), 1e-7);

    for (const g of groups) {
      const zeta = g.root, m = g.multiplicity;
      const m2 = zeta.re * zeta.re + zeta.im * zeta.im;
      const zj = { re: zeta.re / m2, im: -zeta.im / m2 };

      const dR = localPrincipalResidues(num, den, zeta, m);
      const A = branchCoeffsFromResidues(dR, zeta, m);

      // a_j = w₀·exp(r#(z_j)); φ Taylor = w₀·exp(Taylor of r# at z_j).
      const aj = C.mul(w0, cexp(evalRat(zj)));
      const NT = polyTaylorAt(num, zj, m), DT = polyTaylorAt(den, zj, m);
      const rT = T.mul(NT, T.reciprocal(DT, m), m);          // Taylor of r# at z_j
      const expR = T.exp(rT, m);
      const phiFull = expR.map(c => C.mul(w0, c));           // phiFull[0] = a_j
      const phiTilde = T.zero(m + 1);
      for (let i = 1; i <= m; i++) phiTilde[i] = C.clone(phiFull[i]);

      // C_j: invert the LQD (★): A = inverseFaberAtPole( D(C), φ̃ ),
      // D_{s} = a_j·C_s + C_{s+1}. Probe + back-substitute.
      const Mmat = [];
      for (let k = 0; k < m; k++) Mmat.push(new Array(m));
      for (let c = 0; c < m; c++) {
        const e = new Array(m).fill(null).map((_, i) => ({ re: i === c ? 1 : 0, im: 0 }));
        const D = QD.LqdCommon.modifiedResidues({ poles: [{ a: aj, principal: e }] })[0];
        const Aprobe = QD.Faber.inverseFaberAtPole(D, phiTilde);
        for (let k = 0; k < m; k++) Mmat[k][c] = Aprobe[k];
      }
      const Cj = solveUpperTriComplex(Mmat, A, m);

      hPoles.push({ a: aj, principal: Cj });
      poleData.push({ z: zj, w: aj, multiplicity: m });
    }

    return { hData: { poles: hPoles }, poleData, w0, warnings };
  }

  // ===========================================================================
  // SINGULAR weighted forward kernels (0 ∈ Ω). φ carries a Blaschke factor
  // b_{z₀}(z) whose zero z₀ ∈ 𝔻 is the preimage of the origin. z₀ is a FREE input
  // in both (it sets where 0 sits inside Ω); by Theorem 4.3.3 the domain is a
  // weighted QD for any rational kernel with a univalent φ. The two kernels use
  // DIFFERENT routes to the residues:
  //   • powerQD_singular     φ = b_{z₀}·(R#)^{1/α}   — `boundedPowerQDSingular`,
  //       computed from the authoritative forward map Theorem 4.3.5; h = finite
  //       poles + an origin term r₀/w. See that function's header for the formula.
  //   • boundedLQD_singular  φ = γ·b_{z₀}·exp(r#)    — `boundedLogQDSingular`,
  //       residues recovered by INVERTING the family's per-pole (★) builder
  //       (`computeTargets`, block-diagonal per pole, via solveResiduesViaProbe —
  //       the LQD (★) has no clean closed forward form here); h gains q/w.
  // ===========================================================================

  // Recover the per-pole residues C_j by probing a family's (★) builder (used by
  // boundedLogQDSingular). For each pole j, the map C_j → A_j (target branch
  // coeffs) is upper-triangular; build it column-by-column with unit residues and
  // back-substitute. `poles` is [{a, m, A}] (node, multiplicity, target branch coeffs).
  function solveResiduesViaProbe(fam, phi, poles) {
    const zeros = (m) => Array.from({ length: m }, () => ({ re: 0, im: 0 }));
    const out = [];
    for (let j = 0; j < poles.length; j++) {
      const mj = poles[j].m;
      const M = []; for (let k = 0; k < mj; k++) M.push(new Array(mj));
      for (let c = 0; c < mj; c++) {
        const probe = { poles: poles.map((p, i) => ({
          a: p.a,
          principal: (i === j)
            ? Array.from({ length: p.m }, (_, t) => ({ re: t === c ? 1 : 0, im: 0 }))
            : zeros(p.m),
        })) };
        const A = fam.computeTargets(phi, probe).A[j];
        for (let k = 0; k < mj; k++) M[k][c] = A[k];
      }
      out.push(solveUpperTriComplex(M, poles[j].A, mj));
    }
    return out;
  }

  // Shared front-end: parse a rational kernel + validate it's analytic on 𝔻̄
  // (poles strictly outside), returning { num, den, evalRat, groups }. `needNonVanishing`
  // additionally rejects zeros of the kernel inside 𝔻̄ (for the αth root).
  function parseKernel(K, validateTol, needNonVanishing, label) {
    let num, den;
    if (Array.isArray(K)) { num = trimTrailingZeros(K.slice()); den = [{ re: 1, im: 0 }]; }
    else { num = trimTrailingZeros(K.num.slice()); den = trimTrailingZeros(K.den.slice()); }
    if (den.length > 1) for (const r of polynomialRoots(den)) {
      if (Math.hypot(r.re, r.im) <= 1 + validateTol)
        throw new Error(label + ': kernel has a pole at |z| ≤ 1; not analytic on 𝔻̄.');
    }
    if (needNonVanishing && num.length > 1) for (const r of polynomialRoots(num)) {
      if (Math.hypot(r.re, r.im) <= 1 + validateTol)
        throw new Error(label + ': kernel has a zero at |z| ≤ 1; the αth root is not single-valued.');
    }
    const evalRat = (z) => C.div(evalPolyAscending(num, z), evalPolyAscending(den, z));
    const groups = den.length > 1 ? groupRootsByMultiplicity(polynomialRoots(den), 1e-7) : [];
    return { num, den, evalRat, groups };
  }

  function validateZ0(z0, label) {
    if (!z0 || C.abs2(z0) < 1e-18) throw new Error(label + ': z₀ must be nonzero (interior preimage of the origin).');
    if (C.abs(z0) >= 1 - 1e-9) throw new Error(label + ': z₀ must satisfy 0 < |z₀| < 1.');
  }

  // Weighted area t = ∫_Ω |w|^{2(α−1)} dA via the Green's-form boundary integral
  //   t = (1/(αN)) Σ conj(w)·|w|^{2(α−1)}·(φ'·z),  z = e^{iθ}
  // (φ'·z obtained by a radial central difference; φ from QD.evalPhi). Returns the
  // complex value (≈ real for a true PQD); used as the t/w normalization term.
  function weightedAreaPQD(phi, alpha, N) {
    const h = 1e-6; let acc = { re: 0, im: 0 };
    for (let n = 0; n < N; n++) {
      const th = 2 * Math.PI * n / N, z = { re: Math.cos(th), im: Math.sin(th) };
      const w = QD.evalPhi(z, phi);
      const zp = { re: z.re * (1 + h), im: z.im * (1 + h) }, zm = { re: z.re * (1 - h), im: z.im * (1 - h) };
      const wp = QD.evalPhi(zp, phi), wm = QD.evalPhi(zm, phi);
      const pz = { re: (wp.re - wm.re) / (2 * h), im: (wp.im - wm.im) / (2 * h) };   // φ'·z
      const w2 = w.re * w.re + w.im * w.im;
      if (w2 < 1e-30) continue;
      acc = C.add(acc, C.mul(C.scale(C.conj(w), Math.pow(w2, alpha - 1)), pz));
    }
    return C.scale(acc, 1 / (alpha * N));
  }

  // boundedPowerQDSingular: rational R#(z) + α + z₀ → hData (bounded SINGULAR PQD,
  // 0∈Ω), via the AUTHORITATIVE forward map Theorem 4.3.5 (Eq 4.13):
  //   h(w) = (1/(α·w))·Φ_φ( AnalyticIn_{𝔻∁}[ r·r# ] )(w) + t/w,
  // where r#=R#, r(z)=conj(r#(1/conj z)) is its Schwarz reflection (so r·r# is
  // rational with in-𝔻 poles exactly at the node-preimages z_j), Φ_φ is the
  // forward Faber transform (forwardLocalPrincipal), and t = ∫_Ω|w|^{2(α−1)}dA is
  // the weighted area — the t/w NORMALIZATION term (returned separately, NOT a
  // quadrature node). By Theorem 4.3.3 any rational R# with a univalent φ is a PQD,
  // so z₀ is FREE and realizability ⟺ univalence of φ.
  function boundedPowerQDSingular(Rhash, alpha, z0, options) {
    options = options || {};
    const label = 'Direct.boundedPowerQDSingular';
    if (!(alpha > 0) || Math.abs(alpha - 1) < 1e-12) throw new Error(label + ': need α > 0, α ≠ 1.');
    validateZ0(z0, label);
    const { num, den, evalRat, groups } = parseKernel(Rhash, options.validateTol || 1e-6, true, label);
    if (den.length <= 1) return { hData: { poles: [] }, poleData: [], z0: C.clone(z0), t: { re: 0, im: 0 }, warnings: ['R# is a polynomial (no finite poles) — degenerate.'] };
    if ((num.length - 1) > (den.length - 1)) throw new Error(label + ': R# must be a proper rational (deg num ≤ deg den).');

    const absZ0 = C.abs(z0);
    const R0 = evalRat({ re: 0, im: 0 });
    const w0 = C.scale(C.cpow(R0, 1 / alpha), absZ0);            // |z₀|·(R#(0))^{1/α}
    const anchorArg0 = Math.atan2(R0.im, R0.re);                 // = α·arg(w0)
    const evalRHashRaw = (z) => evalRat(z);

    // Family phi (branches from R#'s principal parts) — for evalPhi / univalence /
    // the family verifier on the returned φ.
    const branches = [], geom = [];
    for (const g of groups) {
      const zeta = g.root, m = g.multiplicity;
      const m2 = zeta.re * zeta.re + zeta.im * zeta.im;
      const zj = { re: zeta.re / m2, im: -zeta.im / m2 };
      const A = branchCoeffsFromResidues(localPrincipalResidues(num, den, zeta, m), zeta, m);
      branches.push({ z: zj, A }); geom.push({ zj, m });
    }
    // Guard: z₀ must not coincide with a node-preimage z_j (else a_j = φ(z_j) = 0,
    // i.e. a finite node collides with the origin — degenerate).
    for (const g of geom) {
      if (C.abs(C.sub(z0, g.zj)) < 1e-6)
        throw new Error(label + ': z₀ coincides with a node preimage (node at the origin) — choose a different z₀.');
    }

    const phi = { family: 'powerQD_singular', unbounded: false, alpha, w0, z0: C.clone(z0), branches };

    // REALIZABILITY ⟺ univalence of φ (Thm 4.3.3 ⇒ QD-property is then automatic).
    let univalent = false;
    try { univalent = QD.isBoundaryUnivalent(phi); } catch (e) { univalent = false; }
    if (!univalent) {
      return { hData: { poles: [] }, poleData: [], w0, z0: C.clone(z0), t: { re: 0, im: 0 }, phi, univalent: false,
               warnings: ['Not realizable: φ = b_{z₀}·(R#)^{1/α} is not univalent (boundary self-intersects). Adjust R#, z₀, or α.'] };
    }

    // r·r# (rational). r = reflection of R#=N/D: r = z^{degD−degN}·Ñ/D̃ with
    // Ñ=reverseConjugate(N), D̃=reverseConjugate(D). r·r# = z^{…}·(Ñ·N)/(D̃·D).
    const Ntil = reverseConjugate(num), Dtil = reverseConjugate(den);
    let numRR = mulPolys(Ntil, num);
    if (den.length - num.length > 0) numRR = shiftPolynomialUp(numRR, den.length - num.length);
    const denRR = mulPolys(Dtil, den);

    // Per node: AnalyticIn_{𝔻∁}[r·r#] principal part at z_j → Φ_φ → (1/α)·(w·h)
    // modified residues → h residues C_j.
    const hPoles = [], poleData = [];
    for (const g of geom) {
      const zj = g.zj, m = g.m;
      const d = localPrincipalResidues(numRR, denRR, zj, m);     // residues of r·r# at z_j

      const aj = QD.evalPhi(zj, phi);                            // a_j = φ(z_j)
      // Anchored singular Taylor φ̃_j = b_{z₀}·(R#)^{1/α} at z_j (constant dropped).
      const argZj = QD.PqdCommon.argContAt(null, zj, evalRHashRaw, anchorArg0, { re: 0, im: 0 });
      const RT = T.mul(polyTaylorAt(num, zj, m), T.reciprocal(polyTaylorAt(den, zj, m), m), m);
      const L = T.log(RT, m); L[0] = { re: L[0].re, im: argZj };
      const root = T.exp(scaleTaylorReal(L, 1 / alpha), m);      // (R#)^{1/α} anchored
      const phiFull = T.mul(QD.LqdCommon.blaschkeTaylor(zj, z0, m), root, m);
      const phiTilde = T.zero(m + 1);
      for (let i = 1; i <= m; i++) phiTilde[i] = C.clone(phiFull[i]);

      // Φ_φ(g) principal part at a_j; ÷α gives the (w·h) modified residues D̃_k.
      const Dt = forwardLocalPrincipal(d, phiTilde).map(p => C.scale(p, 1 / alpha));
      // h residues from w·h: D̃_k = a_j·C_k + C_{k+1} ⇒ C_k = (D̃_k − C_{k+1})/a_j.
      const Cj = new Array(m);
      for (let k = m; k >= 1; k--) {
        const above = (k < m) ? Cj[k] : { re: 0, im: 0 };
        Cj[k - 1] = C.div(C.sub(Dt[k - 1], above), aj);
      }
      hPoles.push({ a: aj, principal: Cj });
      poleData.push({ z: zj, w: aj, multiplicity: m });
    }

    const t = weightedAreaPQD(phi, alpha, 2000);                 // weighted area ∫_Ω|w|^{2(α−1)}dA
    // Origin charge (the t/w term of Eq 4.13): the f=1 quadrature identity forces
    //   ∫_Ω|w|^{2(α−1)}dA = Σ_j C_{j,1} + r₀,  so r₀ = t − Σ C_{j,1}.
    // h(w) = Σ_finite C_{j,k}/(w−a_j)^k + r₀/w. r₀ vanishes exactly on the
    // mass-constrained "no-origin-charge" class the inverse solver targets.
    let sumC1 = { re: 0, im: 0 };
    for (const p of hPoles) sumC1 = C.add(sumC1, p.principal[0]);
    const originResidue = C.sub(t, sumC1);

    return { hData: { poles: hPoles }, poleData, w0, z0: C.clone(z0), t, originResidue, phi, univalent: true, warnings: [] };
  }

  // boundedLogQDSingular: rational r#(z) + w₀ + z₀ → hData (bounded LQD, 0∈Ω).
  // h has an ORIGIN pole at w=0 with residue q = ln|γ|² + r#(z₀) + conj(r#(1/conj z₀)),
  // returned separately (the inverse takes it via opts.q).
  function boundedLogQDSingular(rhash, w0, z0, options) {
    options = options || {};
    const label = 'Direct.boundedLogQDSingular';
    if (!w0 || C.abs(w0) < 1e-14) throw new Error(label + ': w₀ (= φ(0)) must be nonzero.');
    validateZ0(z0, label);
    const { num, den, groups } = parseKernel(rhash, options.validateTol || 1e-6, false, label);
    if (den.length <= 1) return { hData: { poles: [] }, poleData: [], w0, z0, q: { re: 0, im: 0 }, warnings: ['r# is a polynomial (no finite poles) — degenerate.'] };

    const absZ0 = C.abs(z0);
    const gamma = C.scale(w0, 1 / absZ0);                  // γ = w₀/|z₀|
    const branches = [], geom = [];
    for (const g of groups) {
      const zeta = g.root, m = g.multiplicity;
      const m2 = zeta.re * zeta.re + zeta.im * zeta.im;
      const zj = { re: zeta.re / m2, im: -zeta.im / m2 };
      const A = branchCoeffsFromResidues(localPrincipalResidues(num, den, zeta, m), zeta, m);
      branches.push({ z: zj, A }); geom.push({ zj, m, A });
    }
    const phi = { family: 'boundedLQD_singular', unbounded: false, w0, z0: C.clone(z0), gamma, branches };
    const poles = geom.map(g => ({ a: QD.evalPhi(g.zj, phi), m: g.m, A: g.A }));
    const Cs = solveResiduesViaProbe(QD.Family.boundedLQD_singular, phi, poles);

    // Origin residue (the (●₀) q-equation), via the family's r̃# (constant 0).
    const rZ0 = QD.LqdCommon.evalRHash(z0, phi);
    const oneOverConjZ0 = C.scale(z0, 1 / C.abs2(z0));
    const rInv = QD.LqdCommon.evalRHash(oneOverConjZ0, phi);
    const q = { re: Math.log(C.abs2(gamma)) + rZ0.re + rInv.re, im: rZ0.im - rInv.im };
    phi.q = q;   // complete the phi for the family identity verifier

    return {
      hData: { poles: poles.map((p, j) => ({ a: p.a, principal: Cs[j] })) },
      poleData: geom.map((g, j) => ({ z: g.zj, w: poles[j].a, multiplicity: g.m })),
      w0, z0: C.clone(z0), gamma, q, phi, warnings: [],
    };
  }

  // ===========================================================================
  // UNBOUNDED weighted forward kernels (∞ ∈ Ω) — Theorem 4.3.7 (Laurent-at-∞).
  // ---------------------------------------------------------------------------
  // φ : 𝔻* → Ω, φ(∞)=∞, φ'(∞)=c>0. The rational class is again the KERNEL r#,
  // now analytic on the closed EXTERIOR |z|≥1 (poles strictly INSIDE 𝔻):
  //   • PQD:  φ = z·(r#(z))^{1/α},   r#(∞) = c^α            (Family.unboundedPQD)
  //   • LQD:  φ = c·z·exp(r̃#(z)),     r#(∞) = 0 (gauge)      (Family.unboundedLQD)
  // r#'s in-disk poles are at ζ_p = 1/conj(z_j); the node-preimages z_j = 1/conj(ζ_p)
  // therefore lie OUTSIDE the disk (|z_j| > 1). A pole of r# AT z = 0 (order N_G)
  // is the Laurent-at-∞ block Σ_l G_l/z^l and is what gives h a polynomial part at
  // ∞ (degree N_G − 1); other in-disk poles ζ_p ≠ 0 are the finite-node branches.
  //
  // The returned h is { finite poles a_j = φ(z_j) } + a polynomial-at-∞ part
  // polyPart (present iff r# has a pole at 0). There is NO origin term: the author-
  // confirmed Theorem 4.3.7 forward map (Eq UPQDDirectProblemSol) reads
  //   h(w) = (1/(α·w))·Φ_φ( AnalyticIn_𝔻[ r·r# ] )(w) − t/w,
  //          t = ∫_{Ω^c} |w|^{2(α−1)} dA(w)   (over the bounded complement K),
  // but the (1/(α·w))·Φ_φ(…) term has a pole at w=0 of residue EXACTLY +t, so the
  // −t/w cancels it: the realized h is analytic at 0 (net origin residue ≡ 0 for
  // every realizable config — verified empirically to ~1e-16 across the §DF cases).
  // This is the unbounded mirror of bounded Eq 4.13 with two differences (the
  // analytic projection swaps to the disk interior — principal parts at the
  // EXTERIOR z_j — and the t/w sign flips), but in the {poles, polyPart}
  // representation shared with the inverse solver / Schwarz / sphere subsystems
  // the origin term is absent.
  //
  // IMPLEMENTATION: rather than evaluate the closed Φ_φ(…) form (whose global
  // w=0 cancellation and ∞-growth are delicate), we INVERT the inverse solver's
  // own tested forward chain — the exact codification of Theorem 4.3.7:
  //   • finite poles  — invert (★)_A (computeTargets(…).A, upper-triangular per
  //     pole) via solveResiduesViaProbe, exactly as the bounded NON-singular
  //     kernels do;
  //   • polyPart      — invert (★)_F (computeTargets(…).F, linear) via
  //     forwardPolyPartAtInfinity.
  // This round-trips against the inverse solver by construction and is validated
  // at machine precision by the family identity verifier (§DF). REALIZABILITY ⟺
  // univalence of φ (Theorem 4.3.3); a non-univalent kernel is reported as
  // not-realizable rather than returning bad data.
  // ===========================================================================

  // Parse a rational kernel that must be analytic (and, if needNonVanishing, also
  // non-vanishing) on the closed EXTERIOR |z| ≥ 1 — i.e. poles (and zeros) strictly
  // inside 𝔻. Exterior sibling of parseKernel. Returns { num, den, evalRat }.
  function parseKernelExterior(K, validateTol, needNonVanishing, label) {
    let num, den;
    if (Array.isArray(K)) { num = trimTrailingZeros(K.slice()); den = [{ re: 1, im: 0 }]; }
    else { num = trimTrailingZeros(K.num.slice()); den = trimTrailingZeros(K.den.slice()); }
    if (den.length > 1) for (const r of polynomialRoots(den)) {
      if (Math.hypot(r.re, r.im) >= 1 - validateTol)
        throw new Error(label + ': r# has a pole at |z| ≥ 1; r# must be analytic on the closed exterior |z| ≥ 1.');
    }
    if (needNonVanishing && num.length > 1) for (const r of polynomialRoots(num)) {
      if (Math.hypot(r.re, r.im) >= 1 - validateTol)
        throw new Error(label + ': r# has a zero at |z| ≥ 1; (r#)^{1/α} is not single-valued on the exterior.');
    }
    const evalRat = (z) => C.div(evalPolyAscending(num, z), evalPolyAscending(den, z));
    return { num, den, evalRat };
  }

  // Split a rational kernel r# = num/den (analytic on |z|≥1) into the inverse
  // solver's phi data. Returns { rInf, polyA, geom } where:
  //   rInf  = r#(∞) (= c^α; leading-coeff ratio; requires deg num = deg den),
  //   polyA = [G_1,…,G_{N_G}] — the principal part of r# at z=0 (the Laurent block;
  //           [] if r# has no pole at 0),
  //   geom  = [{ zeta, zj, m, A }] for each in-disk pole ζ_p ≠ 0 (m = multiplicity),
  //           with z_j = 1/conj(ζ_p) the exterior node-preimage and A its kernel
  //           Möbius-branch coefficients.
  function splitUnboundedKernel(num, den, label) {
    if ((num.length - 1) !== (den.length - 1))
      throw new Error(label + ': r# must satisfy deg(num) = deg(den) so that r#(∞) = c^α is finite and nonzero.');
    const rInf = C.div(num[num.length - 1], den[den.length - 1]);
    const groups = den.length > 1 ? groupRootsByMultiplicity(polynomialRoots(den), 1e-7) : [];
    let polyA = [];
    const geom = [];
    for (const g of groups) {
      const zeta = g.root, m = g.multiplicity;
      if (C.abs(zeta) < 1e-7) {                                   // pole of r# at z=0 ⇒ Laurent block
        polyA = localPrincipalResidues(num, den, { re: 0, im: 0 }, m);
        continue;
      }
      const m2 = zeta.re * zeta.re + zeta.im * zeta.im;
      const zj = { re: zeta.re / m2, im: -zeta.im / m2 };          // 1/conj(ζ_p), exterior
      const A = branchCoeffsFromResidues(localPrincipalResidues(num, den, zeta, m), zeta, m);
      geom.push({ zeta, zj, m, A });
    }
    return { rInf, polyA, geom };
  }

  // The hardwired constant r0Const of the family's r# (= c^α for non-sing PQD,
  // |c·z₀|^α for singular). It is NOT r#(∞): the family's branch basis u^k =
  // z^k/(1−conj(z_j)z)^k has a nonzero value at ∞ (u(∞) = −1/conj(z_j)), so to
  // make the family's r# EQUAL the input r# the constant must absorb those
  // branch-at-∞ contributions:
  //   r0Const = r#(∞) − Σ_j Σ_k conj(A_{j,k})·(−1/conj(z_j))^k.
  // (The Laurent block G_l/z^l vanishes at ∞ and does not contribute.) For a
  // valid unbounded PQD r0Const must be REAL POSITIVE (= c^α). [The LQD families
  // are immune to this because their r# − r#(∞) gauge cancels any constant.]
  function unboundedR0Const(rInf, geom) {
    let s = C.clone(rInf);
    for (const g of geom) {
      const base = C.scale(C.inv(C.conj(g.zj)), -1);          // −1/conj(z_j)
      let pw = C.clone(base);                                  // base^k, k=1…
      for (let k = 1; k <= g.A.length; k++) {
        s = C.sub(s, C.mul(C.conj(g.A[k - 1]), pw));
        pw = C.mul(pw, base);
      }
    }
    return s;
  }

  // Invert the unbounded family's (★)_F block: given phi (hence the kernel's
  // Laurent block G_l), solve the LINEAR system fam.computeTargets(phi, {polyPart}).F
  // = 0 for h's polynomial-at-∞ part polyPart (length n+1). Probes the tested
  // forward residual on unit polyPart vectors and back-substitutes — the system is
  // triangular (polyPart[m] first enters the s^{n−m} residual coefficient, with
  // diagonal c^m ≠ 0).
  function forwardPolyPartAtInfinity(fam, phi, n) {
    if (n < 0) return [];
    const N = n + 1;
    const zeroPP = Array.from({ length: N }, () => ({ re: 0, im: 0 }));
    const b = fam.computeTargets(phi, { poles: [], polyPart: zeroPP }).F;   // firstTerm (polyPart = 0)
    const cols = [];
    for (let m = 0; m < N; m++) {
      const e = zeroPP.map((_, i) => ({ re: i === m ? 1 : 0, im: 0 }));
      const Fe = fam.computeTargets(phi, { poles: [], polyPart: e }).F;
      cols.push(b.map((bi, i) => C.sub(bi, Fe[i])));                        // hTerm(e_m)
    }
    const pp = new Array(N);
    for (let i = 0; i < N; i++) {                                           // row i introduces m = n−i
      const mNew = n - i;
      let acc = C.clone(b[i]);
      for (let mm = mNew + 1; mm < N; mm++) acc = C.sub(acc, C.mul(cols[mm][i], pp[mm]));
      if (C.abs(cols[mNew][i]) < 1e-300) throw new Error('forwardPolyPartAtInfinity: singular (★)_F system');
      pp[mNew] = C.div(acc, cols[mNew][i]);
    }
    return pp;
  }

  // unboundedPowerQD: rational r#(z) + α → hData (unbounded power-weighted QD,
  // 0∉Ω, ∞∈Ω) via Theorem 4.3.7 (see section header). Returns
  //   { hData:{poles, polyPart}, poleData, c, alpha, phi, univalent, warnings }.
  function unboundedPowerQD(rHash, alpha, options) {
    options = options || {};
    const label = 'Direct.unboundedPowerQD';
    if (!(alpha > 0) || Math.abs(alpha - 1) < 1e-12) throw new Error(label + ': need α > 0, α ≠ 1.');
    const { num, den } = parseKernelExterior(rHash, options.validateTol || 1e-6, true, label);
    const { rInf, polyA, geom } = splitUnboundedKernel(num, den, label);
    const r0Const = unboundedR0Const(rInf, geom);
    if (Math.abs(r0Const.im) > 1e-6 || r0Const.re <= 0)
      throw new Error(label + ': r#(∞)−Σbranch(∞) must be real positive (= c^α) for a valid unbounded PQD. Got ' + complexFmt(r0Const) + '.');
    const c = Math.pow(r0Const.re, 1 / alpha);

    const phi = { family: 'unboundedPQD', unbounded: true, alpha, c, polyA, branches: geom.map(g => ({ z: g.zj, A: g.A })) };
    const fam = QD.Family.unboundedPQD;

    // REALIZABILITY ⟺ univalence of φ (Thm 4.3.3); otherwise report not-realizable.
    let univalent = false;
    try { univalent = QD.isBoundaryUnivalent(phi); } catch (e) { univalent = false; }
    if (!univalent) {
      return { hData: { poles: [], polyPart: [] }, poleData: [], c, alpha, phi, univalent: false,
               warnings: ['Not realizable: φ = z·(r#)^{1/α} is not univalent (boundary self-intersects). Adjust r# or α.'] };
    }
    // The NON-singular family requires 0 ∉ Ω (0 ∈ K). If this r# yields 0 ∈ Ω the
    // domain is SINGULAR — h would carry a Blaschke/origin structure the
    // non-singular verifier cannot represent — so report it (use the singular
    // kernel with a z₀ instead).
    if (QD.originInsideOmega(phi)) {
      return { hData: { poles: [], polyPart: [] }, poleData: [], c, alpha, phi, univalent: true, originInside: true,
               warnings: ['0 ∈ Ω for this r#: the domain is SINGULAR. Use unboundedPowerQDSingular (provide z₀).'] };
    }

    // Finite poles: invert (★)_A (computeTargets.A) per exterior node z_j.
    const poles = geom.map(g => ({ a: fam.evalPhi(g.zj, phi), m: g.m, A: g.A }));
    const Cs = solveResiduesViaProbe(fam, phi, poles);
    const hPoles = poles.map((p, j) => ({ a: p.a, principal: Cs[j] }));
    const poleData = geom.map((g, j) => ({ z: g.zj, w: poles[j].a, multiplicity: g.m }));

    // Polynomial-at-∞ part (degree n = N_G − 1, inferred from r#'s Laurent block):
    // invert (★)_F (computeTargets.F, linear in polyPart).
    const polyPart = forwardPolyPartAtInfinity(fam, phi, polyA.length - 1);

    return { hData: { poles: hPoles, polyPart }, poleData, c, alpha, phi, univalent: true, warnings: [] };
  }

  // unboundedPowerQDSingular: rational r#(z) + α → hData (unbounded SINGULAR PQD,
  // 0∈Ω, ∞∈Ω). φ = z·b_{z₀}·(r#)^{1/α}, with z₀ ∈ 𝔻* (|z₀|>1) the origin-preimage.
  //
  // z₀ is NOT a free parameter here (unlike the bounded singular case): the
  // "no pole at 0" structure forces the z₀-closure r(z₀)=0 (thesis Prop 4.6.3),
  // where r = reflection of r#. r(z₀)=0 ⟺ r#(1/conj z₀)=0, i.e. z₀ = 1/conj(ρ) for
  // a ZERO ρ of r# (a root of the numerator, |ρ|<1 ⇒ |z₀|>1). So z₀ is DERIVED
  // from r#'s numerator roots; options.z0 (optional) selects among them (the num
  // root nearest 1/conj(z₀_hint)). h = finite poles a_j + polyPart, NO origin term
  // (0∈Ω ⇒ w=0 is not a legal pole site). Returns
  //   { hData:{poles, polyPart}, poleData, c, alpha, z0, phi, univalent, warnings }.
  function unboundedPowerQDSingular(rHash, alpha, options) {
    options = options || {};
    const label = 'Direct.unboundedPowerQDSingular';
    if (!(alpha > 0) || Math.abs(alpha - 1) < 1e-12) throw new Error(label + ': need α > 0, α ≠ 1.');
    const { num, den } = parseKernelExterior(rHash, options.validateTol || 1e-6, true, label);
    if (num.length <= 1)
      throw new Error(label + ': r# must have a zero inside 𝔻 (non-constant numerator) to place the origin-preimage z₀.');
    const { rInf, polyA, geom } = splitUnboundedKernel(num, den, label);
    const r0Const = unboundedR0Const(rInf, geom);             // = |c·z₀|^α
    if (Math.abs(r0Const.im) > 1e-6 || r0Const.re <= 0)
      throw new Error(label + ': r#(∞)−Σbranch(∞) must be real positive (= |c·z₀|^α). Got ' + complexFmt(r0Const) + '.');

    // z₀ = 1/conj(ρ), ρ a zero of r# inside 𝔻 (root of num). Order candidates by
    // proximity to 1/conj(z₀_hint) when a hint is given.
    const inv = (p) => C.scale(C.conj(p), 1 / C.abs2(p));     // 1/conj(p) = p/|p|²
    let roots = polynomialRoots(num).filter(r => C.abs(r) < 1 - 1e-9);
    if (roots.length === 0) throw new Error(label + ': r# has no zero strictly inside 𝔻; cannot place z₀.');
    if (options.z0) { const tgt = inv(options.z0); roots = roots.slice().sort((p, q) => C.abs(C.sub(p, tgt)) - C.abs(C.sub(q, tgt))); }

    const fam = QD.Family.unboundedPQD_singular;
    let chosen = null, lastWarn = '';
    for (const rho of roots) {
      const z0 = inv(rho);                                    // |z₀| > 1
      const c = Math.pow(r0Const.re, 1 / alpha) / C.abs(z0);  // r0Const = |c·z₀|^α
      const phi = { family: 'unboundedPQD_singular', unbounded: true, alpha, c, z0,
                    polyA, branches: geom.map(g => ({ z: g.zj, A: g.A })) };
      let univalent = false;
      try { univalent = QD.isBoundaryUnivalent(phi); } catch (e) { univalent = false; }
      if (!univalent) { lastWarn = 'φ not univalent'; continue; }
      let originIn = false;
      try { originIn = QD.originInsideOmega(phi); } catch (e) { originIn = false; }
      if (!originIn) { lastWarn = '0 ∉ Ω for this z₀'; continue; }
      chosen = { phi, c, z0 }; break;
    }
    if (!chosen) {
      return { hData: { poles: [], polyPart: [] }, poleData: [], alpha, univalent: false,
               warnings: ['Not realizable: no zero of r# yields a univalent singular φ with 0∈Ω (' + lastWarn + ').'] };
    }
    const { phi, c, z0 } = chosen;

    // Finite poles via (★)_A inversion; polyPart via (★)_F inversion. No origin term.
    const poles = geom.map(g => ({ a: fam.evalPhi(g.zj, phi), m: g.m, A: g.A }));
    const Cs = solveResiduesViaProbe(fam, phi, poles);
    const hPoles = poles.map((p, j) => ({ a: p.a, principal: Cs[j] }));
    const poleData = geom.map((g, j) => ({ z: g.zj, w: poles[j].a, multiplicity: g.m }));
    const polyPart = forwardPolyPartAtInfinity(fam, phi, polyA.length - 1);

    return { hData: { poles: hPoles, polyPart }, poleData, c, alpha, z0, phi, univalent: true, warnings: [] };
  }

  // ---------------------------------------------------------------------------
  // UNBOUNDED LOG-weighted forward kernels. φ = c·z·exp(r̃# + B(1/z)) (non-sing) /
  // c·|z₀|·z·b_{z₀}·exp(r̃#) (sing), r̃# = r# − r#(∞) (gauge ∞-absorption). The
  // input EXPONENT kernel is rational with poles INSIDE 𝔻: a pole at z=0 (order
  // N_β) is the polynomial-h block B(1/z)=Σ β_l/z^l; poles ζ_p≠0 are the r#
  // branches. c (= φ'(∞)) is a free real-positive input (NOT derived). Finite
  // poles come from inverting (★)_A (computeTargets.A) by probe; the poly-h
  // block, when present, from the coupled (★)_F linear system (β = target_F(P̃),
  // P̃ = [Σ C_{j,1}, polyPart…]). NB the LQD-SINGULAR family carries an ORIGIN
  // pole q/w (the author's q-formula; unlike the PQD-singular family which has
  // none) — returned separately as `q`.
  // ---------------------------------------------------------------------------

  // Dense complex linear solve M·x = b (Gaussian elimination, partial pivot).
  // n small (poly-h degree); used for the LQD coupled (★)_F block.
  function solveComplexLinear(M, b, n) {
    const A = M.map(row => row.map(c => C.clone(c)));
    const x = b.map(c => C.clone(c));
    for (let col = 0; col < n; col++) {
      let piv = col, best = C.abs(A[col][col]);
      for (let r = col + 1; r < n; r++) { const v = C.abs(A[r][col]); if (v > best) { best = v; piv = r; } }
      if (best < 1e-300) throw new Error('solveComplexLinear: singular system');
      if (piv !== col) { const tr = A[col]; A[col] = A[piv]; A[piv] = tr; const tx = x[col]; x[col] = x[piv]; x[piv] = tx; }
      const d = A[col][col];
      for (let r = col + 1; r < n; r++) {
        const f = C.div(A[r][col], d);
        for (let cc = col; cc < n; cc++) A[r][cc] = C.sub(A[r][cc], C.mul(f, A[col][cc]));
        x[r] = C.sub(x[r], C.mul(f, x[col]));
      }
    }
    const out = new Array(n);
    for (let i = n - 1; i >= 0; i--) {
      let acc = C.clone(x[i]);
      for (let j = i + 1; j < n; j++) acc = C.sub(acc, C.mul(A[i][j], out[j]));
      out[i] = C.div(acc, A[i][i]);
    }
    return out;
  }

  // Solve h's polyPart (length N) from the LQD (★)_F block, given the known β
  // (phi.lqdBeta) and the finite poles (hPoles, whose Σ C_{j,1} feeds P̃[0]).
  // β = target_F(P̃) is linear in polyPart; probe + dense solve.
  function forwardBetaToPolyPart(fam, phi, hPoles, N) {
    if (N <= 0) return [];
    const zeroPP = Array.from({ length: N }, () => ({ re: 0, im: 0 }));
    const base = fam.computeTargets(phi, { poles: hPoles, polyPart: zeroPP }).F;   // target_F at polyPart=0
    const cols = [];
    for (let m = 0; m < N; m++) {
      const e = zeroPP.map((_, i) => ({ re: i === m ? 1 : 0, im: 0 }));
      const fm = fam.computeTargets(phi, { poles: hPoles, polyPart: e }).F;
      cols.push(fm.map((v, i) => C.sub(v, base[i])));                              // ∂ target_F / ∂ polyPart[m]
    }
    const rows = [];
    for (let i = 0; i < N; i++) rows.push(cols.map(col => col[i]));
    const rhs = phi.lqdBeta.map((b, i) => C.sub(b, base[i]));
    return solveComplexLinear(rows, rhs, N);
  }

  // Shared front-end for the unbounded LQD kernels: parse the exponent kernel
  // (poles inside 𝔻), split into β (pole at 0) + branches (poles ζ_p ≠ 0).
  // Returns { num, den, lqdBeta, geom } (geom: [{zeta, zj, m, A}]).
  function splitUnboundedLogKernel(rHash, validateTol, label) {
    const { num, den } = parseKernelExterior(rHash, validateTol, false, label);
    const groups = den.length > 1 ? groupRootsByMultiplicity(polynomialRoots(den), 1e-7) : [];
    let lqdBeta = [];
    const geom = [];
    for (const g of groups) {
      const zeta = g.root, m = g.multiplicity;
      if (C.abs(zeta) < 1e-7) { lqdBeta = localPrincipalResidues(num, den, { re: 0, im: 0 }, m); continue; }
      const m2 = zeta.re * zeta.re + zeta.im * zeta.im;
      const zj = { re: zeta.re / m2, im: -zeta.im / m2 };
      const A = branchCoeffsFromResidues(localPrincipalResidues(num, den, zeta, m), zeta, m);
      geom.push({ zeta, zj, m, A });
    }
    return { num, den, lqdBeta, geom };
  }

  // unboundedLogQD: rational exponent r#(z) + c → hData (unbounded log-weighted
  // QD, 0∉Ω, ∞∈Ω). φ = c·z·exp(r̃# + B(1/z)). Returns
  //   { hData:{poles, polyPart}, poleData, c, phi, univalent, warnings }.
  function unboundedLogQD(rHash, c, options) {
    options = options || {};
    const label = 'Direct.unboundedLogQD';
    if (!(c > 0)) throw new Error(label + ": c (= φ'(∞)) must be a positive real.");
    const { lqdBeta, geom } = splitUnboundedLogKernel(rHash, options.validateTol || 1e-6, label);

    const phi = { family: 'unboundedLQD', unbounded: true, c, lqdBeta, branches: geom.map(g => ({ z: g.zj, A: g.A })) };
    const fam = QD.Family.unboundedLQD;

    let univalent = false;
    try { univalent = QD.isBoundaryUnivalent(phi); } catch (e) { univalent = false; }
    if (!univalent) {
      return { hData: { poles: [], polyPart: [] }, poleData: [], c, phi, univalent: false,
               warnings: ['Not realizable: φ = c·z·exp(r#) is not univalent (boundary self-intersects). Adjust r# or c.'] };
    }
    if (QD.originInsideOmega(phi)) {
      return { hData: { poles: [], polyPart: [] }, poleData: [], c, phi, univalent: true, originInside: true,
               warnings: ['0 ∈ Ω for this r#: the domain is SINGULAR. Use unboundedLogQDSingular.'] };
    }

    const poles = geom.map(g => ({ a: fam.evalPhi(g.zj, phi), m: g.m, A: g.A }));
    const Cs = solveResiduesViaProbe(fam, phi, poles);
    const hPoles = poles.map((p, j) => ({ a: p.a, principal: Cs[j] }));
    const poleData = geom.map((g, j) => ({ z: g.zj, w: poles[j].a, multiplicity: g.m }));
    const polyPart = forwardBetaToPolyPart(fam, phi, hPoles, lqdBeta.length);

    return { hData: { poles: hPoles, polyPart }, poleData, c, phi, univalent: true, warnings: [] };
  }

  // unboundedLogQDSingular: rational exponent r#(z) + c + z₀ → hData (unbounded
  // SINGULAR LQD, 0∈Ω, ∞∈Ω). φ = c·|z₀|·z·b_{z₀}·exp(r̃# + B(1/z)). UNLIKE the
  // PQD-singular family, h carries an ORIGIN pole q/w with (author's q-formula)
  //   q = ln(c²|z₀|²) + r̃#(z₀) + conj(r̃#(1/conj z₀)) + B(1/z₀) + conj(B(conj z₀)),
  // returned separately (the inverse solver takes it via opts.q). z₀ ∈ 𝔻* is a
  // FREE input here (the q-equation absorbs the origin consistency, as in the
  // bounded LQD-singular case). Returns
  //   { hData:{poles, polyPart}, poleData, c, z0, q, phi, univalent, warnings }.
  function unboundedLogQDSingular(rHash, c, z0, options) {
    options = options || {};
    const label = 'Direct.unboundedLogQDSingular';
    if (!(c > 0)) throw new Error(label + ": c (= φ'(∞)) must be a positive real.");
    if (!z0 || C.abs(z0) <= 1 + 1e-9) throw new Error(label + ': z₀ must be exterior (|z₀| > 1) — the origin-preimage in 𝔻*.');
    const { lqdBeta, geom } = splitUnboundedLogKernel(rHash, options.validateTol || 1e-6, label);

    const phi = { family: 'unboundedLQD_singular', unbounded: true, c, z0: C.clone(z0), lqdBeta,
                  branches: geom.map(g => ({ z: g.zj, A: g.A })) };
    const fam = QD.Family.unboundedLQD_singular;

    let univalent = false;
    try { univalent = QD.isBoundaryUnivalent(phi); } catch (e) { univalent = false; }
    if (!univalent) {
      return { hData: { poles: [], polyPart: [] }, poleData: [], c, z0: C.clone(z0), phi, univalent: false,
               warnings: ['Not realizable: φ = c·|z₀|·z·b_{z₀}·exp(r#) is not univalent. Adjust r#, c, or z₀.'] };
    }

    const poles = geom.map(g => ({ a: fam.evalPhi(g.zj, phi), m: g.m, A: g.A }));
    const Cs = solveResiduesViaProbe(fam, phi, poles);
    const hPoles = poles.map((p, j) => ({ a: p.a, principal: Cs[j] }));
    const poleData = geom.map((g, j) => ({ z: g.zj, w: poles[j].a, multiplicity: g.m }));
    const polyPart = forwardBetaToPolyPart(fam, phi, hPoles, lqdBeta.length);

    // Origin residue q (the (●₀) q-equation; author's full formula):
    //   q = ln(c²|z₀|²) + R(z₀) + conj(R(1/conj z₀)),   R(z) = r̃#(z) + B(1/z),
    //   r̃# = r# − r#(∞). Note R(1/conj z₀) = r̃#(1/conj z₀) + B(conj z₀), and
    //   B(conj z₀) = B(1/z)|_{z = 1/conj z₀} (since B is a function of 1/z).
    const rInf = QD.LqdCommon.rHashAtInfinity(phi);
    const Bz = (z) => QD.LqdCommon.evalB_OverZ(phi, z);                  // B(1/z)
    const Rfull = (z) => C.add(C.sub(QD.LqdCommon.evalRHash(z, phi), rInf), Bz(z));
    const oneOverConjZ0 = C.scale(z0, 1 / C.abs2(z0));                   // 1/conj(z₀) = z₀/|z₀|²
    const q = C.add(C.add({ re: Math.log(c * c * C.abs2(z0)), im: 0 }, Rfull(z0)),
                    C.conj(Rfull(oneOverConjZ0)));
    phi.q = q;

    return { hData: { poles: hPoles, polyPart }, poleData, c, z0: C.clone(z0), q, phi, univalent: true, warnings: [] };
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
  Direct.boundedPowerQD     = boundedPowerQD;
  Direct.boundedLogQD       = boundedLogQD;
  Direct.boundedPowerQDSingular = boundedPowerQDSingular;
  Direct.boundedLogQDSingular   = boundedLogQDSingular;
  Direct.unboundedPowerQD       = unboundedPowerQD;
  Direct.unboundedPowerQDSingular = unboundedPowerQDSingular;
  Direct.unboundedLogQD         = unboundedLogQD;
  Direct.unboundedLogQDSingular = unboundedLogQDSingular;
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
