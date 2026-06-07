// =============================================================================
// faber-analysis.js  —  Faber polynomials of the complement of a UQD.
//
// For a classical UNBOUNDED quadrature domain (family 'unboundedQD'), the
// solved map φ is the EXTERIOR conformal map  {|z|>1} → Ω,  φ(∞)=∞, φ'(∞)=c>0.
// That makes φ exactly the exterior map of the bounded complement K = ℂ\Ω, so
// the Faber polynomials of K are read straight off φ's Laurent expansion at ∞:
//
//     φ(z) = c·z + c₀ + c₁/z + c₂/z² + …            (QD.phiLaurentAtInfinity)
//
// The Faber polynomials F_n(ζ) (degree n in the IMAGE-plane variable ζ — the
// same plane where ∂Ω lives) satisfy the three-term-with-history recurrence:
//
//     F₀(ζ) = 1
//     F₁(ζ) = (ζ − c₀)/c
//     c·F_{n+1}(ζ) = (ζ − c₀)·F_n(ζ) − Σ_{k=1}^{n} c_k·F_{n−k}(ζ) − n·c_n
//
// (Derived from ψ'(z)/(ψ(z)−ζ) = Σ F_n(ζ) z^{−n−1}; verified against the disk
// c_k=0 ⇒ F_n=ζ^n and the Joukowski/interval c₁=1 ⇒ F_n=2·T_n(ζ/2), Chebyshev.)
//
// Faber roots cluster in/around K (the bounded "hole" of the unbounded domain).
//
// SCOPE: classical 'unboundedQD' only. The PQD/LQD families carry an extra
// power/Blaschke weight, so their φ-Laurent is NOT the plain exterior-map
// expansion and this clean identity does not apply — callers must gate on the
// family tag (the UI does).
//
// Pure module (no DOM): loads cleanly in node-test.js. Depends on QD.Poly
// (poly-helpers.js), global Complex (complex.js), and QD.phiLaurentAtInfinity
// (solver-uqd.js) — all present in the worker bundle / shared QD namespace.
//
// API (QD.FaberAnalysis):
//   faberPolynomials(phi, N) → { c, c0, coeffs:[F₀..F_N] }   (ascending Complex[])
//   faberPolynomial(phi, n)  → Complex[]                     (single F_n)
//   polynomialRoots(coeffsAsc, opts) → { roots:[Complex], converged, iterations, degree }
//   formatFaberPoly(Fn, opts) → string                       (readable ζ-expression)
//   faberConvergence(phi, N) → [{ n, converged, residual, roots:[Complex] }]
// =============================================================================

(function (global) {
  'use strict';

  // QD-namespace resolution — same idiom every solver/analysis file uses.
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD) return;

  // Complex is a global in the browser/worker; in the node-test vm it is also
  // installed on the shared context. Fall back to QD.Complex if needed.
  const C = (typeof Complex !== 'undefined') ? Complex : QD.Complex;

  // ---------------------------------------------------------------------------
  // Faber polynomial coefficient lists F₀..F_N (ascending-power Complex[]).
  // ---------------------------------------------------------------------------
  function faberPolynomials(phi, N) {
    if (!phi || !phi.unbounded) {
      throw new Error('faberPolynomials: requires an unbounded conformal map');
    }
    const c = phi.c;
    if (!(typeof c === 'number' && c > 0 && isFinite(c))) {
      throw new Error('faberPolynomials: capacity c = φ\'(∞) must be a positive finite number');
    }
    N = Math.max(0, Math.floor(N || 0));

    // Laurent coeffs c₀..c_{N−1}. For F_{n+1} (n = N−1 max) we need c_n = c_{N−1}.
    const lc = QD.phiLaurentAtInfinity(phi, Math.max(1, N));
    const at = (k) => (k < lc.length && lc[k]) ? lc[k] : { re: 0, im: 0 };
    const c0 = at(0);
    const invC = { re: 1 / c, im: 0 };

    // ζ − c₀  as ascending [ −c₀ , 1 ].
    const zMinusC0 = [{ re: -c0.re, im: -c0.im }, { re: 1, im: 0 }];

    const coeffs = [];
    coeffs[0] = [{ re: 1, im: 0 }];                              // F₀ = 1
    if (N >= 1) {
      coeffs[1] = QD.Poly.scale(zMinusC0.slice(), invC);         // F₁ = (ζ − c₀)/c
    }
    for (let n = 1; n < N; n++) {
      // c·F_{n+1} = (ζ − c₀)·F_n − Σ_{k=1}^{n} c_k·F_{n−k} − n·c_n
      const term1 = QD.Poly.mul(zMinusC0, coeffs[n]);
      let sum = QD.Poly.zero();
      for (let k = 1; k <= n; k++) {
        sum = QD.Poly.add(sum, QD.Poly.scale(coeffs[n - k], at(k)));
      }
      let next = QD.Poly.add(term1, QD.Poly.neg(sum));
      const cn = at(n);
      next[0] = { re: next[0].re - n * cn.re, im: next[0].im - n * cn.im };
      coeffs[n + 1] = QD.Poly.scale(next, invC);
    }
    return { c, c0, coeffs };
  }

  function faberPolynomial(phi, n) {
    return faberPolynomials(phi, n).coeffs[n];
  }

  // ---------------------------------------------------------------------------
  // Complex polynomial root-finder: Durand–Kerner (Weierstrass) + Newton polish.
  //
  // coeffs: ascending-power Complex[] (index i = coeff of ζ^i).
  // opts:   { maxIter=200, tol=1e-12, polish=true }
  //
  // Monomial-basis root-finding is ill-conditioned at high degree, so on a
  // failure to converge we return converged:false (callers surface a warning)
  // rather than emitting garbage silently.
  // ---------------------------------------------------------------------------
  function polynomialRoots(coeffs, opts) {
    opts = opts || {};
    const tol = opts.tol != null ? opts.tol : 1e-12;
    const maxIter = opts.maxIter != null ? opts.maxIter : 200;
    const polish = opts.polish !== false;

    // Strip trailing (highest-degree) near-zero coefficients to find true degree.
    const a = coeffs ? coeffs.slice() : [];
    while (a.length > 1 && Math.hypot(a[a.length - 1].re, a[a.length - 1].im) < 1e-14) a.pop();
    const d = a.length - 1;
    if (d <= 0) return { roots: [], converged: true, iterations: 0, degree: Math.max(0, d) };

    // Monic normalization (improves conditioning of the DK iteration).
    const lead = a[d];
    const mon = a.map(co => C.div(co, lead));                    // ascending; mon[d] = 1

    // Cauchy root bound  R = 1 + max_{k<d} |a_k|  (monic).
    let maxAbs = 0;
    for (let k = 0; k < d; k++) {
      const r = Math.hypot(mon[k].re, mon[k].im);
      if (r > maxAbs) maxAbs = r;
    }
    const R = 1 + maxAbs;

    // Horner evaluation of the monic polynomial at a point.
    const evalP = (pt) => {
      let acc = { re: mon[d].re, im: mon[d].im };
      for (let k = d - 1; k >= 0; k--) acc = C.add(C.mul(acc, pt), mon[k]);
      return acc;
    };
    // Horner evaluation of P'(pt) = Σ_{k≥1} k·mon[k]·ζ^{k−1}.
    const evalDP = (pt) => {
      let der = { re: 0, im: 0 };
      for (let k = d; k >= 1; k--) der = C.add(C.mul(der, pt), C.scale(mon[k], k));
      return der;
    };

    // Initialize on a circle of radius R, angle 2πj/d + 0.4. The phase offset
    // breaks symmetry so no iterate lands exactly on a real root of a highly
    // symmetric polynomial (ζ^n, Chebyshev), which would stall the iteration.
    let z = new Array(d);
    for (let j = 0; j < d; j++) {
      const ang = 2 * Math.PI * j / d + 0.4;
      z[j] = { re: R * Math.cos(ang), im: R * Math.sin(ang) };
    }

    let iter = 0, converged = false;
    for (; iter < maxIter; iter++) {
      let maxDelta = 0;
      const zNew = new Array(d);
      for (let j = 0; j < d; j++) {
        const pj = evalP(z[j]);
        let denom = { re: 1, im: 0 };
        for (let k = 0; k < d; k++) {
          if (k === j) continue;
          denom = C.mul(denom, C.sub(z[j], z[k]));
        }
        let delta;
        if (C.abs2(denom) < 1e-300) delta = { re: 0, im: 0 };
        else delta = C.div(pj, denom);
        zNew[j] = C.sub(z[j], delta);
        const dm = Math.hypot(delta.re, delta.im);
        if (dm > maxDelta) maxDelta = dm;
      }
      z = zNew;
      if (maxDelta < tol) { converged = true; iter++; break; }
    }

    if (polish) {
      for (let j = 0; j < d; j++) {
        for (let s = 0; s < 8; s++) {
          const pj = evalP(z[j]);
          const der = evalDP(z[j]);
          if (C.abs2(der) < 1e-300) break;
          const step = C.div(pj, der);
          z[j] = C.sub(z[j], step);
          if (Math.hypot(step.re, step.im) < 1e-15) break;
        }
      }
    }

    return { roots: z, converged, iterations: iter, degree: d };
  }

  // ---------------------------------------------------------------------------
  // Human-readable expression for a Faber polynomial (descending powers).
  //   opts: { varSym='ζ', digits=4, tol=1e-9 }
  // ---------------------------------------------------------------------------
  function formatFaberPoly(Fn, opts) {
    opts = opts || {};
    const v = opts.varSym || 'ζ';
    const digits = opts.digits != null ? opts.digits : 4;
    const tol = opts.tol != null ? opts.tol : 1e-9;
    if (!Fn || !Fn.length) return '0';

    // Exponents render as Unicode superscripts (ζ², not ζ^2) so the formula
    // matches the coefficient table in the UI.
    const supMap = '⁰¹²³⁴⁵⁶⁷⁸⁹';
    const sup = (k) => String(k).split('').map(d => supMap[+d] || d).join('');

    let out = '';
    let any = false;
    for (let k = Fn.length - 1; k >= 0; k--) {
      const co = Fn[k] || { re: 0, im: 0 };
      if (Math.hypot(co.re, co.im) < tol) continue;
      const powStr = k === 0 ? '' : (k === 1 ? v : v + sup(k));
      const realOnly = Math.abs(co.im) < tol;
      let sign, body;
      if (realOnly) {
        let r = co.re;
        const rr = Math.round(r);
        if (Math.abs(r - rr) < tol) r = rr;
        sign = r < 0 ? '−' : '+';
        const aval = Math.abs(r);
        if (powStr && Math.abs(aval - 1) < tol) {
          body = powStr;                                         // ±1·ζ^k → ζ^k
        } else {
          body = Number(aval.toPrecision(digits)).toString() + powStr;
        }
      } else {
        sign = '+';
        body = '(' + C.format(co, { digits, tol }) + ')' + powStr;
      }
      if (!any) { out = (sign === '−' ? '−' : '') + body; any = true; }
      else { out += ' ' + sign + ' ' + body; }
    }
    return any ? out : '0';
  }

  // ---------------------------------------------------------------------------
  // Workhorse for the UI: build F₁..F_N, then root-find each order. Returns
  // per-order { n, converged, residual (max |F_n(root)|), roots }.
  // ---------------------------------------------------------------------------
  function faberConvergence(phi, N) {
    const { coeffs } = faberPolynomials(phi, N);
    const out = [];
    for (let n = 1; n <= N; n++) {
      const Fn = coeffs[n];
      const r = polynomialRoots(Fn);
      let maxRes = 0;
      for (const root of r.roots) {
        let acc = { re: 0, im: 0 };
        for (let k = Fn.length - 1; k >= 0; k--) acc = C.add(C.mul(acc, root), Fn[k]);
        const m = Math.hypot(acc.re, acc.im);
        if (m > maxRes) maxRes = m;
      }
      out.push({ n, converged: r.converged, residual: maxRes, roots: r.roots });
    }
    return out;
  }

  QD.FaberAnalysis = {
    faberPolynomials,
    faberPolynomial,
    polynomialRoots,
    formatFaberPoly,
    faberConvergence,
  };

})(typeof globalThis !== 'undefined' ? globalThis : this);
