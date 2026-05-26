// =============================================================================
// complex.js -- Complex number arithmetic
//
// Complex numbers are plain objects {re: number, im: number}. We avoid math.js
// in the inner solver loop for performance, but the data format is identical
// to math.complex(re, im) so it interops trivially.
// =============================================================================

const Complex = {
  ZERO: () => ({re: 0, im: 0}),
  ONE:  () => ({re: 1, im: 0}),
  I:    () => ({re: 0, im: 1}),

  c(re, im = 0) { return {re, im}; },
  clone(a) { return {re: a.re, im: a.im}; },

  add(a, b) { return {re: a.re + b.re, im: a.im + b.im}; },
  sub(a, b) { return {re: a.re - b.re, im: a.im - b.im}; },
  neg(a)    { return {re: -a.re, im: -a.im}; },
  mul(a, b) { return {re: a.re*b.re - a.im*b.im, im: a.re*b.im + a.im*b.re}; },
  scale(a, s) { return {re: a.re*s, im: a.im*s}; },
  conj(a)   { return {re: a.re, im: -a.im}; },

  // In-place arithmetic variants. Write the result into a caller-supplied
  // `out` object instead of allocating a fresh {re,im}. Use in tight inner
  // loops (per-pixel solver, branch sums) to remove allocator + GC pressure.
  // SAFE TO ALIAS: `out` may be the same object as `a` or `b`.
  //
  // The functional variants above stay for readability; these are an
  // additive perf-only API — callers opt in.
  mulInto(a, b, out) {
    const re = a.re*b.re - a.im*b.im;
    out.im   = a.re*b.im + a.im*b.re;
    out.re   = re;
    return out;
  },
  addInto(a, b, out) { out.re = a.re + b.re; out.im = a.im + b.im; return out; },
  subInto(a, b, out) { out.re = a.re - b.re; out.im = a.im - b.im; return out; },
  scaleInto(a, s, out) { out.re = a.re*s; out.im = a.im*s; return out; },
  // Accumulator: out += a*b (used for residual / branch-sum kernels).
  addMulInto(a, b, out) {
    out.re += a.re*b.re - a.im*b.im;
    out.im += a.re*b.im + a.im*b.re;
    return out;
  },

  inv(a) {
    const d = a.re*a.re + a.im*a.im;
    if (d === 0) throw new Error("Complex.inv: division by zero");
    return {re: a.re/d, im: -a.im/d};
  },
  div(a, b) {
    const d = b.re*b.re + b.im*b.im;
    if (d === 0) throw new Error("Complex.div: division by zero");
    return {re: (a.re*b.re + a.im*b.im)/d, im: (a.im*b.re - a.re*b.im)/d};
  },

  abs(a)  { return Math.hypot(a.re, a.im); },
  abs2(a) { return a.re*a.re + a.im*a.im; },
  arg(a)  { return Math.atan2(a.im, a.re); },

  // Integer power (positive, zero, or negative)
  pow(a, n) {
    if (n === 0) return {re: 1, im: 0};
    if (n < 0)   return Complex.inv(Complex.pow(a, -n));
    let result = {re: 1, im: 0};
    let base = {re: a.re, im: a.im};
    let e = n;
    while (e > 0) {
      if (e & 1) result = Complex.mul(result, base);
      base = Complex.mul(base, base);
      e >>= 1;
    }
    return result;
  },

  eq(a, b, tol = 1e-12) {
    return Math.abs(a.re - b.re) < tol && Math.abs(a.im - b.im) < tol;
  },

  // Parse a complex number from a string.
  // Accepts: "1+2i", "1-2i", "3", "2i", "-i", "+i", "1.5e-3+2.1e2i", etc.
  parse(str) {
    if (typeof str === 'number') return {re: str, im: 0};
    if (typeof str === 'object' && str !== null && 're' in str) return Complex.clone(str);
    if (typeof str !== 'string') return null;

    let s = str.trim().replace(/\s+/g, '').replace(/I/g, 'i').replace(/\*/g, '');
    if (s === '') return null;

    if (s === 'i'  || s === '+i') return {re: 0, im: 1};
    if (s === '-i')               return {re: 0, im: -1};

    // Tokenize: split on + or - that's not at position 0 and not after e/E
    const tokens = [];
    let current = '';
    let i0 = 0;
    if (s[0] === '+' || s[0] === '-') { current = s[0]; i0 = 1; }
    for (let i = i0; i < s.length; i++) {
      const c = s[i];
      const prev = i > 0 ? s[i-1] : '';
      if ((c === '+' || c === '-') && prev !== 'e' && prev !== 'E') {
        if (current) tokens.push(current);
        current = c;
      } else {
        current += c;
      }
    }
    if (current) tokens.push(current);

    let re = 0, im = 0;
    for (const t of tokens) {
      if (t.endsWith('i')) {
        const numPart = t.slice(0, -1);
        let v;
        if (numPart === '' || numPart === '+') v = 1;
        else if (numPart === '-') v = -1;
        else { v = parseFloat(numPart); if (isNaN(v)) return null; }
        im += v;
      } else {
        const v = parseFloat(t);
        if (isNaN(v)) return null;
        re += v;
      }
    }
    return {re, im};
  },

  // toString: legacy fixed-decimal formatter, used by ui.js call sites that
  // need a specific decimal-place count (e.g. canvas overlay labels).
  // For human-readable expression-style output (KaTeX / paste-roundtrip /
  // display cards) use Complex.format instead.
  toString(a, digits = 4) {
    const r = a.re, i = a.im;
    if (Math.abs(i) < 1e-12) return Number(r.toFixed(digits)).toString();
    if (Math.abs(r) < 1e-12) {
      if (Math.abs(i - 1) < 1e-12) return 'i';
      if (Math.abs(i + 1) < 1e-12) return '-i';
      return Number(i.toFixed(digits)).toString() + 'i';
    }
    const rStr = Number(r.toFixed(digits)).toString();
    const sign = i >= 0 ? '+' : '-';
    const iAbs = Math.abs(i);
    const iStr = Math.abs(iAbs - 1) < 1e-12 ? '' : Number(iAbs.toFixed(digits)).toString();
    return rStr + sign + iStr + 'i';
  },

  // format: unified expression-style complex formatter.
  //
  //   opts.digits   significant figures, default 6 (toPrecision-style)
  //   opts.tol      tolerance for "snap to zero" and "snap to integer",
  //                 default 1e-12
  //
  // Snaps any component within `tol` of an integer to that integer, so
  // numerical drift like 0.99999999999 renders as "1" cleanly. Handles the
  // standard short forms: 0, ±i, ±1±i, pure real, pure imaginary.
  //
  // Replaces direct-common.js's formatComplex and direct-ui.js's
  // coeffToString / complexToString / complexToKatex (which had four
  // near-identical bodies).
  format(a, opts) {
    if (!a) return '0';
    opts = opts || {};
    const digits = opts.digits ?? 6;
    const tol    = opts.tol    ?? 1e-12;

    // Snap a real to the nearest integer if within tol; else within tol of
    // zero, return 0; else pass through.
    const snap = (x) => {
      if (!isFinite(x)) return x;
      if (Math.abs(x) < tol) return 0;
      const r = Math.round(x);
      if (Math.abs(x - r) < tol) return r;
      return x;
    };
    // Format a real (after snap): integers as-is, others via toPrecision.
    const fmt = (x) => {
      if (!isFinite(x)) return String(x);
      if (Number.isInteger(x)) return String(x);
      return Number(x.toPrecision(digits)).toString();
    };

    const r = snap(a.re), i = snap(a.im);
    if (i === 0) return fmt(r);
    if (r === 0) {
      if (i ===  1) return 'i';
      if (i === -1) return '-i';
      return fmt(i) + 'i';
    }
    if (i ===  1) return fmt(r) + '+i';
    if (i === -1) return fmt(r) + '-i';
    const sign = i < 0 ? '' : '+';   // fmt(i) carries its own '-' when i < 0
    return fmt(r) + sign + fmt(i) + 'i';
  }
};

if (typeof module !== 'undefined' && module.exports) module.exports = Complex;
