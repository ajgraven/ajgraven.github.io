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
  }
};

if (typeof module !== 'undefined' && module.exports) module.exports = Complex;
