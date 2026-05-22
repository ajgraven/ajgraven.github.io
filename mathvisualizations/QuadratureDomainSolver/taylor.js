// =============================================================================
// taylor.js -- Truncated Taylor series arithmetic over Complex
//
// A Taylor series is an array of Complex; t[i] is the coefficient of (x-x0)^i.
// All operations truncate to a specified length.
// =============================================================================

const Taylor = {
  // Zero series of length n
  zero(n) {
    const r = new Array(n);
    for (let i = 0; i < n; i++) r[i] = {re: 0, im: 0};
    return r;
  },

  // Constant series c, padded to length n
  constant(c, n) {
    const r = Taylor.zero(n);
    r[0] = Complex.clone(c);
    return r;
  },

  // Add two series (output length = max of inputs)
  add(p, q) {
    const n = Math.max(p.length, q.length);
    const r = Taylor.zero(n);
    for (let i = 0; i < n; i++) {
      if (i < p.length) r[i] = Complex.add(r[i], p[i]);
      if (i < q.length) r[i] = Complex.add(r[i], q[i]);
    }
    return r;
  },

  // p - q
  sub(p, q) {
    return Taylor.add(p, q.map(Complex.neg));
  },

  // Multiply, truncated to length L+1 (default: full product)
  mul(p, q, L) {
    if (L === undefined) L = (p.length - 1) + (q.length - 1);
    const r = Taylor.zero(L + 1);
    for (let i = 0; i <= L && i < p.length; i++) {
      for (let j = 0; j <= L - i && j < q.length; j++) {
        r[i + j] = Complex.add(r[i + j], Complex.mul(p[i], q[j]));
      }
    }
    return r;
  },

  // Multiply each coefficient by a Complex scalar
  scaleComplex(p, c) {
    return p.map(x => Complex.mul(x, c));
  },

  // p^k, truncated to L+1
  pow(p, k, L) {
    if (k === 0) {
      const r = Taylor.zero(L + 1);
      r[0] = {re: 1, im: 0};
      return r;
    }
    let result = Taylor.truncate(p, L);
    for (let i = 1; i < k; i++) {
      result = Taylor.mul(result, p, L);
    }
    return result;
  },

  // Truncate or pad to length L+1
  truncate(p, L) {
    const r = Taylor.zero(L + 1);
    for (let i = 0; i <= L && i < p.length; i++) r[i] = Complex.clone(p[i]);
    return r;
  },

  // -------------------------------------------------------------------------
  // Compositional inverse:
  //   given p = [0, c_1, c_2, ...] with c_1 != 0,
  //   return q = [0, d_1, d_2, ...] of length L+1 such that p(q(t)) = t.
  //
  // Key observation: in [t^n] p(q), the coefficient q_n appears ONLY via the
  // k=1 term (since q_0 = 0 forces every monomial in q^k of degree n with k>=2
  // to use indices that all lie in {1,...,n-k+1} <= n-1). So
  //   [t^n] p(q) = c_1 * q_n + (terms in q_1,...,q_{n-1}),
  // giving a direct recursion q_n = -known / c_1.
  // -------------------------------------------------------------------------
  invert(p, L) {
    if (p.length < 2) throw new Error("Taylor.invert: input has no t^1 term");
    if (Complex.abs2(p[0]) > 1e-20) throw new Error("Taylor.invert: nonzero constant term");
    const c1 = p[1];
    if (Complex.abs2(c1) < 1e-30) throw new Error("Taylor.invert: c_1 is zero");

    const q = Taylor.zero(L + 1);
    q[1] = Complex.inv(c1);

    // We will iteratively compute q^k for k=2,...,n (where n is the order we're
    // solving for). Within the n-loop, q[n] is treated as 0 (still unset).
    for (let n = 2; n <= L; n++) {
      // Build q_no_n = q with q[n] forced to 0 (it already is, since we're
      // working on this index). Length n+1.
      const q_no_n = Taylor.truncate(q, n);

      // Σ_{k=2}^{n} c_k * [t^n] (q_no_n)^k
      let known = {re: 0, im: 0};
      let qPow = q_no_n;  // q^1
      const kMax = Math.min(n, p.length - 1);
      for (let k = 2; k <= kMax; k++) {
        qPow = Taylor.mul(qPow, q_no_n, n);
        known = Complex.add(known, Complex.mul(p[k], qPow[n]));
      }

      // Solve c_1 * q[n] + known = 0
      q[n] = Complex.neg(Complex.div(known, c1));
    }

    return q;
  },

  // -------------------------------------------------------------------------
  // Compositional exp:
  //   given p = [p_0, p_1, p_2, ...], return q = exp(p) truncated to L+1.
  //
  // Standard ODE recursion: q' = p' · q, so reading [t^l] gives
  //   (l+1) q_{l+1} = Σ_{k=0..l} (k+1) p_{k+1} q_{l-k},
  // and q_0 = exp(p_0). Requires no constraint on p_0 (the constant term
  // just scales the output by a factor exp(p_0)).
  //
  // For our LQD use, p_0 is typically nonzero (it's r#(z_0), the value
  // of the rational function at the expansion point), and the constant
  // scale exp(p_0) is just absorbed downstream.
  // -------------------------------------------------------------------------
  exp(p, L) {
    const q = Taylor.zero(L + 1);

    // q_0 = exp(p_0). For complex p_0 = a + bi:
    //   exp(a + bi) = exp(a) · (cos(b) + i sin(b))
    const p0 = (p && p.length > 0) ? p[0] : {re: 0, im: 0};
    const ea = Math.exp(p0.re);
    q[0] = { re: ea * Math.cos(p0.im), im: ea * Math.sin(p0.im) };

    for (let l = 0; l < L; l++) {
      // q_{l+1} = (1/(l+1)) · Σ_{k=0..l} (k+1) p_{k+1} q_{l-k}
      let acc = {re: 0, im: 0};
      for (let k = 0; k <= l; k++) {
        const pk1 = (k + 1 < p.length) ? p[k + 1] : null;
        if (!pk1 || (pk1.re === 0 && pk1.im === 0)) continue;
        const term = Complex.scale(Complex.mul(pk1, q[l - k]), k + 1);
        acc = Complex.add(acc, term);
      }
      q[l + 1] = Complex.scale(acc, 1 / (l + 1));
    }

    return q;
  },

  // -------------------------------------------------------------------------
  // Compositional log:
  //   given p = [p_0, p_1, ...] with p_0 != 0, return q = log(p) truncated
  //   to L+1.
  //
  // Derivation: q = log(p) ⇒ p·q' = p'. Reading [t^{l-1}] for l ≥ 1:
  //   Σ_{k=0..l-1} p_k · (l-k) · q_{l-k} = l · p_l
  // so
  //   q_l = p_l/p_0  −  Σ_{k=1..l-1} ((l-k)/l) · (p_k/p_0) · q_{l-k}.
  //
  // q_0 = log(p_0) is taken on the PRINCIPAL branch (atan2 ∈ (-π, π]).
  // For our LQD-singular use, p is the Taylor of φ/b_{z_0} at z = z_j (a
  // non-vanishing holomorphic function on 𝔻̄ by construction), so p_0 ≠ 0
  // and the branch choice is locally unambiguous on a small disk around z_j.
  // -------------------------------------------------------------------------
  log(p, L) {
    if (!p || p.length === 0) throw new Error("Taylor.log: empty input");
    if (Complex.abs2(p[0]) < 1e-30) throw new Error("Taylor.log: zero constant term");

    const q = Taylor.zero(L + 1);
    // q_0 = log(p_0): complex principal log
    q[0] = { re: 0.5 * Math.log(Complex.abs2(p[0])), im: Math.atan2(p[0].im, p[0].re) };

    const invP0 = Complex.inv(p[0]);

    for (let l = 1; l <= L; l++) {
      // pl/p0
      const pl = (l < p.length) ? p[l] : {re: 0, im: 0};
      let acc = Complex.mul(pl, invP0);
      // − Σ_{k=1..l-1} ((l-k)/l) · (p_k/p_0) · q_{l-k}
      for (let k = 1; k < l; k++) {
        const pk = (k < p.length) ? p[k] : null;
        if (!pk || (pk.re === 0 && pk.im === 0)) continue;
        const term = Complex.scale(
          Complex.mul(Complex.mul(pk, invP0), q[l - k]),
          (l - k) / l
        );
        acc = Complex.sub(acc, term);
      }
      q[l] = acc;
    }

    return q;
  },

  // -------------------------------------------------------------------------
  // Multiplicative reciprocal: 1/p as a Taylor of length L+1.
  //   Requires p[0] ≠ 0. Standard recursion from p · q = 1:
  //     q_0 = 1/p_0,
  //     q_l = -(1/p_0) · Σ_{k=1..l} p_k · q_{l-k},   l ≥ 1.
  // -------------------------------------------------------------------------
  reciprocal(p, L) {
    if (!p || p.length === 0) throw new Error("Taylor.reciprocal: empty input");
    const p0 = p[0];
    if (Complex.abs2(p0) < 1e-30) throw new Error("Taylor.reciprocal: zero constant term");
    const invP0 = Complex.inv(p0);
    const q = Taylor.zero(L + 1);
    q[0] = Complex.clone(invP0);
    for (let l = 1; l <= L; l++) {
      let acc = { re: 0, im: 0 };
      const kMax = Math.min(l, p.length - 1);
      for (let k = 1; k <= kMax; k++) {
        acc = Complex.add(acc, Complex.mul(p[k], q[l - k]));
      }
      q[l] = Complex.neg(Complex.mul(invP0, acc));
    }
    return q;
  },

  // -------------------------------------------------------------------------
  // Compositional q ∘ p:
  //   Given Taylor series q = [q_0, q_1, ...] and p = [p_0, p_1, ...] with
  //   p_0 = 0 (a key requirement for term-by-term convergence), compute
  //   q(p(t)) = Σ_k q_k · p(t)^k truncated to t^L.
  //
  //   p_0 ≠ 0 is rejected because the constant-term issue is ill-defined for
  //   a truncated series (you'd need all of q to evaluate exactly). Callers
  //   should subtract p[0] off before composing.
  // -------------------------------------------------------------------------
  compose(q, p, L) {
    if (!p || p.length === 0) throw new Error("Taylor.compose: empty p");
    if (Complex.abs2(p[0]) > 1e-20) {
      throw new Error("Taylor.compose: requires p[0] = 0 (got " + p[0].re + "+" + p[0].im + "i)");
    }
    const out = Taylor.zero(L + 1);
    out[0] = Complex.clone(q[0] || { re: 0, im: 0 });
    let pPow = Taylor.zero(L + 1);
    pPow[0] = { re: 1, im: 0 };                   // p^0 = 1
    for (let k = 1; k < q.length; k++) {
      pPow = Taylor.mul(pPow, p, L);
      const qk = q[k];
      if (qk.re === 0 && qk.im === 0) continue;
      for (let i = 0; i <= L; i++) {
        out[i] = Complex.add(out[i], Complex.mul(qk, pPow[i]));
      }
    }
    return out;
  }
};

if (typeof module !== 'undefined' && module.exports) module.exports = Taylor;
