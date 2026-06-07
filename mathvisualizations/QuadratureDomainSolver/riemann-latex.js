// =============================================================================
// riemann-latex.js -- Pure LaTeX generation for the Riemann-map card.
//
// Extracted from ui.js so it can be unit-tested in Node (the DOM rendering
// stays in ui.js). `QD.RiemannLatex.build(phi)` returns the three pieces the
// card needs — { symbolic, numeric, params } — as plain strings / data, with
// NO DOM dependency. node-test.js renders `symbolic` and `numeric` through
// KaTeX with throwOnError:true for every solved family, which guards the class
// of bug where a formula silently renders as a KaTeX error node (e.g. the
// `\\[2pt]` row separator that KaTeX 0.16.11 rejects).
//
//   build(phi) → {
//     symbolic: <LaTeX>,                 // the (1) symbolic identity
//     numeric:  <LaTeX>,                 // the (2) value-substituted closed form
//     params:   [ { name:<LaTeX>, value:<complex> }, … ],  // (3) parameter rows
//   }
//
// Per-family pieces live in RIEMANN_FRAGMENTS, keyed by phi.family (with
// `_boundedQD` / `_unboundedQD` fallbacks for legacy untagged φ). Each fragment
// declares either a flat additive form (numericLeader/numericTrailer) or a
// `wrap` (prefix/suffix/bodyLeader/bodyTrailer) for families whose φ wraps the
// branch-sum in a delimiter — exp(…) for LQDs, (…)^{1/α} for PQDs.
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD) return;

  // Format a complex number as LaTeX-safe text (no math-mode commands; pure
  // real / pure imag / general handled naturally).
  function katexCmpx(c) {
    const re = c.re, im = c.im;
    const fmt = (x) => Number(x.toFixed(5)).toString();
    if (Math.abs(im) < 1e-12) return fmt(re);
    if (Math.abs(re) < 1e-12) {
      if (Math.abs(im - 1) < 1e-12) return 'i';
      if (Math.abs(im + 1) < 1e-12) return '-i';
      return fmt(im) + 'i';
    }
    const sign = im >= 0 ? ' + ' : ' - ';
    const iAbs = Math.abs(im);
    const iPart = Math.abs(iAbs - 1) < 1e-12 ? '' : fmt(iAbs);
    return fmt(re) + sign + iPart + 'i';
  }

  // Wrap in \left(...\right) unless the value is a bare non-negative real.
  function katexCmpxParen(c) {
    const s = katexCmpx(c);
    if (/^\d+(\.\d+)?$/.test(s)) return s;
    return String.raw`\left(${s}\right)`;
  }

  // r#(∞) — delegate to QD.LqdCommon.rHashAtInfinity (solver-lqd-common.js).
  function rHashAtInfinityForDisplay(phi) {
    return QD.LqdCommon.rHashAtInfinity(phi);
  }

  const RIEMANN_FRAGMENTS = {
    // Bounded classical QD (legacy untagged phi or family === 'boundedQD'):
    // φ(z) = w_0 + Σ … . No leader transforms; sum is appended directly.
    '_boundedQD': {
      symbolic: (sumBody) => String.raw`\varphi(z) \;=\; w_0 \;+\; ${sumBody}`,
      numericLeader: (phi) => [String.raw`\varphi(z) \;\approx\; ${katexCmpx(phi.w0)}`],
      numericTrailer: () => [],
      extraParameterRows: () => [],
    },
    // Unbounded classical QD:
    // φ(z) = c·z + Σ_l F_l / z^l + Σ … .
    '_unboundedQD': {
      symbolic: (sumBody, unb, mInf) => {
        const polySym = mInf >= 0
          ? String.raw` \;+\; \sum_{l=0}^{m_\infty} \frac{F_{l}}{z^{l}}`
          : '';
        return String.raw`\varphi(z) \;=\; c\, z${polySym} \;+\; ${sumBody}`;
      },
      numericLeader: (phi) => {
        const rows = [String.raw`\varphi(z) \;\approx\; ${Number(phi.c.toFixed(5))}\, z`];
        const polyA = phi.polyA || [];
        for (let l = 0; l < polyA.length; l++) {
          const F = polyA[l];
          if (l === 0) {
            rows.push(String.raw`+\, ${katexCmpxParen(F)}`);
          } else {
            const zl = l === 1 ? 'z' : `z^{${l}}`;
            rows.push(String.raw`+\, \dfrac{${katexCmpxParen(F)}}{${zl}}`);
          }
        }
        return rows;
      },
      numericTrailer: () => [],
      extraParameterRows: () => [],
    },
    // Bounded non-singular LQD: φ(z) = w_0 · exp(Σ …).
    'boundedLQD': {
      symbolic: (sumBody) =>
        String.raw`\varphi(z) \;=\; w_0 \cdot \exp\!\left(${sumBody}\right)`,
      wrap: {
        prefix: (phi) => String.raw`${katexCmpx(phi.w0)} \cdot \exp`,
        suffix: () => '',
        bodyLeader: () => [],
        bodyTrailer: () => [],
      },
      extraParameterRows: () => [],
    },
    // Unbounded non-singular LQD: φ(z) = c·z · exp(r#(z) − r#(∞)) on 𝔻*.
    // The (− r#(∞)) absorbs the ∞-gauge so the displayed leading coefficient
    // really is c. The numerical expansion shows the actual A_{j,k} the solver
    // returned and the additional − r#(∞) term (computed in closed form).
    'unboundedLQD': {
      symbolic: (sumBody) =>
        String.raw`\varphi(z) \;=\; c\, z \cdot \exp\!\left(${sumBody} \;-\; r_\#(\infty)\right)`,
      wrap: {
        prefix: (phi) => String.raw`${Number(phi.c.toFixed(5))}\, z \cdot \exp`,
        suffix: () => '',
        bodyLeader: () => [],
        bodyTrailer: (phi) => [
          String.raw`-\, ${katexCmpxParen(rHashAtInfinityForDisplay(phi))}`,
        ],
      },
      extraParameterRows: (phi) => [
        { name: String.raw`r_\#(\infty)`, value: rHashAtInfinityForDisplay(phi) },
      ],
    },
    // Unbounded singular LQD: φ(z) = c·|z₀|·z·b_{z₀}(z)·exp(r#(z) − r#(∞)).
    'unboundedLQD_singular': {
      // Blaschke definition on its own aligned line (it's too wide to share a
      // line with the φ identity — was forcing a horizontal scrollbar).
      symbolic: (sumBody) =>
        String.raw`\begin{aligned}\varphi(z) &= c\cdot|z_0|\cdot z\cdot b_{z_0}(z) \cdot \exp\!\left(${sumBody} \;-\; r_\#(\infty)\right) \\ b_{z_0}(z) &= -\tfrac{\overline{z_0}}{|z_0|}\cdot\tfrac{z - z_0}{1 - \overline{z_0}\, z}\end{aligned}`,
      wrap: {
        prefix: (phi) => {
          const absZ0 = QD.Complex.abs(phi.z0).toFixed(5);
          const z0Latex = katexCmpxParen(phi.z0);
          return String.raw`${Number(phi.c.toFixed(5))}\cdot ${absZ0}\cdot z\cdot b_{${z0Latex}}(z) \cdot \exp`;
        },
        suffix: () => '',
        bodyLeader: () => [],
        bodyTrailer: (phi) => [
          String.raw`-\, ${katexCmpxParen(rHashAtInfinityForDisplay(phi))}`,
        ],
      },
      extraParameterRows: (phi) => [
        { name: String.raw`z_0`,        value: phi.z0 },
        { name: String.raw`q`,          value: phi.q },
        { name: String.raw`r_\#(\infty)`, value: rHashAtInfinityForDisplay(phi) },
      ],
    },
    // Bounded power-weighted QD (Family.powerQD, α ≥ 2):
    //   φ(z) = (R#(z))^{1/α},   R#(z) = w_0^α + Σ_k Σ_j conj(A_{j,k}) z^k/(1-conj(z_k)z)^k
    // The sumBody is R#'s pole sum; we wrap the whole rational in (· · ·)^{1/α}.
    // The (·)^{1/α} denotes the single continuous branch ANCHORED at φ(0)=w_0
    // (see QD.PqdCommon.phiAnchored), not the principal branch — so the form is
    // valid for off-axis poles (|arg a| > π/α) too.
    'powerQD': {
      symbolic: (sumBody, _unb, _mInf, phi) => {
        const a = (phi && phi.alpha) || 2;
        return String.raw`\varphi(z) \;=\; \left(w_0^{${a}} \;+\; ${sumBody}\right)^{1/${a}}`;
      },
      wrap: {
        prefix: () => '',
        suffix: (phi) => String.raw`^{1/${phi.alpha || 2}}`,
        bodyLeader: (phi) => {
          const a = phi.alpha || 2;
          let r0 = { re: 1, im: 0 };           // r0 = w0^α
          for (let k = 0; k < a; k++) r0 = QD.Complex.mul(r0, phi.w0);
          return [katexCmpx(r0)];
        },
        bodyTrailer: () => [],
      },
      extraParameterRows: (phi) => [
        { name: String.raw`\alpha`, value: { re: phi.alpha || 2, im: 0 } },
      ],
    },
    // Bounded singular LQD: φ(z) = γ · b_{z_0}(z) · exp(Σ …).
    'boundedLQD_singular': {
      // Blaschke definition on its own aligned line (too wide to share a line).
      symbolic: (sumBody) =>
        String.raw`\begin{aligned}\varphi(z) &= \gamma \cdot b_{z_0}(z) \cdot \exp\!\left(${sumBody}\right) \\ b_{z_0}(z) &= -\tfrac{\overline{z_0}}{|z_0|}\cdot\tfrac{z - z_0}{1 - \overline{z_0}\, z}\end{aligned}`,
      wrap: {
        prefix: (phi) => String.raw`${katexCmpxParen(phi.gamma)} \cdot b_{${katexCmpxParen(phi.z0)}}(z) \cdot \exp`,
        suffix: () => '',
        bodyLeader: () => [],
        bodyTrailer: () => [],
      },
      extraParameterRows: (phi) => [
        { name: String.raw`z_0`,    value: phi.z0 },
        { name: String.raw`\gamma`, value: phi.gamma },
        { name: String.raw`q`,      value: phi.q },
      ],
    },
    // Bounded singular PQD (Family.powerQD_singular, 0 ∈ Ω):
    //   φ(z) = b_{z_0}(z) · (R#(z))^{1/α},  R#(z) = w_0^α/|z_0|^α + Σ conj(A) z^k/(1-conj(z_j)z)^k
    'powerQD_singular': {
      symbolic: (sumBody, _unb, _mInf, phi) => {
        const a = (phi && phi.alpha) || 2;
        return String.raw`\varphi(z) \;=\; b_{z_0}(z)\,\left(\tfrac{w_0^{${a}}}{|z_0|^{${a}}} \;+\; ${sumBody}\right)^{1/${a}}`;
      },
      wrap: {
        prefix: (phi) => String.raw`b_{${katexCmpxParen(phi.z0)}}(z)\,`,
        suffix: (phi) => String.raw`^{1/${phi.alpha || 2}}`,
        bodyLeader: (phi) => {
          const a = phi.alpha || 2;
          const z0a = Math.pow(QD.Complex.abs2(phi.z0), 0.5 * a);  // |z₀|^α
          const r0 = QD.Complex.scale(QD.Complex.cpow(phi.w0, a), z0a > 0 ? 1 / z0a : 0);
          return [katexCmpx(r0)];
        },
        bodyTrailer: () => [],
      },
      extraParameterRows: (phi) => [
        { name: String.raw`z_0`,    value: phi.z0 },
        { name: String.raw`\alpha`, value: { re: phi.alpha || 2, im: 0 } },
      ],
    },
    // Unbounded PQD (Family.unboundedPQD, 0 ∉ Ω):
    //   φ(z) = z·(r#(z))^{1/α},  r#(z) = c^α + Σ G_l/z^l + Σ conj(A) z^k/(1-conj(z_j)z)^k
    'unboundedPQD': {
      symbolic: (sumBody, _unb, _mInf, phi) => {
        const a = (phi && phi.alpha) || 2;
        return String.raw`\varphi(z) \;=\; z\,\left(c^{${a}} \;+\; ${sumBody}\right)^{1/${a}}`;
      },
      wrap: {
        prefix: () => String.raw`z\,`,
        suffix: (phi) => String.raw`^{1/${phi.alpha || 2}}`,
        bodyLeader: (phi) => [katexCmpx({ re: Math.pow(phi.c, phi.alpha || 2), im: 0 })],
        bodyTrailer: () => [],
      },
      extraParameterRows: (phi) => [
        { name: String.raw`c`,      value: { re: phi.c, im: 0 } },
        { name: String.raw`\alpha`, value: { re: phi.alpha || 2, im: 0 } },
      ],
    },
    // Unbounded singular PQD (Family.unboundedPQD_singular, 0 ∈ Ω):
    //   φ(z) = z·b_{z₀}(z)·(r#(z))^{1/α},  r#(∞)=|cz₀|^α,  r(z₀)=0 (Prop 4.6.3).
    'unboundedPQD_singular': {
      symbolic: (sumBody, _unb, _mInf, phi) => {
        const a = (phi && phi.alpha) || 2;
        return String.raw`\varphi(z) \;=\; z\,b_{z_0}(z)\,\left(|cz_0|^{${a}} \;+\; ${sumBody}\right)^{1/${a}}`;
      },
      wrap: {
        prefix: (phi) => String.raw`z\,b_{${katexCmpxParen(phi.z0)}}(z)\,`,
        suffix: (phi) => String.raw`^{1/${phi.alpha || 2}}`,
        bodyLeader: (phi) => [katexCmpx({ re: Math.pow(phi.c * QD.Complex.abs(phi.z0), phi.alpha || 2), im: 0 })],
        bodyTrailer: () => [],
      },
      extraParameterRows: (phi) => [
        { name: String.raw`z_0`,    value: phi.z0 },
        { name: String.raw`c`,      value: { re: phi.c, im: 0 } },
        { name: String.raw`\alpha`, value: { re: phi.alpha || 2, im: 0 } },
      ],
    },
  };

  function getRiemannFragment(phi) {
    if (phi.family && RIEMANN_FRAGMENTS[phi.family]) return RIEMANN_FRAGMENTS[phi.family];
    // Legacy: family tag absent → QD/UQD by phi.unbounded.
    return phi.unbounded ? RIEMANN_FRAGMENTS['_unboundedQD'] : RIEMANN_FRAGMENTS['_boundedQD'];
  }

  // build(phi) — pure assembly of the three card pieces (no DOM).
  function build(phi) {
    const maxOrder = phi.branches.reduce((m, b) => Math.max(m, b.A.length), 1);
    const unb = !!phi.unbounded;
    const polyA = phi.polyA || [];
    const m_inf = polyA.length - 1;
    const frag = getRiemannFragment(phi);

    // --- (1) symbolic identity ---
    const sumBody = maxOrder === 1
      ? String.raw`\sum_{j} \overline{A_{j}}\, \frac{z}{1 - \overline{z_{j}}\, z}`
      : String.raw`\sum_{j,k} \overline{A_{j,k}}\, \frac{z^{k}}{\bigl(1 - \overline{z_{j}}\, z\bigr)^{k}}`;
    const symbolic = frag.symbolic(sumBody, unb, m_inf, phi);

    // --- (2) closed-form expression with values substituted ---
    // Build the additive sum-of-branches rows. Each row starts with "+\, ".
    const sumTermsRows = [];
    phi.branches.forEach((br) => {
      const zjC = QD.Complex.conj(br.z);
      br.A.forEach((a, k) => {
        const aC = QD.Complex.conj(a);
        const power = k + 1;
        const numerator = power === 1 ? 'z' : `z^{${power}}`;
        if (QD.Complex.abs2(zjC) < 1e-20) {
          // z_j = 0 ⇒ the denominator (1 − conj(z_j) z)^k collapses to 1, so the
          // term is just conj(A_j) z^k — avoids the ugly "1 − 0·z".
          sumTermsRows.push(String.raw`+\, ${katexCmpxParen(aC)}\, ${numerator}`);
        } else {
          const denomCore = String.raw`1 - ${katexCmpxParen(zjC)}\, z`;
          const denom = power === 1 ? `\\bigl(${denomCore}\\bigr)` : `\\bigl(${denomCore}\\bigr)^{${power}}`;
          sumTermsRows.push(String.raw`+\, ${katexCmpxParen(aC)}\, \dfrac{${numerator}}{${denom}}`);
        }
      });
    });

    // The closed-form expansion is rendered on a SINGLE line (no `aligned`
    // row-stacking — that forced line breaks even when the formula would fit).
    // Width is handled by `.rm-numer { overflow-x: auto }`, so long expressions
    // scroll horizontally rather than wrap at arbitrary points.
    let numeric;
    if (frag.wrap) {
      // Families whose φ wraps the branch-sum in a delimiter — exp(…) for the
      // LQDs, (…)^{1/α} for the PQDs. The body (constant + branch terms + any
      // trailer) goes inside a single \left( … \right) that auto-sizes to the
      // (single-line) content.
      const w = frag.wrap;
      const lead  = w.bodyLeader  ? w.bodyLeader(phi)  : [];
      const trail = w.bodyTrailer ? w.bodyTrailer(phi) : [];
      const body = [];
      for (const r of lead) body.push(r);
      const stripFirst = lead.length === 0;   // exp(Σ…): no constant before the sum
      if (sumTermsRows.length === 0 && stripFirst) {
        body.push('0');                        // empty sum → literal 0
      } else {
        sumTermsRows.forEach((row, i) => {
          body.push((i === 0 && stripFirst) ? row.replace(/^\+\\,\s*/, '') : row);
        });
      }
      for (const r of trail) body.push(r);
      const inner = body.join(' ');            // "const +\, t1 +\, t2 -\, r#(∞)"
      numeric = String.raw`\varphi(z) \;\approx\; ${w.prefix(phi)}{\left(${inner}\right)}${w.suffix(phi)}`;
    } else {
      // Flat additive families (classical bounded / unbounded QD): leader +
      // branch terms + trailer, concatenated inline ("φ(z) ≈ w₀ +\, t1 +\, t2").
      const numerRows = frag.numericLeader(phi);
      for (const row of sumTermsRows) numerRows.push(row);
      for (const row of frag.numericTrailer(phi)) numerRows.push(row);
      numeric = numerRows.join(' ');
    }

    // --- (3) parameter rows (pure data; caller renders the table) ---
    const params = [];
    if (unb) {
      params.push({ name: String.raw`c`, value: { re: phi.c, im: 0 } });
    } else {
      params.push({ name: String.raw`w_0`, value: phi.w0 });
    }
    for (const r of frag.extraParameterRows(phi)) params.push(r);
    phi.branches.forEach((br, j) => {
      params.push({ name: String.raw`z_{${j + 1}}`, value: br.z });
      br.A.forEach((a, k) => {
        const sub = br.A.length === 1 ? `${j + 1}` : `${j + 1},${k + 1}`;
        params.push({ name: String.raw`A_{${sub}}`, value: a });
      });
    });
    if (unb && m_inf >= 0) {
      for (let l = 0; l <= m_inf; l++) {
        params.push({ name: String.raw`F_{${l}}`, value: polyA[l] });
      }
    }

    return { symbolic, numeric, params };
  }

  QD.RiemannLatex = { build, katexCmpx, katexCmpxParen, getRiemannFragment };

})();
