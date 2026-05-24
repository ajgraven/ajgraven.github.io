// =============================================================================
// solver-faber.js -- Shared inverse Faber transform primitives.
//
// Every inverse-problem family in the app (classical bounded/unbounded QDs,
// all four LQDs) eventually solves an inverse Faber transform:
//
//   • At a finite pole a_j of h, with z_j = φ⁻¹(a_j), match a local target
//     principal part (in t = z − z_j) given by a residue list. The result
//     is the coefficient list A_{j,1..m} of the Möbius-style local series
//
//        Σ_k conj(A_{j,k}) · z^k / (1 − conj(z_j) z)^k.
//
//   • At ∞ (unbounded variants only), match the polynomial part of h via
//     the Laurent coefficients f_l of φ around ∞.
//
// Math reference: §3.2 / §4.3 (classical & PQD) and §5.4 (LQD) of the
// thesis. The primitive is family-agnostic; family-specific code is
// confined to (i) pre-processing residues (e.g. LQDs use modified residues
// D_{j,s} = a_j C_{j,s} + C_{j,s+1}) and (ii) building the input phiTilde
// from each family's own Taylor evaluator.
//
// This file must be loaded after solver.js (it depends on QD existing) and
// before any family solver (the families call QD.Faber.* during registration).
// =============================================================================
'use strict';

(function () {

  // QD namespace is created by solver.js — solver.js must load BEFORE this
  // file, so the namespace exists. (We don't depend on Family.)
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD) {
    throw new Error("solver-faber.js: solver.js must be loaded first");
  }

  // ===========================================================================
  // (★)_A  Inverse Faber transform at one finite pole.
  // ---------------------------------------------------------------------------
  // Given:
  //   residues = [r_1, r_2, ..., r_m]                 length m, complex
  //              — the target principal-part coefficients at the pole, with
  //                r_s being the coefficient of (w − a)^{-s} after whatever
  //                preprocessing the caller has done (raw C_{j,s} for QDs,
  //                modified D_{j,s} for LQDs, etc.).
  //   phiTilde = [phiTilde_0, phiTilde_1, ..., phiTilde_m]   length m+1, complex
  //              — Taylor of (φ(z_j + t) − a_j) in t, with phiTilde_0 = 0
  //                by convention (caller fills it). phiTilde_i is the i-th
  //                Taylor coefficient (= φ^{(i)}(z_j) / i!).
  //
  // Returns:
  //   A = [A_1, A_2, ..., A_m]                              length m, complex
  //   where  A_k = Σ_{s ≥ k} (s/k) · r_s · [t^s] ψ̃^k(t)
  //   and ψ̃(t) is the formal inverse of phiTilde (so phiTilde(ψ̃(t)) = t).
  //
  // Requires QD.Taylor (defined in taylor.js) for the series inversion and
  // multiplication.
  // ===========================================================================
  function inverseFaberAtPole(residues, phiTilde) {
    const m = residues.length;
    if (m === 0) return [];
    if (phiTilde.length < m + 1) {
      throw new Error(`inverseFaberAtPole: phiTilde must have length ≥ m+1 (got ${phiTilde.length}, need ${m + 1})`);
    }
    // We use the global Taylor / Complex from the existing shared scripts.
    // ψ̃(t) = phiTilde⁻¹(t), truncated to order m.
    const psiTilde = Taylor.invert(phiTilde, m);
    const psiPow = [null];
    psiPow[1] = Taylor.truncate(psiTilde, m);
    for (let k = 2; k <= m; k++) {
      psiPow[k] = Taylor.mul(psiPow[k - 1], psiTilde, m);
    }
    const A = new Array(m);
    for (let k = 1; k <= m; k++) {
      let acc = { re: 0, im: 0 };
      for (let s = k; s <= m; s++) {
        const ts = psiPow[k][s];
        const term = Complex.scale(Complex.mul(residues[s - 1], ts), s / k);
        acc = Complex.add(acc, term);
      }
      A[k - 1] = acc;
    }
    return A;
  }

  // ===========================================================================
  // (★)_F  Inverse Faber transform at ∞ for the polynomial part of h.
  // ---------------------------------------------------------------------------
  // Used by unbounded families when h has a nontrivial polynomial part
  //   h(w) = Σ_{l=0..m_∞} C_{∞,l} w^l
  // matched against the Laurent expansion of φ around ∞:
  //   φ(z) = c·z + f_0 + f_1/z + f_2/z² + ...
  //
  // Given:
  //   polyPart = [C_{∞,0}, C_{∞,1}, ..., C_{∞,m_∞}]   length m_∞+1, complex
  //   f        = [f_0, f_1, ..., f_{m_∞ - 1}]         length m_∞, complex
  //              (Caller supplies the first m_∞ Laurent coeffs; we don't peek
  //               at the φ struct.)
  //   c        = real positive  (leading coefficient at ∞)
  //
  // Returns:
  //   F = [F_0, F_1, ..., F_{m_∞}]
  //   F_l = Σ_{l' ≥ l} conj(C_{∞,l'} · [u^{l' − l}] g(u)^{l'})
  //   where g(u) = c + Σ_{i ≥ 1} f_{i-1} · u^i.  (g is the unit-disk image of
  //   φ(1/u) · u via u = 1/z, after factoring out 1/u; see solver-uqd.js
  //   comments for the derivation.)
  // ===========================================================================
  function inverseFaberAtInfinity(polyPart, f, c) {
    const m_inf = polyPart.length - 1;
    if (m_inf < 0) return [];
    if (f.length < m_inf) {
      throw new Error(`inverseFaberAtInfinity: f must have length ≥ m_∞ (got ${f.length}, need ${m_inf})`);
    }
    // g(u) = c + Σ_{i ≥ 1} f_{i-1} · u^i, as a Taylor in u up to order m_inf.
    const g = new Array(m_inf + 1);
    g[0] = { re: c, im: 0 };
    for (let l = 1; l <= m_inf; l++) g[l] = Complex.clone(f[l - 1]);

    // g(u)^l' for l' = 0..m_inf
    const gPow = [Taylor.zero(m_inf + 1)];
    gPow[0][0] = { re: 1, im: 0 };
    for (let lp = 1; lp <= m_inf; lp++) {
      gPow[lp] = Taylor.mul(gPow[lp - 1], g, m_inf);
    }

    const F = new Array(m_inf + 1);
    for (let l = 0; l <= m_inf; l++) {
      let acc = { re: 0, im: 0 };
      for (let lp = l; lp <= m_inf; lp++) {
        const C = polyPart[lp];
        const t = gPow[lp][lp - l];
        acc = Complex.add(acc, Complex.conj(Complex.mul(C, t)));
      }
      F[l] = acc;
    }
    return F;
  }

  // ===========================================================================
  // Expose
  // ===========================================================================
  QD.Faber = {
    inverseFaberAtPole,
    inverseFaberAtInfinity,
  };

})();
