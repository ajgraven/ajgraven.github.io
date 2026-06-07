// =============================================================================
// schwarz-analysis.js — Phase S4: Schwarz-function analysis utilities.
//
// Three independent kernels:
//   1. explicitSigmaForm(schwarz)
//        → { family, phiText, phiLatex, fText, fLatex, sigmaText, sigmaLatex }
//      Per-family closed-form expressions with current numerical values plugged
//      in. Used by the "Show σ(w)" panel (E13 / §9.1).
//
//   2. findSigmaSingularities(schwarz)
//        → { poles: [{w, kind, label}], branchPoints: [{w, label, severity}] }
//      σ-poles via F-pole pullback through φ; σ-branch points via
//      QD.findCriticalPoints (zeros of φ' mapped to w-plane). Used by the
//      singularity-analyzer overlay (F3).
//
//   3. computeSigmaLevelCurves(schwarz, opts)
//        → { abs: ContourPolyline[], arg: ContourPolyline[] }
//      Marching-squares contours of |σ(w)| (solid) and arg(σ(w)) (dashed).
//      Used by the level-curve overlay (F12 / §A6).
// =============================================================================

(function () {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD || !QD.Schwarz) {
    throw new Error("schwarz-analysis.js: schwarz-common.js must be loaded first");
  }

  const C = QD.Complex;

  // ---------------------------------------------------------------------------
  // Helpers — number → text and number → LaTeX. Compact forms; small reals get
  // 4-decimal text, complex numbers print as `a + b i` or LaTeX `a + bi`.
  // ---------------------------------------------------------------------------
  function fmtReal(x, digits) {
    digits = digits != null ? digits : 4;
    if (!isFinite(x)) return String(x);
    if (Math.abs(x) < 5 * Math.pow(10, -(digits + 1))) return '0';
    return (+x.toFixed(digits)).toString();
  }
  function fmtCText(c) {
    if (!c) return '0';
    const re = fmtReal(c.re), im = fmtReal(c.im);
    if (im === '0') return re;
    if (re === '0') return im === '1' ? 'i' : (im === '-1' ? '-i' : im + 'i');
    return re + (c.im >= 0 ? '+' : '') + im + 'i';
  }
  // (A LaTeX numeric-substitution formatter was removed here as dead code —
  // current LaTeX strings use abstract symbols. See HANDOFF "deferred" TODOs
  // if a numerical-substitution σ panel is revived.)
  function sub(k) {
    const digits = ['₀','₁','₂','₃','₄','₅','₆','₇','₈','₉'];
    return String(k).split('').map(d => digits[+d] || d).join('');
  }

  // ---------------------------------------------------------------------------
  // 1. explicitSigmaForm — per-family formula expressions.
  // ---------------------------------------------------------------------------
  function explicitSigmaForm(schwarz) {
    if (!schwarz || !schwarz._phi) {
      return { family: '(no φ captured)', phiText: '', phiLatex: '',
               fText: '', fLatex: '', sigmaText: '', sigmaLatex: '' };
    }
    const phi    = schwarz._phi;
    const family = schwarz.family;

    // ------- bounded classical QD -------
    if (family === 'boundedQD') {
      const w0 = phi.w0 || { re: 0, im: 0 };
      const branches = phi.branches || [];
      // φ(z) = w₀ + Σⱼ Σₖ conj(A_{j,k}) z^k / (1 − conj(z_j) z)^k
      const phiTerms = [];
      const fTerms   = [];
      for (let j = 0; j < branches.length; j++) {
        const z_j = branches[j].z;
        for (let k = 1; k <= branches[j].A.length; k++) {
          const A   = branches[j].A[k - 1];
          const Aco = C.conj(A);
          phiTerms.push(
            `(${fmtCText(Aco)}) z^${k}/(1−(${fmtCText(C.conj(z_j))})z)^${k}`
          );
          fTerms.push(`(${fmtCText(A)})/(z−(${fmtCText(z_j)}))^${k}`);
        }
      }
      const phiText = `φ(z) = ${fmtCText(w0)}` +
                     (phiTerms.length ? ` + ${phiTerms.join(' + ')}` : '');
      const fText   = `F(z) = ${fmtCText(C.conj(w0))}` +
                     (fTerms.length ? ` + ${fTerms.join(' + ')}` : '');
      // LaTeX versions
      const phiLatexTerms = branches.flatMap((br, j) => {
        return br.A.map((A, kIdx) => {
          const k = kIdx + 1;
          return `\\overline{A_{${j + 1},${k}}}\\,\\frac{z^{${k}}}{(1-\\overline{z_{${j + 1}}}z)^{${k}}}`;
        });
      });
      const fLatexTerms = branches.flatMap((br, j) => {
        return br.A.map((_, kIdx) => {
          const k = kIdx + 1;
          return `\\frac{A_{${j + 1},${k}}}{(z-z_{${j + 1}})^{${k}}}`;
        });
      });
      const phiLatex = `\\varphi(z) = w_0 ${phiLatexTerms.length ? '+ ' + phiLatexTerms.join(' + ') : ''}`;
      const fLatex   = `F(z) = \\overline{w_0} ${fLatexTerms.length ? '+ ' + fLatexTerms.join(' + ') : ''}`;
      const sigmaText  = 'σ(w) = conj(F(ψ(w)))     where ψ = φ⁻¹';
      const sigmaLatex = '\\sigma(w) = \\overline{F(\\psi(w))}, \\quad \\psi = \\varphi^{-1}';
      return { family, phiText, phiLatex, fText, fLatex, sigmaText, sigmaLatex };
    }

    // ------- unbounded classical QD -------
    if (family === 'unboundedQD') {
      const c = phi.c;
      const F = phi.polyA || phi.F || [];
      const branches = phi.branches || [];
      const phiTerms = [];
      const fTerms   = [];
      for (let l = 0; l < F.length; l++) {
        phiTerms.push(`(${fmtCText(F[l])})/z^${l}`);
        fTerms.push(`(${fmtCText(C.conj(F[l]))})z^${l}`);
      }
      for (let j = 0; j < branches.length; j++) {
        const z_j = branches[j].z;
        for (let k = 1; k <= branches[j].A.length; k++) {
          const A = branches[j].A[k - 1];
          phiTerms.push(`(${fmtCText(C.conj(A))}) z^${k}/(1−(${fmtCText(C.conj(z_j))})z)^${k}`);
          fTerms.push(`(${fmtCText(A)})/(z−(${fmtCText(z_j)}))^${k}`);
        }
      }
      const phiText = `φ(z) = ${fmtReal(c)}·z` +
                     (phiTerms.length ? ` + ${phiTerms.join(' + ')}` : '');
      const fText   = `F(z) = ${fmtReal(c)}/z` +
                     (fTerms.length ? ` + ${fTerms.join(' + ')}` : '');
      const phiLatex = `\\varphi(z) = c\\,z ${F.length ? '+ \\sum_{l=0}^{m_\\infty} \\frac{F_l}{z^l}' : ''} ${branches.length ? '+ \\text{(branches)}' : ''}`;
      const fLatex   = `F(z) = \\frac{c}{z} ${F.length ? '+ \\sum_{l=0}^{m_\\infty} \\overline{F_l}\\,z^l' : ''} ${branches.length ? '+ \\text{(branches)}' : ''}`;
      const sigmaText  = 'σ(w) = conj(F(ψ(w)))     where ψ = φ⁻¹';
      const sigmaLatex = '\\sigma(w) = \\overline{F(\\psi(w))}, \\quad \\psi = \\varphi^{-1}';
      return { family, phiText, phiLatex, fText, fLatex, sigmaText, sigmaLatex };
    }

    // ------- bounded power-weighted QD (Family.powerQD, α ≥ 2) -------
    if (family === 'powerQD') {
      const alpha    = phi.alpha;
      const w0       = phi.w0 || { re: 0, im: 0 };
      const branches = phi.branches || [];
      // R#(z) = w₀^α + Σ conj(A_{j,k}) z^k/(1−conj(z_j)z)^k
      // R(z)  = conj(w₀)^α + Σ A_{j,k}/(z−z_j)^k
      let w0Pow = { re: 1, im: 0 };
      for (let k = 0; k < alpha; k++) w0Pow = C.mul(w0Pow, w0);
      const rhTerms = [];
      const rTerms  = [];
      for (let j = 0; j < branches.length; j++) {
        const z_j = branches[j].z;
        for (let k = 1; k <= branches[j].A.length; k++) {
          const A   = branches[j].A[k - 1];
          const Aco = C.conj(A);
          rhTerms.push(`(${fmtCText(Aco)}) z^${k}/(1−(${fmtCText(C.conj(z_j))})z)^${k}`);
          rTerms.push(`(${fmtCText(A)})/(z−(${fmtCText(z_j)}))^${k}`);
        }
      }
      const phiText = `φ(z) = (R#(z))^{1/${alpha}},   R#(z) = ${fmtCText(w0Pow)}` +
                     (rhTerms.length ? ` + ${rhTerms.join(' + ')}` : '');
      const fText   = `F(z) = (R(z))^{1/${alpha}},   R(z) = ${fmtCText(C.conj(w0Pow))}` +
                     (rTerms.length ? ` + ${rTerms.join(' + ')}` : '');
      const phiLatex = `\\varphi(z) = \\left(w_0^{${alpha}} + \\sum_{j,k}\\overline{A_{j,k}}\\,\\frac{z^{k}}{(1-\\overline{z_j}z)^{k}}\\right)^{1/${alpha}}`;
      const fLatex   = `F(z) = \\left(\\overline{w_0}^{${alpha}} + \\sum_{j,k}\\frac{A_{j,k}}{(z-z_j)^{k}}\\right)^{1/${alpha}}`;
      const sigmaText  = 'σ(w) = conj(F(ψ(w)))     where ψ = φ⁻¹';
      const sigmaLatex = '\\sigma(w) = \\overline{F(\\psi(w))}, \\quad \\psi = \\varphi^{-1}';
      return { family, phiText, phiLatex, fText, fLatex, sigmaText, sigmaLatex };
    }

    // ------- unbounded power-weighted QD (Family.unboundedPQD) -------
    if (family === 'unboundedPQD') {
      const alpha    = phi.alpha;
      const c        = phi.c;
      const branches = phi.branches || [];
      const G        = phi.polyA || [];
      const rhTerms = [], rTerms = [];
      for (let l = 1; l <= G.length; l++) {
        rhTerms.push(`(${fmtCText(G[l - 1])})/z^${l}`);
        rTerms.push(`(${fmtCText(C.conj(G[l - 1]))}) z^${l}`);
      }
      for (let j = 0; j < branches.length; j++) {
        const z_j = branches[j].z;
        for (let k = 1; k <= branches[j].A.length; k++) {
          const A = branches[j].A[k - 1];
          rhTerms.push(`(${fmtCText(C.conj(A))}) z^${k}/(1−(${fmtCText(C.conj(z_j))})z)^${k}`);
          rTerms.push(`(${fmtCText(A)})/(z−(${fmtCText(z_j)}))^${k}`);
        }
      }
      const ca = Math.pow(c, alpha);
      const phiText = `φ(z) = z·(R#(z))^{1/${alpha}},   R#(z) = ${ca.toPrecision(6)}` +
                     (rhTerms.length ? ` + ${rhTerms.join(' + ')}` : '');
      const fText   = `F(z) = (1/z)·(R(z))^{1/${alpha}},   R(z) = ${ca.toPrecision(6)}` +
                     (rTerms.length ? ` + ${rTerms.join(' + ')}` : '');
      const phiLatex = `\\varphi(z) = z\\left(c^{${alpha}} + \\sum_l \\frac{G_l}{z^l} + \\sum_{j,k}\\overline{A_{j,k}}\\,\\frac{z^{k}}{(1-\\overline{z_j}z)^{k}}\\right)^{1/${alpha}}`;
      const fLatex   = `F(z) = \\tfrac1z\\left(c^{${alpha}} + \\sum_l \\overline{G_l} z^l + \\sum_{j,k}\\frac{A_{j,k}}{(z-z_j)^{k}}\\right)^{1/${alpha}}`;
      const sigmaText  = 'σ(w) = conj(F(ψ(w)))     where ψ = φ⁻¹';
      const sigmaLatex = '\\sigma(w) = \\overline{F(\\psi(w))}, \\quad \\psi = \\varphi^{-1}';
      return { family, phiText, phiLatex, fText, fLatex, sigmaText, sigmaLatex };
    }

    // ------- unbounded SINGULAR power-weighted QD (Family.unboundedPQD_singular) -------
    if (family === 'unboundedPQD_singular') {
      const alpha = phi.alpha, c = phi.c, z0 = phi.z0 || { re: 0, im: 0 };
      const cz0a = Math.pow(c * Math.hypot(z0.re, z0.im), alpha);
      const phiText = `φ(z) = z·b_{z₀}(z)·(R#(z))^{1/${alpha}},  R#(∞)=${cz0a.toPrecision(6)} (=|cz₀|^${alpha}),  z₀=${fmtCText(z0)}`;
      const fText   = `F(z) = (1/z)·b#_{z₀}(z)·(R(z))^{1/${alpha}}`;
      const phiLatex = `\\varphi(z) = z\\,b_{z_0}(z)\\left(|cz_0|^{${alpha}} + \\sum_l \\frac{G_l}{z^l} + \\sum_{j,k}\\overline{A_{j,k}}\\,\\frac{z^{k}}{(1-\\overline{z_j}z)^{k}}\\right)^{1/${alpha}}`;
      const fLatex   = `F(z) = \\tfrac1z\\,b^{\\#}_{z_0}(z)\\left(|cz_0|^{${alpha}} + \\sum_l \\overline{G_l} z^l + \\sum_{j,k}\\frac{A_{j,k}}{(z-z_j)^{k}}\\right)^{1/${alpha}}`;
      const sigmaText  = 'σ(w) = conj(F(ψ(w)))     where ψ = φ⁻¹';
      const sigmaLatex = '\\sigma(w) = \\overline{F(\\psi(w))}, \\quad \\psi = \\varphi^{-1}';
      return { family, phiText, phiLatex, fText, fLatex, sigmaText, sigmaLatex };
    }

    // ------- bounded SINGULAR power-weighted QD (Family.powerQD_singular) -------
    if (family === 'powerQD_singular') {
      const alpha    = phi.alpha;
      const w0       = phi.w0 || { re: 0, im: 0 };
      const z0       = phi.z0 || { re: 0, im: 0 };
      const branches = phi.branches || [];
      // R# constant = w₀^α/|z₀|^α (R = its conjugate). cpow handles any α.
      const z0absPowA = Math.pow(C.abs2(z0), 0.5 * alpha);
      const r0     = C.scale(C.cpow(w0, alpha), z0absPowA > 0 ? 1 / z0absPowA : 0);
      const r0Conj = C.conj(r0);
      const rhTerms = [], rTerms = [];
      for (let j = 0; j < branches.length; j++) {
        const z_j = branches[j].z;
        for (let k = 1; k <= branches[j].A.length; k++) {
          const A = branches[j].A[k - 1];
          rhTerms.push(`(${fmtCText(C.conj(A))}) z^${k}/(1−(${fmtCText(C.conj(z_j))})z)^${k}`);
          rTerms.push(`(${fmtCText(A)})/(z−(${fmtCText(z_j)}))^${k}`);
        }
      }
      const phiText = `φ(z) = b_{z₀}(z)·(R#(z))^{1/${alpha}},   R#(z) = ${fmtCText(r0)}` +
                     (rhTerms.length ? ` + ${rhTerms.join(' + ')}` : '') +
                     `   [z₀ = ${fmtCText(z0)}]`;
      const fText   = `F(z) = b#_{z₀}(z)·(R(z))^{1/${alpha}},   R(z) = ${fmtCText(r0Conj)}` +
                     (rTerms.length ? ` + ${rTerms.join(' + ')}` : '');
      const phiLatex = `\\varphi(z) = b_{z_0}(z)\\left(\\frac{\\overline{w_0}^{${alpha}}}{|z_0|^{${alpha}}} + \\sum_{j,k}\\overline{A_{j,k}}\\,\\frac{z^{k}}{(1-\\overline{z_j}z)^{k}}\\right)^{1/${alpha}}`;
      const fLatex   = `F(z) = b^{\\#}_{z_0}(z)\\left(\\frac{w_0^{${alpha}}}{|z_0|^{${alpha}}} + \\sum_{j,k}\\frac{A_{j,k}}{(z-z_j)^{k}}\\right)^{1/${alpha}}`;
      const sigmaText  = 'σ(w) = conj(F(ψ(w)))     where ψ = φ⁻¹';
      const sigmaLatex = '\\sigma(w) = \\overline{F(\\psi(w))}, \\quad \\psi = \\varphi^{-1}';
      return { family, phiText, phiLatex, fText, fLatex, sigmaText, sigmaLatex };
    }

    // ------- LQD families: present the compositional structure -------
    if (family === 'boundedLQD') {
      const w0 = phi.w0 || { re: 0, im: 0 };
      const phiText = `φ(z) = (${fmtCText(w0)}) · exp(r#(z))`;
      const fText   = `F(z) = (${fmtCText(C.conj(w0))}) · exp(R##(z))`;
      const phiLatex = `\\varphi(z) = w_0\\,\\exp(r^\\#(z))`;
      const fLatex   = `F(z) = \\overline{w_0}\\,\\exp(R^{\\#\\#}(z))`;
      return _lqdSigmaResult(family, phiText, phiLatex, fText, fLatex);
    }
    if (family === 'boundedLQD_singular') {
      const g = phi.gamma, z0 = phi.z0;
      const phiText  = `φ(z) = (${fmtCText(g)}) · b_{${fmtCText(z0)}}(z) · exp(r#(z))`;
      const fText    = `F(z) = (${fmtCText(C.conj(g))}) · b#_{${fmtCText(z0)}}(z) · exp(R##(z))`;
      const phiLatex = `\\varphi(z) = \\gamma\\,b_{z_0}(z)\\,\\exp(r^\\#(z))`;
      const fLatex   = `F(z) = \\overline{\\gamma}\\,b^\\#_{z_0}(z)\\,\\exp(R^{\\#\\#}(z))`;
      return _lqdSigmaResult(family, phiText, phiLatex, fText, fLatex);
    }
    if (family === 'unboundedLQD') {
      const c = phi.c;
      const phiText  = `φ(z) = (${fmtReal(c)})·z · exp(r̃#(z) + B(1/z))`;
      const fText    = `F(z) = (${fmtReal(c)})/z · exp(R̃##(z) − conj(r̃#(∞)) + conj(B(z)))`;
      const phiLatex = `\\varphi(z) = c\\,z\\,\\exp\\big(\\widetilde{r}^\\#(z) + B(1/z)\\big)`;
      const fLatex   = `F(z) = \\tfrac{c}{z}\\,\\exp\\big(\\widetilde{R}^{\\#\\#}(z) - \\overline{\\widetilde{r}^\\#(\\infty)} + \\overline{B(z)}\\big)`;
      return _lqdSigmaResult(family, phiText, phiLatex, fText, fLatex);
    }
    if (family === 'unboundedLQD_singular') {
      const c = phi.c, z0 = phi.z0;
      const phiText  = `φ(z) = (${fmtReal(c)})·|z₀|·z · b_{${fmtCText(z0)}}(z) · exp(r̃#(z) + B(1/z))`;
      const fText    = `F(z) = (${fmtReal(c)})·|z₀|/z · b#_{${fmtCText(z0)}}(z) · exp(R̃##(z) − conj(r̃#(∞)) + conj(B(z)))`;
      const phiLatex = `\\varphi(z) = c\\,|z_0|\\,z\\,b_{z_0}(z)\\,\\exp\\big(\\widetilde{r}^\\#(z)+B(1/z)\\big)`;
      const fLatex   = `F(z) = \\tfrac{c\\,|z_0|}{z}\\,b^\\#_{z_0}(z)\\,\\exp\\big(\\widetilde{R}^{\\#\\#}(z) - \\overline{\\widetilde{r}^\\#(\\infty)} + \\overline{B(z)}\\big)`;
      return _lqdSigmaResult(family, phiText, phiLatex, fText, fLatex);
    }

    // Fallback (unknown family): minimal info.
    const phiText = 'φ(z) = (custom; family ' + (family || 'unknown') + ')';
    const fText   = 'F(z) = (per-family Schwarz pullback)';
    return {
      family: family || 'unknown',
      phiText, phiLatex: phiText,
      fText,   fLatex:   fText,
      sigmaText:  'σ(w) = conj(F(ψ(w)))',
      sigmaLatex: '\\sigma(w) = \\overline{F(\\psi(w))}',
    };
  }
  function _lqdSigmaResult(family, phiText, phiLatex, fText, fLatex) {
    return {
      family, phiText, phiLatex, fText, fLatex,
      sigmaText:  'σ(w) = conj(F(ψ(w)))     where ψ = φ⁻¹',
      sigmaLatex: '\\sigma(w) = \\overline{F(\\psi(w))}, \\quad \\psi = \\varphi^{-1}',
    };
  }

  // ---------------------------------------------------------------------------
  // 2. findSigmaSingularities — poles + branch points of σ in Ω.
  //
  // Poles of σ in Ω correspond to poles of F (in z-space) pulled through φ.
  // For each finite pole z_j of F (which sits inside 𝔻 or 𝔻* depending on
  // family), the σ-pole in w-space is φ(z_jReflected) where z_jReflected =
  // 1/conj(z_j) is on the other side of the unit circle. We map z_jReflected
  // through schwarz.evalPhi to get the w-plane location.
  //
  // For LQD families F has essential singularities (exp of a rational), not
  // poles — we still report the locations as "essential" so the user knows
  // σ has bad behaviour there, but no residue.
  //
  // Branch points: zeros of φ' inside the relevant disk. We delegate to
  // QD.findCriticalPoints (the critical-set kernel that powers the
  // inverse-tab overlay) and pick the in-domain ones.
  // ---------------------------------------------------------------------------
  function findSigmaSingularities(schwarz) {
    if (!schwarz || !schwarz._phi) return { poles: [], branchPoints: [] };
    const phi      = schwarz._phi;
    const family   = schwarz.family;
    const branches = phi.branches || [];

    // ---- σ-poles: F-pole pullback ----
    const polesOut = [];
    const isLqd = family && family.indexOf('LQD') >= 0;
    const kind  = isLqd ? 'essential' : 'pole';

    for (let j = 0; j < branches.length; j++) {
      const z_j = branches[j].z;
      const absZj = Math.hypot(z_j.re, z_j.im);
      if (absZj < 1e-12) continue;                      // can't reflect z=0
      // Reflection: z_j → 1/conj(z_j) (the "other side of |z|=1").
      const z_jR = { re: z_j.re / (absZj * absZj), im: z_j.im / (absZj * absZj) };
      let wPole;
      try { wPole = schwarz.evalPhi(z_jR); }
      catch (_) { continue; }
      if (!isFinite(wPole.re) || !isFinite(wPole.im)) continue;
      polesOut.push({
        w: wPole,
        kind,
        order: branches[j].A.length,
        label: 'a' + sub(j + 1),
      });
    }

    // unboundedQD: F has an additional pole at z = 0 → σ-pole at φ(∞) = ∞.
    if (family === 'unboundedQD') {
      // φ(∞) is unbounded by construction — no marker.
    }

    // ---- σ-branch points: zeros of φ' inside the relevant disk ----
    const branchPointsOut = [];
    if (QD.findCriticalPoints) {
      let cs;
      try { cs = QD.findCriticalPoints(phi); }
      catch (_) { cs = null; }
      if (cs && cs.points) {
        for (const p of cs.points) {
          // p has { z, w, absZ, inDomain, severity }. Keep critical-set
          // results that are inside or close to the boundary (severities
          // 'critical' and 'near'); skip 'safe' (irrelevant to σ).
          if (p.severity === 'safe') continue;
          branchPointsOut.push({
            w:        p.w,
            label:    'φ\'=0',
            severity: p.severity || 'near',
            absZ:     p.absZ,
          });
        }
      }
    }

    return { poles: polesOut, branchPoints: branchPointsOut };
  }

  // ---------------------------------------------------------------------------
  // 3. computeSigmaLevelCurves — marching squares for |σ| and arg(σ).
  //
  // Algorithm: sample σ on a regular gridSize × gridSize grid spanning the
  // viewport. For each requested isovalue, march each cell (find the up-to-2
  // edge crossings of the iso-line, emit segments). Output: one
  // ContourPolyline per (level, segment) pair: { value, kind: 'abs'|'arg',
  // segments: [{x0,y0,x1,y1}, ...] }.
  //
  // arg(σ) wraps at ±π; marching squares on the raw arg field would create
  // spurious crossings at the seam. We unwrap by working with the principal
  // branch and rejecting cells whose arg-range exceeds π — a coarse but
  // effective filter for the visualization purpose.
  // ---------------------------------------------------------------------------
  function computeSigmaLevelCurves(schwarz, opts) {
    opts = opts || {};
    const gridSize  = (opts.gridSize  != null) ? opts.gridSize  : 128;
    const viewport  = opts.viewport;             // { reMin, reMax, imMin, imMax }
    if (!schwarz || !schwarz.sigma || !viewport) return { abs: [], arg: [] };
    const absLevels = opts.absLevels || _autoAbsLevels(viewport);
    const argLevels = opts.argLevels || [
      -2.617, -1.571, -0.524, 0, 0.524, 1.571, 2.617,        // -5π/6 to 5π/6 / 6
    ];

    // Build the |σ| + arg(σ) field. Cells with σ undefined (outside Ω) get
    // NaN; marching squares skips any cell containing a NaN corner.
    const N = gridSize;
    const dx = (viewport.reMax - viewport.reMin) / (N - 1);
    const dy = (viewport.imMax - viewport.imMin) / (N - 1);
    const absF = new Float64Array(N * N);
    const argF = new Float64Array(N * N);
    for (let iy = 0; iy < N; iy++) {
      const im = viewport.imMin + iy * dy;
      for (let ix = 0; ix < N; ix++) {
        const re = viewport.reMin + ix * dx;
        const w = { re, im };
        if (!schwarz.isInOmega(w)) {
          absF[iy * N + ix] = NaN;
          argF[iy * N + ix] = NaN;
          continue;
        }
        let sv;
        try { sv = schwarz.sigma(w); }
        catch (_) { absF[iy * N + ix] = NaN; argF[iy * N + ix] = NaN; continue; }
        if (!sv || !isFinite(sv.re) || !isFinite(sv.im)) {
          absF[iy * N + ix] = NaN;
          argF[iy * N + ix] = NaN;
          continue;
        }
        absF[iy * N + ix] = Math.hypot(sv.re, sv.im);
        argF[iy * N + ix] = Math.atan2(sv.im, sv.re);
      }
    }

    const absContours = _marchingSquares(absF, N, viewport, absLevels, /*skipWrap*/ false);
    const argContours = _marchingSquares(argF, N, viewport, argLevels, /*skipWrap*/ true);
    for (const c of absContours) c.kind = 'abs';
    for (const c of argContours) c.kind = 'arg';
    return { abs: absContours, arg: argContours };
  }

  function _autoAbsLevels(viewport) {
    // Default geometric ramp of iso-magnitudes — covers typical σ ranges
    // for both bounded (σ ≤ a few) and unbounded (σ can be large) Ω.
    return [0.25, 0.5, 1, 2, 4, 8];
  }

  function _marchingSquares(field, N, viewport, levels, skipWrap) {
    const out = [];
    const reMin = viewport.reMin, imMin = viewport.imMin;
    const dx = (viewport.reMax - viewport.reMin) / (N - 1);
    const dy = (viewport.imMax - viewport.imMin) / (N - 1);
    for (const lv of levels) {
      const segments = [];
      for (let iy = 0; iy < N - 1; iy++) {
        for (let ix = 0; ix < N - 1; ix++) {
          // Corners (BL, BR, TR, TL).
          const i00 = iy * N + ix;
          const i10 = iy * N + ix + 1;
          const i11 = (iy + 1) * N + ix + 1;
          const i01 = (iy + 1) * N + ix;
          const v00 = field[i00], v10 = field[i10],
                v11 = field[i11], v01 = field[i01];
          if (!isFinite(v00) || !isFinite(v10) || !isFinite(v11) || !isFinite(v01)) continue;
          if (skipWrap) {
            const lo = Math.min(v00, v10, v11, v01);
            const hi = Math.max(v00, v10, v11, v01);
            if (hi - lo > Math.PI) continue;          // arg-seam straddling cell
          }
          // Cell type: 4-bit code where bit k is set if corner k > level.
          let code = 0;
          if (v00 > lv) code |= 1;
          if (v10 > lv) code |= 2;
          if (v11 > lv) code |= 4;
          if (v01 > lv) code |= 8;
          if (code === 0 || code === 15) continue;

          // Edge intersection helpers.
          const x0 = reMin + ix * dx;
          const y0 = imMin + iy * dy;
          function lerp(a, b) { return (lv - a) / (b - a); }
          // Edges: 0 (bottom): (x0,y0)→(x0+dx,y0)  v00,v10
          //        1 (right ): (x0+dx,y0)→(x0+dx,y0+dy)  v10,v11
          //        2 (top   ): (x0+dx,y0+dy)→(x0,y0+dy)  v11,v01
          //        3 (left  ): (x0,y0+dy)→(x0,y0)  v01,v00
          function edgePoint(e) {
            if (e === 0) { const t = lerp(v00, v10); return { x: x0 + t*dx, y: y0 }; }
            if (e === 1) { const t = lerp(v10, v11); return { x: x0 + dx,   y: y0 + t*dy }; }
            if (e === 2) { const t = lerp(v11, v01); return { x: x0 + (1-t)*dx, y: y0 + dy }; }
            return { x: x0, y: y0 + (1 - lerp(v01, v00)) * dy };
          }
          // Lookup: each non-zero code emits 1 or 2 line segments. Sat-edge
          // pairs come from the standard MS table.
          const segs = _MS_TABLE[code];
          for (const seg of segs) {
            const a = edgePoint(seg[0]);
            const b = edgePoint(seg[1]);
            segments.push({ x0: a.x, y0: a.y, x1: b.x, y1: b.y });
          }
        }
      }
      out.push({ value: lv, segments });
    }
    return out;
  }

  // Standard marching-squares case table. For ambiguous codes 5 and 10 we
  // pick a consistent diagonal pairing; the result is visually clean for
  // the level-curve overlay even if it's not the unique mathematical
  // resolution.
  const _MS_TABLE = {
    1:  [[0, 3]],
    2:  [[0, 1]],
    3:  [[1, 3]],
    4:  [[1, 2]],
    5:  [[0, 1], [2, 3]],
    6:  [[0, 2]],
    7:  [[2, 3]],
    8:  [[2, 3]],
    9:  [[0, 2]],
    10: [[1, 2], [0, 3]],
    11: [[1, 2]],
    12: [[1, 3]],
    13: [[0, 1]],
    14: [[0, 3]],
  };

  // ---------------------------------------------------------------------------
  // Wire onto QD.Schwarz.
  // ---------------------------------------------------------------------------
  QD.Schwarz.explicitSigmaForm       = explicitSigmaForm;
  QD.Schwarz.findSigmaSingularities  = findSigmaSingularities;
  QD.Schwarz.computeSigmaLevelCurves = computeSigmaLevelCurves;
})();
