// =============================================================================
// ui-strings.js  —  Single source of truth for editable UI prose (QD.Strings).
//
// THIS IS THE FILE TO EDIT when reworking descriptions / helptext / tooltips /
// example blurbs. Each value below is the exact text shown to the user. The code
// reads from here; nothing else needs to change.
//
// How the text reaches the screen:
//   • JS prose (the "?" help-popovers, Faber/oracle help, example blurbs, the
//     solve-failure guidance) is read directly, e.g. QD.Strings.help.cCard.
//   • Static HTML text (the family hints, card hints, overlay tooltips/notes,
//     coachmark) is injected by QD.Strings.apply() into elements that carry a
//     data-str* attribute in index.html:
//        data-str="path"        -> el.textContent = QD.Strings.path
//        data-str-html="path"   -> el.innerHTML   = QD.Strings.path   (tags render)
//        data-str-title="path"  -> el.title       = QD.Strings.path   (PLAIN text)
//     apply() runs once at page load (inline, just after the script loader).
//
// EDITING RULES
//   • help.* / familyHints.* / hints.* / notes.* / faber.help / oracle.help are
//     HTML — <b> <i> <sub> <sup> <code> <br> <kbd> all work; keep entities
//     (&gt; &amp;) as entities.
//   • tooltips.* are PLAIN TEXT (title attributes) — no HTML tags; use Unicode
//     (Fₙ, φ′, ζ², Dₙ) instead.
//   • Preserve the true minus '−' (U+2212), sub/superscripts, and Greek/blackboard
//     glyphs (φ θ κ λ ζ ∂Ω 𝔻 𝔻* ℂ Ω̄) exactly.
//
// NOT centralized here (edit in the named module; see HELPTEXT.md for the full map):
//   • Pure control labels / buttons / options / tab names (index.html).
//   • Value-interleaved micro-labels and computed status lines (ui-solve.js,
//     ui-faber.js property rows, oracle row names in thesis-examples.js).
//   • Family external-field labels and validation errors (ui-modes.js).
//
// After editing ANY app/ file: run `npm run version:sync` (then lint + test).
// =============================================================================

(function (global) {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' && module.exports ? module.exports : (global.QD || (global.QD = {})));

  const S = {

    // ---- "?" help-popovers (HTML) -------------------------------------------
    help: {
      // App title — "What is this?" intro.
      intro:
        `<b>What is this?</b> A <b>quadrature domain</b> Ω is a region where the area integral of any
         analytic function f equals a linear combination of point evaluations of f and its derivatives,
         encoded by the <b>quadrature function h(w)</b>. This tool solves the <i>inverse</i>
         problem: given h(w) it finds the domain Ω and its conformal map φ(z), then analyzes the
         boundary (cusps, curvature, symmetry, accuracy).<br><br>
         <b>To start:</b> pick a <b>Preset</b> or <b>Thesis example</b>, drag the red poles on the plot
         to reshape Ω, and read the solution properties in the panel. The <b>Schwarz dynamics</b> and
         <b>Parameter slice</b> tabs enable further analysis. Press <b>?</b> for shortcuts.`,

      domainType:
        `<b>Domain type.</b> The quadrature identity Ω must satisfy.
         <b>QD</b>: classical (unweighted). <b>PQD</b>: power weight
         |w|<sup>2(α−1)</sup> (α on the card below). <b>LQD</b>: log weight 1/|w|².
         <b>Bounded</b> = finite Ω; <b>unbounded</b> = Ω reaches ∞.
         <b>Singular</b> = 0 ∈ Ω (φ gains a Blaschke factor); non-singular = 0 ∉ Ω.
         Classical QDs have no singular variant.`,

      hCard:
        `<b>Quadrature function h(w).</b> A rational function encoding the quadrature identity.
         Edit poles + residues structurally below, or enter a rational function
         in the textbox at the top. The inverse solver finds Ω whose
         quadrature identity matches h.`,

      mapParams:
        `<b>Map parameters.</b> The scalar knobs of the Riemann map, shown per family.
         <b>PQD power α</b> (PQD modes): the weight is |w|<sup>2(α−1)</sup>; α = 1 is
         classical. <b>Center φ(0)</b> (bounded families): the image of 0 ∈ 𝔻 — a free
         parameter (Manual, or Auto = pole centroid, recomputed as you drag a pole);
         implicit for unbounded families. <b>Residue q at origin</b> (singular LQDs,
         0 ∈ Ω): the residue of the log-weighted Schwarz function at w=0, linked to
         the finite poles and any polynomial part by a closed-form constraint. 
         The unbounded conformal radius c has its own control beside φ(z) in the
         Domain-type card.`,

      cCard:
        `<b>Conformal radius c = φ′(∞).</b> Unbounded QDs form a
         one-parameter family with respect to c — sweep the slider to explore it; past the
         critical c* the simply-connected QD ceases to exist (its boundary forms cusps,
         then self-intersections). <b>Estimate max c</b> finds c* automatically (bracket +
         bisection on the solver's univalence + identity gate), then caps the slider
         at c* and jumps to it — the extremal domain (its boundary just cusps).`,

      solverSettings:
        `<b>Solver settings.</b> The <i>Aggressiveness</i> preset
         (Quick / Standard / Thorough) balances Newton iterations, identity-check
         samples, and how many alternate branches are sought; fine-tune via
         <i>Search options</i>. Also here: the boundary-sample count, the
         vector-field overlay (Pólya h̄(w), or the family-specific external potential),
         the critical-set overlay, and <i>Auto-switch singular ⇄ non-singular PQD</i> 
         — which re-solves in the correct family when a PQD boundary crosses the origin.`,

      searchOptions:
        `<b>Search options.</b> Each phase is a distinct strategy for finding a φ
         consistent with h(w). Direct = single Newton from the initial guess;
         continuation = parameter-homotopy from a related solved scenario;
         multistart = many random seeds; diverse + deflation = explicit
         branch-finding.`,

      status:
        `<b>Status.</b> Live readout of the solver: convergence diagnostics,
         identity residual, univalence, and which branches succeeded.`,

      geom:
        `<b>Geometric properties.</b> Special univalence classes of the solved Ω,
         checked asynchronously after each solve. <i>Star-like</i>: every ray from
         the center (w₀ for bounded; ∞ for unbounded) stays in Ω
         — Re(z·φ′/(φ−w₀)) &gt; 0. <i>Convex</i> (bounded): Re(1 + z·φ″/φ′) &gt; 0.
         <i>Spiral-like</i>: a log-spiral generalization of star-like; λ is the
         optimal spiral angle. The hierarchy is convex ⟹ star-like ⟹ spiral-like,
         all ⟹ univalent.`,

      cusps:
        `<b>Boundary singularities.</b> Cusps of ∂Ω, found asynchronously after each
         solve. A cusp sits where the Riemann map's derivative vanishes on the unit
         circle, φ′(e<sup>iθ</sup>) = 0; the order m of that zero fixes the local
         <i>type</i> (p,q) = (m+1, m+2): m = 1 is the ordinary 3⁄2-power (2,3) cusp.
         A filled ● / magenta triangle marks an actual cusp; a hollow ○ marks an
         <i>incipient</i> one — a φ′-zero near but not yet on ∂𝔻, shown with its
         distance d (a "how close to a cusp" gauge). The (p,q) type is read exactly
         from φ's Taylor coefficients and cross-checked numerically.`,

      alternates:
        `<b>Alternate solutions.</b> When more than one φ satisfies the same h
         (multiple branches), the solver lists them here. Click an alternate to
         promote it to the primary.`,
    },

    // ---- Family hint blocks under the Domain-type card (HTML) ---------------
    familyHints: {
      pqdBounded:
        `<b>Power-weighted QD.</b> Bounded Ω, weight <code>|w|<sup>2(α−1)</sup></code>; <code>φ = (R<sup>#</sup>)<sup>1/α</sup></code>.
        <details><summary>Details</summary>
        <code>∫<sub>Ω</sub> f(w)·|w|<sup>2(α−1)</sup> dA = ∮<sub>∂Ω</sub> f·h dw</code>
        for real <code>α &gt; 0, α ≠ 1</code> (α=1 is the classical case). Riemann map
        factors as <code>φ(z) = (R<sup>#</sup>(z))<sup>1/α</sup></code>, where
        <code>R</code> is rational with poles at the preimages <code>z<sub>j</sub></code> of <code>p<sub>j</sub></code>.
        Realizability constraint: <code>C<sub>j,1</sub></code> must be large enough that
        the preimages stay inside <code>𝔻</code>.
        </details>`,

      pqdBoundedSingular:
        `<b>Bounded singular PQD.</b> <code>0 ∈ Ω</code>; φ gains a Blaschke factor pinned by the mass constraint.
        <details><summary>Details</summary>
        The Riemann map gains a
        Blaschke factor, <code>φ(z) = b<sub>z₀</sub>(z)·(R<sup>#</sup>(z))<sup>1/α</sup></code>,
        where <code>φ(z<sub>0</sub>) = 0</code> (the origin's preimage). No point
        charge: the <code>|w|<sup>2(α−1)</sup></code> weight makes the quadrature
        data unique, and the <b>mass/area constraint</b> (the <code>f=1</code> case
        of the identity) pins <code>|z<sub>0</sub>|</code>. <code>w₀ = φ(0) ≠ 0</code>
        is a non-origin interior point.
        </details>`,

      pqdUnbounded:
        `<b>Unbounded PQD.</b> Ω reaches ∞ (<code>0 ∉ Ω</code>); <code>φ = z·(r<sup>#</sup>)<sup>1/α</sup></code> on <code>𝔻*</code>, conformal radius c.
        <details><summary>Details</summary>
        Ω is unbounded with bounded complement K
        (<code>0 ∉ Ω</code>). The test class is <code>A₀(Ω)</code> (analytic,
        vanishing at ∞), so the weighted integral converges for all α &gt; 0.
        Riemann map <code>φ(z) = z·(r<sup>#</sup>(z))<sup>1/α</sup></code> on
        <code>𝔻*</code> with <code>r<sup>#</sup>(∞) = c<sup>α</sup></code>
        (conformal radius <code>c</code>, a user input). Supports finite poles
        and a polynomial part of any degree (pole at ∞), e.g. monomial
        <code>h = α·k·w<sup>k−1</sup></code> (Thm 4.5.3).
        </details>`,

      pqdUnboundedSingular:
        `<b>Unbounded singular PQD.</b> Ω unbounded with <code>0 ∈ Ω</code>; φ gains a Blaschke factor on <code>𝔻*</code>.
        <details><summary>Details</summary>
        The Riemann map gains a Blaschke factor,
        <code>φ(z) = z·b<sub>z₀</sub>(z)·(r<sup>#</sup>(z))<sup>1/α</sup></code> on
        <code>𝔻*</code>, where <code>z₀ ∈ 𝔻*</code> is the origin's preimage
        (<code>φ(z₀)=0</code>) and <code>r<sup>#</sup>(∞)=|cz₀|<sup>α</sup></code>.
        No point charge; z₀ is pinned by <code>r(z₀)=0</code> (Prop 4.6.3, h
        analytic at 0). c is the conformal radius.
        </details>`,

      lqdBounded:
        `<b>Log-weighted QD.</b> Weight <code>1/|w|²</code>, non-singular (<code>0 ∉ Ω̄</code>).
        <details><summary>Details</summary>
        <code>∫<sub>Ω</sub> f(w)/|w|² dA = ∮<sub>∂Ω</sub> f·h dw</code>.
        Non-singular only: requires <code>0 ∉ Ω̄</code>, so <code>w₀ = φ(0)</code> must be nonzero.
        </details>`,

      lqdSingular:
        `<b>Singular log-weighted QD.</b> <code>0 ∈ Ω</code>; residue q at the origin; <code>φ = γ·b<sub>z₀</sub>·exp(r#)</code>.
        <details><summary>Details</summary>
        Test functions
        <code>f ∈ L¹<sub>a</sub>(Ω; ρ₀)</code> automatically vanish at 0. Riemann map
        factors as <code>φ(z) = γ · b<sub>z₀</sub>(z) · exp(r#(z))</code>. The pole of
        <code>h</code> at the origin contributes residue <code>q</code> (complex);
        <code>q = 0</code> is the degenerate case with no log-charge.
        </details>`,

      lqdUnbounded:
        `<b>Unbounded log-weighted QD.</b> <code>0 ∉ Ω̄</code>, <code>∞ ∈ Ω</code>; <code>φ = c·z·exp(r#−r#(∞))</code> on <code>𝔻*</code>.
        <details><summary>Details</summary>
        Test functions <code>f ∈ L¹<sub>a</sub>(Ω; ρ₀)</code> vanish at <code>∞</code>.
        Riemann map: <code>φ(z) = c·z · exp(r#(z) − r#(∞))</code> on <code>𝔻*</code>,
        with conformal radius <code>c = φ′(∞)</code>. The <code>− r#(∞)</code>
        subtraction pins the leading coefficient at ∞ to exactly <code>c</code>.
        </details>`,

      lqdUnboundedSingular:
        `<b>Unbounded singular log-weighted QD.</b> <code>0 ∈ Ω</code> and <code>∞ ∈ Ω</code>; Blaschke factor + residue q.
        <details><summary>Details</summary>
        Test functions <code>f</code> vanish at BOTH 0 and ∞ (e.g.
        <code>w/(w−b)<sup>k</sup></code> for <code>k ≥ 2</code>, <code>b ∈ K</code>).
        Riemann map: <code>φ(z) = c·|z₀|·z·b<sub>z₀</sub>(z) · exp(r#(z) − r#(∞))</code>
        with <code>z₀ ∈ 𝔻*</code> the preimage of 0, and <code>q</code> (complex)
        the residue of <code>h</code> at the origin. <code>q = 0</code> is the
        degenerate case with no log-charge.
        </details>`,
    },

    // ---- Card hints (HTML) --------------------------------------------------
    hints: {
      hSum:
        `h(w) = Σ<sub>j</sub> Σ<sub>s=1..m<sub>j</sub></sub> C<sub>j,s</sub> / (w − a<sub>j</sub>)<sup>s</sup>.
         Enter complex values as e.g. <code>1+2i</code>.`,

      newHere:
        `New here? Pick a <strong>Preset</strong> below to start, drag poles
         directly on the plot to reposition them, and press <kbd>?</kbd> for
         keyboard shortcuts.`,

      polyPart:
        `Polynomial part of h(w) = Σ<sub>l=0..m<sub>∞</sub></sub> C<sub>∞,l</sub> w<sup>l</sup>.
         Quadrature nodes at infinity.`,

      alpha:
        `Weight is <code>|w|<sup>2(α−1)</sup></code> for any real <code>α &gt; 0</code>
         (<code>α ≠ 1</code>). <code>α &gt; 1</code>: weight vanishes at 0;
         <code>α &lt; 1</code>: weight blows up at 0 (the LQD-limit regime,
         <code>→ |w|<sup>−2</sup></code> as <code>α → 0⁺</code>). <code>α = 1</code>
         is classical QD — use the Classical group above.`,

      qResidue:
        `Complex residue of <code>h</code> at the simple pole at <code>w = 0</code>
         (the singular pole inside Ω). Dial via text or the |q| / arg(q) sliders.`,

      faberCard:
        `Faber polynomials F<sub>n</sub>(ζ) of the bounded complement K = ℂ∖Ω, read off φ's
         Laurent expansion at ∞. Their roots cluster inside K.`,

      searchOptions:
        `Overrides for the aggressiveness preset. Leave a field blank to use
         the preset's default. Phase checkboxes are enabled by default.`,

      directProblem:
        `Given a Riemann map <code>φ : 𝔻 → Ω</code>, compute the quadrature
         function <code>h</code> for which <code>Ω ∈ QD(h)</code>. Supports
         bounded classical QD (polynomial and rational <code>φ</code>),
         unbounded classical QD (Laurent at ∞), and a numerical-fallback
         mode for arbitrary expressions.`,

      schwarzDynamics:
        `For a Riemann map <code>φ : 𝔻 → Ω</code>, the Schwarz reflection
         <code>σ(w) = conj(F(ψ(w)))</code> is anti-meromorphic in Ω and fixes
         <code>∂Ω</code>. This tab colors each <code>w ∈ Ω</code> by the
         smallest <code>n</code> for which <code>σⁿ(w) ∈ Ω^c</code> (the
         "fundamental tile"). Points whose orbit diverges to ∞ before
         re-entering <code>Ω^c</code> are flagged as the escaping set.`,

      coachmark:
        `<strong>On the plot:</strong> drag poles to reshape · double-click to add a pole · scroll to zoom · drag empty space to pan`,
    },

    // ---- Tooltips (title="" — PLAIN TEXT, no HTML) --------------------------
    tooltips: {
      weightQD: `Classical quadrature domain — unweighted area integral ∫_Ω f dA.`,
      weightPQD: `Power-weighted QD — weight |w|^(2(α−1)) in the area integral.`,
      weightLQD: `Log-weighted QD — weight 1/|w|² in the area integral.`,
      riemannSym: `Show / hide the symbolic form`,
      thesisSelect: `Curated canonical quadrature domains, each with an analytic oracle. Loading one frames the view, turns on the annotated-phenomena overlay, and shows a computed-vs-expected oracle card.`,
      estimateMaxC: `Estimate the critical conformal radius c* — the largest c with a valid unbounded QD. Sets the slider to c* (the extremal, just-cusped domain).`,
      criticalPoints: `Plot the w-plane images of {z : φ'(z) = 0}. Critical points inside the relevant disk (𝔻 for bounded, 𝔻* for unbounded) predict univalence loss; those near |z|=1 predict imminent degeneracy.`,
      curvature: `Color ∂Ω by its curvature |κ| (cool → hot). The hottest stretch marks the sharpest bend; κ → ∞ at a cusp. Reads the same data as the Geometry & accuracy panel.`,
      phenomena: `Annotate the phenomena the cusp / critical-set overlays don't: the harmonic-measure hot spot (the tip, where ρ = 1/(2π|φ′|) peaks), the maximum-curvature point on ∂Ω, and the domain's symmetry axes (dashed) with its dihedral (Dₙ) or cyclic (Zₙ) symmetry group.`,
      faberRoots: `Plot the roots of the Faber polynomials Fₙ of the bounded complement K (classical unbounded QD only). Roots cluster inside K. Teal circles = union of all roots up to N; violet diamonds = the single selected Fₙ. Drive N / the mode from the Faber polynomials card.`,
      fit: `Fit view to data`,
      reset: `Reset view`,
      dock: `Dock the panel into the sidebar (clear the plot)`,
      collapse: `Collapse panel`,
    },

    // ---- Overlay legend notes (.ov-note — HTML) ----------------------------
    notes: {
      criticalPoints: `zeros of φ′ — red inside, orange near ∂Ω`,
      curvature: `|κ| on ∂Ω, blue → red (→ Geometry &amp; accuracy)`,
      phenomena: `tip · max-κ · symmetry axes`,
      faberRoots: `roots of Fₙ — cluster in K (UQD only)`,
    },

    // ---- Faber-polynomials card (HTML help + whole-sentence messages) -------
    faber: {
      help:
        `Faber polynomials F_n(ζ) of the bounded complement K = ℂ∖Ω. For a classical ` +
        `unbounded QD the solved map φ is the EXTERIOR map of K, so the F_n are read off ` +
        `φ(z) = c·z + c₀ + c₁/z + … directly. Their roots cluster inside K ` +
        `— tick "Plot roots on domain" to overlay them. The capacity ` +
        `cap(K) = c = φ′(∞). High orders are ill-conditioned; non-convergence is flagged.`,
      pending: `solving… Faber analysis pending`,
      unavailablePrefix: `Faber analysis unavailable: `,
    },

    // ---- Analytic-oracle card (HTML help + whole-sentence messages) --------
    oracle: {
      help:
        `Curated canonical quadrature domains, each with an ANALYTIC ORACLE — the ` +
        `closed-form quantities a correct solve must reproduce (area, symmetry, cusps, ` +
        `c*, accuracy). Rows show computed vs expected: ✓ pass, ⚠ marginal, ✗ off. The ` +
        `c* row is verified on demand (it runs the conformal-radius estimator).`,
      pickExample: `— pick an example —`,
      pending: `solving… oracle pending`,
      noValid: `no valid solution to check against the oracle`,
      matches: `✓ matches the analytic oracle`,
      someDiffer: `⚠ some rows differ from the oracle`,
      verifyCmax: `Verify c* (slow)`,
      estimating: `estimating c*…`,
    },

    // ---- Thesis-example blurbs (HTML) --------------------------------------
    blurbs: {
      disk:
        `The simplest QD: h = 1/w gives the unit disk. Area π, constant curvature, ` +
        `uniform harmonic measure — an exact analytic baseline.`,
      twoPointSym:
        `h = 1.5/(w−1) + 1.5/(w+1): two equal poles at ±1. The domain is a symmetric ` +
        `oval with D₂ symmetry (a half-turn + two mirror axes).`,
      triangle:
        `Three equal poles at the cube roots of unity → a rounded triangle with D₃ ` +
        `symmetry (3-fold rotation + 3 mirror axes).`,
      square4pole:
        `Four equal poles at ±1, ±i → a rounded square with D₄ symmetry.`,
      cardioidUnbounded:
        `h = 1.5/w + 0.5/w² on the exterior. As the conformal radius c grows the ` +
        `boundary cusps at c* ≈ 1.46 (a Hele-Shaw / Polubarinova–Galin blow-up).`,
      deltoidUnbounded:
        `h = w² on the exterior → the deltoid (3-cusp hypocycloid). The QD branch ends ` +
        `at c* ≈ 0.5, where three Z₃-symmetric cusps form simultaneously.`,
      singlePoleUnbounded:
        `h = 1/(w−2): one simple pole off the origin. A smooth, mirror-symmetric ` +
        `exterior QD — the canonical one-pole unbounded example.`,
    },

    // ---- Solve-failure guidance (plain text) -------------------------------
    guidance: {
      noSolutionPrefix: `No quadrature domain found. Suggestions: `,
      tryHarder: `try the “Try harder (exhaustive search)” button or raise Aggressiveness`,
      lqd: `this h may have no log-weighted QD — try smaller residues, or adjust c`,
      pqd: `PQDs need a large-enough residue and an interior w₀ — try a bigger |C| or move w₀`,
      poles: `move poles away from each other and the boundary, or adjust residue magnitudes`,
    },
  };

  // Resolve a dotted path ("help.cCard") against S.
  function get(path) {
    if (!path) return undefined;
    return path.split('.').reduce((o, k) => (o == null ? undefined : o[k]), S);
  }

  // Populate every [data-str*] element from S. Safe to call multiple times.
  function apply(root) {
    if (typeof document === 'undefined') return;
    root = root || document;
    const each = (attr, set) => {
      const els = root.querySelectorAll('[' + attr + ']');
      for (let i = 0; i < els.length; i++) {
        const key = els[i].getAttribute(attr);
        const val = get(key);
        if (val == null) {
          if (typeof console !== 'undefined') console.warn('[ui-strings] missing key: ' + key);
          continue;
        }
        set(els[i], val);
      }
    };
    each('data-str',       (el, v) => { el.textContent = v; });
    each('data-str-html',  (el, v) => { el.innerHTML = v; });
    each('data-str-title', (el, v) => { el.setAttribute('title', v); });
  }

  S.get = get;
  S.apply = apply;
  QD.Strings = S;

})(typeof globalThis !== 'undefined' ? globalThis : this);
