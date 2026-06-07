// =============================================================================
// direct-ui.js -- Direct-problem tab UI.
//
// A compact segmented "Domain type" control (mirroring the inverse tab's) sets
// the family along three axes — Weight {classical | power(PQD) | log(LQD)} ×
// Domain {bounded | unbounded | numerical} × singular(0∈Ω) — stored UNIFIED on
// directState.{weight, mode, singular}. The three φ-input cards (one visible per
// Domain) hold the φ / kernel coefficients and the weight PARAMETER inputs
// (α / w₀ / z₀ / c). applyDirectMode() is the single canonical refresh: it sets
// φ-card + param-row visibility and syncs the segmented control (INVARIANT: any
// directState.mode/weight/singular change must go through it).
//
// The three Domains share the same plot canvas:
//   1. Bounded — polynomial or rational φ (or, weighted, the rational KERNEL
//      R#/r#); structured coefficient fields + a live-parsing "paste expression".
//   2. Unbounded — Laurent-at-∞ φ = c·z + Σ F_l/z^l (or the weighted KERNEL r#);
//      structured (c, F_l) coefficient fields / a kernel paste field.
//   3. Numerical — free-form math.js expression in z (classical only); DFT-
//      extracted polynomial truncation produces an approximate h with an
//      analyticity diagnostic.
//
// Each mode pushes:
//   • h to QD.Direct._sendHToInverseTab     (pre-fill QD tab inverse view and switch)
//   • ∂Ω points to QD.Direct._setPlotBoundary   (live boundary preview)
// Both hooks are installed by ui.js after DOMContentLoaded. If the hooks
// aren't installed yet (race), we no-op and re-try on tab swap.
//
// State is local to this module (directState). The shared canvas + the
// DomainPlot renderer are owned by ui.js.
// =============================================================================
'use strict';

(function () {

  // ---------------------------------------------------------------------------
  // State
  // ---------------------------------------------------------------------------
  const directState = {
    // 'bounded' or 'unbounded'.
    mode: 'bounded',

    // Bounded mode:
    //   'polynomial' kind  → φ(z) = Σ c_k z^k          (uses coeffs)
    //   'rational' kind    → φ(z) = N(z) / D(z)        (uses coeffsNum, coeffsDen)
    coeffsKind: 'polynomial',
    coeffs: ['0', '1'],                                   // polynomial path
    coeffsNum: ['0', '1'],                                // rational numerator (default = z)
    coeffsDen: ['1'],                                     // rational denominator (default = 1)

    // WEIGHT × SINGULAR — UNIFIED across bounded + unbounded (one axis each, like
    // the inverse tab). Live in the "Domain type" card's segmented control, not in
    // the φ cards. weight ∈ {'classical','power','log'}; singular = 0 ∈ Ω (Blaschke
    // b_{z₀}; power/log only). Forward-kernel meaning of weight:
    //   bounded:    power → φ=(R#)^{1/α}; log → φ=w₀·exp(r#)
    //   unbounded:  power → φ=z·(r#)^{1/α} (c derived from r#); log → φ=c·z·exp(r#)
    // weight≠classical forces coeffsKind='rational' (the input is the KERNEL R#/r#).
    weight: 'classical',
    singular: false,
    // Bounded weight PARAMETERS (stay in the bounded φ card):
    alpha: '2',                    // PQD exponent (power weight |w|^{2(α−1)}), α>0, ≠1
    logW0: '2',                    // LQD gauge w₀ = φ(0)
    z0: '0.3',                     // preimage of the origin (|z₀| < 1), singular only

    // Unbounded mode: φ(z) = c·z + Σ_l F_l/z^l. Strings.
    cValue: '1',                   // conformal radius (positive real)
    Fcoeffs: [],                   // [F_0, F_1, ..., F_{m-1}], strings; empty ⇒ φ = c·z

    // Unbounded weight PARAMETERS (stay in the unbounded φ card; ∞∈Ω, Thm 4.3.7):
    unsAlpha: '2',                 // PQD exponent
    unsZ0: '1.3',                  // origin preimage |z₀|>1 (LQD-singular: free; PQD-singular: hint)
    unsKernelNum: ['0.81', '-1.725'], // r# numerator (strings; default = a valid non-sing PQD kernel, c=0.9)
    unsKernelDen: ['1', '-2.5'],   // r# denominator (strings)
    unsKernelExpr: '',             // last kernel paste expression (rational in z)
    lastWeight: 'classical',       // weight that produced lastH (for Send/Verify dispatch)

    // Numerical mode: any math.js expression in z.
    numExpr: 'z + 0.2*sin(z)',     // default: a non-polynomial example
    numMaxOrder: 12,

    // Last successfully computed h, and the c that produced it.
    lastH: null,
    lastC: 1,
    // Last bounded-weighted solve metadata (set by recomputeBounded; used by the
    // Verify button and the origin-term display):
    lastPhi: null,                 // the built φ (for the family identity verifier)
    lastSingular: false,           // whether the last solve was a singular (0∈Ω) weighted QD
    lastQ: null,                   // LQD-singular origin residue q  (h has + q/w)
    lastOriginRes: null,           // PQD-singular origin residue r₀ = ∫|w|^{2(α−1)}dA − ΣC (h has + r₀/w)

    // True when the user's most recent action was typing in the paste field
    // (so we don't clobber what they're typing with auto-regenerated form).
    expressionInput: false,
  };

  // ---------------------------------------------------------------------------
  // Presets
  // ---------------------------------------------------------------------------
  const PHI_PRESETS_BOUNDED = [
    // Polynomial
    { id: 'unit-disk',     label: 'Unit disk:  φ = z',                            kind: 'polynomial', coeffs: ['0', '1'] },
    { id: 'shifted-disk',  label: 'Shifted disk:  φ = (1+i) + 2z',                kind: 'polynomial', coeffs: ['1+i', '2'] },
    { id: 'tilted-disk',   label: 'Tilted disk:  φ = (1+i)·z',                    kind: 'polynomial', coeffs: ['0', '1+i'] },
    { id: 'quadratic',     label: 'Quadratic:  φ = z + 0.1·z²',                   kind: 'polynomial', coeffs: ['0', '1', '0.1'] },
    { id: 'cubic',         label: 'Cubic:  φ = z + 0.1·z² − 0.05·z³',             kind: 'polynomial', coeffs: ['0', '1', '0.1', '-0.05'] },
    { id: 'cassini-ish',   label: 'Smoothed Cassini:  φ = z + 0.2·z³',            kind: 'polynomial', coeffs: ['0', '1', '0', '0.2'] },
    // Rational
    { id: 'mobius',        label: 'Möbius:  φ = z / (1 − 0.3z)',                  kind: 'rational',
      num: ['0', '1'], den: ['1', '-0.3'] },
    { id: 'mobius-2',      label: 'Two-pole:  φ = z / ((1 − 0.3z)(1 − 0.4z))',    kind: 'rational',
      num: ['0', '1'], den: ['1', '-0.7', '0.12'] },
    { id: 'shifted-rat',   label: 'Shifted rational:  φ = (z + 0.5i) / (1 − 0.4z)', kind: 'rational',
      num: ['0.5i', '1'], den: ['1', '-0.4'] },
    { id: 'repeated',      label: 'Repeated pole:  φ = z / (1 − 0.3z)²',          kind: 'rational',
      num: ['0', '1'], den: ['1', '-0.6', '0.09'] },
  ];

  const PHI_PRESETS_UNBOUNDED = [
    { id: 'ext-unit',      label: 'Exterior of unit disk:  φ = z',                c: '1',   F: [] },
    { id: 'ext-r2',        label: 'Exterior of disk r=2:  φ = 2z',                c: '2',   F: [] },
    { id: 'ext-shifted',   label: 'Ext of disk r=1.5 at 1+i:  φ = 1.5z + (1+i)',  c: '1.5', F: ['1+i'] },
    { id: 'ext-tilted',    label: 'Ext of disk r=0.5 at -2−i:  φ = 0.5z + (−2−i)', c: '0.5', F: ['-2-i'] },
    // Higher-Laurent example (non-QD generically; for exploration only):
    { id: 'ellipse-like',  label: '(Not a classical QD) φ = z + 0.3/z',           c: '1',   F: ['0', '0.3'] },
  ];

  // ---------------------------------------------------------------------------
  // Mount API (HANDOFF #30): the Direct UI is no longer a stand-alone tab.
  // It's mounted into a sub-container of the QD tab and activated by ui.js's
  // setViewMode('direct'). The host calls QD.Direct._mountUI() on the first
  // switch to direct view (idempotent) and QD.Direct._activate() each time
  // the user toggles back to direct view (re-pushes φ boundary to the canvas).
  // ---------------------------------------------------------------------------
  let mounted = false;
  // Forward bindings for the Phase-3 extracted modules (assigned by the install
  // after makeOutputCard below; called by name from card handlers + _activate).
  let recomputeAndRender, runVerify;
  function _mountUI() {
    if (mounted) return;
    mountDirectSidebar();
    mounted = true;
  }
  function _activate() {
    if (!mounted) return;
    recomputeAndRender();
  }
  // Expose on the QD.Direct namespace (created by direct-common.js, populated
  // here at script-evaluation time so ui.js can call without a load-order race).
  if (typeof QD !== 'undefined') {
    QD.Direct = QD.Direct || {};
    QD.Direct._mountUI   = _mountUI;
    QD.Direct._activate  = _activate;
  }

  function mountDirectSidebar() {
    const root = document.getElementById('controls-direct');
    if (!root) return;
    root.innerHTML = '';
    root.appendChild(makeDomainTypeCard());
    root.appendChild(makePhiCardBounded());      // initially visible
    root.appendChild(makePhiCardUnbounded());    // initially hidden
    root.appendChild(makePhiCardNumerical());    // initially hidden
    root.appendChild(makeOutputCard());
    applyDirectMode();
    attachDirectHelp();      // HANDOFF #33
  }

  function attachDirectHelp() {
    if (!window.QD || !window.QD.QoL || !window.QD.QoL.attachHelp) return;
    const H = window.QD.QoL.attachHelp;
    const root = document.getElementById('controls-direct');
    if (!root) return;
    const cards = root.querySelectorAll('section.card');
    if (cards[0]) H(cards[0].querySelector('h2'),
      `<b>Domain type.</b> Bounded: φ(z) is a polynomial mapping 𝔻 → Ω.
       Unbounded: φ(z) = c·z + lower-order terms mapping 𝔻* → Ω. Numerical:
       paste any analytic-in-𝔻 expression; the kernel infers an order-N
       polynomial via DFT and reports a non-analyticity diagnostic.`);
    // Three φ cards (only one visible at a time) all get the same help.
    for (let i = 1; i <= 3; i++) {
      if (cards[i]) H(cards[i].querySelector('h2'),
        `<b>Riemann map φ(z).</b> Edit coefficients (or paste a math.js
         expression). The Direct kernel computes the explicit
         quadrature data h(w) such that ∮<sub>∂Ω</sub> f h dw = ∮<sub>∂Ω</sub>
         f dψ for analytic test functions f.`);
    }
    if (cards[4]) H(cards[4].querySelector('h2'),
      `<b>Output h(w).</b> The computed quadrature data, with a copy
       button and a "Send to inverse" affordance that pushes the data
       into the QD tab and runs the inverse solver — verifies the
       direct kernel by round-tripping.`);
  }

  // CANONICAL Direct-tab mode refresh (mirror of ui.js applyModeVisuals). INVARIANT:
  // any code that changes directState.mode / weight / singular MUST call this — it
  // is the single source of truth for φ-card visibility, the weight-dependent
  // param-row visibility inside the active φ card, and the Domain-type segmented
  // control + singular-checkbox sync.
  function applyDirectMode() {
    const root = document.getElementById('controls-direct');
    if (!root) return;
    const mode = directState.mode;                 // 'bounded' | 'unbounded' | 'numerical'
    const numerical = (mode === 'numerical');
    // Numerical is classical free-form only — it carries no weight/singular.
    const weight = numerical ? 'classical' : directState.weight;
    const singular = directState.singular && (weight === 'power' || weight === 'log');

    const show = (sel, on) => { const el = root.querySelector(sel); if (el) el.style.display = on ? '' : 'none'; };
    // φ-card visibility.
    show('.dir-phi-card-bounded',   mode === 'bounded');
    show('.dir-phi-card-unbounded', mode === 'unbounded');
    show('.dir-phi-card-numerical', numerical);

    // Domain-type segmented control + singular checkbox sync.
    root.querySelectorAll('#dir-dm-weight .seg-btn').forEach(b => {
      b.classList.toggle('active', b.dataset.weight === weight);
      b.disabled = numerical;                       // weight is meaningless for numerical
    });
    root.querySelectorAll('#dir-dm-domain .seg-btn').forEach(b => b.classList.toggle('active', b.dataset.domain === mode));
    const sing = root.querySelector('#dir-dm-singular');
    if (sing) {
      const allowed = !numerical && weight !== 'classical';
      sing.checked = singular;
      sing.disabled = !allowed;
      if (sing.parentElement) sing.parentElement.style.opacity = allowed ? '' : '0.45';
    }

    // Weight-dependent param rows inside the (bounded / unbounded) φ cards.
    applyBoundedWeightRows(root, weight, singular);
    applyUnboundedWeightRows(root, weight, singular);
  }

  // Bounded φ card: show the α (power) / w₀ (log) / z₀ (singular) param rows + the
  // kernel-input hint, driven by the unified weight/singular.
  function applyBoundedWeightRows(root, weight, singular) {
    const show = (sel, on) => { const el = root.querySelector(sel); if (el) el.style.display = on ? '' : 'none'; };
    show('.dir-phi-alpha-row', weight === 'power');
    show('.dir-phi-logw0-row', weight === 'log');
    show('.dir-phi-z0-row', (weight === 'power' || weight === 'log') && singular);
    const hintB = root.querySelector('.dir-phi-weight-hint');
    if (hintB) {
      if (weight === 'power') {
        hintB.style.display = '';
        hintB.innerHTML = singular
          ? 'Enter the rational kernel <strong>R#(z)</strong> + the origin preimage <strong>z₀</strong> (|z₀|&lt;1); the domain is the SINGULAR power-weighted QD (0 ∈ Ω) with φ = b<sub>z₀</sub>·(R#)<sup>1/α</sup> (realizable ⟺ φ univalent). h gains an origin term r₀/w (shown below).'
          : 'Enter the rational kernel <strong>R#(z)</strong>; the domain is the power-weighted QD with φ = (R#)<sup>1/α</sup> (weight |w|<sup>2(α−1)</sup>). R# must be analytic and non-vanishing on 𝔻̄.';
      } else if (weight === 'log') {
        hintB.style.display = '';
        hintB.innerHTML = singular
          ? 'Enter the rational kernel <strong>r#(z)</strong> + the origin preimage <strong>z₀</strong> (|z₀|&lt;1); the domain is the SINGULAR log-weighted QD (0 ∈ Ω) with φ = γ·b<sub>z₀</sub>·exp(r#), γ = w₀/|z₀|. h gains an origin pole q/w (shown below).'
          : 'Enter the rational kernel <strong>r#(z)</strong>; the domain is the log-weighted QD with φ = w₀·exp(r#). r# must be analytic on 𝔻̄.';
      } else {
        hintB.style.display = 'none';
      }
    }
  }

  // Unbounded φ card: classical-Laurent vs weighted-kernel blocks, α/c/z₀ rows + hint.
  function applyUnboundedWeightRows(root, weight, singular) {
    const weighted = (weight === 'power' || weight === 'log');
    const show = (sel, on) => { const el = root.querySelector(sel); if (el) el.style.display = on ? '' : 'none'; };
    show('.dir-phi-uns-classical', !weighted);
    show('.dir-phi-uns-kernel', weighted);
    show('.dir-phi-uns-alpha-row', weight === 'power');
    // c is a user input for classical + log; DERIVED for power (φ'(∞) from r#).
    show('.dir-phi-uns-c-row', weight !== 'power');
    show('.dir-phi-uns-z0-row', weighted && singular);
    const hintB = root.querySelector('.dir-phi-uns-weight-hint');
    if (hintB) {
      if (weight === 'power') hintB.innerHTML = singular
        ? 'SINGULAR unbounded power QD (0∈Ω): φ = z·b<sub>z₀</sub>·(r#)<sup>1/α</sup>. z₀ is pinned by r(z₀)=0 (= 1/conj of a zero of r#); the z₀ field selects among zeros. No origin term.'
        : 'Unbounded power QD: φ = z·(r#)<sup>1/α</sup>, φ′(∞)=c derived from r#. h = finite poles + polynomial-at-∞.';
      else if (weight === 'log') hintB.innerHTML = singular
        ? 'SINGULAR unbounded log QD (0∈Ω): φ = c·|z₀|·z·b<sub>z₀</sub>·exp(r#). z₀ (|z₀|>1) is a free input; h gains an origin pole q/w (shown below).'
        : 'Unbounded log QD: φ = c·z·exp(r#), c = φ′(∞) (input). h = finite poles + polynomial-at-∞.';
    }
  }

  // ---------------------------------------------------------------------------
  // Domain-type card — compact segmented Weight × Domain × singular control,
  // styled like the inverse tab's #dm-* control (shared .segmented / .seg-btn CSS).
  // ---------------------------------------------------------------------------
  function makeDomainTypeCard() {
    const card = section('Domain type', `
      <div class="row" style="margin-bottom:6px; align-items:center;">
        <span class="domain-mode-group-label" style="width:52px; margin:0;">Weight</span>
        <div class="segmented" id="dir-dm-weight" role="group" aria-label="Weight class">
          <button type="button" class="seg-btn" data-weight="classical">QD</button>
          <button type="button" class="seg-btn" data-weight="power">PQD</button>
          <button type="button" class="seg-btn" data-weight="log">LQD</button>
        </div>
      </div>
      <div class="row" style="margin-bottom:6px; align-items:center;">
        <span class="domain-mode-group-label" style="width:52px; margin:0;">Domain</span>
        <div class="segmented" id="dir-dm-domain" role="group" aria-label="Domain extent">
          <button type="button" class="seg-btn" data-domain="bounded">Bounded</button>
          <button type="button" class="seg-btn" data-domain="unbounded">Unbounded</button>
          <button type="button" class="seg-btn" data-domain="numerical">Numerical</button>
        </div>
      </div>
      <label class="row" style="margin:0 0 2px; gap:6px;">
        <input type="checkbox" id="dir-dm-singular"> singular (0 ∈ Ω)
      </label>
      <div class="hint" style="margin-top: 6px;">
        <strong>Weight</strong>: QD (classical) · PQD (power |w|<sup>2(α−1)</sup>) · LQD (log).
        <strong>Domain</strong>: Bounded (𝔻→Ω) · Unbounded (∞∈Ω) · Numerical (free-form φ, classical only).
        <strong>singular</strong>: 0 ∈ Ω (Blaschke factor; PQD/LQD).
      </div>
    `);
    card.querySelectorAll('#dir-dm-weight .seg-btn').forEach(b => b.addEventListener('click', () => {
      if (directState.mode === 'numerical') return;          // weight locked for numerical
      directState.weight = b.dataset.weight;
      // Weighted kernels need a rational kernel: force rational kind on the bounded card.
      if (directState.weight !== 'classical' && directState.coeffsKind !== 'rational') {
        directState.coeffsKind = 'rational';
        directState.expressionInput = false;
        const bcard = document.querySelector('.dir-phi-card-bounded');
        if (bcard) { renderCoeffFields(bcard); setPasteFromCoeffs(bcard); }
      }
      applyDirectMode();
      recomputeAndRender();
    }));
    card.querySelectorAll('#dir-dm-domain .seg-btn').forEach(b => b.addEventListener('click', () => {
      directState.mode = b.dataset.domain;
      directState.expressionInput = false;
      applyDirectMode();
      recomputeAndRender();
    }));
    const sing = card.querySelector('#dir-dm-singular');
    if (sing) sing.addEventListener('change', () => {
      directState.singular = sing.checked;
      applyDirectMode();
      recomputeAndRender();
    });
    return card;
  }

  // ---------------------------------------------------------------------------
  // φ-input card
  // ---------------------------------------------------------------------------
  // Two equivalent input paths kept in sync:
  //
  //   • Structured: cₖ text fields and ± degree buttons (one row per coeff).
  //   • Paste: a single math.js expression in z; parsed in real time
  //     (debounced) on every keystroke; success → fields populate.
  //
  // The "expressionInput" flag remembers which side last accepted user input.
  // When structured fields change, we regenerate the canonical paste string
  // (unless the user is actively typing in the paste field). When the paste
  // field successfully parses, we populate the structured fields.
  // ---------------------------------------------------------------------------
  function makePhiCardBounded() {
    // Marker class so applyDirectMode() can show/hide this card.
    const card = section('Riemann map φ(z) — bounded', `
      <div class="hint">
        Polynomial φ(z) = Σ<sub>k=0..n</sub> c<sub>k</sub> z<sup>k</sup>, with
        c<sub>0</sub> = φ(0) = w₀ and c<sub>1</sub> ≠ 0. Complex literals like
        <code>1+2i</code>, <code>i</code>, <code>-0.5i</code> are accepted.
      </div>
      <div class="row" style="margin-bottom: 8px;">
        <label>Preset:
          <select class="dir-phi-preset" style="width: 280px;">
            <option value="">— custom —</option>
            ${PHI_PRESETS_BOUNDED.map(p => `<option value="${p.id}">${p.label}</option>`).join('')}
          </select>
        </label>
      </div>

      <!-- Weight PARAMETER inputs (the Weight/singular selectors live in the
           Domain-type card; visibility driven by applyBoundedWeightRows). -->
      <div class="row" style="margin-bottom: 6px; gap: 10px; align-items: center;">
        <label class="dir-phi-alpha-row" style="display:none;">α =
          <input type="text" class="dir-phi-alpha" value="2" style="width: 56px; font-family: ui-monospace, monospace;">
        </label>
        <label class="dir-phi-logw0-row" style="display:none;">w₀ =
          <input type="text" class="dir-phi-logw0" value="2" style="width: 70px; font-family: ui-monospace, monospace;">
        </label>
        <label class="dir-phi-z0-row" style="display:none;">z₀ =
          <input type="text" class="dir-phi-z0" value="0.3" style="width: 70px; font-family: ui-monospace, monospace;">
        </label>
      </div>
      <div class="dir-phi-weight-hint hint" style="display:none; margin-bottom: 6px;"></div>

      <!-- Expression input -->
      <div class="row" style="margin-bottom: 4px; align-items: stretch;">
        <label style="flex: 1 1 auto; display: flex; align-items: center; gap: 6px;">
          <span style="white-space: nowrap;">φ(z) =</span>
          <input type="text" class="dir-phi-paste"
                 placeholder="e.g. z + 0.1*z^2 - 0.05i*z^3  |  (z+1)^2 - 1  |  (1+i)*z"
                 style="flex: 1 1 auto; min-width: 200px; font-family: ui-monospace, monospace;">
        </label>
        <span class="dir-phi-status" style="margin-left: 6px; font-size: 16px; line-height: 1.6;" aria-live="polite"></span>
      </div>
      <div class="dir-phi-paste-msg" style="font-size: 11px; min-height: 1.2em; margin-bottom: 6px;"></div>
      <details style="margin-bottom: 8px;">
        <summary style="cursor: pointer; font-size: 11px; color: #5677a8;">Supported grammar</summary>
        <div class="hint" style="margin-top: 4px;">
          <code>z</code>, <code>i</code>, real / complex literals, operators
          <code>+ − * /</code>, parentheses, <code>^</code> with an integer
          exponent, and function calls (e.g. <code>exp(1+i)</code>) whose
          arguments are pure constants. <strong>Rational</strong> expressions
          like <code>z/(1-0.3z)</code> or <code>z/2 + 1/(z+2)</code> are
          auto-reduced to P(z)/Q(z) form. The denominator Q must have no
          zeros in 𝔻̄.
        </div>
      </details>

      <!-- Structured coefficient fields -->
      <div class="hint">Coefficient fields (kept in sync with the expression):</div>
      <div class="dir-phi-coeffs"></div>
      <div class="row" style="margin-top: 6px;">
        <button class="small dir-phi-add">+ Increase degree</button>
        <button class="small dir-phi-rm" style="margin-left: 4px;">− Decrease degree</button>
      </div>

      <div class="dir-phi-warnings" style="margin-top: 8px; font-size: 11px; color: #b8860b;"></div>
    `);

    // Initial render of structured fields
    renderCoeffFields(card);
    // Initial paste-field content from default state
    setPasteFromCoeffs(card);

    // Preset dropdown
    card.querySelector('.dir-phi-preset').addEventListener('change', e => {
      const p = PHI_PRESETS_BOUNDED.find(p => p.id === e.target.value);
      if (!p) return;
      const kind = p.kind || 'polynomial';
      directState.coeffsKind = kind;
      if (kind === 'rational') {
        directState.coeffsNum = (p.num || ['0', '1']).slice();
        directState.coeffsDen = (p.den || ['1']).slice();
      } else {
        directState.coeffs = (p.coeffs || ['0', '1']).slice();
      }
      directState.expressionInput = false;
      renderCoeffFields(card);
      setPasteFromCoeffs(card);
      setStatus(card, 'ok', '');
      recomputeAndRender();
    });

    // Degree adjustment (polynomial mode only — rational uses paste field
    // or per-side controls inside renderRationalCoeffPanel).
    function adjustPolyDegree(delta) {
      if (directState.coeffsKind !== 'polynomial') return;
      if (delta > 0) {
        if (directState.coeffs.length >= 30) return;
        directState.coeffs.push('0');
      } else {
        if (directState.coeffs.length <= 2) return;
        directState.coeffs.pop();
      }
      directState.expressionInput = false;
      renderCoeffFields(card);
      setPasteFromCoeffs(card);
      recomputeAndRender();
    }
    card.querySelector('.dir-phi-add').addEventListener('click', () => adjustPolyDegree(+1));
    card.querySelector('.dir-phi-rm').addEventListener('click', () => adjustPolyDegree(-1));

    // Paste-expression: debounced real-time parsing
    const pasteInput = card.querySelector('.dir-phi-paste');
    let pasteTimer = null;
    pasteInput.addEventListener('input', () => {
      directState.expressionInput = true;
      // Immediate "checking" visual; debounce the actual parse.
      setStatus(card, 'pending', '');
      if (pasteTimer) clearTimeout(pasteTimer);
      pasteTimer = setTimeout(() => { tryParsePaste(card); }, 150);
    });
    pasteInput.addEventListener('keydown', e => {
      if (e.key === 'Enter') {
        e.preventDefault();
        if (pasteTimer) clearTimeout(pasteTimer);
        tryParsePaste(card);
      }
    });

    // Weight PARAMETER inputs (the Weight/singular selectors are in the Domain-type
    // card; row visibility is handled centrally by applyBoundedWeightRows).
    const alphaInp = card.querySelector('.dir-phi-alpha');
    if (alphaInp) alphaInp.addEventListener('input', () => { directState.alpha = alphaInp.value; recomputeAndRender(); });
    const logw0Inp = card.querySelector('.dir-phi-logw0');
    if (logw0Inp) logw0Inp.addEventListener('input', () => { directState.logW0 = logw0Inp.value; recomputeAndRender(); });
    const z0Inp = card.querySelector('.dir-phi-z0');
    if (z0Inp) z0Inp.addEventListener('input', () => { directState.z0 = z0Inp.value; recomputeAndRender(); });

    card.classList.add('dir-phi-card-bounded');
    return card;
  }

  // ---------------------------------------------------------------------------
  // φ-input card (UNBOUNDED mode)
  // ---------------------------------------------------------------------------
  // φ(z) = c·z + F_0 + F_1/z + ... + F_{m-1}/z^{m-1}
  //
  // Single c text input (real positive). m+1 (well, m: F_0..F_{m-1}) complex
  // text inputs for F_l. Add/remove buttons for m. Preset dropdown. No
  // paste-expression for unbounded mode (it's awkward to write Laurents
  // unambiguously) — structured fields only.
  // ---------------------------------------------------------------------------
  function makePhiCardUnbounded() {
    const card = section('Riemann map φ(z) — unbounded', `
      <div class="hint">
        Laurent at ∞: φ(z) = c·z + F<sub>0</sub> + F<sub>1</sub>/z + F<sub>2</sub>/z² + … where c &gt; 0.
        For a classical QD, only the case F<sub>l</sub>=0 for l≥1 gives a rational
        h with finite poles (exterior of a disk). Higher Laurent terms can be
        explored but produce h's polynomial part only.
      </div>
      <div class="row" style="margin-bottom: 8px;">
        <label>Preset:
          <select class="dir-phi-uns-preset" style="width: 280px;">
            <option value="">— custom —</option>
            ${PHI_PRESETS_UNBOUNDED.map(p => `<option value="${p.id}">${p.label}</option>`).join('')}
          </select>
        </label>
      </div>

      <!-- Weight PARAMETER inputs (Weight/singular live in the Domain-type card;
           row visibility handled centrally by applyUnboundedWeightRows). -->
      <div class="row" style="margin-bottom: 6px; gap: 10px; align-items: center;">
        <label class="dir-phi-uns-alpha-row" style="display:none;">α =
          <input type="text" class="dir-phi-uns-alpha" value="${escapeAttr(directState.unsAlpha)}" style="width: 56px; font-family: ui-monospace, monospace;">
        </label>
        <label class="dir-phi-uns-z0-row" style="display:none;">z₀ =
          <input type="text" class="dir-phi-uns-z0" value="${escapeAttr(directState.unsZ0)}" style="width: 70px; font-family: ui-monospace, monospace;">
        </label>
      </div>

      <div class="row">
        <label class="dir-phi-uns-c-row">c (φ′(∞)) =
          <input type="text" class="cnum dir-phi-uns-c" value="${directState.cValue}" style="width: 80px;">
        </label>
      </div>

      <!-- Classical Laurent inputs (hidden when a weight is selected). -->
      <div class="dir-phi-uns-classical">
        <div class="hint">Laurent coefficients F<sub>l</sub> (l = 0, 1, …):</div>
        <div class="dir-phi-uns-Fcoeffs"></div>
        <div class="row" style="margin-top: 6px;">
          <button class="small dir-phi-uns-add">+ Add F<sub>l</sub></button>
          <button class="small dir-phi-uns-rm" style="margin-left: 4px;">− Remove last</button>
        </div>
      </div>

      <!-- Weighted kernel input (shown when weight = power/log). -->
      <div class="dir-phi-uns-kernel" style="display:none;">
        <div class="dir-phi-uns-weight-hint hint" style="margin-bottom: 6px;"></div>
        <div class="row" style="align-items: stretch;">
          <label style="flex: 1 1 auto; display: flex; align-items: center; gap: 6px;">
            <span style="white-space: nowrap;">r#(z) =</span>
            <input type="text" class="dir-phi-uns-kexpr"
                   placeholder="e.g. (0.81 - 1.725z)/(1 - 2.5z)"
                   style="flex: 1 1 auto; min-width: 200px; font-family: ui-monospace, monospace;">
          </label>
          <span class="dir-phi-uns-kstatus" style="margin-left: 6px; font-size: 16px; line-height: 1.6;" aria-live="polite"></span>
        </div>
        <div class="hint" style="margin-top: 4px;">
          The rational kernel must be analytic on the closed exterior |z|≥1 (poles
          strictly inside 𝔻); a pole at z=0 gives h a polynomial part at ∞.
          Realizable ⟺ φ univalent. (Power: c is derived from r#. Singular PQD:
          z₀ is derived from a zero of r#; the z₀ field is an optional hint.)
        </div>
      </div>

      <div class="dir-phi-uns-warnings" style="margin-top: 8px; font-size: 11px; color: #b8860b;"></div>
    `);

    renderUnboundedFcoeffs(card);

    // ---- Weighted-unbounded controls (Weight/singular live in the Domain-type
    // card; row visibility handled centrally by applyUnboundedWeightRows). ----
    const setKStatus = (state, msg) => { const el = card.querySelector('.dir-phi-uns-kstatus'); if (!el) return; el.textContent = state === 'ok' ? '✓' : state === 'err' ? '✗' : ''; el.title = msg || ''; };
    const setKExprFromState = () => {
      const inp = card.querySelector('.dir-phi-uns-kexpr'); if (!inp) return;
      try {
        const P = directState.unsKernelNum.map(parseComplex), Q = directState.unsKernelDen.map(parseComplex);
        const pStr = QD.Direct.polynomialToString(P), qStr = QD.Direct.polynomialToString(Q);
        const wrap = s => /[ +\-]/.test(s.trim()) ? '(' + s + ')' : s;
        inp.value = wrap(pStr) + ' / ' + wrap(qStr);
      } catch (e) { /* leave as-is */ }
    };
    setKExprFromState();

    card.querySelector('.dir-phi-uns-alpha').addEventListener('input', e => { directState.unsAlpha = e.target.value; recomputeAndRender(); });
    card.querySelector('.dir-phi-uns-z0').addEventListener('input', e => { directState.unsZ0 = e.target.value; recomputeAndRender(); });
    const kexpr = card.querySelector('.dir-phi-uns-kexpr');
    let kTimer = null;
    const parseKernelExpr = () => {
      const expr = kexpr.value.trim(); if (!expr) { setKStatus('idle', ''); return; }
      const mathLib = (typeof window !== 'undefined') ? window.math : null;
      if (!mathLib) { setKStatus('err', 'math.js not loaded'); return; }
      let parsed; try { parsed = QD.Direct.parseRationalInZ(expr, mathLib); }
      catch (err) { setKStatus('err', err.message || String(err)); return; }
      const num = Array.isArray(parsed) ? parsed : parsed.num;
      const den = Array.isArray(parsed) ? [{ re: 1, im: 0 }] : parsed.den;
      directState.unsKernelNum = num.map(coeffToString);
      directState.unsKernelDen = den.map(coeffToString);
      setKStatus('ok', '');
      const pre = card.querySelector('.dir-phi-uns-preset'); if (pre) pre.value = '';
      recomputeAndRender();
    };
    kexpr.addEventListener('input', () => { if (kTimer) clearTimeout(kTimer); kTimer = setTimeout(parseKernelExpr, 150); });
    kexpr.addEventListener('keydown', ev => { if (ev.key === 'Enter') { ev.preventDefault(); if (kTimer) clearTimeout(kTimer); parseKernelExpr(); } });

    card.querySelector('.dir-phi-uns-c').addEventListener('input', e => {
      directState.cValue = e.target.value;
      recomputeAndRender();
    });

    card.querySelector('.dir-phi-uns-preset').addEventListener('change', e => {
      const p = PHI_PRESETS_UNBOUNDED.find(p => p.id === e.target.value);
      if (!p) return;
      directState.cValue = p.c;
      directState.Fcoeffs = p.F.slice();
      card.querySelector('.dir-phi-uns-c').value = p.c;
      renderUnboundedFcoeffs(card);
      recomputeAndRender();
    });

    card.querySelector('.dir-phi-uns-add').addEventListener('click', () => {
      if (directState.Fcoeffs.length >= 8) return;
      directState.Fcoeffs.push('0');
      renderUnboundedFcoeffs(card);
      recomputeAndRender();
    });
    card.querySelector('.dir-phi-uns-rm').addEventListener('click', () => {
      if (directState.Fcoeffs.length === 0) return;
      directState.Fcoeffs.pop();
      renderUnboundedFcoeffs(card);
      recomputeAndRender();
    });

    card.classList.add('dir-phi-card-unbounded');
    return card;
  }

  function renderUnboundedFcoeffs(card) {
    const box = card.querySelector('.dir-phi-uns-Fcoeffs');
    box.innerHTML = '';
    for (let l = 0; l < directState.Fcoeffs.length; l++) {
      const row = document.createElement('div');
      row.className = 'row';
      const label = document.createElement('label');
      label.innerHTML = `F<sub>${l}</sub> = `;
      const inp = document.createElement('input');
      inp.type = 'text';
      inp.className = 'cnum';
      inp.value = directState.Fcoeffs[l];
      inp.addEventListener('input', () => {
        directState.Fcoeffs[l] = inp.value;
        recomputeAndRender();
      });
      label.appendChild(inp);
      row.appendChild(label);
      box.appendChild(row);
    }
  }

  // ---------------------------------------------------------------------------
  // φ-input card (NUMERICAL mode)
  // ---------------------------------------------------------------------------
  // Free-form math.js expression in z. Live-parses + samples on |z|=1 +
  // DFT-extracts Taylor coefficients + calls boundedQD on the truncated
  // polynomial approximation.
  // ---------------------------------------------------------------------------
  const NUM_PRESETS = [
    { id: 'identity',  label: 'φ = z',                              expr: 'z' },
    { id: 'cubic',     label: 'φ = z + 0.1z² − 0.05z³',             expr: 'z + 0.1*z^2 - 0.05*z^3' },
    { id: 'exp',       label: 'φ = z·exp(z/4)',                     expr: 'z * exp(z/4)' },
    { id: 'sin',       label: 'φ = z + 0.2·sin(z)',                 expr: 'z + 0.2*sin(z)' },
    { id: 'rational',  label: 'φ = z/(1 − 0.3z)',                   expr: 'z / (1 - 0.3*z)' },
    { id: 'log',       label: 'φ = log(1+z)  (slow Taylor decay)',  expr: 'log(1+z)' },
    { id: 'nonana',    label: '(Non-analytic) φ = conj(z)',         expr: 'conj(z)' },
  ];

  function makePhiCardNumerical() {
    const card = section('Riemann map φ(z) — numerical (free-form)', `
      <div class="hint">
        Type any math.js expression in <code>z</code> (e.g. <code>z + 0.2*sin(z)</code>,
        <code>z·exp(z/4)</code>, <code>z/(1 − 0.3z)</code>). The app samples on
        |z|=1, infers the polynomial Taylor approximation of φ at z=0 via DFT, and
        computes h. Non-analytic φ (e.g. <code>conj(z)</code>) is flagged.
      </div>
      <div class="row" style="margin-bottom: 6px;">
        <label>Preset:
          <select class="dir-phi-num-preset" style="width: 280px;">
            <option value="">— custom —</option>
            ${NUM_PRESETS.map(p => `<option value="${p.id}">${p.label}</option>`).join('')}
          </select>
        </label>
      </div>
      <div class="row" style="align-items: stretch;">
        <label style="flex: 1 1 auto; display: flex; align-items: center; gap: 6px;">
          <span>φ(z) =</span>
          <input type="text" class="dir-phi-num-expr"
                 placeholder="e.g. z + 0.2*sin(z)"
                 value="${escapeAttr(directState.numExpr)}"
                 style="flex: 1 1 auto; font-family: ui-monospace, monospace;">
        </label>
        <span class="dir-phi-num-status" style="margin-left: 6px; font-size: 16px;"></span>
      </div>
      <div class="dir-phi-num-msg" style="font-size: 11px; min-height: 1.2em; margin-bottom: 6px;"></div>

      <details style="margin-bottom: 4px;">
        <summary style="cursor:pointer; user-select:none; font-size:12px; color:#555;">Advanced</summary>
        <div class="row" style="margin-top:6px;">
          <label>Truncation degree (DFT cap):
            <input type="number" class="dir-phi-num-maxorder" min="1" max="32"
                   value="${directState.numMaxOrder}" style="width: 64px;">
          </label>
        </div>
      </details>

      <div class="dir-phi-num-diag" style="margin-top: 8px; font-size: 11px; color: #5677a8; font-family: ui-monospace, monospace;"></div>
      <div class="dir-phi-num-warnings" style="margin-top: 4px; font-size: 11px; color: #b8860b;"></div>
    `);

    // Expression input (debounced re-evaluation)
    const expr  = card.querySelector('.dir-phi-num-expr');
    const order = card.querySelector('.dir-phi-num-maxorder');
    let timer = null;
    function schedule() {
      if (timer) clearTimeout(timer);
      timer = setTimeout(() => recomputeAndRender(), 200);
    }
    expr.addEventListener('input', () => {
      directState.numExpr = expr.value;
      schedule();
    });
    order.addEventListener('input', () => {
      const v = parseInt(order.value, 10);
      if (!Number.isNaN(v) && v >= 1 && v <= 32) directState.numMaxOrder = v;
      schedule();
    });

    // Preset dropdown
    card.querySelector('.dir-phi-num-preset').addEventListener('change', e => {
      const p = NUM_PRESETS.find(p => p.id === e.target.value);
      if (!p) return;
      directState.numExpr = p.expr;
      expr.value = p.expr;
      recomputeAndRender();
    });

    card.classList.add('dir-phi-card-numerical');
    return card;
  }

  function escapeAttr(s) {
    return String(s).replace(/&/g,'&amp;').replace(/"/g,'&quot;').replace(/</g,'&lt;');
  }

  // Try parsing the paste field. Accepts polynomial OR rational expressions in z.
  // - Polynomial result → coeffsKind='polynomial', populate coeffs.
  // - Rational result   → coeffsKind='rational', populate coeffsNum + coeffsDen.
  function tryParsePaste(card) {
    const pasteInput = card.querySelector('.dir-phi-paste');
    const expr = pasteInput.value.trim();
    if (!expr) { setStatus(card, 'idle', ''); return; }
    const mathLib = (typeof window !== 'undefined') ? window.math : null;
    if (!mathLib) { setStatus(card, 'err', 'math.js not loaded'); return; }
    let parsed;
    try { parsed = QD.Direct.parseRationalInZ(expr, mathLib); }
    catch (err) {
      setStatus(card, 'err', err.message || String(err));
      return;
    }
    if (Array.isArray(parsed)) {
      // Polynomial result.
      directState.coeffsKind = 'polynomial';
      directState.coeffs = parsed.map(coeffToString);
      // c_1 ≠ 0 sanity: required for φ to be locally univalent at z=0.
      const c1 = parsed[1] || {re:0, im:0};
      if (Math.hypot(c1.re, c1.im) < 1e-14) {
        setStatus(card, 'err', 'c₁ = 0; φ not locally univalent at 0');
        return;
      }
      setStatus(card, 'ok', `polynomial, degree ${parsed.length - 1}`);
    } else {
      // Rational result {num, den}.
      directState.coeffsKind = 'rational';
      directState.coeffsNum = parsed.num.map(coeffToString);
      directState.coeffsDen = parsed.den.map(coeffToString);
      setStatus(card, 'ok',
        `rational, deg(num)=${parsed.num.length - 1}, deg(den)=${parsed.den.length - 1}`);
    }
    renderCoeffFields(card);
    const preset = card.querySelector('.dir-phi-preset');
    if (preset) preset.value = '';
    recomputeAndRender();
  }

  // Populate the paste field with the canonical form of the current coeffs.
  // Handles both polynomial and rational kinds. Skipped while the user is
  // typing in the paste field (directState.expressionInput).
  function setPasteFromCoeffs(card) {
    if (directState.expressionInput) return;
    const pasteInput = card.querySelector('.dir-phi-paste');
    if (!pasteInput) return;
    if (directState.coeffsKind === 'rational') {
      let P, Q;
      try {
        P = directState.coeffsNum.map(parseComplex);
        Q = directState.coeffsDen.map(parseComplex);
      } catch (e) { return; }
      const pStr = QD.Direct.polynomialToString(P);
      const qStr = QD.Direct.polynomialToString(Q);
      // Wrap each side in parens when it has multiple terms; '/' otherwise.
      const wrapNeed = s => /[ +\-]/.test(s.trim());
      const lhs = wrapNeed(pStr) ? '(' + pStr + ')' : pStr;
      const rhs = wrapNeed(qStr) ? '(' + qStr + ')' : qStr;
      pasteInput.value = lhs + ' / ' + rhs;
      setStatus(card, 'ok', `rational, deg(num)=${P.length - 1}, deg(den)=${Q.length - 1}`);
      return;
    }
    let coeffsC;
    try { coeffsC = directState.coeffs.map(parseComplex); }
    catch (e) { return; }
    if (QD.Direct.polynomialToString) {
      pasteInput.value = QD.Direct.polynomialToString(coeffsC);
    }
    setStatus(card, 'ok', `polynomial, degree ${coeffsC.length - 1}`);
  }

  // Status indicator: 'ok' (green ✓), 'err' (red ✗), 'pending' (gray spinner),
  // 'idle' (cleared).
  function setStatus(card, kind, detail) {
    const status = card.querySelector('.dir-phi-status');
    const msg    = card.querySelector('.dir-phi-paste-msg');
    if (!status || !msg) return;
    switch (kind) {
      case 'ok':
        status.textContent = '✓'; status.style.color = '#2a8f2a';
        msg.style.color = '#2a8f2a'; msg.textContent = detail || '';
        break;
      case 'err':
        status.textContent = '✗'; status.style.color = '#b53030';
        msg.style.color = '#b53030'; msg.textContent = detail || '';
        break;
      case 'pending':
        status.textContent = '…'; status.style.color = '#888';
        msg.style.color = '#888'; msg.textContent = 'parsing…';
        break;
      case 'idle':
      default:
        status.textContent = ''; msg.textContent = '';
        break;
    }
  }

  // Render the structured coefficient fields. Dispatches on coeffsKind:
  //   • polynomial → one coefficient list (c_0 .. c_n) + ± degree row visible
  //   • rational   → two lists (P and Q) + ± degree row hidden (each has its own)
  function renderCoeffFields(card) {
    const box = card.querySelector('.dir-phi-coeffs');
    box.innerHTML = '';
    // Show/hide the polynomial-mode degree-adjustment buttons.
    const addBtn = card.querySelector('.dir-phi-add');
    const polyDegRow = addBtn && addBtn.parentElement;
    if (polyDegRow) polyDegRow.style.display = (directState.coeffsKind === 'rational') ? 'none' : '';
    if (directState.coeffsKind === 'rational') {
      renderRationalCoeffPanel(box, card);
    } else {
      renderPolynomialCoeffPanel(box, card);
    }
  }

  function renderPolynomialCoeffPanel(box, card) {
    for (let k = 0; k < directState.coeffs.length; k++) {
      const row = document.createElement('div');
      row.className = 'row';
      const label = document.createElement('label');
      label.innerHTML = `c<sub>${k}</sub>${k === 0 ? ' = w₀' : ''} = `;
      const inp = document.createElement('input');
      inp.type = 'text';
      inp.className = 'cnum';
      inp.value = directState.coeffs[k];
      inp.addEventListener('input', () => {
        directState.coeffs[k] = inp.value;
        directState.expressionInput = false;
        setPasteFromCoeffs(card);
        recomputeAndRender();
      });
      if (k === 1) {
        const checkValidity = () => {
          let c;
          try { c = parseComplex(inp.value); inp.classList.remove('invalid'); }
          catch (e) { inp.classList.add('invalid'); return; }
          if (Math.hypot(c.re, c.im) < 1e-14) inp.classList.add('invalid');
          else inp.classList.remove('invalid');
        };
        inp.addEventListener('input', checkValidity);
        checkValidity();
      }
      label.appendChild(inp);
      row.appendChild(label);
      box.appendChild(row);
    }
  }

  function renderRationalCoeffPanel(box, card) {
    const wrap = document.createElement('div');
    wrap.innerHTML = `
      <div class="hint" style="margin-bottom:4px;">
        Numerator P(z) = Σ p<sub>k</sub> z<sup>k</sup>
      </div>
      <div class="dir-phi-num-coeffs"></div>
      <div class="row" style="margin-top: 4px;">
        <button class="small dir-num-add">+ deg(P)</button>
        <button class="small dir-num-rm" style="margin-left:4px;">− deg(P)</button>
      </div>
      <div class="hint" style="margin-top:10px; margin-bottom:4px;">
        Denominator Q(z) = Σ q<sub>k</sub> z<sup>k</sup>
        <span style="color:#888;">(must have no zeros in 𝔻̄)</span>
      </div>
      <div class="dir-phi-den-coeffs"></div>
      <div class="row" style="margin-top: 4px;">
        <button class="small dir-den-add">+ deg(Q)</button>
        <button class="small dir-den-rm" style="margin-left:4px;">− deg(Q)</button>
      </div>
    `;
    box.appendChild(wrap);
    const onAny = () => {
      directState.expressionInput = false;
      setPasteFromCoeffs(card);
      recomputeAndRender();
    };
    fillCoeffList(wrap.querySelector('.dir-phi-num-coeffs'), 'p', directState.coeffsNum, onAny);
    fillCoeffList(wrap.querySelector('.dir-phi-den-coeffs'), 'q', directState.coeffsDen, onAny);

    function adjustDegree(arr, delta, minLen) {
      if (delta > 0) {
        if (arr.length >= 30) return;
        arr.push('0');
      } else {
        if (arr.length <= minLen) return;
        arr.pop();
      }
      renderCoeffFields(card);
      onAny();
    }
    wrap.querySelector('.dir-num-add').addEventListener('click', () => adjustDegree(directState.coeffsNum, +1, 1));
    wrap.querySelector('.dir-num-rm').addEventListener('click', () => adjustDegree(directState.coeffsNum, -1, 2));
    wrap.querySelector('.dir-den-add').addEventListener('click', () => adjustDegree(directState.coeffsDen, +1, 1));
    wrap.querySelector('.dir-den-rm').addEventListener('click', () => adjustDegree(directState.coeffsDen, -1, 1));
  }

  function fillCoeffList(box, sym, arr, onChange) {
    box.innerHTML = '';
    for (let k = 0; k < arr.length; k++) {
      const row = document.createElement('div');
      row.className = 'row';
      const label = document.createElement('label');
      label.innerHTML = `${sym}<sub>${k}</sub> = `;
      const inp = document.createElement('input');
      inp.type = 'text';
      inp.className = 'cnum';
      inp.value = arr[k];
      inp.addEventListener('input', () => { arr[k] = inp.value; onChange(); });
      label.appendChild(inp);
      row.appendChild(label);
      box.appendChild(row);
    }
  }

  // ---------------------------------------------------------------------------
  // Output card: computed h + Send-to-inverse button
  // ---------------------------------------------------------------------------
  function makeOutputCard() {
    const card = section('Computed quadrature function h(w)', `
      <div class="dir-h-display" style="font-family: ui-monospace, monospace; font-size: 12px; line-height: 1.6;"></div>
      <div class="dir-h-katex rm-sym" style="margin-top: 8px;"></div>
      <div class="row" style="margin-top: 10px;">
        <button class="primary dir-send-btn">Send to inverse mode →</button>
        <button class="dir-verify-btn" style="margin-left: 6px;">Verify ↻</button>
        <span class="dir-send-msg" style="margin-left: 8px; font-size: 11px;"></span>
      </div>
      <div class="dir-verify-result" style="margin-top: 6px; font-size: 12px; font-family: ui-monospace, monospace;"></div>
      <div class="dir-error" style="margin-top: 8px; font-size: 12px; color: #b53030; font-family: ui-monospace, monospace;"></div>
    `);
    card.querySelector('.dir-send-btn').addEventListener('click', () => {
      const msg = card.querySelector('.dir-send-msg');
      if (!directState.lastH) {
        msg.style.color = '#b53030';
        msg.textContent = 'Compute a valid h first.';
        return;
      }
      const hook = window.QD && window.QD.Direct && window.QD.Direct._sendHToInverseTab;
      if (!hook) {
        msg.style.color = '#b53030';
        msg.textContent = 'Send hook not installed yet (try again in a moment).';
        return;
      }
      // Tag the FULL family identity so the inverse tab lands in the matching
      // mode (the receiving _sendHToInverseTab composes the mode from these).
      const unbounded = directState.mode === 'unbounded';
      const opts = { unbounded };
      if (unbounded) opts.c = directState.lastC;
      if (directState.lastWeight === 'power') {
        opts.alpha = parseFloat(unbounded ? directState.unsAlpha : directState.alpha);
      } else if (directState.lastWeight === 'log') {
        opts.lqd = true;
      }
      if (directState.lastWeight !== 'classical' && directState.lastSingular) opts.singular = true;
      if (directState.lastQ) opts.q = directState.lastQ;
      hook(directState.lastH, opts);
      msg.style.color = '#2a8f2a';
      msg.textContent = 'Sent. Switched to inverse view.';
    });

    // ---- Verify button: round-trip via the inverse solver. ----
    card.querySelector('.dir-verify-btn').addEventListener('click', () => {
      runVerify(card);
    });

    return card;
  }

  // ===========================================================================
  // Recompute pipeline + Verify -> direct-recompute.js / direct-verify.js
  // (Phase-3 item E). One shared dCtx; the card-builder handlers + _activate
  // call the captured names. parseComplex / coeffToString / section stay here
  // (host card builders use them too).
  // ===========================================================================
  const dCtx = { directState, parseComplex, isMounted: () => mounted };
  ({ recomputeAndRender } = window.QD_UI.installDirectRecompute(dCtx));
  ({ runVerify } = window.QD_UI.installDirectVerify(dCtx));

  // ---------------------------------------------------------------------------
  // String <-> Complex helpers (parser lives in QD.Direct.parseRationalInZ)
  // ---------------------------------------------------------------------------
  function parseComplex(s) {
    if (typeof s !== 'string') return { re: Number(s) || 0, im: 0 };
    s = s.trim();
    if (s === '') return { re: 0, im: 0 };
    // Use math.js to parse "1+2i" robustly.
    if (typeof math !== 'undefined' && math.complex) {
      try {
        const v = math.complex(s);
        return { re: Number(v.re), im: Number(v.im) };
      } catch (e) { /* fall through */ }
    }
    // Fallback: try as a number.
    const n = Number(s);
    if (!Number.isNaN(n)) return { re: n, im: 0 };
    throw new Error("Can't parse complex value: " + s);
  }

  // Thin wrappers around Complex.format (which handles the ±i short forms,
  // integer snap, and zero detection in one place). The three names below
  // are kept as separate functions for readability at the call sites —
  // they all map to the same primitive.
  function coeffToString(c)    { return QD.Complex.format(c); }

  // ---------------------------------------------------------------------------
  // Section helper (matches existing card markup)
  // ---------------------------------------------------------------------------
  function section(title, innerHTML) {
    const sec = document.createElement('section');
    sec.className = 'card';
    sec.innerHTML = `<h2>${title}</h2>${innerHTML}`;
    return sec;
  }

}());
