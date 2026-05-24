// =============================================================================
// direct-ui.js -- Direct-problem tab UI.
//
// Three modes share the same plot canvas:
//   1. Bounded — polynomial or rational φ; structured coefficient fields
//      side-by-side with a live-parsing "paste expression" text input.
//   2. Unbounded — Laurent-at-∞ φ = c·z + Σ F_l/z^l; structured (c, F_l)
//      coefficient fields.
//   3. Numerical — free-form math.js expression in z; DFT-extracted
//      polynomial truncation produces an approximate h with an analyticity
//      diagnostic.
//
// Each mode pushes:
//   • h to QD.Direct._sendHToInverseTab     (pre-fill QD/LQD tab and switch)
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

    // Unbounded mode: φ(z) = c·z + Σ_l F_l/z^l. Strings.
    cValue: '1',                   // conformal radius (positive real)
    Fcoeffs: [],                   // [F_0, F_1, ..., F_{m-1}], strings; empty ⇒ φ = c·z

    // Numerical mode: any math.js expression in z.
    numExpr: 'z + 0.2*sin(z)',     // default: a non-polynomial example
    numMaxOrder: 12,

    // Last successfully computed h, and the c that produced it.
    lastH: null,
    lastC: 1,

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
  // Lazy mount on first tab activation
  // ---------------------------------------------------------------------------
  let mounted = false;
  document.addEventListener('tab-changed', function (e) {
    if (e.detail && e.detail.tab === 'direct') {
      if (!mounted) { mountDirectSidebar(); mounted = true; }
      // Push the current ∂Ω to the canvas every time we switch in.
      recomputeAndRender();
    }
  });

  function mountDirectSidebar() {
    const root = document.getElementById('controls-direct');
    if (!root) return;
    root.innerHTML = '';
    root.appendChild(makeDomainTypeCard());
    root.appendChild(makePhiCardBounded());      // initially visible
    root.appendChild(makePhiCardUnbounded());    // initially hidden
    root.appendChild(makePhiCardNumerical());    // initially hidden
    root.appendChild(makeOutputCard());
    applyModeVisibility();
  }

  function applyModeVisibility() {
    const root = document.getElementById('controls-direct');
    if (!root) return;
    const bounded   = root.querySelector('.dir-phi-card-bounded');
    const unbounded = root.querySelector('.dir-phi-card-unbounded');
    const numerical = root.querySelector('.dir-phi-card-numerical');
    if (bounded)   bounded.style.display   = directState.mode === 'bounded'   ? '' : 'none';
    if (unbounded) unbounded.style.display = directState.mode === 'unbounded' ? '' : 'none';
    if (numerical) numerical.style.display = directState.mode === 'numerical' ? '' : 'none';
  }

  // ---------------------------------------------------------------------------
  // Domain-type card
  // ---------------------------------------------------------------------------
  function makeDomainTypeCard() {
    const card = section('Domain type', `
      <div class="row">
        <label><input type="radio" name="dir-domain-mode" value="bounded" ${directState.mode==='bounded'?'checked':''}> Bounded</label>
        <label style="margin-left:14px;"><input type="radio" name="dir-domain-mode" value="unbounded" ${directState.mode==='unbounded'?'checked':''}> Unbounded</label>
      </div>
      <div class="row">
        <label><input type="radio" name="dir-domain-mode" value="numerical" ${directState.mode==='numerical'?'checked':''}> Numerical (any expression)</label>
      </div>
      <div class="hint" style="margin-top: 6px;">
        <strong>Bounded</strong>: φ(z) = Σ c<sub>k</sub> z<sup>k</sup> (polynomial).
        <strong>Unbounded</strong>: φ(z) = c·z + Σ F<sub>l</sub>/z<sup>l</sup>.
        <strong>Numerical</strong>: free-form math.js expression in z; we infer the
        bounded-QD polynomial by DFT and report a non-analyticity diagnostic.
      </div>
    `);
    card.querySelectorAll('input[name="dir-domain-mode"]').forEach(r => {
      r.addEventListener('change', () => {
        directState.mode = r.value;
        directState.expressionInput = false;
        applyModeVisibility();
        recomputeAndRender();
      });
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
    // Marker class so applyModeVisibility() can show/hide this card.
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

      <div class="row">
        <label>c (φ′(∞)) =
          <input type="text" class="cnum dir-phi-uns-c" value="${directState.cValue}" style="width: 80px;">
        </label>
      </div>

      <div class="hint">Laurent coefficients F<sub>l</sub> (l = 0, 1, …):</div>
      <div class="dir-phi-uns-Fcoeffs"></div>
      <div class="row" style="margin-top: 6px;">
        <button class="small dir-phi-uns-add">+ Add F<sub>l</sub></button>
        <button class="small dir-phi-uns-rm" style="margin-left: 4px;">− Remove last</button>
      </div>

      <div class="dir-phi-uns-warnings" style="margin-top: 8px; font-size: 11px; color: #b8860b;"></div>
    `);

    renderUnboundedFcoeffs(card);

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

      <div class="row">
        <label>Truncation degree (DFT cap):
          <input type="number" class="dir-phi-num-maxorder" min="1" max="32"
                 value="${directState.numMaxOrder}" style="width: 64px;">
        </label>
      </div>

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
    const rmBtn  = card.querySelector('.dir-phi-rm');
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
      const opts = (directState.mode === 'unbounded')
        ? { unbounded: true, c: directState.lastC }
        : { unbounded: false };
      hook(directState.lastH, opts);
      msg.style.color = '#2a8f2a';
      msg.textContent = 'Sent. Switched to QD/LQD tab.';
    });

    // ---- Verify button: round-trip via the inverse solver. ----
    card.querySelector('.dir-verify-btn').addEventListener('click', () => {
      runVerify(card);
    });

    return card;
  }

  // ===========================================================================
  // Verify: check the BOUNDARY IDENTITY directly.
  // ---------------------------------------------------------------------------
  //   For any classical QD (bounded or unbounded), the Schwarz function
  //   satisfies σ(w) = w̄ on ∂Ω, and h is the meromorphic representative
  //   of σ. So at every boundary sample w_n = φ(e^{iθ_n}) the identity
  //
  //       h(w_n)  =  conj(w_n)
  //
  //   must hold. This is exact for any φ that genuinely defines a QD (the
  //   symbolic kernel produces h such that the identity is satisfied to
  //   machine precision); for the numerical kernel applied to a polynomial-
  //   truncated φ, the residual is the truncation error; for non-QD φ
  //   (e.g. higher-Laurent unbounded shapes) the residual is large and
  //   confirms that Ω is not actually a QD.
  //
  //   No inverse solver involved — h is explicit from φ via Faber, so the
  //   verification is purely a forward evaluation.
  // ===========================================================================
  function runVerify(card) {
    const resBox = card.querySelector('.dir-verify-result');
    const overlayHook = window.QD && window.QD.Direct && window.QD.Direct._setPlotOverlay;
    if (overlayHook) overlayHook(null);                  // clear any stale overlay

    if (!directState.lastH) {
      resBox.style.color = '#b53030';
      resBox.textContent = 'Compute a valid h first.';
      return;
    }

    // Sample φ at N points on |z|=1.
    const N = 500;
    let phiPts;
    try { phiPts = sampleAnalyticPhi(N); }
    catch (e) {
      resBox.style.color = '#b53030';
      resBox.textContent = 'Could not sample φ: ' + (e.message || e);
      return;
    }

    // The correct identity is: h(φ(z)) − conj(φ(z)) is analytic in 𝔻 (after
    // composing with φ), so its Fourier expansion on |z|=1 has only non-
    // negative-frequency terms. We measure the negative-frequency Fourier
    // mass — this should be ≈ 0 for any valid classical QD.
    const v = QD.Direct.verifyBoundaryIdentity(directState.lastH, phiPts);

    // Relative score: negMass normalised by the boundary-data scale.
    const relNeg = v.scale > 0 ? v.negMass / v.scale : v.negMass;
    let color;
    if      (relNeg < 1e-8) color = '#2a8f2a';
    else if (relNeg < 1e-2) color = '#b8860b';
    else                    color = '#b53030';

    resBox.style.color = color;
    resBox.innerHTML =
      'Fourier diagnostic (' + N + ' samples on |z|=1):<br>' +
      '&nbsp;&nbsp;negative-freq mass = <strong>' + v.negMass.toExponential(2) +
      '</strong> (relative ' + relNeg.toExponential(2) + ')<br>' +
      '&nbsp;&nbsp;<span style="color:#888">zero-mode mass = ' + v.zeroMass.toExponential(2) +
      ', positive-freq mass = ' + v.posMass.toExponential(2) + '</span>';
  }

  // Sample the user's input φ at N uniform θ on |z|=1, in the mode-appropriate way.
  function sampleAnalyticPhi(N) {
    if (directState.mode === 'bounded') {
      if (directState.coeffsKind === 'rational') {
        const P = directState.coeffsNum.map(parseComplex);
        const Q = directState.coeffsDen.map(parseComplex);
        const pts = new Array(N);
        for (let n = 0; n < N; n++) {
          const t = 2 * Math.PI * n / N;
          const z = { re: Math.cos(t), im: Math.sin(t) };
          const pv = QD.Direct.evalPolyAscending(P, z);
          const qv = QD.Direct.evalPolyAscending(Q, z);
          const d2 = qv.re*qv.re + qv.im*qv.im;
          pts[n] = { re: (pv.re*qv.re + pv.im*qv.im) / d2,
                     im: (pv.im*qv.re - pv.re*qv.im) / d2 };
        }
        return pts;
      }
      const cs = directState.coeffs.map(parseComplex);
      return QD.Direct.sampleBoundaryPolynomial(cs, N);
    } else if (directState.mode === 'unbounded') {
      const c = Number(directState.cValue);
      const F = directState.Fcoeffs.map(parseComplex);
      return QD.Direct.sampleBoundaryLaurent(c, F, N);
    } else {
      // Numerical: re-evaluate the user's expression.
      if (!window.math) throw new Error('math.js not loaded');
      const compiled = window.math.parse(directState.numExpr).compile();
      const pts = new Array(N);
      for (let n = 0; n < N; n++) {
        const theta = 2 * Math.PI * n / N;
        const v = compiled.evaluate({ z: window.math.complex(Math.cos(theta), Math.sin(theta)) });
        if (typeof v === 'number') pts[n] = { re: v, im: 0 };
        else pts[n] = { re: v.re, im: v.im };
      }
      return pts;
    }
  }

  // ---------------------------------------------------------------------------
  // Recompute h and redraw ∂Ω. Dispatches on directState.mode.
  // ---------------------------------------------------------------------------
  function recomputeAndRender() {
    if (!mounted) return;
    const root = document.getElementById('controls-direct');
    if (!root) return;

    const hDisp   = root.querySelector('.dir-h-display');
    const hKatex  = root.querySelector('.dir-h-katex');
    const errBox  = root.querySelector('.dir-error');
    if (hDisp) hDisp.textContent = '';
    if (hKatex) hKatex.innerHTML = '';
    if (errBox) errBox.textContent = '';

    if (directState.mode === 'bounded') {
      recomputeBounded(root, hDisp, hKatex, errBox);
    } else if (directState.mode === 'unbounded') {
      recomputeUnbounded(root, hDisp, hKatex, errBox);
    } else {
      recomputeNumerical(root, hDisp, hKatex, errBox);
    }
  }

  function recomputeNumerical(root, hDisp, hKatex, errBox) {
    const card    = root.querySelector('.dir-phi-card-numerical');
    const status  = card && card.querySelector('.dir-phi-num-status');
    const msg     = card && card.querySelector('.dir-phi-num-msg');
    const diag    = card && card.querySelector('.dir-phi-num-diag');
    const warnBox = card && card.querySelector('.dir-phi-num-warnings');
    if (status)  status.textContent = '';
    if (msg)     msg.textContent = '';
    if (diag)    diag.textContent = '';
    if (warnBox) warnBox.textContent = '';

    const exprStr = directState.numExpr.trim();
    if (!exprStr) {
      if (msg) { msg.style.color = '#888'; msg.textContent = 'enter an expression'; }
      return;
    }
    if (typeof window === 'undefined' || !window.math || !window.math.parse) {
      if (errBox) errBox.textContent = 'math.js not loaded';
      return;
    }

    // Parse expression once; build a numeric phiFn.
    let node;
    try { node = window.math.parse(exprStr); }
    catch (err) {
      if (status) { status.textContent = '✗'; status.style.color = '#b53030'; }
      if (msg) { msg.style.color = '#b53030'; msg.textContent = 'parse: ' + (err.message || err); }
      return;
    }
    let compiled;
    try { compiled = node.compile(); }
    catch (err) {
      if (status) { status.textContent = '✗'; status.style.color = '#b53030'; }
      if (msg) { msg.style.color = '#b53030'; msg.textContent = 'compile: ' + (err.message || err); }
      return;
    }
    const phiFn = z => {
      const v = compiled.evaluate({ z: window.math.complex(z.re, z.im) });
      if (typeof v === 'number') return { re: v, im: 0 };
      if (v && typeof v.re === 'number' && typeof v.im === 'number') return { re: v.re, im: v.im };
      throw new Error('expression did not evaluate to a complex/number');
    };

    let result;
    try {
      result = QD.Direct.numericalBoundedQD(phiFn, {
        numSamples: 256,
        maxOrder: directState.numMaxOrder,
      });
    } catch (err) {
      if (status) { status.textContent = '✗'; status.style.color = '#b53030'; }
      if (msg) { msg.style.color = '#b53030'; msg.textContent = err.message || String(err); }
      return;
    }

    directState.lastH = result.hData;
    directState.lastC = 0;                          // bounded mode (numerical reduces to bounded)

    if (status) { status.textContent = '✓'; status.style.color = '#2a8f2a'; }
    if (msg) {
      msg.style.color = '#2a8f2a';
      msg.textContent = 'truncated at degree ' + result.truncationOrder
                      + ' (analyticity score = ' + result.analyticityScore.toExponential(2) + ')';
    }
    if (diag) {
      const lines = ['Recovered Taylor coefficients of φ at z=0:'];
      for (let k = 0; k <= Math.min(result.truncationOrder, 6); k++) {
        const c = result.taylorCoeffs[k];
        lines.push('  c_' + k + ' = ' + formatNumLocal(c.re) + (c.im >= 0 ? '+' : '') + formatNumLocal(c.im) + 'i');
      }
      if (result.truncationOrder > 6) lines.push('  …');
      diag.textContent = lines.join('\n');
    }
    if (warnBox && result.warnings.length) warnBox.textContent = '⚠ ' + result.warnings.join('; ');

    displayH(hDisp, hKatex, result.hData, /*isUnbounded=*/false);

    // Live ∂Ω preview: sample using the user's phiFn directly.
    try {
      const N = 400;
      const pts = new Array(N);
      for (let n = 0; n < N; n++) {
        const theta = 2 * Math.PI * n / N;
        pts[n] = phiFn({ re: Math.cos(theta), im: Math.sin(theta) });
      }
      pushBoundaryToPlot(pts, false);
    } catch (e) { /* preview is best-effort */ }
  }

  function formatNumLocal(x) {
    if (!isFinite(x)) return String(x);
    if (Math.abs(x) < 1e-12) return '0';
    if (Math.abs(x - Math.round(x)) < 1e-10) return String(Math.round(x));
    return Number(x.toPrecision(6)).toString();
  }

  function recomputeBounded(root, hDisp, hKatex, errBox) {
    const warnBox = root.querySelector('.dir-phi-warnings');
    if (warnBox) warnBox.textContent = '';

    if (directState.coeffsKind === 'rational') {
      let P, Q;
      try {
        P = directState.coeffsNum.map(parseComplex);
        Q = directState.coeffsDen.map(parseComplex);
      } catch (err) {
        if (errBox) errBox.textContent = 'Coefficient parse error: ' + err.message;
        return;
      }
      let result;
      try { result = QD.Direct.boundedQDRational(P, Q); }
      catch (err) {
        if (errBox) errBox.textContent = err.message;
        directState.lastH = null;
        return;
      }
      directState.lastH = result.hData;
      directState.lastC = 0;
      if (warnBox && result.warnings.length) warnBox.textContent = '⚠ ' + result.warnings.join('; ');

      displayH(hDisp, hKatex, result.hData);
      const N = 400;
      const pts = new Array(N);
      for (let n = 0; n < N; n++) {
        const t = 2 * Math.PI * n / N;
        const z = { re: Math.cos(t), im: Math.sin(t) };
        const pv = QD.Direct.evalPolyAscending(P, z);
        const qv = QD.Direct.evalPolyAscending(Q, z);
        const d2 = qv.re*qv.re + qv.im*qv.im;
        pts[n] = { re: (pv.re*qv.re + pv.im*qv.im) / d2,
                   im: (pv.im*qv.re - pv.re*qv.im) / d2 };
      }
      pushBoundaryToPlot(pts, false);
      return;
    }

    // Polynomial path.
    let coeffs;
    try { coeffs = directState.coeffs.map(parseComplex); }
    catch (err) {
      if (errBox) errBox.textContent = 'Coefficient parse error: ' + err.message;
      return;
    }
    let result;
    try { result = QD.Direct.boundedQD(coeffs); }
    catch (err) {
      if (errBox) errBox.textContent = err.message;
      directState.lastH = null;
      return;
    }
    directState.lastH = result.hData;
    directState.lastC = 0;
    if (warnBox && result.warnings.length) warnBox.textContent = '⚠ ' + result.warnings.join('; ');

    displayH(hDisp, hKatex, result.hData);
    pushBoundaryToPlot(QD.Direct.sampleBoundaryPolynomial(coeffs, 400), false);
  }

  function recomputeUnbounded(root, hDisp, hKatex, errBox) {
    const warnBox = root.querySelector('.dir-phi-uns-warnings');
    if (warnBox) warnBox.textContent = '';

    let c;
    try {
      const parsed = parseComplex(directState.cValue);
      if (Math.abs(parsed.im) > 1e-12 || parsed.re <= 0 || !isFinite(parsed.re)) {
        throw new Error("c must be a positive real number");
      }
      c = parsed.re;
    } catch (err) {
      if (errBox) errBox.textContent = 'c parse error: ' + err.message;
      return;
    }

    let F;
    try { F = directState.Fcoeffs.map(parseComplex); }
    catch (err) {
      if (errBox) errBox.textContent = 'F coefficient parse error: ' + err.message;
      return;
    }

    let result;
    try { result = QD.Direct.unboundedQD(c, F); }
    catch (err) {
      if (errBox) errBox.textContent = err.message;
      directState.lastH = null;
      return;
    }
    directState.lastH = result.hData;
    directState.lastC = c;
    if (warnBox && result.warnings.length) warnBox.textContent = '⚠ ' + result.warnings.join('; ');

    displayH(hDisp, hKatex, result.hData, /*isUnbounded=*/true, c);
    pushBoundaryToPlot(QD.Direct.sampleBoundaryLaurent(c, F, 400), /*unbounded=*/true);
  }

  // ---------------------------------------------------------------------------
  // Display the computed h in both text and KaTeX forms.
  // ---------------------------------------------------------------------------
  function displayH(hDisp, hKatex, hData, isUnbounded, cValue) {
    const lines = ['h(w) = '];
    const polyPart = hData.polyPart || [];
    if (polyPart.length > 0) {
      const polyTerms = [];
      for (let l = 0; l < polyPart.length; l++) {
        const c = polyPart[l];
        if (Math.abs(c.re) < 1e-14 && Math.abs(c.im) < 1e-14) continue;
        polyTerms.push('  ' + complexToString(c) + (l === 0 ? '' : ' · w' + (l === 1 ? '' : '^' + l)));
      }
      if (polyTerms.length) lines.push.apply(lines, polyTerms);
    }
    for (const pole of hData.poles) {
      for (let k = 0; k < pole.principal.length; k++) {
        const C = pole.principal[k];
        const denPow = (k === 0) ? '' : '^' + (k + 1);
        lines.push('  ' + complexToString(C) + ' / (w − ' + complexToString(pole.a) + ')' + denPow);
      }
    }
    if (lines.length === 1) lines.push('  0');
    if (hDisp) hDisp.textContent = lines.join('\n');

    if (hKatex && window.katex) {
      try {
        let body = '';
        let first = true;
        for (let l = 0; l < polyPart.length; l++) {
          const c = polyPart[l];
          if (Math.abs(c.re) < 1e-14 && Math.abs(c.im) < 1e-14) continue;
          const cstr = complexToKatex(c);
          body += (first ? '' : ' + ') + (l === 0 ? cstr
                  : (l === 1 ? cstr + '\\,w' : cstr + '\\,w^{' + l + '}'));
          first = false;
        }
        for (const pole of hData.poles) {
          for (let k = 0; k < pole.principal.length; k++) {
            const C = pole.principal[k];
            const cstr = complexToKatex(C);
            const den = (k === 0)
              ? `(w - (${complexToKatex(pole.a)}))`
              : `(w - (${complexToKatex(pole.a)}))^{${k + 1}}`;
            body += (first ? '' : ' + ') + '\\frac{' + cstr + '}{' + den + '}';
            first = false;
          }
        }
        if (first) body = '0';
        window.katex.render('h(w) = ' + body, hKatex, { throwOnError: false });
      } catch (e) { /* leave text fallback */ }
    }
  }

  function pushBoundaryToPlot(pts, unbounded) {
    const setBdy = window.QD && window.QD.Direct && window.QD.Direct._setPlotBoundary;
    if (setBdy) setBdy(pts, { unbounded: !!unbounded });
  }

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
  function complexToString(c)  { return QD.Complex.format(c); }
  function complexToKatex(c)   { return QD.Complex.format(c); }

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
