// =============================================================================
// aqd-ui.js -- Frontend for the AQD tab.
//
// Mirrors the layout of the QD/LQD sidebar (Domain type / Riemann-map center /
// R input / h input / status / formula) but operates on AQD-specific state.
// The plot canvas is shared with the QD/LQD tab via the 'tab-changed' event.
//
// Stage 1 status: lazy-mounts the sidebar on first AQD-tab activation. R and
// h inputs use a generic 'rational function' widget (rfn-widget) that reads
// and writes the AQD state. No solver yet — that comes in Stage 2.
// =============================================================================
'use strict';

(function () {

  // ---------------------------------------------------------------------------
  // State (single source of truth for the AQD tab)
  // ---------------------------------------------------------------------------
  const aqdState = {
    domainMode: 'bounded',
    // R(w) = polyPart polynomial + Σ poles. Default: R = w  (classical QD).
    R: {
      polyPart: [{re: 0, im: 0}, {re: 1, im: 0}],
      poles: [],
    },
    // Quadrature function h: same shape as R (and as hData in the QD tab).
    h: {
      polyPart: [],
      poles: [ { a: '0', principal: ['1'] } ],
    },
    // Riemann-map center: 'auto' = centroid of h-poles, 'manual' = w0Manual.
    w0Mode: 'auto',
    w0Manual: '0',
  };

  // ---------------------------------------------------------------------------
  // Presets
  // ---------------------------------------------------------------------------
  const R_PRESETS = [
    { id: 'identity', label: 'R = w  (classical QD; α = 1)',
      R: { polyPart: ['0', '1'], poles: [] } },
    { id: 'w2-half',  label: 'R = w² / 2  (PQD ν=2; α = w)',
      R: { polyPart: ['0', '0', '0.5'], poles: [] } },
    { id: 'w3-third', label: 'R = w³ / 3  (PQD ν=3; α = w²)',
      R: { polyPart: ['0', '0', '0', '0.33333333333333333'], poles: [] } },
    { id: 'neg-1-w',  label: 'R = -1/w  (α = 1/w²; singular weight at 0)',
      R: { polyPart: [], poles: [ { a: '0', principal: ['-1'] } ] } },
    { id: 'rational-mix', label: 'R = w + 1/(w-2)  (α = 1 - 1/(w-2)²)',
      R: { polyPart: ['0', '1'],
           poles: [ { a: '2', principal: ['1'] } ] } },
  ];

  // ---------------------------------------------------------------------------
  // Lazy mount
  // ---------------------------------------------------------------------------
  let mounted = false;
  document.addEventListener('tab-changed', function (e) {
    if (e.detail && e.detail.tab === 'aqd' && !mounted) {
      mountAqdSidebar();
      mounted = true;
    }
  });

  function mountAqdSidebar() {
    const root = document.getElementById('controls-aqd');
    if (!root) return;
    root.innerHTML = '';   // clear the Stage-0 placeholder

    root.appendChild(makeDomainTypeCard());
    root.appendChild(makeRCard());
    root.appendChild(makeHCard());
    root.appendChild(makeW0Card());
    root.appendChild(makeStatusCard());
  }

  // ---------------------------------------------------------------------------
  // Card: Domain type
  // ---------------------------------------------------------------------------
  function makeDomainTypeCard() {
    const card = section('Domain type', `
      <div class="domain-mode-group">
        <div class="domain-mode-group-label">Algebraic QD (AQD)</div>
        <div class="row">
          <label><input type="radio" name="aqd-domain-mode" value="bounded" checked> Bounded</label>
          <label style="margin-left:14px;"><input type="radio" name="aqd-domain-mode" value="unbounded" disabled>
            Unbounded <span style="color:#888; font-size:11px;">(Stage 3)</span></label>
        </div>
      </div>
      <div class="hint" style="margin-top: 6px;">
        Weight <code>ρ = |α|²</code> with <code>α = R′</code> and
        <code>R</code> rational (LQDs excluded). The inverse problem solves
        for the rational <code>R∘φ</code> via Thm 6.4.1, then inverts
        <code>R</code> at sample points to recover <code>φ</code>.
      </div>
    `);
    card.querySelectorAll('input[name="aqd-domain-mode"]').forEach(r => {
      r.addEventListener('change', () => { aqdState.domainMode = r.value; });
    });
    return card;
  }

  // ---------------------------------------------------------------------------
  // Card: primitive R
  // ---------------------------------------------------------------------------
  function makeRCard() {
    const card = section('Primitive R(w)', `
      <div class="hint">
        R(w) = Σ<sub>k</sub> c<sub>k</sub> w<sup>k</sup> + Σ<sub>j</sub> Σ<sub>s</sub> D<sub>j,s</sub>/(w − p<sub>j</sub>)<sup>s</sup>.
        Then α = R′. (Pole shape automatically guarantees α has zero residues.)
      </div>
      <div class="row" style="margin-bottom: 8px;">
        <label>Preset:
          <select class="aqd-R-preset" style="width: 260px;">
            <option value="">— custom —</option>
            ${R_PRESETS.map(p => `<option value="${p.id}">${p.label}</option>`).join('')}
          </select>
        </label>
      </div>
      <div class="aqd-R-widget"></div>
    `);
    const widget = card.querySelector('.aqd-R-widget');
    renderRfnWidget(widget, aqdState.R, { name: 'R', allowPoly: true, polyVarLabel: 'w' });
    const preset = card.querySelector('.aqd-R-preset');
    preset.addEventListener('change', e => {
      const p = R_PRESETS.find(p => p.id === e.target.value);
      if (!p) return;
      aqdState.R = cloneRfnFromStrings(p.R);
      renderRfnWidget(widget, aqdState.R, { name: 'R', allowPoly: true, polyVarLabel: 'w' });
    });
    return card;
  }

  // ---------------------------------------------------------------------------
  // Card: quadrature function h
  // ---------------------------------------------------------------------------
  function makeHCard() {
    const card = section('Quadrature function h(w)', `
      <div class="hint">
        h(w) = Σ<sub>j</sub> Σ<sub>s</sub> C<sub>j,s</sub> / (w − a<sub>j</sub>)<sup>s</sup>.
        Use <code>1+2i</code> etc. for complex values.
      </div>
      <div class="aqd-h-widget"></div>
    `);
    const widget = card.querySelector('.aqd-h-widget');
    renderRfnWidget(widget, aqdState.h, { name: 'h', allowPoly: false });
    return card;
  }

  // ---------------------------------------------------------------------------
  // Card: Riemann-map center φ(0)
  // ---------------------------------------------------------------------------
  function makeW0Card() {
    const card = section('Riemann map center φ(0)', `
      <div class="row">
        <label><input type="radio" name="aqd-w0mode" value="auto" checked> Centroid of h-poles</label>
      </div>
      <div class="row">
        <label>
          <input type="radio" name="aqd-w0mode" value="manual"> Manual:
          <input type="text" class="cnum aqd-w0-manual" value="0" disabled>
        </label>
      </div>
    `);
    const manualInput = card.querySelector('.aqd-w0-manual');
    card.querySelectorAll('input[name="aqd-w0mode"]').forEach(r => {
      r.addEventListener('change', () => {
        aqdState.w0Mode = r.value;
        manualInput.disabled = (r.value !== 'manual');
      });
    });
    manualInput.addEventListener('input', () => {
      aqdState.w0Manual = manualInput.value;
    });
    return card;
  }

  // ---------------------------------------------------------------------------
  // Card: status (Stage 1 placeholder; Stage 2 will hook the solver here)
  // ---------------------------------------------------------------------------
  function makeStatusCard() {
    const card = section('Status', `
      <div id="aqd-status">Idle. Solver lands in Stage 2.</div>
    `);
    return card;
  }

  // ===========================================================================
  // Generic rational-function widget
  // ---------------------------------------------------------------------------
  // Renders a list of (pole, multiplicity, residues) plus an optional
  // polynomial-part editor. Writes back into `state` in place.
  //
  // state shape (strings for editability; converted to Complex on demand):
  //   { polyPart: [string, ...], poles: [{ a: string, principal: [string, ...] }, ...] }
  // ===========================================================================
  function renderRfnWidget(container, state, opts) {
    container.innerHTML = '';
    const polesList = el('div', 'aqd-rfn-poles');
    container.appendChild(polesList);
    const addBtn = el('button', 'small');
    addBtn.textContent = '+ Add pole';
    addBtn.addEventListener('click', () => {
      state.poles.push({ a: '0', principal: ['1'] });
      renderRfnWidget(container, state, opts);
    });
    container.appendChild(addBtn);

    state.poles.forEach((pole, idx) => {
      polesList.appendChild(renderRfnPoleBlock(state, idx, opts, container));
    });

    if (opts.allowPoly) {
      const polySection = el('div');
      polySection.style.marginTop = '12px';
      polySection.style.paddingTop = '10px';
      polySection.style.borderTop = '1px dashed #d8dbe0';
      polySection.innerHTML = `
        <div class="hint">
          Polynomial part: Σ<sub>k=0..m</sub> c<sub>k</sub> ${opts.polyVarLabel || 'w'}<sup>k</sup>.
        </div>
        <div class="row">
          <label>Degree m:
            <input type="number" class="aqd-poly-deg" min="-1" max="6"
                   value="${state.polyPart.length - 1}" style="width: 56px;">
          </label>
          <span style="font-size:11px; color:#777;">(−1 = no polynomial part)</span>
        </div>
        <div class="aqd-poly-coefs"></div>
      `;
      container.appendChild(polySection);

      const degInput = polySection.querySelector('.aqd-poly-deg');
      degInput.addEventListener('input', () => {
        let deg = parseInt(degInput.value, 10);
        if (Number.isNaN(deg) || deg < -1) deg = -1;
        const newLen = deg + 1;
        if (state.polyPart.length < newLen) {
          while (state.polyPart.length < newLen) state.polyPart.push('0');
        } else if (state.polyPart.length > newLen) {
          state.polyPart.length = newLen;
        }
        renderRfnWidget(container, state, opts);
      });

      const coefBox = polySection.querySelector('.aqd-poly-coefs');
      state.polyPart.forEach((c, k) => {
        const row = el('div', 'row');
        const lab = document.createElement('label');
        lab.innerHTML = `c<sub>${k}</sub> = `;
        const inp = el('input', 'cnum');
        inp.type = 'text';
        inp.value = typeof c === 'string' ? c : (c.re + (c.im < 0 ? '' : '+') + c.im + 'i');
        inp.addEventListener('input', () => { state.polyPart[k] = inp.value; });
        lab.appendChild(inp);
        row.appendChild(lab);
        coefBox.appendChild(row);
      });
    }
  }

  function renderRfnPoleBlock(state, idx, opts, container) {
    const pole = state.poles[idx];
    const block = el('div', 'pole');
    block.innerHTML = `
      <div class="pole-header">
        <span class="pole-num">${opts.name}-pole #${idx + 1}</span>
        <div>
          <button class="small aqd-pole-up" title="Move up">↑</button>
          <button class="small aqd-pole-dn" title="Move down">↓</button>
          <button class="small danger aqd-pole-rm">Remove</button>
        </div>
      </div>
      <div class="row">
        <label>Location:
          <input type="text" class="cnum aqd-pole-a" value="${pole.a}">
        </label>
        <label>Order m:
          <input type="number" class="aqd-pole-order" min="1" max="6" value="${pole.principal.length}" style="width: 56px;">
        </label>
      </div>
      <div class="residues"></div>
    `;
    const residues = block.querySelector('.residues');
    pole.principal.forEach((C, s) => {
      const row = el('div', 'residue-row');
      const lab = el('span', 'label-fixed');
      lab.innerHTML = `${opts.name === 'R' ? 'D' : 'C'}<sub>${idx+1},${s+1}</sub>`;
      const inp = el('input', 'cnum');
      inp.type = 'text';
      inp.value = typeof C === 'string' ? C : (C.re + (C.im < 0 ? '' : '+') + C.im + 'i');
      inp.addEventListener('input', () => { pole.principal[s] = inp.value; });
      row.appendChild(lab); row.appendChild(inp);
      residues.appendChild(row);
    });
    block.querySelector('.aqd-pole-a').addEventListener('input', e => {
      pole.a = e.target.value;
    });
    block.querySelector('.aqd-pole-order').addEventListener('input', e => {
      let m = parseInt(e.target.value, 10);
      if (Number.isNaN(m) || m < 1) m = 1;
      while (pole.principal.length < m) pole.principal.push('0');
      pole.principal.length = m;
      renderRfnWidget(container, state, opts);
    });
    block.querySelector('.aqd-pole-rm').addEventListener('click', () => {
      state.poles.splice(idx, 1);
      renderRfnWidget(container, state, opts);
    });
    block.querySelector('.aqd-pole-up').addEventListener('click', () => {
      if (idx > 0) { const t = state.poles[idx-1]; state.poles[idx-1] = pole; state.poles[idx] = t; renderRfnWidget(container, state, opts); }
    });
    block.querySelector('.aqd-pole-dn').addEventListener('click', () => {
      if (idx < state.poles.length - 1) { const t = state.poles[idx+1]; state.poles[idx+1] = pole; state.poles[idx] = t; renderRfnWidget(container, state, opts); }
    });
    return block;
  }

  // ---------------------------------------------------------------------------
  // Small helpers
  // ---------------------------------------------------------------------------
  function section(title, innerHTML) {
    const sec = document.createElement('section');
    sec.className = 'card';
    sec.innerHTML = `<h2>${title}</h2>${innerHTML}`;
    return sec;
  }
  function el(tag, cls) {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    return e;
  }
  // Clone a preset (where residues/polyPart entries are string forms) into a
  // fresh state object. Keeps strings so the UI shows them verbatim.
  function cloneRfnFromStrings(R) {
    return {
      polyPart: (R.polyPart || []).map(c => String(c)),
      poles: (R.poles || []).map(p => ({
        a: String(p.a),
        principal: p.principal.map(c => String(c)),
      })),
    };
  }

  // Expose a tiny inspection hook so Stage 2 can read aqdState.
  if (typeof window !== 'undefined') {
    window.QD = window.QD || {};
    window.QD.Aqd = window.QD.Aqd || {};
    window.QD.Aqd._uiState = aqdState;       // intentionally underscored
  }

}());
