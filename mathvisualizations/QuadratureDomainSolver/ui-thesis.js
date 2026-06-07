// =============================================================================
// ui-thesis.js  —  Thesis-example gallery + analytic-oracle card (#8).
//
// Installs via QD_UI.installThesis(uiCtx). Populates the #thesis-select gallery
// from QD.ThesisExamples; selecting one calls uiCtx.loadThesisExample(ex) (the
// ui.js-internal half: switch family, apply config, frame the view, enable the
// annotated-phenomena overlay, solve). On every solve while an example is active
// it runs QD.checkOracle and renders the #oracle-card (computed vs expected, with
// pass/warn/fail). The c* row is opt-in behind a button (heavy + off-thread via
// the primary-solver worker). A manual edit or a family-preset load fires the
// `qd-customized` event, which clears the card. UI-only; no solver math.
// =============================================================================

(function () {
  'use strict';

  function esc(s) {
    return String(s).replace(/[&<>"]/g, c =>
      ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  }
  function fmt(x) {
    if (typeof x !== 'number') return esc(String(x));
    return (Math.abs(x) >= 1e4 || (x !== 0 && Math.abs(x) < 1e-3))
      ? x.toExponential(2) : String(+(+x).toPrecision(5));
  }
  function icon(status) {
    if (status === 'pass') return '<span class="ok">✓</span>';
    return status === 'warn' ? '<span class="warn">⚠</span>' : '<span class="err">✗</span>';
  }

  function installThesis(ctx) {
    const QD = window.QD;
    const $ = ctx.$;
    const examples = (QD && QD.ThesisExamples) || [];

    const sel     = $('#thesis-select');
    const card    = $('#oracle-card');
    const content = $('#oracle-content');
    if (!sel || !card || !content) return {};

    // Oracle-card prose lives in QD.Strings.oracle (ui-strings.js).
    const STR = (QD.Strings && QD.Strings.oracle) || {};

    let active = null;       // the example currently on exhibit (null ⇒ card hidden)
    let suppressClear = false; // true while we're applying a thesis load

    // Populate the gallery.
    sel.innerHTML = '<option value="">' + esc(STR.pickExample || '— pick an example —') + '</option>' +
      examples.map(e => `<option value="${esc(e.id)}">${esc(e.label)}</option>`).join('');

    sel.addEventListener('change', (e) => {
      const ex = examples.find(x => x.id === e.target.value);
      if (!ex) { clearCard(); return; }
      // loadThesisExample applies config (setC/setQ/…) which fires `qd-customized`
      // synchronously; suppress the clear so this load doesn't cancel itself, then
      // (re-)establish the active example afterward.
      suppressClear = true;
      ctx.loadThesisExample(ex);     // family + config + view + overlay + solve
      suppressClear = false;
      active = ex;
      sel.value = ex.id;
      card.classList.remove('hidden');
      content.innerHTML = `<div class="hint">${esc(ex.blurb || '')}</div>` +
        `<div class="key" style="margin-top:6px;">${esc(STR.pending || 'solving… oracle pending')}</div>`;
    });

    // Leaving the example (manual edit / family preset) clears the card.
    document.addEventListener('qd-customized', clearCard);
    function clearCard() {
      if (suppressClear) return;
      active = null;
      sel.value = '';
      card.classList.add('hidden');
    }

    // Help popover on the card header.
    if (QD.QoL && QD.QoL.attachHelp && STR.help) {
      const h = card.querySelector('h2');
      if (h) QD.QoL.attachHelp(h, STR.help);
    }

    // Re-check the oracle whenever a fresh primary solution lands for the active
    // example. Cheap rows auto-run; c* stays behind the button.
    if (QD && QD.PrimarySolution && QD.PrimarySolution.subscribe) {
      QD.PrimarySolution.subscribe((env) => {
        if (!active || !active.oracle) return;
        if (!env || !env.success || !env.primary || !env.primary.phi || !env.hData) {
          content.innerHTML = `<div class="hint">${esc(active.blurb || '')}</div>` +
            `<div class="warn">${esc(STR.noValid || 'no valid solution to check against the oracle')}</div>`;
          return;
        }
        const ex = active;
        Promise.resolve(QD.checkOracle(env.primary.phi, env.hData, ex.oracle, { includeCmax: false }))
          .then((res) => { if (active === ex) renderRows(ex, res, env); })
          .catch(() => { /* leave the pending text */ });
      });
    }

    function renderRows(ex, res, env) {
      const rows = res.rows.map(r =>
        `<div class="geom-row">${icon(r.status)} <span class="key">${esc(r.name)}:</span> ` +
        `${fmt(r.computed)} <span class="key">(exp ${fmt(r.expected)})</span></div>`).join('');
      const summary = res.allPass
        ? '<span class="ok">' + esc(STR.matches || '✓ matches the analytic oracle') + '</span>'
        : '<span class="warn">' + esc(STR.someDiffer || '⚠ some rows differ from the oracle') + '</span>';
      const cmaxBtn = (ex.oracle.cMax != null)
        ? `<button id="oracle-cmax-btn" class="oracle-cmax" style="margin-top:8px;">${esc(STR.verifyCmax || 'Verify c* (slow)')}</button>`
        : '';
      content.innerHTML =
        `<div class="hint">${esc(ex.blurb || '')}</div>` +
        `<div style="margin:6px 0;">${summary}</div>` + rows + cmaxBtn;
      const btn = $('#oracle-cmax-btn');
      if (btn) btn.addEventListener('click', () => verifyCmax(ex, env, btn));
    }

    async function verifyCmax(ex, env, btn) {
      btn.disabled = true;
      btn.textContent = STR.estimating || 'estimating c*…';
      // Off-thread via the primary-solver worker, exactly like the "Estimate max c"
      // button — so the (heavy) bracket search doesn't block the UI.
      const PSW = QD.PrimarySolverWorker;
      const solveFn = (PSW && typeof PSW.solve === 'function')
        ? (h, o) => PSW.solve(h, o)
        : (h, o) => QD.solveInverseQD(h, o);
      try {
        const res = await QD.checkOracle(env.primary.phi, env.hData, ex.oracle,
          { includeCmax: true, solveFn });
        if (active === ex) renderRows(ex, res, env);   // now includes the c* rows
      } catch (e) {
        btn.disabled = false;
        btn.textContent = STR.verifyCmax || 'Verify c* (slow)';
      }
    }

    return {};
  }

  window.QD_UI = window.QD_UI || {};
  window.QD_UI.installThesis = installThesis;
})();
