// =============================================================================
// ui-url-state.js -- URL/hash state serialization (B1) for the Inverse tab.
//
// Extracted from ui.js by the Phase-3 UI modularization (item E). Pulled out as
// a factory on the shared QD_UI namespace so it receives its ui.js closures
// (state, MODES, PRESETS, $, and the four apply-helpers) via dependency
// injection rather than relying on shared script-scope variables (which don't
// cross <script> tags). See ui-domain-plot.js for the original template.
//
// ui.js installs it via:
//   const { writeUrlState, applyUrlState } = window.QD_UI.installUrlState(uiCtx);
//
// The two returned functions are then used exactly as before: writeUrlState is
// called after each solve / tab switch (rAF-coalesced), applyUrlState restores
// a shared/bookmarked config on load. Bodies are identical to what previously
// lived in ui.js — only the dependency source changed (module scope → `ui`).
// =============================================================================

(function (global) {
  'use strict';
  global.QD_UI = global.QD_UI || {};

  global.QD_UI.installUrlState = function installUrlState(ui) {
    // Stable deps destructured to their original ui.js names, so the bodies
    // below are verbatim. All exist on `ui` by the time this factory runs
    // (installed near the end of ui.js, after every dependency is defined).
    const state             = ui.state;
    const MODES             = ui.MODES;
    const PRESETS           = ui.PRESETS;
    const $                 = ui.$;
    const applyModeVisuals  = ui.applyModeVisuals;
    const setC              = ui.setC;
    const setQ              = ui.setQ;
    const parseAndApplyHText = ui.parseAndApplyHText;

    // -----------------------------------------------------------------------
    // URL/hash state (B1) — shareable, bookmarkable, reload-restorable config.
    //
    // We serialize the user-meaningful inputs (mode, the h(w) text, the
    // normalization gauges w₀/c/α/q, aggressiveness, and the active tab) into
    // location.hash. The h-text round-trips both poles AND the polynomial part
    // (refreshHText → formatH), so it alone captures the full quadrature data;
    // parseAndApplyHText rebuilds the structured grid from it on restore. We use
    // history.replaceState (not assignment to location.hash) so writing the URL
    // never pushes a back-button entry or re-navigates.
    // -----------------------------------------------------------------------
    function _activeTabId() {
      const el = document.querySelector('.tab-btn.active');
      return (el && el.dataset.tab) || 'qd';
    }

    let _writeUrlScheduled = false;
    function writeUrlState() {
      // Coalesce bursts (a slider drag fires many solves) into one history write
      // per frame.
      if (_writeUrlScheduled) return;
      _writeUrlScheduled = true;
      const raf = (typeof requestAnimationFrame === 'function')
        ? requestAnimationFrame
        : (fn) => setTimeout(() => fn(), 16);
      raf(() => {
        _writeUrlScheduled = false;
        try {
          const p = new URLSearchParams();
          p.set('mode', state.mode);
          const hText = ($('#h-text') && $('#h-text').value || '').trim();
          if (hText) p.set('h', hText);
          if (state.w0Mode) p.set('w0m', state.w0Mode);
          if (state.w0Manual) p.set('w0', state.w0Manual);
          if (state.c != null) p.set('c', String(state.c));
          if (state.alpha != null && state.alpha !== 1) p.set('a', String(state.alpha));
          if (state.q && state.q !== '0') p.set('q', state.q);
          if (state.aggressiveness) p.set('agg', state.aggressiveness);
          const tab = _activeTabId();
          if (tab && tab !== 'qd') p.set('tab', tab);
          const hash = '#' + p.toString();
          // Avoid redundant history churn when nothing changed.
          if (hash !== location.hash) {
            history.replaceState(null, '', location.pathname + location.search + hash);
          }
        } catch (e) { /* never let URL bookkeeping break the app */ }
      });
    }

    // Restore state from location.hash on load. Returns true if a hash was applied
    // (so the caller can skip the default-config solve). Sets mode + gauges FIRST,
    // then the h-text, then parses it (which schedules the solve), then the tab.
    function applyUrlState() {
      let hash = (location.hash || '').replace(/^#/, '');
      if (!hash) return false;
      let p;
      try { p = new URLSearchParams(hash); } catch (e) { return false; }
      if (![...p.keys()].length) return false;

      // 1. Mode (drives card visibility + which gauges matter). applyModeVisuals
      //    forces α back to 1 for non-PQD modes, so it must run BEFORE we set α.
      const mode = p.get('mode');
      if (mode && MODES[mode]) {
        state.mode = mode;
        applyModeVisuals();   // also syncs the compact domain-type control
      }
      // 2. Gauges.
      if (p.has('a')) {
        const a = +p.get('a');
        if (a > 0 && a !== 1) { state.alpha = a; const inp = $('#alpha-input'); if (inp) inp.value = String(a); }
      }
      if (p.has('c')) { const c = +p.get('c'); if (c > 0) setC(c); }
      if (p.has('w0m')) {
        const m = p.get('w0m');
        if (m === 'auto' || m === 'manual') {
          state.w0Mode = m;
          const r = document.querySelector(`input[name="w0mode"][value="${m}"]`);
          if (r) r.checked = true;
          const wManual = $('#w0-manual');
          if (wManual) wManual.disabled = (m !== 'manual');
        }
      }
      if (p.has('w0')) {
        state.w0Manual = p.get('w0');
        const wManual = $('#w0-manual');
        if (wManual) wManual.value = state.w0Manual;
      }
      if (p.has('q')) setQ(p.get('q'));
      if (p.has('agg') && PRESETS[p.get('agg')]) {
        state.aggressiveness = p.get('agg');
        const aggSel = $('#aggressiveness');
        if (aggSel) aggSel.value = state.aggressiveness;
      }
      // 3. h(w): set the text and parse it (rebuilds the pole grid + poly + solves).
      if (p.has('h')) {
        const inp = $('#h-text');
        if (inp) { inp.value = p.get('h'); parseAndApplyHText(); }
      }
      // 4. Active tab (deferred a tick so the QD solve kicks off first).
      const tab = p.get('tab');
      if (tab && tab !== 'qd') {
        const btn = document.querySelector(`.tab-btn[data-tab="${tab}"]`);
        if (btn) setTimeout(() => btn.click(), 0);
      }
      return true;
    }

    return { writeUrlState, applyUrlState };
  };
})(typeof window !== 'undefined' ? window : globalThis);
