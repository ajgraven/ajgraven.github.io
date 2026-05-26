// =============================================================================
// qol.js -- Shared quality-of-life primitives (HANDOFF #33).
//
// Three small DOM helpers, exposed on QD.QoL, used everywhere across the UI:
//
//   attachHelp(headerEl, htmlOrFn, opts)
//     Adds a "?" button to a card's <h2>. Click toggles a popover with the
//     text (string or () => string). At most one popover open across the
//     page. Click-outside + Esc closes.
//
//   attachHoverTooltip(canvasEl, formatter, opts)
//     Registers a throttled (rAF-coalesced) mousemove handler on canvasEl.
//     formatter(cssX, cssY) returns null (hide) or an HTML string.
//     A single page-level tooltip div is reused across canvases. Hides on
//     mouseleave / mouseout. Position auto-flips to stay in viewport.
//
//   copyButton(getText, opts)
//     Returns a small button element. Click → navigator.clipboard.writeText
//     + brief toast confirmation. Falls back gracefully if clipboard API
//     unavailable (selects and shows fallback hint).
//
// Plus:
//   openShortcutsOverlay() / closeShortcutsOverlay()
//     Page-level shortcut cheatsheet, anchored to the page (not a card).
//   wireGlobalKeyboardShortcuts()
//     Esc → close all popovers/tooltips; '?' → toggle shortcuts overlay.
//
// No external dependencies. Safe to load before / after other QD modules.
// =============================================================================

(function (global) {
  'use strict';

  const QD = global.QD = global.QD || {};
  const QoL = QD.QoL = {};

  // ---------------------------------------------------------------------------
  // Singletons
  // ---------------------------------------------------------------------------
  let _activePopover = null;          // currently-open help popover element
  let _tooltipEl     = null;          // shared hover-tooltip div
  let _shortcutsEl   = null;          // shortcuts overlay element

  function _ensureTooltipEl() {
    if (_tooltipEl) return _tooltipEl;
    const d = document.createElement('div');
    d.className = 'hover-tooltip';
    d.style.display = 'none';
    d.setAttribute('role', 'tooltip');
    document.body.appendChild(d);
    _tooltipEl = d;
    return d;
  }

  // ---------------------------------------------------------------------------
  // attachHelp
  // ---------------------------------------------------------------------------
  function attachHelp(headerEl, htmlOrFn, opts) {
    if (!headerEl) return null;
    // Avoid double-attaching.
    if (headerEl.querySelector('.help-btn')) return headerEl.querySelector('.help-btn');
    opts = opts || {};
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'help-btn';
    btn.textContent = '?';
    btn.title = 'Show help for this section';
    btn.setAttribute('aria-label', 'Help');
    btn.setAttribute('aria-haspopup', 'dialog');
    btn.addEventListener('click', (e) => {
      e.stopPropagation();
      _toggleHelpPopover(btn, htmlOrFn, opts);
    });
    headerEl.appendChild(btn);
    return btn;
  }

  function _toggleHelpPopover(anchorEl, htmlOrFn, opts) {
    if (_activePopover && _activePopover._anchor === anchorEl) {
      _closeActivePopover();
      return;
    }
    _closeActivePopover();
    const pop = document.createElement('div');
    pop.className = 'help-popover';
    pop.setAttribute('role', 'dialog');
    let html = '';
    try {
      html = (typeof htmlOrFn === 'function') ? htmlOrFn() : String(htmlOrFn || '');
    } catch (e) { html = '<em>(help text error)</em>'; }
    pop.innerHTML = html;
    pop._anchor = anchorEl;
    document.body.appendChild(pop);
    _positionPopoverNear(pop, anchorEl);
    _activePopover = pop;
    // Click outside closes.
    setTimeout(() => {
      document.addEventListener('click', _maybeCloseOnOutsideClick, true);
    }, 0);
  }

  function _maybeCloseOnOutsideClick(e) {
    if (!_activePopover) return;
    if (_activePopover.contains(e.target)) return;
    if (_activePopover._anchor === e.target) return;
    _closeActivePopover();
  }

  function _closeActivePopover() {
    if (!_activePopover) return;
    _activePopover.remove();
    _activePopover = null;
    document.removeEventListener('click', _maybeCloseOnOutsideClick, true);
  }

  function _positionPopoverNear(pop, anchorEl) {
    const ar = anchorEl.getBoundingClientRect();
    const pr = pop.getBoundingClientRect();
    const vw = window.innerWidth, vh = window.innerHeight;
    // Default: just below the button, left edge aligned with button.
    let left = ar.left;
    let top  = ar.bottom + 6;
    // Flip horizontally if overflow.
    if (left + pr.width > vw - 8) left = Math.max(8, vw - pr.width - 8);
    // Flip vertically if overflow.
    if (top + pr.height > vh - 8) top = Math.max(8, ar.top - pr.height - 6);
    pop.style.left = (left + window.scrollX) + 'px';
    pop.style.top  = (top  + window.scrollY) + 'px';
  }

  // ---------------------------------------------------------------------------
  // attachHoverTooltip
  // ---------------------------------------------------------------------------
  function attachHoverTooltip(canvasEl, formatter, opts) {
    if (!canvasEl || typeof formatter !== 'function') return;
    opts = opts || {};
    let pendingFrame = 0;
    let lastEvt = null;

    function onMove(e) {
      lastEvt = e;
      if (pendingFrame) return;
      pendingFrame = requestAnimationFrame(() => {
        pendingFrame = 0;
        if (!lastEvt) return;
        const rect = canvasEl.getBoundingClientRect();
        const x = lastEvt.clientX - rect.left;
        const y = lastEvt.clientY - rect.top;
        let result;
        try { result = formatter(x, y, canvasEl); }
        catch (err) { result = null; }
        if (!result) { _hideTooltip(); return; }
        const tip = _ensureTooltipEl();
        tip.innerHTML = (typeof result === 'string') ? result : (result.html || '');
        tip.style.display = 'block';
        _positionTooltipNear(tip, lastEvt.clientX, lastEvt.clientY);
      });
    }

    function onLeave() {
      _hideTooltip();
    }

    canvasEl.addEventListener('mousemove', onMove);
    canvasEl.addEventListener('mouseleave', onLeave);
    canvasEl.addEventListener('mouseout', onLeave);
    // Return a detach handle for completeness.
    return function detach() {
      canvasEl.removeEventListener('mousemove', onMove);
      canvasEl.removeEventListener('mouseleave', onLeave);
      canvasEl.removeEventListener('mouseout', onLeave);
    };
  }

  function _hideTooltip() {
    if (_tooltipEl) _tooltipEl.style.display = 'none';
  }

  function _positionTooltipNear(tip, cx, cy) {
    const r = tip.getBoundingClientRect();
    const vw = window.innerWidth, vh = window.innerHeight;
    let left = cx + 14;
    let top  = cy + 14;
    if (left + r.width  > vw - 8) left = Math.max(8, cx - 14 - r.width);
    if (top  + r.height > vh - 8) top  = Math.max(8, cy - 14 - r.height);
    tip.style.left = (left + window.scrollX) + 'px';
    tip.style.top  = (top  + window.scrollY) + 'px';
  }

  // ---------------------------------------------------------------------------
  // copyButton
  // ---------------------------------------------------------------------------
  function copyButton(getText, opts) {
    opts = opts || {};
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'copy-btn';
    btn.title = opts.title || 'Copy to clipboard';
    btn.setAttribute('aria-label', btn.title);
    btn.textContent = opts.label || '⧉';
    btn.addEventListener('click', async (e) => {
      e.stopPropagation();
      let text;
      try { text = (typeof getText === 'function') ? getText() : String(getText || ''); }
      catch (err) { text = ''; }
      if (!text) { _showToast('Nothing to copy', btn); return; }
      let ok = false;
      try {
        if (navigator.clipboard && navigator.clipboard.writeText) {
          await navigator.clipboard.writeText(text);
          ok = true;
        }
      } catch (err) { /* fall through to fallback */ }
      if (!ok) {
        // Fallback for non-secure contexts: use a temporary textarea.
        try {
          const ta = document.createElement('textarea');
          ta.value = text;
          ta.style.position = 'fixed'; ta.style.opacity = '0';
          document.body.appendChild(ta);
          ta.select();
          ok = document.execCommand && document.execCommand('copy');
          document.body.removeChild(ta);
        } catch (err) { ok = false; }
      }
      _showToast(ok ? 'Copied' : 'Copy failed', btn);
    });
    return btn;
  }

  function _showToast(msg, anchorEl) {
    const t = document.createElement('div');
    t.className = 'copy-toast';
    t.textContent = msg;
    document.body.appendChild(t);
    const r = anchorEl.getBoundingClientRect();
    t.style.left = (r.left + window.scrollX) + 'px';
    t.style.top  = (r.bottom + 6 + window.scrollY) + 'px';
    setTimeout(() => {
      t.classList.add('fade');
      setTimeout(() => t.remove(), 350);
    }, 750);
  }

  // ---------------------------------------------------------------------------
  // Shortcuts overlay
  // ---------------------------------------------------------------------------
  function openShortcutsOverlay(items) {
    closeShortcutsOverlay();
    const list = items || _defaultShortcuts();
    const html = `
      <div class="help-popover" style="max-width:340px;">
        <h3 style="margin:0 0 6px 0; font-size:13px;">Keyboard shortcuts</h3>
        <table class="shortcuts-table">
          ${list.map(it => `<tr><td><kbd>${it.key}</kbd></td><td>${it.desc}</td></tr>`).join('')}
        </table>
        <div class="hint" style="margin-top:6px;">Press Esc or ? to dismiss.</div>
      </div>
    `;
    const wrap = document.createElement('div');
    wrap.className = 'shortcuts-overlay';
    wrap.innerHTML = html;
    document.body.appendChild(wrap);
    _shortcutsEl = wrap;
  }
  function closeShortcutsOverlay() {
    if (_shortcutsEl) { _shortcutsEl.remove(); _shortcutsEl = null; }
  }
  function _defaultShortcuts() {
    return [
      { key: '?',   desc: 'Show / hide this shortcut list' },
      { key: 'Esc', desc: 'Close help popovers and tooltips' },
      { key: 'Enter', desc: 'In a Param-slice axis field: render slice' },
    ];
  }

  function wireGlobalKeyboardShortcuts() {
    document.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') {
        let handled = false;
        if (_activePopover) { _closeActivePopover(); handled = true; }
        if (_shortcutsEl)   { closeShortcutsOverlay(); handled = true; }
        if (_tooltipEl && _tooltipEl.style.display !== 'none') {
          _hideTooltip(); handled = true;
        }
        if (handled) e.preventDefault();
        return;
      }
      if (e.key === '?' && !_isTypingTarget(e.target)) {
        if (_shortcutsEl) closeShortcutsOverlay();
        else              openShortcutsOverlay();
        e.preventDefault();
      }
    });
  }

  function _isTypingTarget(el) {
    if (!el) return false;
    const tag = (el.tagName || '').toUpperCase();
    return tag === 'INPUT' || tag === 'TEXTAREA' || el.isContentEditable;
  }

  // ---------------------------------------------------------------------------
  // Exports
  // ---------------------------------------------------------------------------
  // Shared HTML escaper — consolidated from previously-duplicated copies
  // in ui.js and param-slice-ui.js (HANDOFF #35). Escapes the full
  // attribute-safe set so it's correct in both text-content and
  // attribute-value positions.
  function escapeHTML(s) {
    return String(s).replace(/[&<>"']/g, c => (
      { '&':'&amp;', '<':'&lt;', '>':'&gt;', '"':'&quot;', "'":'&#39;' }[c]));
  }

  QoL.attachHelp                = attachHelp;
  QoL.escapeHTML                = escapeHTML;
  QoL.attachHoverTooltip        = attachHoverTooltip;
  QoL.copyButton                = copyButton;
  QoL.openShortcutsOverlay      = openShortcutsOverlay;
  QoL.closeShortcutsOverlay     = closeShortcutsOverlay;
  QoL.wireGlobalKeyboardShortcuts = wireGlobalKeyboardShortcuts;

  // Auto-wire on DOM ready unless the loader explicitly opts out.
  if (typeof document !== 'undefined') {
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', wireGlobalKeyboardShortcuts);
    } else {
      wireGlobalKeyboardShortcuts();
    }
  }

  if (typeof module !== 'undefined' && module.exports) module.exports = QoL;
})(typeof window !== 'undefined' ? window : globalThis);
