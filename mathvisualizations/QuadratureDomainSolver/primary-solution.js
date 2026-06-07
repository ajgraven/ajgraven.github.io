// @ts-check
// =============================================================================
// QD.PrimarySolution — typed envelope for the current inverse-solver result,
// with a small subscribe/publish shim so the Schwarz / Sphere / Param-slice
// tabs (and any future readers) don't have to reach into ui.js's internal
// state.current directly.
//
// The envelope shape is exactly what ui.js used to publish onto
// state.current: { success, primary, alternates, hData, w0Used, cUsed,
// unbounded, attempts, selectedIdx }. The ONLY behavior change is that
// readers go through QD.PrimarySolution.get() and subscribers are notified
// when ui.js republishes after a fresh solve or after the alternates list
// grows.
//
// Ownership: ui.js is still the canonical writer (it owns the form data
// and triggers solves); this module is a published mirror. A future
// refactor can flip ownership without changing the reader API.
// =============================================================================

/**
 * @typedef {{ re: number, im: number }} ComplexC
 *   Plain {re, im} complex value — the universal scalar type across the codebase.
 *
 * @typedef {Object} PrimaryEnvelope
 *   Shape published by ui.js's solveAndRender() and read by cross-tab subscribers.
 *   This is exactly what state.current used to hold, lifted into a typed contract.
 * @property {boolean} success                 — true iff a valid solve completed
 * @property {Object} [primary]                — { phi, identity, identityOK, univalent, … }
 * @property {Object[]} [alternates]           — additional valid φ's found in the alt-search
 * @property {Object} [hData]                  — the quadrature data the solver was asked about
 * @property {ComplexC} [w0Used]               — φ(0) used by this solve
 * @property {number} [cUsed]                  — conformal radius for unbounded families
 * @property {boolean} [unbounded]
 * @property {Object[]} [attempts]
 * @property {string} [error]                  — populated when success === false
 * @property {Object} [criticalSet]            — lazily-populated cache of zeros of φ'
 * @property {Object} [geomProps]              — async geometric classification:
 *                                               { convex, starLike, spiralLike, … }
 *                                               from QD.classifyUnivalence (null until computed)
 * @property {Object} [cuspProps]              — async boundary-singularity analysis:
 *                                               { cusps: [{thetaDeg, type:[p,q], orderM,
 *                                               dist, isCusp, …}], … } from
 *                                               QD.classifyCusps (null until computed)
 * @property {boolean} [regimeSwitched]        — §23: PQD auto singular⇄non-singular switch
 *
 * @callback PrimarySubscriber
 * @param {PrimaryEnvelope|null} envelope
 * @returns {void}
 */

(function (global) {
  'use strict';

  // Same QD-namespace resolution idiom every solver file uses.
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD) return;

  /** @type {PrimaryEnvelope|null} */
  let _envelope = null;       // The current envelope, or null if nothing solved yet.
  /** @type {PrimarySubscriber[]} */
  let _subs = [];             // Array of subscriber callbacks.
  let _publishToken = 0;      // Monotonic counter for diagnostics.

  function _notify() {
    // Defensive copy so an unsubscribe-during-notify doesn't perturb iteration.
    // Errors in one subscriber don't break the others.
    const list = _subs.slice();
    for (let i = 0; i < list.length; i++) {
      try { list[i](_envelope); }
      catch (e) {
        if (typeof console !== 'undefined') console.error('PrimarySolution subscriber error:', e);
      }
    }
  }

  QD.PrimarySolution = {
    /**
     * Read the current envelope. Returns null if no solve has happened yet.
     * @returns {PrimaryEnvelope|null}
     */
    get() { return _envelope; },

    /**
     * True iff there is an envelope AND it represents a successful solve.
     * @returns {boolean}
     */
    hasSolution() { return !!(_envelope && _envelope.success); },

    /**
     * Subscribe to publish events. The handler receives the new envelope
     * (or null on clear). Returns an unsubscribe function.
     * @param {PrimarySubscriber} handler
     * @returns {() => void}
     */
    subscribe(handler) {
      if (typeof handler !== 'function') throw new TypeError('subscribe handler must be a function');
      _subs.push(handler);
      return function unsubscribe() {
        const idx = _subs.indexOf(handler);
        if (idx >= 0) _subs.splice(idx, 1);
      };
    },

    /**
     * Replace the envelope wholesale and notify subscribers.
     * @param {PrimaryEnvelope|null} envelope
     */
    publish(envelope) {
      _envelope = envelope || null;
      _publishToken++;
      _notify();
    },

    /**
     * Patch fields on the current envelope in-place (used when ui.js
     * appends a newly-found alternate). Notifies subscribers. No-op if
     * no envelope is currently published.
     * @param {Partial<PrimaryEnvelope>} patch
     */
    update(patch) {
      if (!_envelope || !patch) return;
      Object.assign(_envelope, patch);
      _publishToken++;
      _notify();
    },

    /** Drop the envelope (e.g. on full reset). Notifies subscribers with null. */
    clear() {
      _envelope = null;
      _publishToken++;
      _notify();
    },

    /** Diagnostic — monotonic counter of publish events. Useful in tests. */
    _publishToken() { return _publishToken; },
  };

})(typeof globalThis !== 'undefined' ? globalThis : this);
