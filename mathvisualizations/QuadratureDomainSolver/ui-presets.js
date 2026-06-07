// =============================================================================
// ui-presets.js -- Quadrature-function preset library
//
// Extracted from ui.js by the P0 split (PRIMARY-SOL handoff). Each preset is
// a quadrature function h(w) the user can load with one click from the
// preset dropdown. Selecting a preset replaces the pole data (and, for
// unbounded presets, c). Any manual edit afterward reverts the dropdown to
// "— custom —".
//
// Preset shape:
//   {
//     id:         string,    // unique
//     label:      string,    // shown in dropdown
//     poles: [
//       { a: string, order: int, residues: [string,...] },
//       ...
//     ],
//     polyCoeffs: [string,...] | undefined,   // unbounded only
//     c:          number      | undefined,    // unbounded only
//     w0:         string      | undefined,    // LQD only
//   }
//
// To ADD a preset: copy any entry, change the values. To REORDER: just
// reorder the array — the dropdown is built in array order.
//
// Each preset list is exposed BOTH as a top-level global (the original
// const name in ui.js, kept for zero-churn back-compat) AND under
// window.QD_UI.Presets.* for namespaced future readers.
// =============================================================================

(function (global) {
  'use strict';

  var QD_PRESETS_BOUNDED = [
    { id: 'unit-disk',     label: 'Unit disk:  h = 1/w',
      poles: [ { a: '0', order: 1, residues: ['1'] } ] },

    { id: 'cardioid',      label: 'Cardioid:  h = 1.5/w + 0.5/w²',
      poles: [ { a: '0', order: 2, residues: ['1.5', '0.5'] } ] },

    { id: 'two-point-sym', label: 'Two-point symmetric:  1.5/(w−1) + 1.5/(w+1)',
      poles: [
        { a:  '1', order: 1, residues: ['1.5'] },
        { a: '-1', order: 1, residues: ['1.5'] },
      ] },

    { id: 'triangle', label: 'Equilateral 3-point on unit circle',
      poles: [
        { a:  '1',              order: 1, residues: ['1'] },
        { a: '-0.5+0.8660254i', order: 1, residues: ['1'] },
        { a: '-0.5-0.8660254i', order: 1, residues: ['1'] },
      ] },
  ];

  // -----------------------------------------------------------------
  // Power-weighted QDs (Family.powerQD, α ≥ 2). Bounded PQDs satisfy
  //   ∫_Ω f(w) |w|^{2(α−1)} dA = ∮_∂Ω f h dw,
  // with Riemann-map characterization φ(z) = (R#(z))^{1/α}. Q1 v1 covers
  // α = 2 (single + 2nd-order pole). Realizability constraint for the
  // single-pole case: residue C must satisfy C > (p^α − w₀^α)²/α².
  // -----------------------------------------------------------------
  var QD_PRESETS_BOUNDED_PQD = [
    { id: 'pqd-1pt-a2',
      label: 'α=2 one-pole:  h = 3/(w − 2),  w₀ = 1',
      poles: [ { a: '2', order: 1, residues: ['3'] } ],
      w0: '1',
      alpha: 2 },

    { id: 'pqd-1pt-a2-strong',
      label: 'α=2 one-pole (closer):  h = 4/(w − 1.5),  w₀ = 1',
      poles: [ { a: '1.5', order: 1, residues: ['4'] } ],
      w0: '1',
      alpha: 2 },

    // Q1.2: 2nd-order pole exercising the higher-order (★) formula.
    { id: 'pqd-1pt-a2-m2',
      label: 'α=2 second-order:  h = 0.5/(w−1.2) + 0.05/(w−1.2)²,  w₀ = 1',
      poles: [ { a: '1.2', order: 2, residues: ['0.5', '0.05'] } ],
      w0: '1',
      alpha: 2 },

    // QA: non-integer α. The (★) closed form is α-general; all power ops use
    // Complex.cpow. α<1 is the LQD-limit regime (weight singular at 0, but
    // 0 ∉ Ω̄ so the integral is fine).
    { id: 'pqd-1pt-a1p5',
      label: 'α=1.5 one-pole:  h = 3/(w − 2),  w₀ = 1',
      poles: [ { a: '2', order: 1, residues: ['3'] } ],
      w0: '1',
      alpha: 1.5 },

    { id: 'pqd-1pt-a0p5',
      label: 'α=0.5 one-pole (LQD-limit):  h = 3/(w − 2),  w₀ = 1',
      poles: [ { a: '2', order: 1, residues: ['3'] } ],
      w0: '1',
      alpha: 0.5 },
  ];

  // Bounded SINGULAR PQDs (Family.powerQD_singular, 0 ∈ Ω). φ = b_{z₀}·(R#)^{1/α}.
  // The canonical example h=(63/32)/(w−1), α=2, w₀=1 lands at z₀ = 2/3. The
  // mass/area constraint pins |z₀|; w₀ = φ(0) is a non-origin interior point.
  var QD_PRESETS_BOUNDED_PQD_SINGULAR = [
    { id: 'pqds-1pt-a2',
      label: 'α=2 one-pole:  h = (63/32)/(w − 1),  w₀ = 1   (z₀ = 2/3)',
      poles: [ { a: '1', order: 1, residues: ['1.96875'] } ],
      w0: '1',
      alpha: 2 },

    { id: 'pqds-1pt-a2-b',
      label: 'α=2 one-pole:  h = 3/(w − 1.2),  w₀ = 1.1',
      poles: [ { a: '1.2', order: 1, residues: ['3'] } ],
      w0: '1.1',
      alpha: 2 },

    { id: 'pqds-1pt-a1p5',
      label: 'α=1.5 one-pole:  h = 2.2/(w − 1),  w₀ = 1',
      poles: [ { a: '1', order: 1, residues: ['2.2'] } ],
      w0: '1',
      alpha: 1.5 },
  ];

  // Unbounded PQDs (Family.unboundedPQD, 0 ∉ Ω). φ(z) = z·(r#)^{1/α} on 𝔻*,
  // r#(∞) = c^α. c is the conformal radius (user input). The constant-h case
  // is the thesis Example 4.3.1: φ = c·z·(1 − γ/z)^{1/α}, γ = −α·h₀/c^{2α−1}.
  var QD_PRESETS_UNBOUNDED_PQD = [
    { id: 'upqd-const-a2',
      label: 'α=2 constant:  h = 0.3,  c = 1   (Example 4.3.1)',
      poles: [],
      polyCoeffs: ['0.3'],
      c: 1,
      alpha: 2 },

    { id: 'upqd-1pt-a2',
      label: 'α=2 one-pole:  h = 1/(w − 2.5),  c = 2',
      poles: [ { a: '2.5', order: 1, residues: ['1'] } ],
      c: 2,
      alpha: 2 },

    { id: 'upqd-1pt-a1p5',
      label: 'α=1.5 one-pole:  h = 0.5/(w − 2.5),  c = 2',
      poles: [ { a: '2.5', order: 1, residues: ['0.5'] } ],
      c: 2,
      alpha: 1.5 },

    // Polynomial-h (pole at ∞). Monomial: Ω ∈ QD_a(α·k·w^{k-1}), φ = c·z·(1−γ_k/z^k)^{1/α}.
    { id: 'upqd-mono-a2',
      label: 'α=2 monomial:  h = w,  c = 2   (pole at ∞, Thm 4.5.3)',
      poles: [],
      polyCoeffs: ['0', '1'],
      c: 2,
      alpha: 2 },

    { id: 'upqd-mono2-a2',
      label: 'α=2 monomial:  h = 0.9·w²,  c = 4   (degree-2 pole at ∞)',
      poles: [],
      polyCoeffs: ['0', '0', '0.9'],
      c: 4,
      alpha: 2 },

    { id: 'upqd-poly-pole-a2',
      label: 'α=2 poly + pole:  h = 0.2·w + 0.5/(w − 3),  c = 2',
      poles: [ { a: '3', order: 1, residues: ['0.5'] } ],
      polyCoeffs: ['0', '0.2'],
      c: 2,
      alpha: 2 },
  ];

  // Unbounded SINGULAR PQDs (Family.unboundedPQD_singular, 0 ∈ Ω). φ = z·b_{z₀}·
  // (r#)^{1/α}; z₀ pinned by r(z₀)=0 (Prop 4.6.3). Monomial h=1 → z₀=(−2)^{1/3}.
  var QD_PRESETS_UNBOUNDED_PQD_SINGULAR = [
    { id: 'upqds-mono-a2',
      label: 'α=2 monomial:  h = 1,  c = 1   (z₀ = (−2)^{1/3}, Thm 4.5.2)',
      poles: [],
      polyCoeffs: ['1'],
      c: 1,
      alpha: 2 },

    { id: 'upqds-1pt-a2',
      label: 'α=2 one-pole:  h = 0.5/(w − 2),  c = 1',
      poles: [ { a: '2', order: 1, residues: ['0.5'] } ],
      c: 1,
      alpha: 2 },

    { id: 'upqds-1pt-a1p5',
      label: 'α=1.5 one-pole:  h = 1/(w − 2),  c = 1',
      poles: [ { a: '2', order: 1, residues: ['1'] } ],
      c: 1,
      alpha: 1.5 },
  ];

  var QD_PRESETS_UNBOUNDED = [
    { id: 'unb-1pt-pos',  label: 'One-point positive charge:  h = 1/(w − 2),  c = 2',
      poles: [ { a: '2', order: 1, residues: ['1'] } ],
      c: 2 },

    { id: 'unb-1pt-neg',  label: 'One-point negative charge:  h = −0.5/(w − 2),  c = 0.7',
      poles: [ { a: '2', order: 1, residues: ['-0.5'] } ],
      c: 0.7 },

    { id: 'unb-1pt-imag', label: 'One-point imaginary charge:  h = i/(w − 2),  c = 0.8',
      poles: [ { a: '2', order: 1, residues: ['i'] } ],
      c: 0.8 },

    { id: 'unb-deltoid',  label: 'Deltoid:  h = w²,  c = 0.5',
      poles: [],
      polyCoeffs: ['0', '0', '1'],
      c: 0.5 },

    { id: 'unb-2pt-nonuniq', label: 'Two-point non-uniqueness:  1/(w−1) + 1/(w+1),  c = 0.4',
      poles: [
        { a:  '1', order: 1, residues: ['1'] },
        { a: '-1', order: 1, residues: ['1'] },
      ],
      c: 0.4 },
  ];

  // LQD presets — BOUNDED non-singular (Theorem 5.3.2 single-pole family
  // at the closed-form one-point cases). Ω = {|ln(w/w₀)|² < α} for h =
  // α/(w−w₀), with double-point at α = π².
  var LQD_PRESETS_BOUNDED = [
    { id: 'lqd-1pt-small', label: 'One-pt: α = 0.5 / (w − 1),  w₀ = 1',
      poles: [ { a: '1', order: 1, residues: ['0.5'] } ],
      w0: '1' },

    { id: 'lqd-1pt-medium', label: 'One-pt: α = 2 / (w − 1),  w₀ = 1',
      poles: [ { a: '1', order: 1, residues: ['2'] } ],
      w0: '1' },

    { id: 'lqd-1pt-large', label: 'One-pt: α = 9 / (w − 1),  w₀ = 1  (near critical α = π²)',
      poles: [ { a: '1', order: 1, residues: ['9'] } ],
      w0: '1' },

    { id: 'lqd-1pt-shifted', label: 'Shifted one-pt: α = 0.4 / (w − 2),  w₀ = 2',
      poles: [ { a: '2', order: 1, residues: ['0.4'] } ],
      w0: '2' },

    { id: 'lqd-1pt-complex', label: 'Complex w₀: α = 0.5 / (w − (1+i)),  w₀ = 1+i',
      poles: [ { a: '1+i', order: 1, residues: ['0.5'] } ],
      w0: '1+i' },

    { id: 'lqd-3pt-equi', label: 'Equilateral 3-pt around w₀ = 3 (existence depends on residues)',
      poles: [
        { a:  '3.5',                     order: 1, residues: ['0.2'] },
        { a:  '2.75+0.4330127i',         order: 1, residues: ['0.2'] },
        { a:  '2.75-0.4330127i',         order: 1, residues: ['0.2'] },
      ],
      w0: '3' },
  ];

  // LQD presets — BOUNDED SINGULAR. q is dialed separately via slider.
  var LQD_PRESETS_BOUNDED_SINGULAR = [
    { id: 'lqd-s-thm-562',
      label: 'Thm 5.6.2 family: 0.5/(w − 2),  w₀ = 1  (dial q with slider)',
      poles: [ { a: '2', order: 1, residues: ['0.5'] } ],
      w0: '1' },

    { id: 'lqd-s-shifted',
      label: 'Shifted Thm 5.6.2: 0.3/(w − 1.5),  w₀ = 0.6  (dial q)',
      poles: [ { a: '1.5', order: 1, residues: ['0.3'] } ],
      w0: '0.6' },

    { id: 'lqd-s-2pt-sym',
      label: 'Two-pt symmetric: 0.4/(w−1) + 0.4/(w+1),  w₀ = 0.5+0.5i  (dial q)',
      poles: [
        { a:  '1', order: 1, residues: ['0.4'] },
        { a: '-1', order: 1, residues: ['0.4'] },
      ],
      w0: '0.5+0.5i' },
  ];

  // LQD presets — UNBOUNDED non-singular. Ω unbounded with 0 ∉ Ω̄.
  // Conformal radius c = φ'(∞) > 0. h ≡ 0 → exterior of disk of radius c.
  var LQD_PRESETS_UNBOUNDED = [
    { id: 'lqd-u-trivial', label: 'Trivial:  h = 0,  c = 0.5  (Ω = ext. disk radius c)',
      poles: [], c: 0.5 },

    { id: 'lqd-u-1pt', label: 'One-pt:  h = 1/(w − 2),  c = 0.6',
      poles: [ { a: '2', order: 1, residues: ['1'] } ],
      c: 0.6 },

    { id: 'lqd-u-1pt-small-c', label: 'One-pt, small c:  h = 1/(w − 2),  c = 0.3',
      poles: [ { a: '2', order: 1, residues: ['1'] } ],
      c: 0.3 },

    { id: 'lqd-u-2pt-sym', label: 'Two-pt symmetric:  1/(w−2) + 0.6/(w+1.5),  c = 0.4',
      poles: [
        { a:  '2',    order: 1, residues: ['1'] },
        { a: '-1.5',  order: 1, residues: ['0.6'] },
      ],
      c: 0.4 },
  ];

  // LQD presets — UNBOUNDED SINGULAR. q is dialable (default 0).
  var LQD_PRESETS_UNBOUNDED_SINGULAR = [
    { id: 'lqd-us-1pt', label: 'One-pt:  h = q/w + 1/(w−2),  c = 0.6  (dial q)',
      poles: [ { a: '2', order: 1, residues: ['1'] } ],
      c: 0.6 },

    { id: 'lqd-us-2pt-sym', label: 'Two-pt:  q/w + 1/(w−2) + 0.6/(w+1.5),  c = 0.4',
      poles: [
        { a:  '2',    order: 1, residues: ['1'] },
        { a: '-1.5',  order: 1, residues: ['0.6'] },
      ],
      c: 0.4 },
  ];

  // Expose under both the original top-level names (for zero-churn references
  // from ui.js) AND under window.QD_UI.Presets.* for any future namespaced reader.
  global.QD_PRESETS_BOUNDED            = QD_PRESETS_BOUNDED;
  global.QD_PRESETS_BOUNDED_PQD        = QD_PRESETS_BOUNDED_PQD;
  global.QD_PRESETS_BOUNDED_PQD_SINGULAR = QD_PRESETS_BOUNDED_PQD_SINGULAR;
  global.QD_PRESETS_UNBOUNDED          = QD_PRESETS_UNBOUNDED;
  global.QD_PRESETS_UNBOUNDED_PQD      = QD_PRESETS_UNBOUNDED_PQD;
  global.QD_PRESETS_UNBOUNDED_PQD_SINGULAR = QD_PRESETS_UNBOUNDED_PQD_SINGULAR;
  global.LQD_PRESETS_BOUNDED           = LQD_PRESETS_BOUNDED;
  global.LQD_PRESETS_BOUNDED_SINGULAR  = LQD_PRESETS_BOUNDED_SINGULAR;
  global.LQD_PRESETS_UNBOUNDED         = LQD_PRESETS_UNBOUNDED;
  global.LQD_PRESETS_UNBOUNDED_SINGULAR = LQD_PRESETS_UNBOUNDED_SINGULAR;

  global.QD_UI = global.QD_UI || {};
  global.QD_UI.Presets = {
    bounded:           QD_PRESETS_BOUNDED,
    boundedPQD:        QD_PRESETS_BOUNDED_PQD,
    boundedPQDSing:    QD_PRESETS_BOUNDED_PQD_SINGULAR,
    unbounded:         QD_PRESETS_UNBOUNDED,
    unboundedPQD:      QD_PRESETS_UNBOUNDED_PQD,
    unboundedPQDSing:  QD_PRESETS_UNBOUNDED_PQD_SINGULAR,
    lqdBounded:        LQD_PRESETS_BOUNDED,
    lqdBoundedSing:    LQD_PRESETS_BOUNDED_SINGULAR,
    lqdUnbounded:      LQD_PRESETS_UNBOUNDED,
    lqdUnboundedSing:  LQD_PRESETS_UNBOUNDED_SINGULAR,
  };

})(typeof window !== 'undefined' ? window
                                 : (typeof globalThis !== 'undefined' ? globalThis : this));
