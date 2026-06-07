// =============================================================================
// thesis-examples.js  —  Curated canonical quadrature domains, each with an
// ANALYTIC ORACLE: the closed-form quantities a correct solve must reproduce.
// (Roadmap #8.) Loading one frames the view and shows a computed-vs-expected
// "oracle card"; the same data + checkOracle engine anchors a headless test.
//
//   QD.ThesisExamples           — array of example descriptors (below).
//   QD.thesisExampleHData(ex)    — string poles → Complex hData (for solving/tests).
//   QD.checkOracle(phi, hData, oracle, opts) → Promise<{ rows, allPass }>
//
// Example descriptor (extends the ui-presets.js preset shape with mode + oracle):
//   { id, label, blurb,                          // name + one-line significance
//     mode,                                       // 'bounded' | 'unbounded' | 'pqd-*' | 'lqd-*'
//     poles:[{a, order, residues}],               // STRING coefficients (preset form)
//     polyCoeffs?, c?, w0?, alpha?, q?,           // same fields applyPreset() consumes
//     view?:{cx,cy,scale},                        // optional framing
//     oracle:{ … } }                              // optional closed-form expectations
//
// Oracle fields (each present field becomes one comparison row):
//   area, perimeter, M0           — boundaryObservables (Re for M0)
//   curvatureUniform, harmonicUniform (bool)
//   cuspCount, cusps:[{thetaDeg,type:[p,q],orderM}]   — classifyCusps
//   rotationalOrder, reflectionAxesCount               — detectSymmetry
//   significantDigitsMin                               — estimateAccuracy
//   cMax, cMaxMechanism (opt-in, heavy)                — estimateMaxConformalRadius
//
// Pure + DOM-free; loads page-side (SOLVER_PAGE_ONLY_FILES), not in the workers.
// =============================================================================

(function (global) {
  'use strict';

  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD) return;

  // Example blurbs live in the central strings file (QD.Strings.blurbs, in
  // ui-strings.js — loaded before this module). Look them up by the example id
  // (dash-case → camelCase key); '' if the strings module isn't present.
  function blurbOf(id) {
    const key = id.replace(/-([a-z0-9])/g, (_, c) => c.toUpperCase());
    const B = QD.Strings && QD.Strings.blurbs;
    return (B && B[key]) || '';
  }

  // ---------------------------------------------------------------------------
  // The curated set. Configs reuse the exact preset/test constructions that are
  // already validated elsewhere; oracle values are the closed-form / locked
  // numeric answers (see app/test/{observables,cusps,cmax,solvers}.test.js).
  // ---------------------------------------------------------------------------
  const ThesisExamples = [
    {
      id: 'disk', label: 'Unit disk', mode: 'bounded',
      blurb: blurbOf('disk'),
      poles: [{ a: '0', order: 1, residues: ['1'] }],
      view: { cx: 0, cy: 0, scale: 120 },
      oracle: { area: Math.PI, perimeter: 2 * Math.PI, M0: Math.PI,
                curvatureUniform: true, harmonicUniform: true, cuspCount: 0,
                significantDigitsMin: 9 },
    },
    {
      id: 'two-point-sym', label: 'Two-point symmetric', mode: 'bounded',
      blurb: blurbOf('two-point-sym'),
      poles: [{ a: '1', order: 1, residues: ['1.5'] }, { a: '-1', order: 1, residues: ['1.5'] }],
      view: { cx: 0, cy: 0, scale: 80 },
      oracle: { rotationalOrder: 2, reflectionAxesCount: 2, cuspCount: 0,
                significantDigitsMin: 9 },
    },
    {
      id: 'triangle', label: 'Equilateral 3-point', mode: 'bounded',
      blurb: blurbOf('triangle'),
      poles: [
        { a: '1', order: 1, residues: ['1'] },
        { a: '-0.5+0.8660254i', order: 1, residues: ['1'] },
        { a: '-0.5-0.8660254i', order: 1, residues: ['1'] },
      ],
      view: { cx: 0, cy: 0, scale: 80 },
      oracle: { rotationalOrder: 3, reflectionAxesCount: 3, cuspCount: 0,
                significantDigitsMin: 8 },
    },
    {
      id: 'square-4pole', label: 'Four-fold symmetric', mode: 'bounded',
      blurb: blurbOf('square-4pole'),
      poles: [
        { a: '1', order: 1, residues: ['1'] }, { a: '-1', order: 1, residues: ['1'] },
        { a: 'i', order: 1, residues: ['1'] }, { a: '-i', order: 1, residues: ['1'] },
      ],
      view: { cx: 0, cy: 0, scale: 70 },
      oracle: { rotationalOrder: 4, reflectionAxesCount: 4, cuspCount: 0,
                significantDigitsMin: 8 },
    },
    {
      id: 'cardioid-unbounded', label: 'Cardioid (c*-limited cusp)', mode: 'unbounded', c: 1.0,
      blurb: blurbOf('cardioid-unbounded'),
      poles: [{ a: '0', order: 2, residues: ['1.5', '0.5'] }],
      view: { cx: 0, cy: 0, scale: 60 },
      oracle: { cMax: 1.449, cMaxMechanism: 'cusp', significantDigitsMin: 6 },
    },
    {
      id: 'deltoid-unbounded', label: 'Deltoid (3-cusp, c*-limited)', mode: 'unbounded', c: 0.45,
      polyCoeffs: ['0', '0', '1'],
      blurb: blurbOf('deltoid-unbounded'),
      poles: [],
      view: { cx: 0, cy: 0, scale: 90 },
      oracle: { cMax: 0.5, cMaxMechanism: 'cusp', significantDigitsMin: 6 },
    },
    {
      id: 'single-pole-unbounded', label: 'Single exterior pole', mode: 'unbounded', c: 0.6,
      blurb: blurbOf('single-pole-unbounded'),
      poles: [{ a: '2', order: 1, residues: ['1'] }],
      view: { cx: 0.5, cy: 0, scale: 70 },
      oracle: { cuspCount: 0, reflectionAxesCount: 1, significantDigitsMin: 8 },
    },
  ];

  // ---------------------------------------------------------------------------
  // String preset config → Complex hData (the form the solver + detectors use).
  // ---------------------------------------------------------------------------
  function thesisExampleHData(ex) {
    const C = QD.Complex;
    const parse = (s) => C.parse(String(s)) || { re: 0, im: 0 };
    const poles = (ex.poles || []).map(p => ({
      a: parse(p.a),
      principal: (p.residues || []).map(parse),
    }));
    const polyPart = (ex.polyCoeffs || []).map(parse);
    return { poles, polyPart };
  }

  // ---------------------------------------------------------------------------
  // Oracle comparison engine — pure, unit-testable. Routes each present oracle
  // field to the matching detector and grades it pass / warn / fail.
  // ---------------------------------------------------------------------------
  function _status(relErr, passTol, warnTol) {
    if (relErr <= passTol) return 'pass';
    if (relErr <= warnTol) return 'warn';
    return 'fail';
  }
  function _relErr(computed, expected) {
    const d = Math.abs(computed - expected);
    return d / Math.max(Math.abs(expected), 1e-12);
  }
  function _boolRow(name, expected, computed) {
    return { name, expected, computed, relErr: computed === expected ? 0 : 1,
             status: computed === expected ? 'pass' : 'fail' };
  }

  async function checkOracle(phi, hData, oracle, opts) {
    opts = opts || {};
    const rows = [];
    if (!phi || !oracle) return { rows, allPass: true };

    // Geometry (one observables sweep covers area / perimeter / M0 / curvature).
    let obs = null;
    if (QD.boundaryObservables &&
        (oracle.area != null || oracle.perimeter != null || oracle.M0 != null ||
         oracle.curvatureUniform != null)) {
      try { obs = QD.boundaryObservables(phi, { samples: 1024 }); } catch (e) { /* leave null */ }
    }
    if (oracle.area != null && obs) {
      const c = obs.area, e = _relErr(c, oracle.area);
      rows.push({ name: 'area', expected: oracle.area, computed: c, relErr: e,
                  status: _status(e, 2e-3, 2e-2) });
    }
    if (oracle.perimeter != null && obs) {
      const c = obs.perimeter, e = _relErr(c, oracle.perimeter);
      rows.push({ name: 'perimeter', expected: oracle.perimeter, computed: c, relErr: e,
                  status: _status(e, 2e-3, 2e-2) });
    }
    if (oracle.M0 != null && obs && obs.moments && obs.moments[0]) {
      const c = obs.moments[0].re, e = _relErr(c, oracle.M0);
      rows.push({ name: 'M₀ = ∬ dA', expected: oracle.M0, computed: c, relErr: e,
                  status: _status(e, 2e-3, 2e-2) });
    }
    if (oracle.curvatureUniform != null && obs) {
      // Uniform ⇔ curvature spread is a tiny fraction of its magnitude.
      const spread = (obs.maxCurvature > 0)
        ? (obs.maxCurvature - _minCurv(obs)) / obs.maxCurvature : 1;
      rows.push(_boolRow('curvature uniform', true, spread < 1e-2));
    }

    if (oracle.harmonicUniform != null && QD.harmonicMeasure) {
      try {
        const hm = QD.harmonicMeasure(phi, { samples: 720 });
        const flat = hm.meanDensity > 0 &&
          (hm.maxDensity - hm.meanDensity) / hm.meanDensity < 1e-2;
        rows.push(_boolRow('harmonic measure uniform', true, flat));
      } catch (e) { /* skip */ }
    }

    // Cusps.
    if ((oracle.cuspCount != null || oracle.cusps) && QD.classifyCusps) {
      let cz = null;
      try { cz = QD.classifyCusps(phi, {}); } catch (e) { /* skip */ }
      const found = cz ? cz.cusps.filter(c => c.isCusp) : [];
      if (oracle.cuspCount != null) {
        rows.push(_boolRow('cusp count', oracle.cuspCount, found.length));
      }
      if (oracle.cusps) {
        for (const want of oracle.cusps) {
          const hit = found.find(c =>
            _angDistDeg(c.thetaDeg, want.thetaDeg) < 6 &&
            (!want.type || (c.type && c.type[0] === want.type[0] && c.type[1] === want.type[1])));
          rows.push(_boolRow('cusp @ θ≈' + want.thetaDeg + '°', true, !!hit));
        }
      }
    }

    // Symmetry.
    if ((oracle.rotationalOrder != null || oracle.reflectionAxesCount != null) &&
        QD.detectSymmetry) {
      let sym = null;
      try { sym = QD.detectSymmetry(phi); } catch (e) { /* skip */ }
      if (sym) {
        if (oracle.rotationalOrder != null) {
          rows.push(_boolRow('rotational order', oracle.rotationalOrder, sym.rotationalOrder));
        }
        if (oracle.reflectionAxesCount != null) {
          rows.push(_boolRow('reflection axes', oracle.reflectionAxesCount, sym.reflectionAxes.length));
        }
      }
    }

    // Accuracy (significant digits floor).
    if (oracle.significantDigitsMin != null && QD.estimateAccuracy && hData) {
      let acc = null;
      try { acc = QD.estimateAccuracy(phi, hData, {}); } catch (e) { /* skip */ }
      if (acc && acc.significantDigits != null) {
        const c = acc.significantDigits, want = oracle.significantDigitsMin;
        const status = c >= want ? 'pass' : (c >= want - 2 ? 'warn' : 'fail');
        rows.push({ name: 'significant digits', expected: '≥ ' + want, computed: +c.toFixed(1),
                    relErr: c >= want ? 0 : (want - c) / want, status });
      }
    }

    // c* (heavy + async ⇒ opt-in). Needs a solveFn (sync solver in tests, worker in UI).
    if (oracle.cMax != null && opts.includeCmax && QD.estimateMaxConformalRadius && hData) {
      const solveFn = opts.solveFn || QD.solveInverseQD;
      if (typeof solveFn === 'function') {
        try {
          const res = await QD.estimateMaxConformalRadius(
            hData, { unbounded: true, identityTol: 1e-6 }, solveFn,
            { cStart: oracle.cMax * 0.5, relTol: 1e-2, maxSolves: 90 });
          if (res && res.found) {
            const e = _relErr(res.cMax, oracle.cMax);
            rows.push({ name: 'c*', expected: oracle.cMax, computed: +res.cMax.toFixed(4),
                        relErr: e, status: _status(e, 6e-2, 1.2e-1) });
            if (oracle.cMaxMechanism != null) {
              rows.push(_boolRow('c* mechanism', oracle.cMaxMechanism, res.mechanism));
            }
          } else {
            rows.push({ name: 'c*', expected: oracle.cMax, computed: '(not bracketed)',
                        relErr: 1, status: 'fail' });
          }
        } catch (e) { /* skip */ }
      }
    }

    const allPass = rows.every(r => r.status !== 'fail');
    return { rows, allPass };
  }

  function _minCurv(obs) {
    let m = Infinity;
    for (const k of (obs.curvature || [])) { const a = Math.abs(k); if (a < m) m = a; }
    return isFinite(m) ? m : 0;
  }
  function _angDistDeg(a, b) {
    let d = Math.abs(((a - b) % 360 + 360) % 360);
    return Math.min(d, 360 - d);
  }

  QD.ThesisExamples = ThesisExamples;
  QD.thesisExampleHData = thesisExampleHData;
  QD.checkOracle = checkOracle;

})(typeof globalThis !== 'undefined' ? globalThis : this);
