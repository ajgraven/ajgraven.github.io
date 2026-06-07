// =============================================================================
// observables.js  —  Foundational analysis primitives for a solved QD map φ.
//
// These are the cheap, reusable quantities the dynamics/analysis features build
// on (Hele-Shaw growth, harmonic measure / capacity, DLA, the geometry card).
// Everything is computed from a φ that is ALREADY solved — no new solver math,
// no new families — by sweeping the unit circle and reading φ, φ′, φ″ off the
// family-agnostic Taylor evaluator. Pure + DOM-free; loaded page-side only (like
// critical-set.js / cusps.js / univalence.js), never bundled into the Workers.
//
// API
//   QD.boundaryObservables(phi, { samples = 1024 }) → {
//       samples, w[], curvature[], maxCurvature, argMaxCurvatureTheta,
//       perimeter, area, signedArea, centroid:{re,im}, moments:[{re,im}…k=0..4],
//       minAbsPhiPrime, minAbsPhiPrimeTheta }   (w[] aligned with curvature[])
//   QD.harmonicMeasure(phi, { samples = 720 }) → {
//       samples, theta[], w[], arclength[], density[], perimeter, integral,
//       maxDensity, meanDensity, peakTheta }
//   QD.estimateAccuracy(phi, hData, { samples }) → {
//       residual, relN, rel2N, samplesN, underResolved,
//       significantDigits, conditionEst }
//   QD.minAbsPhiPrimeAngle(phi, samples = 1024) → { theta, value }
//
// MATH NOTES
//  • Boundary curve γ(θ) = φ(e^{iθ}), z = e^{iθ}. With ż = iz, z̈ = −z:
//        γ′ = φ′·(iz),   γ″ = −φ″·z² − φ′·z,
//        κ  = Im(conj(γ′)·γ″) / |γ′|³        (signed curvature).
//    κ → ∞ exactly where φ′ → 0 (a cusp), so maxCurvature flags the cusp.
//  • Area / centroid via the shoelace polygon formulas on the boundary samples.
//    For BOUNDED φ the unit circle traces ∂Ω counter-clockwise (signedArea > 0);
//    for UNBOUNDED φ (exterior map) it traces ∂Ω clockwise around the complement
//    K, so signedArea < 0 — `area` is reported as |signedArea| (the area of K).
//  • Complex area moments  M_k = ∬_Ω w^k dA  via Stokes,
//        M_k = ½ Σ_n w_n^k · conj(w_n) · φ′(z_n) · z_n · Δθ,   Δθ = 2π/N,
//    so Re(M_0) ≈ signedArea (a built-in cross-check) and the M_k are the
//    Hele-Shaw / QD "harmonic moments" — a cheap byproduct of the same sweep.
//  • Harmonic measure of ∂Ω seen from φ(0) (bounded) / ∞ (unbounded) has density
//    w.r.t. arclength  ρ(θ) = 1/(2π·|φ′(e^{iθ})|)  (since |γ′| = |φ′| on |z|=1).
//    ∮ ρ ds = 1 by construction; we report Σ ρ_n·|Δw_n| as a discretization-error
//    self-check (≈ 1 for a well-resolved boundary).
//  • Accuracy: the quadrature-identity verifier's maxRelDiff at N vs 2N nodes
//    separates genuine identity error (relN ≈ rel2N) from contour-quadrature
//    under-resolution (relN ≫ rel2N). significantDigits ≈ −log10(rel2N ∨ residual).
//
// Reuses QD.phiTaylorAt (exact φ, φ′, φ″; family-agnostic), QD.evalPhi,
// QD.Complex, and — for estimateAccuracy — QD.residual / residualNorm / packPhi /
// unpackPhi / numericalJacobian / houseQR and the family identity verifier.
// =============================================================================

(function (global) {
  'use strict';

  // QD-namespace resolution — same idiom every solver/analysis file uses.
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' ? module.exports : null);
  if (!QD) return;

  const C = QD.Complex;
  const TWO_PI = 2 * Math.PI;
  const EPS_MACH = 2.220446049250313e-16;

  // ---------------------------------------------------------------------------
  // Single unit-circle sweep shared by boundaryObservables + harmonicMeasure.
  // Returns parallel arrays of length N: theta, z, w (=φ), phiP (=φ′), phiPP
  // (=φ″), plus the running min|φ′| and its angle. φ, φ′, φ″ come from one
  // QD.phiTaylorAt(z, phi, 2) call per node (a[1]=φ′, a[2]=φ″/2! ⇒ φ″=2·a[2]).
  // ---------------------------------------------------------------------------
  function _sweep(phi, N) {
    const theta = new Array(N), z = new Array(N), w = new Array(N),
          phiP = new Array(N), phiPP = new Array(N);
    let minAbs = Infinity, minTheta = 0;
    for (let n = 0; n < N; n++) {
      const th = TWO_PI * n / N;
      const zn = { re: Math.cos(th), im: Math.sin(th) };
      const t = QD.phiTaylorAt(zn, phi, 2);   // [φ, φ′, φ″/2]
      theta[n] = th;
      z[n] = zn;
      w[n] = t[0];
      phiP[n] = t[1];
      phiPP[n] = C.scale(t[2], 2);            // φ″ = 2·(coeff of t²)
      const ap = C.abs(t[1]);
      if (ap < minAbs) { minAbs = ap; minTheta = th; }
    }
    return { N, theta, z, w, phiP, phiPP, minAbs, minTheta };
  }

  // ---------------------------------------------------------------------------
  // boundaryObservables — curvature, area, perimeter, centroid, moments.
  // ---------------------------------------------------------------------------
  function boundaryObservables(phi, opts) {
    opts = opts || {};
    const N = Math.max(64, opts.samples || 1024);
    const s = _sweep(phi, N);
    const { w, z, phiP, phiPP } = s;
    const dTheta = TWO_PI / N;

    // --- signed curvature κ(θ) of γ(θ)=φ(e^{iθ}) ----------------------------
    const curvature = new Array(N);
    let maxCurvature = 0, argMaxTheta = s.theta[0];
    for (let n = 0; n < N; n++) {
      const zn = z[n];
      const iz = { re: -zn.im, im: zn.re };          // i·z
      const gP = C.mul(phiP[n], iz);                 // γ′ = φ′·iz
      const z2 = C.mul(zn, zn);                       // z²
      const gPP = C.sub(C.neg(C.mul(phiPP[n], z2)),   // −φ″·z² − φ′·z
                        C.mul(phiP[n], zn));
      const gpAbs = C.abs(gP);
      let kappa;
      if (gpAbs < 1e-300) {
        kappa = Infinity;
      } else {
        const num = C.mul(C.conj(gP), gPP).im;       // Im(conj(γ′)·γ″)
        kappa = num / (gpAbs * gpAbs * gpAbs);
      }
      curvature[n] = kappa;
      const ak = Math.abs(kappa);
      if (isFinite(ak) && ak > maxCurvature) { maxCurvature = ak; argMaxTheta = s.theta[n]; }
    }

    // --- perimeter, shoelace signed area, polygon centroid -------------------
    let perimeter = 0, signedArea = 0, cxAcc = 0, cyAcc = 0;
    for (let n = 0; n < N; n++) {
      const a = w[n], b = w[(n + 1) % N];
      perimeter += Math.hypot(b.re - a.re, b.im - a.im);
      const cross = a.re * b.im - b.re * a.im;        // x_n y_{n+1} − x_{n+1} y_n
      signedArea += cross;
      cxAcc += (a.re + b.re) * cross;
      cyAcc += (a.im + b.im) * cross;
    }
    signedArea *= 0.5;
    const area = Math.abs(signedArea);
    let centroid;
    if (Math.abs(signedArea) > 1e-300) {
      centroid = { re: cxAcc / (6 * signedArea), im: cyAcc / (6 * signedArea) };
    } else {
      let mx = 0, my = 0;                              // degenerate fallback: mean
      for (let n = 0; n < N; n++) { mx += w[n].re; my += w[n].im; }
      centroid = { re: mx / N, im: my / N };
    }

    // --- complex area moments M_k = ∬_Ω w^k dA via Stokes (k = 0..4) ---------
    const KMAX = 4;
    const moments = [];
    for (let k = 0; k <= KMAX; k++) {
      let acc = { re: 0, im: 0 };
      for (let n = 0; n < N; n++) {
        const wk = C.pow(w[n], k);                    // w^k
        let term = C.mul(wk, C.conj(w[n]));           // w^k·conj(w)
        term = C.mul(term, phiP[n]);                  // ·φ′
        term = C.mul(term, z[n]);                     // ·z
        acc = C.add(acc, term);
      }
      moments.push(C.scale(acc, 0.5 * dTheta));       // ½·Σ(…)·Δθ
    }

    return {
      samples: N,
      w,                                   // boundary points φ(e^{iθ_n}), aligned with curvature[]
      curvature, maxCurvature, argMaxCurvatureTheta: argMaxTheta,
      perimeter, area, signedArea, centroid, moments,
      minAbsPhiPrime: s.minAbs, minAbsPhiPrimeTheta: s.minTheta,
    };
  }

  // ---------------------------------------------------------------------------
  // harmonicMeasure — density on ∂Ω (∝ 1/|φ′|) w.r.t. arclength.
  // ---------------------------------------------------------------------------
  function harmonicMeasure(phi, opts) {
    opts = opts || {};
    const N = Math.max(64, opts.samples || 720);
    const s = _sweep(phi, N);
    const { theta, w, phiP } = s;

    const density = new Array(N), arclength = new Array(N);
    let perimeter = 0, integral = 0, maxDensity = 0, peakTheta = theta[0], sumDensity = 0;
    let acc = 0;
    for (let n = 0; n < N; n++) {
      const ap = C.abs(phiP[n]);
      const rho = ap < 1e-300 ? Infinity : 1 / (TWO_PI * ap);
      density[n] = rho;
      arclength[n] = acc;
      const next = w[(n + 1) % N];
      const ds = Math.hypot(next.re - w[n].re, next.im - w[n].im);
      acc += ds;
      perimeter += ds;
      if (isFinite(rho)) {
        integral += rho * ds;                          // Σ ρ·|Δw| ≈ 1
        sumDensity += rho;
        if (rho > maxDensity) { maxDensity = rho; peakTheta = theta[n]; }
      }
    }
    const meanDensity = sumDensity / N;

    return {
      samples: N, theta, w, arclength, density,
      perimeter, integral, maxDensity, meanDensity, peakTheta,
    };
  }

  // ---------------------------------------------------------------------------
  // minAbsPhiPrimeAngle — the boundary angle of the smallest |φ′| (near-cusp).
  // Standalone (shared with the adaptive-quadrature work); cheap one-pass scan.
  // ---------------------------------------------------------------------------
  function minAbsPhiPrimeAngle(phi, samples) {
    const N = Math.max(64, samples || 1024);
    let minAbs = Infinity, minTheta = 0;
    for (let n = 0; n < N; n++) {
      const th = TWO_PI * n / N;
      const t = QD.phiTaylorAt({ re: Math.cos(th), im: Math.sin(th) }, phi, 1);
      const ap = C.abs(t[1]);
      if (ap < minAbs) { minAbs = ap; minTheta = th; }
    }
    return { theta: minTheta, value: minAbs };
  }

  // ---------------------------------------------------------------------------
  // estimateAccuracy — confidence numbers for a solved φ against its h-data.
  //   relN vs rel2N separate genuine identity error from quadrature under-
  //   resolution; significantDigits ≈ −log10(rel2N ∨ residual); conditionEst is
  //   the QR condition number of the residual Jacobian (best-effort, may be null).
  // ---------------------------------------------------------------------------
  function _familyOf(phi) {
    const name = phi.family || (phi.unbounded ? 'unboundedQD' : 'boundedQD');
    return (QD.Family && QD.Family[name]) || null;
  }

  function estimateAccuracy(phi, hData, opts) {
    opts = opts || {};
    const out = {
      residual: null, relN: null, rel2N: null, samplesN: null,
      underResolved: false, significantDigits: null, conditionEst: null,
      // Near-cusp regime (#11): nearCusp flags that a φ′ zero is within
      // NEAR_CUSP_BAND of |z| = 1 (a forming cusp); cuspDistance is that radial
      // gap (smaller ⇒ closer to c*); trustedSignal names which check to believe.
      // In the near-cusp regime the quadrature-identity verifier becomes
      // unreliable (the hole thins so interior test points can't clear ∂Ω) and
      // the GEOMETRIC criterion (univalence + critical modulus) governs validity.
      nearCusp: false, cuspDistance: null, trustedSignal: 'identity',
    };
    if (!phi || !hData) return out;

    // Near-cusp detection via the critical set (page-only, like this module). The
    // closest φ′ zero to |z| = 1 sets cuspDistance; within the band ⇒ near-cusp.
    const NEAR_CUSP_BAND = 0.05;
    try {
      if (typeof QD.findCriticalPoints === 'function') {
        const cs = QD.findCriticalPoints(phi);
        let minDist = Infinity;
        for (const pt of (cs && cs.points) || []) {
          const d = Math.abs(pt.absZ - 1);
          if (d < minDist) minDist = d;
        }
        if (isFinite(minDist)) {
          out.cuspDistance = minDist;
          out.nearCusp = minDist <= NEAR_CUSP_BAND;
        }
      }
    } catch (e) { /* leave nearCusp=false */ }
    out.trustedSignal = out.nearCusp ? 'geometry' : 'identity';

    // Solve residual ‖F(φ)‖.
    try { out.residual = QD.residualNorm(QD.residual(phi, hData)); } catch (e) { /* leave null */ }

    // Multi-resolution identity error.
    const fam = _familyOf(phi);
    if (fam && typeof fam.verifyQuadratureIdentity === 'function') {
      const N = Math.max(256, opts.samples || 600);
      out.samplesN = N;
      try {
        const vN  = fam.verifyQuadratureIdentity(phi, hData, { numSamples: N });
        const v2N = fam.verifyQuadratureIdentity(phi, hData, { numSamples: 2 * N });
        out.relN  = vN  ? vN.maxRelDiff  : null;
        out.rel2N = v2N ? v2N.maxRelDiff : null;
        // Under-resolution: the coarse check is ≫ the fine one (and not already tiny).
        if (out.relN != null && out.rel2N != null &&
            isFinite(out.relN) && isFinite(out.rel2N) &&
            out.relN > 1e-9 && out.relN > 16 * out.rel2N) {
          out.underResolved = true;
        }
      } catch (e) { /* leave nulls */ }
    }

    // Significant digits from the best available error floor.
    const floor = Math.max(
      (out.rel2N != null && isFinite(out.rel2N)) ? out.rel2N : 0,
      (out.residual != null && isFinite(out.residual)) ? out.residual : 0,
      EPS_MACH);
    out.significantDigits = Math.max(0, Math.min(16, -Math.log10(floor)));

    // Best-effort conditioning of the residual Jacobian at φ.
    try {
      const v0 = QD.packPhi(phi);
      const J = QD.numericalJacobian(v0, (v) => QD.residual(QD.unpackPhi(v, phi), hData));
      const qr = QD.houseQR(J);
      if (qr && isFinite(qr.condEst)) out.conditionEst = qr.condEst;
    } catch (e) { /* conditionEst stays null */ }

    return out;
  }

  QD.Observables = { boundaryObservables, harmonicMeasure, estimateAccuracy, minAbsPhiPrimeAngle };
  QD.boundaryObservables = boundaryObservables;
  QD.harmonicMeasure = harmonicMeasure;
  QD.estimateAccuracy = estimateAccuracy;
  QD.minAbsPhiPrimeAngle = minAbsPhiPrimeAngle;

})(typeof globalThis !== 'undefined' ? globalThis : this);
