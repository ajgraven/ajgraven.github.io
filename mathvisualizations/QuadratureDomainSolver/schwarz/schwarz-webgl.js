// =============================================================================
// schwarz-webgl.js -- GPU renderer for Schwarz-reflection dynamics.
//
// WebGL 2 fragment-shader implementation of the per-pixel escape-time loop
// used by schwarz-ui.js. The CPU implementation in schwarz-common.js is kept
// as the canonical reference and as a fallback when WebGL 2 isn't available
// (very old browsers, context-lost, --disable-webgl, etc.).
//
// Architecture
// ------------
//   • Full-screen triangle. Fragment shader computes one pixel.
//   • Phi data uploaded as uniforms with fixed maximum sizes:
//        MAX_BRANCHES = 12,  MAX_K = 8,  MAX_LAURENT = 12
//     Covers every shippable preset with ~3× headroom.
//   • All six inverse families dispatched on `u_family` (0..5) inside
//     evalPhi / evalPhiDeriv / evalSchwarz. Shared inverse-Faber branch
//     sums are computed once via branchPhi / branchSchwarz.
//   • In-Ω test via an off-screen MASK TEXTURE: the inverse-tab's sampled
//     ∂Ω polygon is rendered to a 2048² off-screen 2D canvas (filled, no
//     anti-aliasing), uploaded with NEAREST sampling and a transparent
//     border. Fragment shader does texture(mask, uv).r > 0.5; the meaning
//     inverts for unbounded Ω (mask = K-interior, "in Ω" means OUTSIDE
//     the mask).
//   • Colormap as a 256×1 RGBA texture; scale mode (smooth/discrete/log/
//     sqrt/modulo) applied in the shader as a t-transformation before
//     the texture lookup.
//
// Precision
//   • Float32 throughout. Sufficient for iteration depths ≤ ~200 at moderate
//     zoom levels. For deeper zooms with strong magnification, banding can
//     appear; emulated-double-precision is a future option.
//
// Robustness
//   • Newton tolerances tuned for float32 (CONVERGE_SQ=1e-14, FINAL_SQ=1e-10).
//   • 4-seed retry ladder in sigma() so a single bad starting basin doesn't
//     ruin a pixel; reproduces the spatial-coherence advantage the CPU path
//     gets for free via raster-order warm starts.
//
// API
//   createGPURenderer(canvas) → null on failure, else
//     {
//       available: true,
//       setPhi(phi, {boundaryPts, escapeR})    // rebuilds mask texture, packs uniforms
//       setColormap(name)                      // rebuilds colormap texture
//       render(view, {maxIter, scaleMode, modK}) // immediate, synchronous
//       capacityError() → string | null         // last setPhi failure reason
//       destroy()                              // frees all GL resources
//     }
// =============================================================================

(function (global) {
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' && module.exports ? module.exports : (global.QD || (global.QD = {})));

  const Schwarz = QD.Schwarz || (QD.Schwarz = {});

  // Compile-time caps (also baked into the shader source via constants).
  const MAX_BRANCHES = 12;
  const MAX_K        = 8;
  const MAX_LAURENT  = 12;
  // polynomial-h β-correction for unbounded LQDs (HANDOFF #22 / #26). Typical
  // polyPart degree is ≤ 3 in practice; 16 gives ~5× headroom. If a user
  // supplies a higher-degree polyPart, the CPU path takes over via the
  // capacity-error fallback in setPhi.
  const MAX_BETA     = 16;
  // The mask is a binary classifier (inside/outside polygon); we want sharp
  // edges, not smooth interpolation. 2048² gives sub-pixel boundary fidelity
  // at typical viewport sizes and combined with NEAREST sampling (set below)
  // eliminates the "boundary speckle" artifact that linear filtering caused.
  const MASK_SIZE    = 2048;

  // ===========================================================================
  // Shader sources.
  // ===========================================================================
  const VERT_SRC = `#version 300 es
in vec2 a_pos;
out vec2 v_uv;
void main() {
  v_uv = (a_pos + 1.0) * 0.5;
  gl_Position = vec4(a_pos, 0.0, 1.0);
}
`;

  // Fragment shader. Variable-bound loops use the WebGL 2 idiom
  //    for (int i = 0; i < CAP; ++i) { if (i >= u_count) break; ... }
  // so the compiler can pre-allocate registers.
  const FRAG_SRC = `#version 300 es
precision highp float;
precision highp int;

const int MAX_BRANCHES = ${MAX_BRANCHES};
const int MAX_K        = ${MAX_K};
const int MAX_LAURENT  = ${MAX_LAURENT};
const int MAX_BETA     = ${MAX_BETA};
const float EPS_DIV    = 1e-30;

uniform vec2  u_viewCenter;
uniform float u_pxPerUnit;
uniform vec2  u_canvasSize;
uniform int   u_unbounded;             // 0 = bounded Ω, 1 = unbounded Ω
uniform vec2  u_w0;                    // bounded only; unbounded ignores
uniform float u_c;                     // unbounded only
uniform vec2  u_polyA[MAX_LAURENT];
uniform int   u_polyALen;
// Polynomial-h β-correction for unbounded LQDs (HANDOFF #22 / #26):
// φ(z) gains an extra exp(B(1/z)) factor with B(1/z) = Σ_l β_l/z^l.
uniform vec2  u_lqdBeta[MAX_BETA];
uniform int   u_lqdBetaLen;
uniform vec2  u_branchZ[MAX_BRANCHES];
uniform vec2  u_branchA[MAX_BRANCHES * MAX_K];
uniform int   u_branchACount[MAX_BRANCHES];
uniform int   u_nBranches;
uniform int   u_maxIter;
uniform float u_escapeR;
uniform int   u_scaleMode;             // 0=smooth, 1=discrete, 2=log, 3=sqrt, 4=modulo
uniform int   u_modK;                  // for modulo

// Family enum:
//   0 = boundedQD                 (φ = w₀ + Σ branches; F = conj(w₀) + R##)
//   1 = unboundedQD               (φ = c·z + Σ polyA/z^l + Σ branches; F = c/z + Σ conj(polyA)·z^l + R##)
//   2 = boundedLQD                (φ = w₀ · exp(r#); F = conj(w₀) · exp(R##))
//   3 = boundedLQD_singular       (φ = γ · b · exp(r#); F = conj(γ) · b# · exp(R##))
//   4 = unboundedLQD              (φ = c·z · exp(r#−r#∞); F = (c/z) · exp(R## − conj r#∞))
//   5 = unboundedLQD_singular     (φ = c·|z₀|·z·b·exp(r#−r#∞); F = (c·|z₀|/z)·b#·exp(R##−conj r#∞))
uniform int   u_family;
uniform vec2  u_gamma;                 // singular bounded LQD
uniform vec2  u_z0;                    // singular LQD (z₀ ∈ 𝔻 bounded, ∈ 𝔻* unbounded)
uniform float u_absZ0;                 // |z₀|
uniform vec2  u_rInfConj;              // conj(r#(∞)) for unbounded LQDs

uniform sampler2D u_mask;
uniform vec2      u_maskCenter;
uniform vec2      u_maskHalfExtent;    // world half-extent (x,y) of the mask
uniform sampler2D u_colormap;

in  vec2 v_uv;
out vec4 outColor;

// -- Complex arithmetic (vec2 representation, x=re, y=im) ---------------------
vec2 cmul(vec2 a, vec2 b)  { return vec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); }
vec2 cinv(vec2 a)          { float d = dot(a, a); return vec2(a.x, -a.y) / max(d, EPS_DIV); }
vec2 cdiv(vec2 a, vec2 b)  { return cmul(a, cinv(b)); }
vec2 cconj(vec2 a)         { return vec2(a.x, -a.y); }
float cabs(vec2 a)         { return length(a); }
vec2 cexp(vec2 a)          { return exp(a.x) * vec2(cos(a.y), sin(a.y)); }

// Blaschke factor b_{z₀}(z) and its Schwarz reflection b#_{z₀}(z) = 1/b_{z₀}.
// On |z|=1, conj(b) = 1/b, so b# extends as a rational function with a pole
// at z=z₀ and a zero at z=1/conj(z₀).
vec2 blaschke(vec2 z, vec2 z0) {
  vec2 phase = -cconj(z0) / max(u_absZ0, 1e-30);
  vec2 num   = z - z0;
  vec2 den   = vec2(1.0, 0.0) - cmul(cconj(z0), z);
  return cmul(phase, cdiv(num, den));
}
vec2 blaschkeSchwarz(vec2 z, vec2 z0) {
  // b#_{z₀}(z) = -(z₀/|z₀|) · (1 − conj(z₀)z) / (z − z₀)
  vec2 phase = -z0 / max(u_absZ0, 1e-30);
  vec2 num   = vec2(1.0, 0.0) - cmul(cconj(z0), z);
  vec2 den   = z - z0;
  return cmul(phase, cdiv(num, den));
}
vec2 blaschkeLogDeriv(vec2 z, vec2 z0) {
  // d/dz log b = 1/(z−z₀) + conj(z₀)/(1 − conj(z₀)z).
  vec2 z0c = cconj(z0);
  vec2 t1  = cinv(z - z0);
  vec2 t2  = cdiv(z0c, vec2(1.0, 0.0) - cmul(z0c, z));
  return t1 + t2;
}

// -- In-Ω test via mask texture. Returns true if pixel is in Ω. --------------
bool inOmega(vec2 w) {
  vec2 uv = (w - u_maskCenter) / (2.0 * u_maskHalfExtent) + 0.5;
  // CLAMP_TO_EDGE: out-of-range uv reads the edge texel. We render the mask
  // with a transparent (0) border so the polygon never reaches the edge, so
  // out-of-bounds always reads 0.
  if (uv.x < 0.0 || uv.x > 1.0 || uv.y < 0.0 || uv.y > 1.0) {
    // Definitely outside the polygon's bounding region.
    // For bounded: outside Ω → return false.
    // For unbounded: outside K → return true (still in Ω).
    return u_unbounded == 1;
  }
  float v = texture(u_mask, uv).r;
  // mask = 1 inside the polygon, 0 outside.
  // For bounded: polygon = Ω, so inΩ ⇔ v > 0.5.
  // For unbounded: polygon = K, so inΩ ⇔ v < 0.5.
  if (u_unbounded == 1) return v < 0.5;
  return v > 0.5;
}

// -- φ(z) and φ'(z) using inverse-Faber form --------------------------------
// Branch contribution: Σ_j Σ_k conj(A_{j,k}) · u_j(z)^k, u_j = z/(1 − conj(z_j)·z).
void branchPhi(vec2 z, out vec2 sum, out vec2 sumD) {
  sum  = vec2(0.0);
  sumD = vec2(0.0);
  for (int j = 0; j < MAX_BRANCHES; ++j) {
    if (j >= u_nBranches) break;
    vec2 zjC = cconj(u_branchZ[j]);
    vec2 denom = vec2(1.0, 0.0) - cmul(zjC, z);
    vec2 u = cdiv(z, denom);
    vec2 denom2 = cmul(denom, denom);
    // For both sum and sumD: roll k = 1..count loops together.
    int count = u_branchACount[j];
    vec2 uPow  = vec2(1.0, 0.0);                       // u^0
    vec2 uPowKm1 = vec2(1.0, 0.0);                     // u^{k-1}, start at k=1
    vec2 inner = vec2(0.0);                            // Σ_k k · conj(A_{j,k}) · u^{k-1}
    for (int k = 0; k < MAX_K; ++k) {
      if (k >= count) break;
      uPow = cmul(uPow, u);                            // → u^{k+1}
      vec2 Aconj = cconj(u_branchA[j * MAX_K + k]);
      sum = sum + cmul(Aconj, uPow);
      inner = inner + cmul(Aconj, uPowKm1) * float(k + 1);
      uPowKm1 = cmul(uPowKm1, u);
    }
    sumD = sumD + cdiv(inner, denom2);
  }
}

// Schwarz-extension branch contribution: Σ_j Σ_k A_{j,k} / (z − z_j)^k.
vec2 branchSchwarz(vec2 z) {
  vec2 sum = vec2(0.0);
  for (int j = 0; j < MAX_BRANCHES; ++j) {
    if (j >= u_nBranches) break;
    vec2 d = z - u_branchZ[j];
    if (dot(d, d) < EPS_DIV) continue;
    vec2 dInv = cinv(d);
    vec2 dInvPow = vec2(1.0, 0.0);
    int count = u_branchACount[j];
    for (int k = 0; k < MAX_K; ++k) {
      if (k >= count) break;
      dInvPow = cmul(dInvPow, dInv);
      sum = sum + cmul(u_branchA[j * MAX_K + k], dInvPow);
    }
  }
  return sum;
}

// Polynomial-h β-correction helpers (UQDL / UQDLS, HANDOFF #22 / #26).
// Mirror evalBOverZ / evalBOverZDeriv / evalBConjOfZ in schwarz-common.js.
vec2 evalBOverZ(vec2 z) {
  if (u_lqdBetaLen == 0) return vec2(0.0);
  vec2 zInv = cinv(z);
  vec2 pow  = zInv;
  vec2 acc  = vec2(0.0);
  for (int l = 0; l < MAX_BETA; ++l) {
    if (l >= u_lqdBetaLen) break;
    acc = acc + cmul(u_lqdBeta[l], pow);
    pow = cmul(pow, zInv);
  }
  return acc;
}
// d/dz B(1/z) = -Σ_l l · β_l / z^{l+1}.
vec2 evalBOverZDeriv(vec2 z) {
  if (u_lqdBetaLen == 0) return vec2(0.0);
  vec2 zInv = cinv(z);
  vec2 pow  = cmul(zInv, zInv);   // 1/z^2 (l=1 term)
  vec2 acc  = vec2(0.0);
  for (int l = 1; l <= MAX_BETA; ++l) {
    if (l > u_lqdBetaLen) break;
    acc = acc - cmul(u_lqdBeta[l - 1], pow) * float(l);
    pow = cmul(pow, zInv);
  }
  return acc;
}
// Schwarz reflection of B(1/z): conj(B(1/conj(z))) = Σ_l conj(β_l) · z^l.
// Polynomial in z (no constant term), evaluated Horner-style from highest
// degree down.
vec2 evalBConjOfZ(vec2 z) {
  if (u_lqdBetaLen == 0) return vec2(0.0);
  vec2 acc = vec2(0.0);
  for (int l = MAX_BETA; l >= 1; --l) {
    if (l > u_lqdBetaLen) continue;
    acc = cmul(acc, z) + cconj(u_lqdBeta[l - 1]);
  }
  return cmul(acc, z);
}

// φ(z) — dispatches on u_family. Each branch composes the existing
// branchPhi(z) contribution with the family-specific wrapper.
vec2 evalPhi(vec2 z) {
  vec2 sum, sumD;
  branchPhi(z, sum, sumD);
  if (u_family == 0) {
    // boundedQD: w₀ + Σ branches.
    return u_w0 + sum;
  }
  if (u_family == 1) {
    // unboundedQD: c·z + Σ polyA/z^l + Σ branches.
    vec2 acc = u_c * z;
    if (u_polyALen > 0) {
      vec2 zInv = cinv(z);
      vec2 zInvPow = vec2(1.0, 0.0);
      for (int l = 0; l < MAX_LAURENT; ++l) {
        if (l >= u_polyALen) break;
        acc = acc + cmul(u_polyA[l], zInvPow);
        zInvPow = cmul(zInvPow, zInv);
      }
    }
    return acc + sum;
  }
  if (u_family == 2) {
    // boundedLQD: w₀ · exp(r#).
    return cmul(u_w0, cexp(sum));
  }
  if (u_family == 3) {
    // boundedLQD_singular: γ · b · exp(r#).
    return cmul(cmul(u_gamma, blaschke(z, u_z0)), cexp(sum));
  }
  if (u_family == 4) {
    // unboundedLQD: c·z · exp(r̃#(z) + B(1/z)).
    // r̃#(z) = r#(z) - r#(∞); B(1/z) = Σ_l β_l/z^l (HANDOFF #22).
    return cmul(u_c * z, cexp(sum - cconj(u_rInfConj) + evalBOverZ(z)));
  }
  // u_family == 5: unboundedLQD_singular: c·|z₀|·z·b·exp(r̃#(z) + B(1/z)).
  // The synthetic γ-branch (HANDOFF #24 / #26) is merged into u_branchZ /
  // u_branchA on the CPU side, so it is already absorbed into 'sum'.
  vec2 scale = (u_c * u_absZ0) * z;
  return cmul(cmul(scale, blaschke(z, u_z0)),
              cexp(sum - cconj(u_rInfConj) + evalBOverZ(z)));
}

// φ'(z). For QD/UQD: linear combination. For LQDs: log-derivative trick
// φ' = φ · (Σ piece log-derivs).
vec2 evalPhiDeriv(vec2 z) {
  vec2 sum, sumD;
  branchPhi(z, sum, sumD);
  if (u_family == 0) {
    // boundedQD: just sumD.
    return sumD;
  }
  if (u_family == 1) {
    // unboundedQD: c − Σ l·polyA/z^{l+1} + sumD.
    vec2 acc = vec2(u_c, 0.0);
    if (u_polyALen > 1) {
      vec2 zInv = cinv(z);
      vec2 zInvPow = cmul(zInv, zInv);
      for (int l = 1; l < MAX_LAURENT; ++l) {
        if (l >= u_polyALen) break;
        acc = acc - cmul(u_polyA[l], zInvPow) * float(l);
        zInvPow = cmul(zInvPow, zInv);
      }
    }
    return acc + sumD;
  }
  // LQDs: φ' = φ · (log-deriv sum).
  vec2 phiV = evalPhi(z);
  vec2 logDeriv;
  if (u_family == 2) {
    // bounded LQD: r#'.
    logDeriv = sumD;
  } else if (u_family == 3) {
    // bounded singular LQD: b'/b + r#'.
    logDeriv = blaschkeLogDeriv(z, u_z0) + sumD;
  } else if (u_family == 4) {
    // unbounded LQD: 1/z + r#'(z) + B(1/z)'.
    logDeriv = cinv(z) + sumD + evalBOverZDeriv(z);
  } else {
    // unbounded singular LQD: 1/z + b'/b + r#'(z) + B(1/z)'.
    // γ-branch contributions are inside sumD via the merged branches.
    logDeriv = cinv(z) + blaschkeLogDeriv(z, u_z0) + sumD + evalBOverZDeriv(z);
  }
  return cmul(phiV, logDeriv);
}

// Schwarz extension F(z), per family.
vec2 evalSchwarz(vec2 z) {
  vec2 branchPart = branchSchwarz(z);
  if (u_family == 0) {
    return cconj(u_w0) + branchPart;
  }
  if (u_family == 1) {
    vec2 zInv = cinv(z);
    vec2 acc = u_c * zInv;
    vec2 zPow = vec2(1.0, 0.0);
    for (int l = 0; l < MAX_LAURENT; ++l) {
      if (l >= u_polyALen) break;
      acc = acc + cmul(cconj(u_polyA[l]), zPow);
      zPow = cmul(zPow, z);
    }
    return acc + branchPart;
  }
  if (u_family == 2) {
    return cmul(cconj(u_w0), cexp(branchPart));
  }
  if (u_family == 3) {
    return cmul(cmul(cconj(u_gamma), blaschkeSchwarz(z, u_z0)), cexp(branchPart));
  }
  if (u_family == 4) {
    // F(z) = (c/z) · exp(R##(z) − conj(r#(∞)) + conj(B(z))).
    return cmul((u_c) * cinv(z),
                cexp(branchPart - u_rInfConj + evalBConjOfZ(z)));
  }
  // u_family == 5: unbounded singular LQD. γ-branch contributions are inside
  // branchPart via the merged branches; β-correction is in evalBConjOfZ.
  vec2 cOverZ = (u_c * u_absZ0) * cinv(z);
  return cmul(cmul(cOverZ, blaschkeSchwarz(z, u_z0)),
              cexp(branchPart - u_rInfConj + evalBConjOfZ(z)));
}

// Newton in z for φ(z) = w. Up to NEWTON_MAX steps, break on convergence.
//
// Tolerances are sized to float32 reality:
//   ε_float32 ≈ 1.2e-7, so the tightest useful |fz| is ~1e-7.
//   - CONVERGE_SQ = 1e-14  → succeed when |fz| < 1e-7  (machine-best)
//   - FINAL_SQ    = 1e-10  → strict-validation when the loop runs out;
//                            anything looser would let junk pixels through
//                            and is the root cause of the speckle pattern.
const int   NEWTON_MAX  = 40;
const float CONVERGE_SQ = 1e-14;
const float FINAL_SQ    = 1e-10;
const float DIVERGE_SQ  = 1e8;     // |z| > 1e4 → divergence

vec2 invertPhi(vec2 w, vec2 zSeed, out bool ok) {
  vec2 z = zSeed;
  for (int it = 0; it < NEWTON_MAX; ++it) {
    vec2 fz = evalPhi(z) - w;
    if (dot(fz, fz) < CONVERGE_SQ) { ok = true; return z; }
    vec2 dfz = evalPhiDeriv(z);
    if (dot(dfz, dfz) < EPS_DIV) { ok = false; return z; }
    z = z - cdiv(fz, dfz);
    if (any(isnan(z)) || any(isinf(z)) || dot(z, z) > DIVERGE_SQ) { ok = false; return z; }
  }
  // Strict-validate Newton's output. Without this the GPU happily accepts
  // a z that didn't actually converge and produces a speckled pixel.
  vec2 fz = evalPhi(z) - w;
  ok = (dot(fz, fz) < FINAL_SQ);
  return z;
}

// Seed for ψ-Newton. Bounded: z near 0 (cheap linearization). Unbounded:
// z ≈ w/c (good when |w| is large), or radial pull-out when w is near K.
vec2 newtonSeedFresh(vec2 w) {
  if (u_unbounded == 1) {
    vec2 cand = w / u_c;
    float r = length(cand);
    if (r > 1.05) return cand;
    if (r < 1e-12) return vec2(1.1, 0.0);
    return cand * (1.1 / r);
  }
  // bounded: linearize at z=0. φ(0) = w_0; φ'(0) = Σ conj(A_{j,1}).
  vec2 dphi0 = vec2(0.0);
  for (int j = 0; j < MAX_BRANCHES; ++j) {
    if (j >= u_nBranches) break;
    if (u_branchACount[j] >= 1) {
      dphi0 = dphi0 + cconj(u_branchA[j * MAX_K]);
    }
  }
  if (dot(dphi0, dphi0) < EPS_DIV) return vec2(0.0);
  vec2 cand = cdiv(w - u_w0, dphi0);
  float r = length(cand);
  if (r < 0.95) return cand;
  return cand * (0.9 / r);
}

// Acceptance test for a ψ-Newton result. Bounded: z must lie in 𝔻. Unbounded:
// z must lie in 𝔻*. Threshold loosened from 1e-7 (float32 ε) to 1e-4 so
// near-boundary points don't bounce between accepted and rejected on noise.
bool acceptZ(vec2 z) {
  float r = length(z);
  return (u_unbounded == 1) ? (r > 1.0 + 1e-4) : (r < 1.0 - 1e-4);
}

// σ(w) = conj(F(ψ(w))). seedHint speeds Newton (previous iterate's z).
// On GPU each pixel runs Newton independently — no cross-pixel warm start.
// To compensate, when the warm seed fails or converges to the wrong side
// of |z|=1, we retry with a small ladder of alternative seeds.
vec2 sigma(vec2 w, inout vec2 zSeed, out bool ok) {
  vec2 z = invertPhi(w, zSeed, ok);
  if (ok && acceptZ(z)) {
    // Fast path: warm seed worked.
  } else {
    // Retry ladder. Each seed perturbs the fresh-linearization estimate by
    // a different rotation/scale; one of them lands in the correct basin.
    vec2 fresh = newtonSeedFresh(w);
    float fr = max(length(fresh), 1e-20);
    vec2 fhat = fresh / fr;          // unit-direction
    // Seeds: fresh, fresh×0.6 (pulled in/out by 0.4), fresh rotated +90°,
    // fresh rotated −90°. For unbounded we replace the scale variant with
    // a pushed-out version (×1.6) so we stay in 𝔻*.
    vec2 s1 = fresh;
    vec2 s2 = (u_unbounded == 1) ? (fresh * 1.6) : (fresh * 0.6);
    vec2 s3 = vec2(-fhat.y, fhat.x) * fr;       // +90° rotation
    vec2 s4 = vec2( fhat.y,-fhat.x) * fr;       // −90° rotation
    bool found = false;
    vec2 zTry;
    zTry = invertPhi(w, s1, ok); if (ok && acceptZ(zTry)) { z = zTry; found = true; }
    if (!found) { zTry = invertPhi(w, s2, ok); if (ok && acceptZ(zTry)) { z = zTry; found = true; } }
    if (!found) { zTry = invertPhi(w, s3, ok); if (ok && acceptZ(zTry)) { z = zTry; found = true; } }
    if (!found) { zTry = invertPhi(w, s4, ok); if (ok && acceptZ(zTry)) { z = zTry; found = true; } }
    if (!found) { ok = false; return vec2(0.0); }
    ok = true;
  }
  // F has a pole at z=0 in both modes (bounded F = conj(w_0) + Σ A/(z−z_j)^k
  // has its singularity at the z_j, but unbounded G also has a 1/z term).
  // Guard at |z| < 1e-4 — float32-realistic threshold.
  if (dot(z, z) < 1e-8) { ok = false; return vec2(0.0); }
  zSeed = z;
  vec2 Sv = evalSchwarz(z);
  if (any(isnan(Sv)) || any(isinf(Sv))) { ok = false; return vec2(0.0); }
  return cconj(Sv);
}

// Map escape time n to t ∈ [0,1] under the user-selected scale mode.
float computeT(int n) {
  float fn   = float(n);
  float fmax = float(u_maxIter);
  if (u_scaleMode == 2) {
    // Log: emphasize early escapes.
    return clamp(log(fn + 1.0) / log(fmax + 1.0), 0.0, 1.0);
  }
  if (u_scaleMode == 3) {
    return clamp(sqrt(fn / fmax), 0.0, 1.0);
  }
  if (u_scaleMode == 4) {
    // Cyclic: ((n-1) mod K) / K. Reveals fine-grained periodicity.
    int k  = max(u_modK, 1);
    int rm = (n - 1) - (n - 1) / k * k;
    return float(rm) / float(k);
  }
  // smooth (default) and discrete share the linear base; discrete then
  // quantizes to integer-escape-time buckets so each n picks a distinct
  // colormap stop (with LINEAR colormap sampling).
  float t = float(n - 1) / max(fmax - 1.0, 1.0);
  if (u_scaleMode == 1) t = (floor(t * fmax) + 0.5) / fmax;
  return clamp(t, 0.0, 1.0);
}

// Color lookup. kind in {0=fund, 1=esc, 2=int, 3=invalid, 4=outside}.
vec4 kindToColor(int kind, int n) {
  if (kind == 4) return vec4(245.0/255.0, 245.0/255.0, 248.0/255.0, 1.0);
  if (kind == 2) return vec4(28.0/255.0, 28.0/255.0, 36.0/255.0, 1.0);
  if (kind == 1) return vec4(80.0/255.0, 80.0/255.0, 90.0/255.0, 1.0);
  if (kind == 3) return vec4(180.0/255.0, 90.0/255.0, 90.0/255.0, 1.0);
  return texture(u_colormap, vec2(computeT(n), 0.5));
}

void main() {
  // gl_FragCoord origin is bottom-left. Match the screen y-flip used by JS:
  // screen-y grows downward, but the canvas's CSS coordinate system is also
  // top-left. We compute world coords assuming y-up in world space.
  vec2 frag = gl_FragCoord.xy;        // (x: left→right, y: bottom→top)
  // World coord at this fragment center, with y flipped so screen-up is +Im.
  vec2 w;
  w.x = u_viewCenter.x + (frag.x - 0.5 * u_canvasSize.x) / u_pxPerUnit;
  w.y = u_viewCenter.y + (frag.y - 0.5 * u_canvasSize.y) / u_pxPerUnit;

  if (!inOmega(w)) {
    outColor = kindToColor(4, 0);
    return;
  }

  vec2 zSeed = newtonSeedFresh(w);
  bool ok = true;
  for (int n = 1; n <= 256; ++n) {
    if (n > u_maxIter) break;
    vec2 wNext = sigma(w, zSeed, ok);
    if (!ok) {
      outColor = kindToColor(3, n - 1);
      return;
    }
    if (any(isnan(wNext)) || any(isinf(wNext))) {
      outColor = kindToColor(1, n);
      return;
    }
    if (length(wNext) > u_escapeR) {
      outColor = kindToColor(1, n);
      return;
    }
    w = wNext;
    if (!inOmega(w)) {
      outColor = kindToColor(0, n);
      return;
    }
  }
  outColor = kindToColor(2, u_maxIter);
}
`;

  // ===========================================================================
  // GL helpers.
  // ===========================================================================
  function compile(gl, type, src) {
    const sh = gl.createShader(type);
    gl.shaderSource(sh, src);
    gl.compileShader(sh);
    if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
      const log = gl.getShaderInfoLog(sh);
      gl.deleteShader(sh);
      throw new Error('shader compile error:\n' + log);
    }
    return sh;
  }
  function link(gl, vs, fs) {
    const prog = gl.createProgram();
    gl.attachShader(prog, vs); gl.attachShader(prog, fs);
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
      const log = gl.getProgramInfoLog(prog);
      gl.deleteProgram(prog);
      throw new Error('program link error:\n' + log);
    }
    return prog;
  }

  function buildColormapTexture(gl, stops) {
    const N = 256;
    const data = new Uint8Array(N * 4);
    for (let i = 0; i < N; i++) {
      const t = i / (N - 1);
      const n = stops.length - 1;
      const f = t * n;
      const k = Math.min(n - 1, Math.floor(f));
      const u = f - k;
      const a = stops[k], b = stops[k + 1];
      data[i*4]   = Math.round(a[0] + (b[0] - a[0]) * u);
      data[i*4+1] = Math.round(a[1] + (b[1] - a[1]) * u);
      data[i*4+2] = Math.round(a[2] + (b[2] - a[2]) * u);
      data[i*4+3] = 255;
    }
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, N, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, data);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    return tex;
  }
  // Colormap stops (each row = RGB in 0–255). All ~8–10 stops per palette;
  // linear interpolation gives a smooth gradient. Sources: matplotlib
  // (magma/viridis/inferno/plasma/cividis), Google's Turbo, and ad-hoc
  // perceptual palettes (rainbow, iceandfire, twotone).
  const MAGMA = [
    [0,0,4],[28,16,68],[79,18,123],[129,37,129],[181,54,122],
    [229,80,100],[251,135,97],[254,194,135],[252,253,191],
  ];
  const INFERNO = [
    [0,0,4],[31,12,72],[85,15,109],[136,34,106],[186,54,85],
    [227,89,51],[249,140,10],[249,201,50],[252,255,164],
  ];
  const PLASMA = [
    [13,8,135],[75,3,161],[125,3,168],[168,34,150],[203,70,121],
    [229,107,93],[248,148,65],[253,195,40],[240,249,33],
  ];
  const VIRIDIS = [
    [68,1,84],[72,40,120],[62,73,137],[49,104,142],[38,130,142],
    [31,158,137],[53,183,121],[109,205,89],[180,222,44],[253,231,37],
  ];
  const CIVIDIS = [
    [0,32,76],[0,52,110],[40,75,124],[80,100,128],[120,127,128],
    [161,156,124],[197,187,108],[230,219,84],[253,253,51],
  ];
  // Turbo approximation (Google "Turbo: improved rainbow", perceptually uniform).
  const TURBO = [
    [48,18,59],[71,118,238],[26,196,231],[26,231,153],[97,239,71],
    [202,231,33],[255,184,33],[255,113,33],[224,40,9],[122,4,2],
  ];
  const GRAYSCALE = [
    [0,0,0],[64,64,64],[128,128,128],[192,192,192],[255,255,255],
  ];
  const RAINBOW = [
    [148,0,211],[75,0,130],[0,0,255],[0,255,0],[255,255,0],
    [255,127,0],[255,0,0],
  ];
  const ICEANDFIRE = [
    [10,40,100],[60,120,200],[160,210,240],[245,245,245],
    [250,210,90],[235,120,40],[170,30,30],
  ];
  const TWOTONE = [
    [245,245,248],[120,130,200],[40,50,110],[20,30,70],
  ];
  // r#(∞) for unbounded LQDs, computed in JS once per setPhi.
  //   r#(∞) = Σⱼ Σₖ conj(A_{j,k}) · (-1)^k / conj(z_j)^k
  function jsRHashAtInfinity(phi) {
    let accRe = 0, accIm = 0;
    for (const br of phi.branches || []) {
      if (!br.A || !br.A.length) continue;
      // 1 / conj(z_j) = z_j / |z_j|².
      const zd = br.z.re * br.z.re + br.z.im * br.z.im;
      if (zd < 1e-300) continue;
      const cinvRe = br.z.re / zd;       // 1/conj(z_j) = conj(z_j)? no. conj(z) = (re,-im), 1/conj(z) = (re,+im)/|z|²
      const cinvIm = br.z.im / zd;
      // Iterate k = 1..A.length, accumulating conj(A_{j,k}) · (-1)^k / conj(z_j)^k.
      let powRe = 1, powIm = 0;
      for (let k = 1; k <= br.A.length; k++) {
        // power *= 1/conj(z_j)
        const nr = powRe * cinvRe - powIm * cinvIm;
        const ni = powRe * cinvIm + powIm * cinvRe;
        powRe = nr; powIm = ni;
        // conj(A) = (A.re, -A.im)
        const A = br.A[k - 1];
        const cR = A.re,  cI = -A.im;
        // term = conj(A) * power
        const tR = cR * powRe - cI * powIm;
        const tI = cR * powIm + cI * powRe;
        const sign = (k % 2 === 0) ? 1 : -1;
        accRe += sign * tR;
        accIm += sign * tI;
      }
    }
    return { re: accRe, im: accIm };
  }

  const SCALE_MODE_ID = {
    smooth:   0,
    discrete: 1,
    log:      2,
    sqrt:     3,
    modulo:   4,
  };
  function pickColormap(name) {
    switch (name) {
      case 'magma':      return MAGMA;
      case 'inferno':    return INFERNO;
      case 'plasma':     return PLASMA;
      case 'viridis':    return VIRIDIS;
      case 'cividis':    return CIVIDIS;
      case 'turbo':      return TURBO;
      case 'grayscale':  return GRAYSCALE;
      case 'rainbow':    return RAINBOW;
      case 'iceandfire': return ICEANDFIRE;
      case 'twotone':    return TWOTONE;
      case 'cyclic':     return MAGMA.concat(MAGMA.slice().reverse(), MAGMA);
      default:           return MAGMA;
    }
  }

  // Build a 2048² (MASK_SIZE-square) RED-channel mask texture by drawing the
  // polygon into an off-screen 2D canvas, then uploading. Mask covers the
  // polygon's bbox with a configurable padding factor (set per family in
  // setPhi: 2.4× for bounded, 5× for unbounded). NEAREST sampling + no
  // anti-aliasing on the 2D fill = clean binary classifier with no
  // boundary-speckle artifacts.
  function buildMaskTexture(gl, polyPts, padFactor) {
    const off = document.createElement('canvas');
    off.width = MASK_SIZE; off.height = MASK_SIZE;
    const ctx = off.getContext('2d');
    // Compute polygon bbox.
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (const p of polyPts) {
      if (p.re < minX) minX = p.re;
      if (p.re > maxX) maxX = p.re;
      if (p.im < minY) minY = p.im;
      if (p.im > maxY) maxY = p.im;
    }
    const w = maxX - minX, h = maxY - minY;
    const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2;
    const half = Math.max(w, h) / 2 * padFactor;       // square mask centered on polygon
    const maskCenter = { re: cx, im: cy };
    const maskHalfExtent = { x: half, y: half };
    // The mask is a hard binary classifier — turn off anti-aliasing in the
    // 2D context so the polygon edge isn't soft-feathered, then sample with
    // NEAREST below. Together this kills the "ring of speckle on ∂Ω"
    // artifact caused by intermediate alpha values across the boundary.
    ctx.imageSmoothingEnabled = false;
    // Clear to BLACK (mask = 0 outside polygon).
    ctx.fillStyle = '#000';
    ctx.fillRect(0, 0, MASK_SIZE, MASK_SIZE);
    // World → mask pixel: px = (world.re - cx) / half * (MASK_SIZE/2) + MASK_SIZE/2.
    const s = MASK_SIZE / 2 / half;
    ctx.fillStyle = '#fff';                            // mask = 1 inside polygon
    ctx.beginPath();
    const p0 = polyPts[0];
    ctx.moveTo((p0.re - cx) * s + MASK_SIZE / 2,
               // Flip y so canvas-y-down matches our world-y-down indexing in the shader.
               // The shader does uv.y from gl_FragCoord.y (bottom-up). We pass the mask
               // unflipped here; the shader maps w.y in world coords to uv.y via the
               // straight transform. To make canvas-y-down match world-y-up, flip vertically.
               (MASK_SIZE / 2) - (p0.im - cy) * s);
    for (let i = 1; i < polyPts.length; i++) {
      const p = polyPts[i];
      ctx.lineTo((p.re - cx) * s + MASK_SIZE / 2,
                 (MASK_SIZE / 2) - (p.im - cy) * s);
    }
    ctx.closePath();
    ctx.fill();
    // Re-flip vertically when uploading so the shader's uv.y (which goes up
    // with world.y) reads the correct pixel.
    // Easier: just upload as-is and account for the flip in the shader's uv
    // formula. We'll set gl.pixelStorei(UNPACK_FLIP_Y_WEBGL, true) instead.
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
    // R8 internal format: we only read .r in the shader, so the other channels
    // are dead weight. Uses 1 byte/texel (vs 4 for RGBA), so the 2048² mask
    // takes 4 MB instead of 16 MB. Requires WebGL 2 (sized internal format).
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, gl.RED, gl.UNSIGNED_BYTE, off);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 4);
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
    // NEAREST sampling on the mask: it's a binary classifier, not a smooth
    // texture. Linear sampling produced intermediate values around ∂Ω that
    // flipped the > 0.5 threshold on sub-pixel position — exactly the
    // boundary-speckle pattern. With NEAREST + a 2048² mask, every fragment
    // gets a clean 0 or 1 with sub-pixel boundary precision.
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    return { tex, maskCenter, maskHalfExtent };
  }

  // ===========================================================================
  // Public factory.
  // ===========================================================================
  function createGPURenderer(canvas) {
    let gl;
    try { gl = canvas.getContext('webgl2', { antialias: false, preserveDrawingBuffer: false }); }
    catch (e) { gl = null; }
    if (!gl) return null;

    let prog, vs, fs, vbo;
    try {
      vs = compile(gl, gl.VERTEX_SHADER,   VERT_SRC);
      fs = compile(gl, gl.FRAGMENT_SHADER, FRAG_SRC);
      prog = link(gl, vs, fs);
    } catch (e) {
      console.error('schwarz-webgl: shader build failed:', e);
      if (vs) gl.deleteShader(vs);
      if (fs) gl.deleteShader(fs);
      if (prog) gl.deleteProgram(prog);
      return null;
    }

    // Full-screen triangle (3 verts) — cheaper than a quad and avoids the
    // diagonal seam artifact.
    vbo = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vbo);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1,-1, 3,-1, -1, 3]), gl.STATIC_DRAW);
    const a_pos = gl.getAttribLocation(prog, 'a_pos');

    // Uniform locations (look up once).
    const U = {};
    [
      'viewCenter','pxPerUnit','canvasSize','unbounded','w0','c',
      'polyA','polyALen','branchZ','branchA','branchACount','nBranches',
      'maxIter','escapeR','mask','maskCenter','maskHalfExtent','colormap',
      'scaleMode','modK',
      // LQD-specific uniforms (zero/unset for classical families).
      'family','gamma','z0','absZ0','rInfConj',
      // Polynomial-h β-correction (HANDOFF #22 / #26): unbounded LQDs only.
      'lqdBeta','lqdBetaLen',
    ].forEach(name => { U[name] = gl.getUniformLocation(prog, 'u_' + name); });

    // Phi state we'll keep packed and ready.
    const phiState = {
      unbounded: false,
      familyId:  0,
      w0:        new Float32Array(2),
      c:         0,
      gamma:     new Float32Array(2),
      z0:        new Float32Array(2),
      absZ0:     0,
      rInfConj:  new Float32Array(2),
      polyA:     new Float32Array(MAX_LAURENT * 2),
      polyALen:  0,
      branchZ:   new Float32Array(MAX_BRANCHES * 2),
      branchA:   new Float32Array(MAX_BRANCHES * MAX_K * 2),
      branchACount: new Int32Array(MAX_BRANCHES),
      nBranches: 0,
      lqdBeta:   new Float32Array(MAX_BETA * 2),
      lqdBetaLen: 0,
      mask:      null,
      maskCenter: new Float32Array(2),
      maskHalfExtent: new Float32Array(2),
      escapeR:   1e10,
      capacityError: null,
    };
    let colormapTex = buildColormapTexture(gl, MAGMA);
    let activeColormap = 'magma';

    function setPhi(phi, opts) {
      opts = opts || {};
      const polyPts = opts.boundaryPts || [];
      if (!phi) throw new Error('setPhi: null phi');

      // For unboundedLQD_singular, the γ-branch (HANDOFF #24 / #26) is fed
      // to the shader by appending a synthetic branch {z: phi.z0, A: γ} to
      // the uploaded branches list — mirrors withSyntheticBranch in
      // schwarz-common.js and _phiWithSyntheticBranch in the singular
      // solver. Compute the effective branch list now (before capacity
      // checks) so the cap applies to the merged total.
      const isSingularLQD = (phi.family === 'unboundedLQD_singular' ||
                              phi.family === 'boundedLQD_singular');
      const gamma = phi.lqdGamma || [];
      const effBranches =
        (phi.family === 'unboundedLQD_singular' && gamma.length > 0 && phi.z0)
          ? (phi.branches || []).concat([{ z: phi.z0, A: gamma }])
          : (phi.branches || []);

      const nb = effBranches.length;
      if (nb > MAX_BRANCHES) {
        phiState.capacityError = `Too many branches (${nb} > ${MAX_BRANCHES}); falling back to CPU.`;
        return false;
      }
      const polyA = phi.polyA || phi.F || [];
      if (polyA.length > MAX_LAURENT) {
        phiState.capacityError = `Laurent length too large (${polyA.length} > ${MAX_LAURENT}); falling back to CPU.`;
        return false;
      }
      const beta = phi.lqdBeta || [];
      if (beta.length > MAX_BETA) {
        phiState.capacityError = `lqdBeta length too large (${beta.length} > ${MAX_BETA}); falling back to CPU.`;
        return false;
      }
      for (let j = 0; j < nb; j++) {
        if (effBranches[j].A.length > MAX_K) {
          phiState.capacityError = `Branch ${j} A-length too large (${effBranches[j].A.length} > ${MAX_K}); falling back to CPU.`;
          return false;
        }
      }
      phiState.capacityError = null;
      // Mark isSingularLQD as touched (for future use) to silence linters.
      void isSingularLQD;
      phiState.unbounded = !!phi.unbounded;
      // Family ID determines which evalPhi/evalSchwarz branch the shader uses.
      // Classical bounded/unbounded QDs don't set phi.family (gotcha #1) so
      // we infer them from the unbounded flag.
      switch (phi.family) {
        case 'boundedLQD':            phiState.familyId = 2; break;
        case 'boundedLQD_singular':   phiState.familyId = 3; break;
        case 'unboundedLQD':          phiState.familyId = 4; break;
        case 'unboundedLQD_singular': phiState.familyId = 5; break;
        default: phiState.familyId = phi.unbounded ? 1 : 0;
      }
      if (phi.w0) { phiState.w0[0] = phi.w0.re; phiState.w0[1] = phi.w0.im; }
      else        { phiState.w0[0] = 0; phiState.w0[1] = 0; }
      if (phi.gamma) { phiState.gamma[0] = phi.gamma.re; phiState.gamma[1] = phi.gamma.im; }
      else           { phiState.gamma[0] = 0; phiState.gamma[1] = 0; }
      if (phi.z0) {
        phiState.z0[0] = phi.z0.re; phiState.z0[1] = phi.z0.im;
        phiState.absZ0 = Math.hypot(phi.z0.re, phi.z0.im);
      } else {
        phiState.z0[0] = 0; phiState.z0[1] = 0;
        phiState.absZ0 = 0;
      }
      phiState.c = (phi.c != null) ? phi.c : 0;
      phiState.polyALen = polyA.length;
      phiState.polyA.fill(0);
      for (let l = 0; l < polyA.length; l++) {
        phiState.polyA[2*l]   = polyA[l].re;
        phiState.polyA[2*l+1] = polyA[l].im;
      }
      // Pre-compute conj(r#(∞)) for unbounded LQDs. Zero for all other
      // families. NOTE: deferred to AFTER the branch-upload loop below so
      // that the synthetic γ-branch (HANDOFF #24 / #26) is included in the
      // ∞-baseline.
      phiState.rInfConj[0] = 0;
      phiState.rInfConj[1] = 0;
      phiState.nBranches = nb;
      phiState.branchZ.fill(0);
      phiState.branchA.fill(0);
      phiState.branchACount.fill(0);
      for (let j = 0; j < nb; j++) {
        const br = effBranches[j];
        phiState.branchZ[2*j]   = br.z.re;
        phiState.branchZ[2*j+1] = br.z.im;
        phiState.branchACount[j] = br.A.length;
        for (let k = 0; k < br.A.length; k++) {
          phiState.branchA[2 * (j * MAX_K + k)]     = br.A[k].re;
          phiState.branchA[2 * (j * MAX_K + k) + 1] = br.A[k].im;
        }
      }
      // Polynomial-h β-correction (HANDOFF #22 / #26): copy into the
      // GPU-side Float32Array. The shader sees B(1/z) = Σ_l β_l/z^l in
      // φ's exp argument for families 4 (unboundedLQD) and 5
      // (unboundedLQD_singular).
      phiState.lqdBeta.fill(0);
      phiState.lqdBetaLen = beta.length;
      for (let l = 0; l < beta.length; l++) {
        phiState.lqdBeta[2*l]   = beta[l].re;
        phiState.lqdBeta[2*l+1] = beta[l].im;
      }
      // For UQDLS with γ, rInfConj must reflect the MERGED branches so the
      // shader's `sum - cconj(u_rInfConj)` term sees the right ∞-baseline.
      // jsRHashAtInfinity operates on phi.branches, so re-run on the
      // synthetic-branch–merged phi just for the rInf computation. (β does
      // NOT affect r#(∞) — it's a separate exp-argument term.)
      if (phiState.familyId === 4 || phiState.familyId === 5) {
        const rInf = jsRHashAtInfinity({ branches: effBranches });
        phiState.rInfConj[0] = rInf.re;
        phiState.rInfConj[1] = -rInf.im;
      }
      // Build mask from polygon. Pad bounded modestly; unbounded needs more
      // headroom because iterates can wander.
      const padFactor = phiState.unbounded ? 5.0 : 2.4;
      if (phiState.mask) gl.deleteTexture(phiState.mask);
      if (!polyPts.length) {
        phiState.mask = null;
        phiState.maskCenter[0] = 0; phiState.maskCenter[1] = 0;
        phiState.maskHalfExtent[0] = 1; phiState.maskHalfExtent[1] = 1;
      } else {
        const m = buildMaskTexture(gl, polyPts, padFactor);
        phiState.mask = m.tex;
        phiState.maskCenter[0] = m.maskCenter.re;
        phiState.maskCenter[1] = m.maskCenter.im;
        phiState.maskHalfExtent[0] = m.maskHalfExtent.x;
        phiState.maskHalfExtent[1] = m.maskHalfExtent.y;
      }
      phiState.escapeR = opts.escapeR || (phiState.unbounded ? phiState.maskHalfExtent[0] * 6.0 : 1e10);
      return true;
    }

    function setColormap(name) {
      if (name === activeColormap) return;
      activeColormap = name;
      gl.deleteTexture(colormapTex);
      colormapTex = buildColormapTexture(gl, pickColormap(name));
    }

    function render(view, opts) {
      opts = opts || {};
      const W = Math.max(1, Math.floor(view.cssW * (window.devicePixelRatio || 1)));
      const H = Math.max(1, Math.floor(view.cssH * (window.devicePixelRatio || 1)));
      if (canvas.width !== W || canvas.height !== H) {
        canvas.width = W; canvas.height = H;
      }
      gl.viewport(0, 0, W, H);
      gl.useProgram(prog);

      // Bind quad.
      gl.bindBuffer(gl.ARRAY_BUFFER, vbo);
      gl.enableVertexAttribArray(a_pos);
      gl.vertexAttribPointer(a_pos, 2, gl.FLOAT, false, 0, 0);

      // Uniforms.
      gl.uniform2f(U.viewCenter, view.cx, view.cy);
      gl.uniform1f(U.pxPerUnit, view.scale * (window.devicePixelRatio || 1));
      gl.uniform2f(U.canvasSize, W, H);
      gl.uniform1i(U.unbounded, phiState.unbounded ? 1 : 0);
      gl.uniform1i(U.family,    phiState.familyId);
      gl.uniform2fv(U.w0,       phiState.w0);
      gl.uniform2fv(U.gamma,    phiState.gamma);
      gl.uniform2fv(U.z0,       phiState.z0);
      gl.uniform1f(U.absZ0,     phiState.absZ0);
      gl.uniform2fv(U.rInfConj, phiState.rInfConj);
      gl.uniform1f(U.c, phiState.c);
      gl.uniform2fv(U.polyA, phiState.polyA);
      gl.uniform1i(U.polyALen, phiState.polyALen);
      gl.uniform2fv(U.lqdBeta, phiState.lqdBeta);
      gl.uniform1i(U.lqdBetaLen, phiState.lqdBetaLen);
      gl.uniform2fv(U.branchZ, phiState.branchZ);
      gl.uniform2fv(U.branchA, phiState.branchA);
      gl.uniform1iv(U.branchACount, phiState.branchACount);
      gl.uniform1i(U.nBranches, phiState.nBranches);
      gl.uniform1i(U.maxIter, opts.maxIter || 64);
      gl.uniform1f(U.escapeR, phiState.escapeR);
      // Escape-time scale mode (smooth/discrete/log/sqrt/modulo).
      const scaleModeId = SCALE_MODE_ID[opts.scaleMode || 'smooth'] | 0;
      gl.uniform1i(U.scaleMode, scaleModeId);
      gl.uniform1i(U.modK, Math.max(2, (opts.modK | 0) || 8));

      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, phiState.mask);
      gl.uniform1i(U.mask, 0);
      gl.uniform2fv(U.maskCenter, phiState.maskCenter);
      gl.uniform2fv(U.maskHalfExtent, phiState.maskHalfExtent);

      gl.activeTexture(gl.TEXTURE1);
      gl.bindTexture(gl.TEXTURE_2D, colormapTex);
      gl.uniform1i(U.colormap, 1);

      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function destroy() {
      if (phiState.mask) gl.deleteTexture(phiState.mask);
      if (colormapTex) gl.deleteTexture(colormapTex);
      gl.deleteBuffer(vbo);
      gl.deleteProgram(prog);
      gl.deleteShader(vs); gl.deleteShader(fs);
    }

    function capacityError() { return phiState.capacityError; }

    return {
      available: true,
      setPhi,
      setColormap,
      render,
      destroy,
      capacityError,
    };
  }

  Schwarz.createGPURenderer = createGPURenderer;
  Schwarz._gpuCaps = { MAX_BRANCHES, MAX_K, MAX_LAURENT, MAX_BETA };

  // Internal API consumed by sphere-webgl.js to reuse the fractal shader and
  // related helpers without duplicating ~600 lines of GLSL + JS.
  // Treat these as package-private — not part of the public surface.
  Schwarz._shaders    = { vert: VERT_SRC, frag: FRAG_SRC };
  Schwarz._glHelpers  = {
    compile,
    link,
    buildColormapTexture,
    buildMaskTexture,
    pickColormap,
    SCALE_MODE_ID,
    jsRHashAtInfinity,
  };

})(typeof globalThis !== 'undefined' ? globalThis : this);
