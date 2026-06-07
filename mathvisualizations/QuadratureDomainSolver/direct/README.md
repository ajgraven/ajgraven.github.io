# Direct problem (`app/direct/`)

Given a Riemann map φ, compute the quadrature function h. This is the
inverse direction of the Inverse-tab solver — the user supplies φ (or,
for the weighted families, the rational KERNEL that defines φ) and the
kernel returns h's poles + polyPart (+ an origin term for some singular
families).

Coverage: **classical** QD (bounded polynomial / rational, unbounded
Laurent-at-∞, numerical) **and** the **weighted** families — power (PQD,
weight |w|^{2(α−1)}) and log (LQD), each bounded / unbounded and
non-singular / singular (0 ∈ Ω, Blaschke factor).

## Files

| File | Role |
| --- | --- |
| `direct-common.js` | Math kernel: classical (polynomial / rational / Laurent / numerical) **and** weighted (bounded + unbounded, power + log, non-singular + singular) φ → h. Durand–Kerner polynomial root finder. Fourier boundary verifier. Sample-the-boundary primitives. |
| `direct-ui.js` | Direct-mode UI hub inside the QD tab (HANDOFF #30 merged the standalone Direct tab into a `inverse \| direct` view-toggle). A compact segmented **Domain-type** control (`#dir-dm-weight` QD/PQD/LQD × `#dir-dm-domain` Bounded/Unbounded/Numerical × `#dir-dm-singular`, inverse-tab style) selects the family; the weight PARAMETER inputs (α/w₀/z₀/c/kernel) live in the φ-input cards. `applyDirectMode()` is the single canonical mode-refresh. Owns `directState`, the card builders, and the `dCtx` injection + the two module installs below. |
| `direct-recompute.js` | `QD_UI.installDirectRecompute(dCtx)` — the recompute→render pipeline: `recomputeAndRender` dispatches on `directState.mode` to `recomputeBounded`/`recomputeUnbounded`/`recomputeNumerical` (each builds h via the `QD.Direct.*` kernels), renders it (`displayH`), and pushes the weight-honoring boundary to the shared plot (`sampleBoundedPhi`/`pushBoundaryToPlot`). Phase-3 (item E) split. |
| `direct-verify.js` | `QD_UI.installDirectVerify(dCtx)` — the **Verify** button: `runVerify` dispatches per family (unbounded/bounded weighted → family identity verifier or inverse round-trip; classical → Fourier boundary-identity diagnostic on `sampleAnalyticPhi` samples). Phase-3 (item E) split. |

## Public surface (`QD.Direct.*`)

| Function | Use |
| --- | --- |
| `boundedQD(coeffs)` | h from polynomial φ. φ = Σ coeffs[k] · z^k. |
| `boundedQDRational(P, Q)` | h from rational `φ = P(z) / Q(z)`. Validates Q ≠ 0 on 𝔻̄. |
| `unboundedQD(c, F)` | h from Laurent-at-∞ φ = `c·z + F_0 + F_1/z + F_2/z² + …`. |
| `numericalBoundedQD(phiFn, opts)` | h from arbitrary analytic φ. Samples on `|z|=1`, extracts Taylor coefficients via DFT, truncates, falls through to the polynomial path. |
| `boundedPowerQD(R#, α)` | bounded power-weighted QD (0∉Ω); φ = (R#)^{1/α}, R# rational analytic + non-vanishing on 𝔻̄. |
| `boundedLogQD(r#, w₀)` | bounded log-weighted QD (0∉Ω); φ = w₀·exp(r#). |
| `boundedPowerQDSingular(R#, α, z₀)` | bounded singular PQD (0∈Ω); φ = b_{z₀}·(R#)^{1/α} (Thm 4.3.5). h += origin term r₀/w (`originResidue`). |
| `boundedLogQDSingular(r#, w₀, z₀)` | bounded singular LQD (0∈Ω); φ = γ·b_{z₀}·exp(r#). h += origin pole q/w (`q`). |
| `unboundedPowerQD(r#, α)` | unbounded PQD (0∉Ω, ∞∈Ω); φ = z·(r#)^{1/α} (Thm 4.3.7). c derived from r#. |
| `unboundedPowerQDSingular(r#, α, {z0?})` | unbounded singular PQD (0∈Ω); φ = z·b_{z₀}·(r#)^{1/α}. z₀ DERIVED from a zero of r#; no origin term. |
| `unboundedLogQD(r#, c)` | unbounded LQD (0∉Ω, ∞∈Ω); φ = c·z·exp(r̃#). |
| `unboundedLogQDSingular(r#, c, z₀)` | unbounded singular LQD (0∈Ω); φ = c·\|z₀\|·z·b_{z₀}·exp(r̃#). h += origin pole q/w (`q`). |
| `polynomialRoots(coeffs)` | Durand–Kerner complex polynomial root finder. |
| `parseRationalInZ(expr, math)` | math.js AST walker. Returns a polynomial array OR `{ num, den }`. |
| `parsePolynomialInZ(expr, math)` | Thin wrapper that rejects rational results. |
| `polynomialToString(coeffs)` | Canonical-form printer. |
| `evalH(hData, w)` | Evaluate h at a complex point. |
| `verifyBoundaryIdentity(hData, pts)` | Fourier negative-frequency-mass diagnostic on `h(φ(z)) − conj(φ(z))`. |
| `sampleBoundaryPolynomial(coeffs, N)` | ∂Ω samples for polynomial φ. |
| `sampleBoundaryLaurent(c, F, N)` | ∂Ω samples for Laurent φ. |

## The four φ-shape modes

| Mode | Closed-form? | How h is computed |
| --- | --- | --- |
| **Polynomial** | Yes | φ has a single pole of order `n` at `w₀ = φ(0)`. Principal parts computed via the forward Faber formula from the Taylor coefficients of φ. |
| **Rational P/Q** | Yes | σ extends with one pole at φ(0) (if deg P > deg Q) plus one pole per root r_i of Q at `φ(1/conj(r_i))`. Roots found via Durand–Kerner; per-pole principal parts via forward Faber. |
| **Laurent (unbounded)** | Partial | `φ = c·z + F_0 + F_1/z + …`. Polynomial-at-∞ part computed via the dual of `inverseFaberAtInfinity` (triangular system). Finite-pole part handled only for the trivial `c·z + F_0` case (exterior of a disk). |
| **Numerical** | Approximate | User-supplied function `φ : ℂ → ℂ`. Sample on `|z|=1`, DFT → Taylor coefficients, truncate, fall through to polynomial mode. Non-analytic φ (e.g. `conj(z)`) is flagged. |

## Weighted families (PQD / LQD)

The weighted kernels take the rational **kernel** R#/r# (not φ): for PQD
`φ = (R#)^{1/α}` (bounded) / `z·(R#)^{1/α}` (unbounded); for LQD
`φ = w₀·exp(r#)` / `c·z·exp(r̃#)`. Singular variants (0 ∈ Ω) multiply by a
Blaschke factor b_{z₀}. h is read off by inverting the **inverse solver's
own tested (★) chain** — `solveResiduesViaProbe` for the finite poles, the
`forwardPolyPartAtInfinity` / `forwardBetaToPolyPart` (★)_F inversion for
the polynomial-at-∞ part — so the forward map round-trips against the
inverse solver by construction. **Realizability ⟺ univalence of φ**
(Thm 4.3.3): a non-univalent kernel returns `{ univalent:false, warnings }`
rather than bad data; a non-singular kernel that yields 0 ∈ Ω is flagged
`{ originInside:true }` ("use the singular kernel"). Origin terms are
returned **separately** from the finite nodes: bounded-PQD r₀, bounded- /
unbounded-LQD q; unbounded-PQD-singular has none (Thm 4.3.7's −t/w cancels
the Φ_φ pole at 0). See the section banners in `direct-common.js`
("SINGULAR weighted forward kernels", "UNBOUNDED weighted forward kernels
(Theorem 4.3.7)") for the full math + conventions.

## Verifier

For **classical** QD the **Verify** button calls `verifyBoundaryIdentity`,
which checks the Fourier negative-frequency mass of `h(φ(z)) − conj(φ(z))`
on `|z|=1`: ≈ 0 for any valid classical QD, since σ − h is analytic in Ω
and its pullback to 𝔻 has only non-negative frequencies.

For the **weighted** families the classical Fourier check does not apply
(the quadrature identity carries the |w|^{2(α−1)} / log weight). Verify
instead uses the **family identity verifier**
(`QD.Family.<fam>.verifyQuadratureIdentity`) on the built φ + full h — the
trusted oracle the inverse solver uses — at the appropriate test class;
the singular boundary floors near ~1e-6 (vs ~1e-15 non-singular), so the
UI uses a looser strong-pass threshold there. The non-singular weighted
forward also round-trips via `Send to inverse` → `solveInverseQD`.

## Send-to-inverse round-trip

The Direct view's **Send to inverse** button pre-fills the inverse-view
fields with the computed h, switches the view toggle back to inverse,
and triggers `solveAndRender()`. With the P0.2 worker, this runs
without freezing the UI. The visualized boundary should match the
direct-mode ∂Ω; mismatches indicate a numerical issue.

## Known limitations

- **Classical unbounded *rational* φ** is not implemented. `unboundedQD`
  computes the polynomial-at-∞ part for any Laurent φ but recovers the
  finite poles only for the trivial `c·z + F_0` case; a general Laurent
  tail emits a "not yet implemented" warning and returns no finite poles.
  (The *weighted* unbounded families DO handle their full ∞-structure —
  this gap is specific to the classical unbounded-rational forward map.)
- **Higher-order pole of h at the origin** in the singular families is not
  handled (Prop 4.6.3 general case); the singular kernels assume h is
  analytic at 0 apart from the single origin term (rare otherwise).
- **Numerical mode** degrades any analytic-in-𝔻̄ φ to its polynomial
  truncation. Non-analytic φ is flagged with a clear error.

## Where it's called from

| Caller | What it uses |
| --- | --- |
| `direct-ui.js` | the entire kernel |
| `ui.js` (Direct mode toggle) | `QD.Direct._mountUI`, `QD.Direct._activate` (HANDOFF #30) |
| `node-test.js` | All classical + weighted kernels (§DF battery) + the polynomial root finder + parser tests |
