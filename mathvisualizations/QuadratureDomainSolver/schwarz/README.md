# Schwarz dynamics (`app/schwarz/`)

Visualises the dynamics of σ(w) — the Schwarz reflection extended
meromorphically into Ω. Iterating σ from any w₀ partitions the plane
into the "tiling set" (orbits stay bounded forever) and its complement
(orbits escape after `n` steps). Colouring by `n` produces the classic
escape-time fractal.

## Files

| File | Role |
| --- | --- |
| `schwarz-common.js` | Pure math kernel + per-family CPU adapters (`adaptBounded`, `adaptUnbounded`, `adaptBoundedLQD`, …). |
| `schwarz-webgl.js` | WebGL 2 fragment-shader renderer; same σ-iteration on the GPU. |
| `schwarz-cpu-worker.js` | `QD.SchwarzCpuWorker` — dedicated Web Worker that computes the CPU escape-time field off the main thread (rebuilds the Schwarz handle from the serializable φ + boundary samples, streams a transferable field snapshot per pyramid pass). Falls back to the in-page renderer on file:// / no-Worker. |
| `schwarz-ui.js` | Schwarz-tab UI hub: source-φ capture, card builders, `setMode` / view-toggle to Sphere mode, coordinate transforms, and the `sCtx` injection + the four module installs below. |
| `schwarz-paint.js` | 2D-canvas output layer: field / boundary / orbit / preimage-tree / limit-set painters + colormaps. `QD_UI.installSchwarzPaint(sCtx)`. |
| `schwarz-render.js` | Progressive escape-time renderer: debounced `requestRecompute` + GPU one-frame path + CPU 4×4→2×2→1×1 pyramid. `QD_UI.installSchwarzRender(sCtx)`. |
| `schwarz-features.js` | Per-feature compute routines for the analysis / limit-set / forward-dynamics cards: domain-coloring, preimage-tree rebuild, limit-set chaos game, σ level curves, critical orbits, cycle finder, orbit sweep, z-panel pullback, PNG export. `QD_UI.installSchwarzFeatures(sCtx)`. |
| `schwarz-interaction.js` | Canvas hover / wheel / click / dblclick / pin handlers + `attachCanvasHandlers`. `QD_UI.installSchwarzInteraction(sCtx)`. |

The last four are the Phase-3 (item E) factory-module split of the former
2477-line `schwarz-ui.js`; see [ARCHITECTURE.md](../../ARCHITECTURE.md) for the
`installSchwarzX(sCtx)` pattern and install order.

## Public surface (`QD.Schwarz.*`)

| Function | Use |
| --- | --- |
| `buildSchwarzFromPhi(phi, hData, boundaryPts)` | Build `{ sigma, psi, evalPhi, evalF, isInOmega, escapeR, family, unbounded, adapter }` from an inverse-solver φ (any of the ten families). |
| `buildSchwarzFromRational(phi, boundaryPts)` | Same builder for a Direct-tab rational `φ = P/Q`. |
| `escapeTime(w₀, schwarz, {maxIter, escapeR})` | Iterate σ from w₀ until it lands in Ω^c, diverges, or hits `maxIter`. Returns `{kind, n, lastW, firstZ}`. `kind ∈ {'fundamental' (σⁿ left Ω), 'escaped' (diverged / |σⁿ|>escapeR), 'invalid' (σ undefined), 'interior' (still in Ω after maxIter)}`. |
| `makeOrbit(w₀, schwarz, {maxIter})` | Orbit polyline for the click-to-orbit overlay. |
| `pointInPolygon(pt, polyPts)` | Even-odd boundary-curve in-Ω test. |
| `createGPURenderer(canvas)` | WebGL 2 renderer. Returns `{ setPhi, setColormap, render, destroy, capacityError }` or `null` if WebGL 2 is unavailable / shader compile fails. |
| `_gpuCaps` | Object with `{ MAX_BRANCHES: 12, MAX_K: 8, MAX_LAURENT: 12, MAX_BETA: 16 }` — caps the shader uniform allocations. φ exceeding any cap falls back to the CPU path. |

## CPU vs GPU

| Path | When used | Notes |
| --- | --- | --- |
| **GPU** (`createGPURenderer`) | Default when WebGL 2 + caps OK. | Full 1024² frame at `maxIter=128` in ~150 ms. Float32 precision; banding at zoom > 1e6. |
| **CPU** (`escapeTime` per pixel) | Fallback when no WebGL 2, caps exceeded, or explicitly chosen in the UI. | Progressive 4×4 → 2×2 → 1×1. Computed off-thread in `QD.SchwarzCpuWorker` when available (one transferable field snapshot per pass), else in-page chunked across `requestAnimationFrame` ticks. Always available. |

CPU↔GPU parity (for the six classical/LQD families; the four PQD families
are CPU-only): both adapters consume the same `phi` shape (with the
`lqdBeta`/`lqdGamma` fields carried through for unbounded LQDs).
HANDOFF #26 added 5 round-trip tests asserting σ(w) ≈ w on ∂Ω at
3e-13 for the previously-broken unbounded-LQD polyPart case; those
live in [`app/node-test.js`](../node-test.js).

## Source-φ capture (P0.1a integration)

Schwarz pulls its source φ from the Inverse tab via `QD.PrimarySolution`:

```js
const envelope = QD.PrimarySolution.get();
if (envelope && envelope.success) {
  const phi = clonePhi(envelope.primary.phi);
  const hData = envelope.hData;
  // ...
}
```

See [`schwarz-ui.js`](schwarz-ui.js) `captureFromInverseTab` for the
canonical reader. The legacy `state.current` path is kept as a
fallback for the rare case where `QD.PrimarySolution` is unavailable.

## Family coverage

All ten inverse families are supported. The six classical/LQD families
ship CPU + GPU support:

- `boundedQD` / `unboundedQD` (classical)
- `boundedLQD` / `boundedLQD_singular`
- `unboundedLQD` / `unboundedLQD_singular`

The four power-weighted families are **CPU-only** — `createGPURenderer`'s
`setPhi` refuses them and the UI falls back to the CPU `escapeTime` path
(the GPU shader has no non-integer αth-root power):

- `powerQD` / `powerQD_singular`
- `unboundedPQD` / `unboundedPQD_singular`

Bounded-rational direct φ is supported via `buildSchwarzFromRational`
(no Direct-tab "Send to Schwarz" wiring yet — pipe `phi` in manually).

## View modes (plane / z-disk / sphere)

A three-way segmented control in the sidebar switches between three views of the
same σ-iteration:

- **plane** — the w-plane escape-time tiling on Ω (default).
- **z-disk** — the SAME tiling uniformized onto 𝔻 (or 𝔻* for unbounded Ω): each
  pixel takes `z`, lifts `w = φ(z)` (the shader already evaluates `evalPhi`), and
  runs the existing σ escape-time. Renders on the GPU for the six classical/LQD
  families via the `u_viewMode` shader branch (PQDs fall back to CPU, like the
  plane view); the 2-D overlay canvas stays transparent so the GL field shows
  through. See HANDOFF #59 (CPU view) / #60 (GPU). `paintZView(overlayOnly)`
  derives overlay-only mode from `activeRenderer()`.
- **sphere** — the iteration textured onto a Riemann sphere; that adapter lives
  in [`app/sphere/`](../sphere/README.md) and is lazy-mounted on first toggle.

Captured φ + render params (maxIter / colormap / scale / modK) are shared across
all three views, so toggling never re-captures.

## Where it's called from

| Caller | What it uses |
| --- | --- |
| `schwarz-ui.js` (tab activation) | the full module (capture, render, controls) |
| `sphere/sphere-ui.js` | `_gpuCaps`, the shared GPU shader source via `_shaders` / `_glHelpers` |
| `node-test.js` | `buildSchwarzFromPhi`, `escapeTime`, `makeOrbit` (CPU round-trip tests) |
