# Parameter slice (`app/param-slice/`)

A 1-D or 2-D adaptive-mesh sweep over any sweepable parameter of the
current scenario (residue Re/Im, pole position, polynomial coefficient,
`c`, `q`, `w₀`). Each pixel solves the inverse problem at that
parameter value and is classified into one of six categories. Click a
pixel to load that φ back into the QD tab.

## Files

| File | Role |
| --- | --- |
| `param-slice-common.js` | Pure math kernel: `ParamRef` descriptors, `listAvailableParams`, `applyParam`, `classifyResult`, colour LUT. |
| `param-slice-pool.js` | Web Worker pool — runtime Blob bundle of the solver source (no build step). |
| `param-slice-render.js` | `QD_UI.installParamSliceRender(psCtx)` — the adaptive 2-D render engine `runAdaptive2D` (progressive quadtree sweep + warm-hint spatial index + coverage fill). Phase-3 (item E) split out of `param-slice-ui.js`. |
| `param-slice-ui.js` | Tab UI hub (lazy mount). Axes pickers, range fields, run orchestration (`startRun`), click-to-load, mini-preview card; captures `runAdaptive2D` from the render module via a forward-`let`. |

## Public surface

`window.ParamSlice.*` (pure math kernel):

| Function | Use |
| --- | --- |
| `listAvailableParams(scenario)` | Enumerate sweepable parameters for the current scenario. |
| `applyParam(scenario, paramRef, value)` | Return a new scenario with the parameter overridden (immutable). |
| `applyParamInPlace(scenario, paramRef, value)` | Mutate `scenario` in place — used by the worker's per-tile scratch scenario for allocation-free sweeps. |
| `cloneScenario(scenario)` | Deep clone for the worker boundary. |
| `solveOnePoint(scenario, sweep, hint, expectedFamilyTag)` | Solve at one point with an optional warm-start `hint` (a serialized φ from a neighbour). |
| `solveOnePointWithScratch(scratch, sweep, hint, expectedFamilyTag)` | Same as above but reuses a scratch scenario; preferred in tight inner loops. |
| `classifyResult(result, scenario, mode)` | Map `{ success, primary, error }` to one of six classes (see below). |
| `MODE_FAMILY_TAG[mode]` | Expected `phi.family` tag per mode — gates warm-start applicability. |

`window.ParamSlicePool.*` (worker pool):

| Function | Use |
| --- | --- |
| `create(opts)` | Returns a `Pool` (or `MainThreadPool` fallback). Defaults `maxWorkers` to `navigator.hardwareConcurrency` (capped at 16). |
| `createWorkerOnly(opts)` | No fallback — for tests that need to assert the worker path is exercised. |

`Pool` / `MainThreadPool` instances expose:

| Method | Use |
| --- | --- |
| `runSweep({scenario, mode, axes, onTile, onError})` | Linear (non-adaptive) sweep — one tile per pixel row. Returns `{ cancel, done }`. |
| `solveBatch(scenario, mode, points, warmHints)` | Solve a batch of pixel coordinates (used by the adaptive renderer). |
| `cancel()` | Abort pending tiles. |
| `terminate()` | Cancel + terminate the worker processes. |

## Cross-tab contract (with the QD tab)

`window.QD_UI` is a small named-export bag set up by `ui.js` for
subsystems to read/write the QD tab state:

| Hook | Direction | Use |
| --- | --- | --- |
| `window.QD_UI.snapshotScenario()` | QD → param-slice | Returns `{ hData, norm, mode }` — the current scenario the slice tab should sweep. |
| `window.QD_UI.loadScenarioIntoQdTab(scenario, mode)` | param-slice → QD | Reflects a clicked scenario back into the QD tab, switches tabs, and runs the solver. |

The first is called when the slice tab activates (or on the
"refresh from QD tab" button). The second is wired to pixel-click in
the rendered slice.

## Classification

Each pixel falls into exactly one class:

| Class | Pixel colour | Condition |
| --- | --- | --- |
| Valid QD | green; brightness ∝ 1/iter | Newton converged, univalent, identity passes |
| Identity fails | yellow | Newton converged but `verifyQuadratureIdentity` exceeds tol |
| Boundary self-intersects | orange | `isBoundaryUnivalent` returns false |
| Newton diverged | red | Newton hit max-iter without convergence |
| No algebraic root | gray | `solveInverseQD` returned `success: false` after all stages |
| Capability refused | slate | Family doesn't yet support this parameter mix (e.g. polynomial-h bounded LQDs, deferred) |

## Adaptive mesh

`runAdaptive2D` (in `param-slice-render.js`) implements a quadtree
refinement. Coarse pass solves a sparse grid; cells whose 4 corners
disagree on classification OR have iter-spread > `REFINE_ITER_DELTA = 8`
are split. Cross-cell warm-start hints come from a 16×16 spatial
index, so refined sub-pixels typically converge in 1-5 Newton iters
after a cold seed. See HANDOFF #31 / #37 for the design and the
coverage-fill fix.

## Worker bundle

`param-slice-pool.js` builds the worker bundle at runtime:

1. `fetch` each file in `SOLVER_SRC_FILES`.
2. Prepend `var window = self;` so the `(typeof window !==
   'undefined' && window.QD)` namespace idiom resolves inside the
   worker scope.
3. Append a worker-side message handler (`WORKER_HANDLER`).
4. Wrap as a Blob and create N workers from the resulting URL.

Falls back to `MainThreadPool` if the bundle fetch fails (typically
when the page is loaded via `file://`). The fallback maintains the
same async API so callers don't special-case.

## Tests

- `node-test.js` exercises `solveOnePoint` / `solveOnePointWithScratch`
  against synthetic scenarios.
- Capability-refused predicates have wiring tests added in HANDOFF #32.
- Adaptive-mesh predicates (`cellIsHomogeneous`, etc.) have unit tests
  in HANDOFF #31.

## Where it's called from

| Caller | What it uses |
| --- | --- |
| `param-slice-ui.js` (tab activation) | the full module |
| `node-test.js` | `ParamSlice.*` kernels (no Worker in Node — exercises `MainThreadPool` path) |
