# Riemann sphere view (`app/sphere/`)

Stereographic projection of the Schwarz fractal onto a sphere with
∞ → north pole. Lazy-mounted by the Schwarz tab on first toggle to
sphere view (HANDOFF #29).

## Files

| File | Role |
| --- | --- |
| `sphere-common.js` | Pure math kernel: stereographic projection + sphere-mesh builder + Float64 mat4 helpers. |
| `sphere-webgl.js` | Three-pass WebGL 2 renderer: opaque sphere base, fractal-textured pass, glow overlay. |
| `sphere-ui.js` | `QD.SphereView.mount(opts) → handle`: orbit camera, drag / wheel zoom, ResizeObserver, hover tooltip. |

## Public surface (`QD.Sphere.*`)

| Function | Use |
| --- | --- |
| `projectToSphere(w)` | Stereographic ℂ → S² (north-pole convention). |
| `unprojectFromSphere(pt)` | Inverse projection; returns `null` at / near the north pole. |
| `buildSphereMesh(divisions)` | UV-sphere mesh (Float32 vertices + Uint16 indices). |
| `mat4lookAt(eye, target, up)` | Float64 mat4 helper (camera). |
| `mat4perspective(fovy, aspect, znear, zfar)` | Float64 mat4 helper. |
| `mat4invertRigid(m)` | Cheap inverse for rigid-body transforms. |
| `createRenderer(canvas)` | WebGL 2 renderer factory. |

## Adapter contract

The Schwarz tab calls `QD.SphereView.mount({ container, getPhiSnapshot,
getBoundary, getRenderParams, isActive })` on first toggle. Returns
`{ destroy, refresh, resize, ... }`. The captured φ + boundary polygon
+ render params (`maxIter`, colormap, scale, modK) are shared with the
plane view — toggling between plane and sphere never requires a
re-capture.

## Rendering pipeline

1. **Fractal → FBO.** Same WebGL 2 shader source as the plane Schwarz
   view (imported via `QD.Schwarz._shaders` / `_glHelpers`), rendered
   into an offscreen framebuffer attached as a texture.
2. **Sphere pass.** Textured UV-sphere mesh drawn from the orbit
   camera; the fractal FBO is sampled in the fragment shader via the
   stereographic projection.
3. **Overlay pass.** ∂Ω as a closed polyline on the sphere surface;
   finite poles and the north-pole ∞-marker as billboard markers.

## Float-precision boundary

- `sphere-common.js` is Float64 throughout — round-trip
  `projectToSphere` ∘ `unprojectFromSphere` test passes at < 1e-12.
- `sphere-webgl.js` is Float32 (GPU constraint) — the small precision
  loss in the texture lookup is invisible at typical zoom.

## When this view is most useful

- Unbounded Ω where orbits wander to infinity — the spherical wrap
  bounds the picture.
- Families with ∞ ∈ Ω (the unbounded singular LQDs); ∞ becomes the
  north pole, visually centred.

Bounded Ω still works mathematically but is visually uninformative —
the sphere shows a spherical cap of the fractal with the rest being
the all-the-same "outside Ω" color. A soft warning chip appears in
the sphere sidebar in that case.

## P0/P1 integration

- **P0.1a (PrimarySolution).** Sphere snapshots φ + hData via the
  same envelope path as Schwarz; see `sphere-ui.js _clonePhi` which
  carries `lqdBeta` / `lqdGamma` / `q` through (HANDOFF #28).
- **P0.2 (worker).** No direct interaction; the sphere reads the
  already-solved φ.

## Where it's called from

| Caller | What it uses |
| --- | --- |
| `schwarz/schwarz-ui.js` | `QD.SphereView.mount` on first toggle to sphere view |
| `node-test.js` | `projectToSphere`, `unprojectFromSphere`, mat4 helpers (round-trip tests; no WebGL in Node) |
