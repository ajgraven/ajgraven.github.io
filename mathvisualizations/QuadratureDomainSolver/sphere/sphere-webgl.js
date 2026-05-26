// =============================================================================
// sphere-webgl.js  — WebGL 2 renderer for the Riemann-sphere visualization.
//
// Three-pass render per frame:
//
//   Pass 1 (cached): Fractal-to-FBO
//     Renders the Schwarz σ-iteration escape-time fractal to an off-screen
//     RGBA8 texture of size texSize × texSize.  The texture covers the
//     w-plane square [cx−R, cx+R]².  Uses the exact same GLSL shader as
//     schwarz-webgl.js (accessed via QD.Schwarz._shaders).  Cached — only
//     re-rendered when φ, maxIter, or colormap changes.
//
//   Pass 2: Textured sphere
//     Renders a UV-sphere mesh.  The fragment shader does inverse stereographic
//     projection from each surface-normal to a w-plane coordinate, then samples
//     the fractal texture.  Near the north pole (∞) the far-field color is used.
//
//   Pass 3: Overlay
//     Boundary polyline, finite-pole markers, and north-pole ✸ marker, all
//     drawn as 3D line geometry (gl.LINES) projected onto the sphere surface.
//
// API
//   createSphereRenderer(canvas) → null on failure, else:
//     {
//       available: true,
//       setPhi(phi, {boundaryPts, escapeR})
//       setRenderParams({maxIter, colormap, scaleMode, modK, texSize})
//       setDisplayParams({rimDarken, showBoundary, showPoles,
//                         showNorthPole, boundaryColor, poleColor})
//       render(camera)   // camera = {azimuth, elevation, distance}
//       destroy()
//     }
//
// Dependencies:
//   QD.Schwarz._shaders     (fractal GLSL source strings)
//   QD.Schwarz._glHelpers   (compile, link, buildMaskTexture, buildColormapTexture,
//                             pickColormap, SCALE_MODE_ID, jsRHashAtInfinity)
//   QD.Schwarz._gpuCaps     (MAX_BRANCHES, MAX_K, MAX_LAURENT)
//   window.SphereCommon     (buildSphereMesh, mat4lookAt, mat4perspective,
//                             mat4multiply)
// =============================================================================

(function (global) {
  'use strict';
  const QD = (typeof window !== 'undefined' && window.QD)
    ? window.QD
    : (typeof module !== 'undefined' && module.exports ? module.exports : (global.QD || (global.QD = {})));

  const Sphere = QD.Sphere || (QD.Sphere = {});

  // Shader sources for Pass 2 (textured sphere) and Pass 3 (overlay).
  // The Pass 1 shader (fractal) is reused from QD.Schwarz._shaders.

  const SPHERE_VERT = `#version 300 es
in vec3 a_pos;
uniform mat4 u_mvp;
out vec3 v_norm;
void main() {
  v_norm = a_pos;
  gl_Position = u_mvp * vec4(a_pos, 1.0);
}`;

  // Inverse stereographic in the fragment shader.
  // u_fractalCenter: (cx, cy) — center of the fractal view in w-space.
  // u_maskHalfExtent: R — half-width of the fractal texture coverage in world units.
  // u_farFieldColor:  color for the north pole (∞) cap.
  // u_rimDarken:      0..0.5; attenuates rim pixels for 3-D shading cue.
  // u_eyeDir:         normalised direction from sphere center toward camera.
  const SPHERE_FRAG = `#version 300 es
precision highp float;
in vec3 v_norm;
uniform sampler2D u_fractal;
uniform vec2  u_fractalCenter;
uniform float u_maskHalfExtent;
uniform vec4  u_farFieldColor;
uniform float u_rimDarken;
uniform vec3  u_eyeDir;
out vec4 fragColor;
void main() {
  vec3 n = normalize(v_norm);
  float dz = 1.0 - n.z;
  // North-pole cap — inverse stereographic is singular here (w → ∞).
  if (dz < 1e-4) { fragColor = u_farFieldColor; return; }
  float invDZ = 1.0 / dz;
  // Inverse stereographic: w = (x/(1−z), y/(1−z)).
  float u = n.x * invDZ;
  float v = n.y * invDZ;
  // Map w ∈ [cx−R, cx+R]² → [0,1]² texture coordinates.
  vec2 uv = (vec2(u, v) - u_fractalCenter) / (2.0 * u_maskHalfExtent) + 0.5;
  // Out-of-range UV samples the edge texel via CLAMP_TO_EDGE — those regions
  // correspond to w values outside the fractal bounding square, which are
  // already colored as "escaped" by the Schwarz iteration.
  fragColor = texture(u_fractal, uv);
  // Rim shading: darken pixels where the surface normal is nearly perpendicular
  // to the view direction (silhouette edges of the sphere).
  float facing = max(0.0, dot(n, u_eyeDir));
  fragColor.rgb *= 1.0 - u_rimDarken * (1.0 - facing * facing);
}`;

  // Overlay geometry: all line-drawn elements (boundary polyline, pole markers,
  // north-pole ✸). Positions are pushed outward 0.3% from the sphere surface
  // to avoid z-fighting without needing polygon-offset for gl.LINES.
  const OVERLAY_VERT = `#version 300 es
in vec3 a_pos;
uniform mat4 u_mvp;
void main() {
  gl_Position = u_mvp * vec4(a_pos * 1.003, 1.0);
}`;

  const OVERLAY_FRAG = `#version 300 es
precision mediump float;
uniform vec4 u_color;
out vec4 fragColor;
void main() { fragColor = u_color; }`;

  // ---------------------------------------------------------------------------
  // createSphereRenderer  — public factory
  // ---------------------------------------------------------------------------
  function createSphereRenderer(canvas) {
    // Acquire WebGL 2 context.
    let gl;
    try { gl = canvas.getContext('webgl2', { antialias: true, preserveDrawingBuffer: false }); }
    catch (e) { gl = null; }
    if (!gl) return null;

    // Pull shared helpers from schwarz-webgl.js (must be loaded first).
    const SW = QD.Schwarz;
    if (!SW || !SW._shaders || !SW._glHelpers || !SW._gpuCaps) {
      console.error('sphere-webgl: QD.Schwarz internal API not available — load schwarz-webgl.js first.');
      return null;
    }
    const { MAX_BRANCHES, MAX_K, MAX_LAURENT, MAX_BETA } = SW._gpuCaps;
    const { vert: FRACTAL_VERT, frag: FRACTAL_FRAG } = SW._shaders;
    const H = SW._glHelpers;  // {compile, link, buildMaskTexture, buildColormapTexture,
                               //  pickColormap, SCALE_MODE_ID, jsRHashAtInfinity}

    // Compile all three programs.
    let fractalProg, sphereProg, overlayProg;
    let fracVS, fracFS, sphVS, sphFS, ovVS, ovFS;
    try {
      fracVS = H.compile(gl, gl.VERTEX_SHADER,   FRACTAL_VERT);
      fracFS = H.compile(gl, gl.FRAGMENT_SHADER, FRACTAL_FRAG);
      fractalProg = H.link(gl, fracVS, fracFS);

      sphVS = H.compile(gl, gl.VERTEX_SHADER,   SPHERE_VERT);
      sphFS = H.compile(gl, gl.FRAGMENT_SHADER, SPHERE_FRAG);
      sphereProg = H.link(gl, sphVS, sphFS);

      ovVS = H.compile(gl, gl.VERTEX_SHADER,   OVERLAY_VERT);
      ovFS = H.compile(gl, gl.FRAGMENT_SHADER, OVERLAY_FRAG);
      overlayProg = H.link(gl, ovVS, ovFS);
    } catch (e) {
      console.error('sphere-webgl: shader compile failed:', e);
      return null;
    }

    // ---- Uniform locations (fractal program — mirrors schwarz-webgl.js) ----
    const UF = {};
    [
      'viewCenter','pxPerUnit','canvasSize',
      'unbounded','family','w0','gamma','z0','absZ0','rInfConj',
      'c','polyA','polyALen','branchZ','branchA','branchACount','nBranches',
      'maxIter','escapeR',
      'mask','maskCenter','maskHalfExtent',
      'colormap','scaleMode','modK',
      // Polynomial-h β-correction (HANDOFF #22 / #26 / #28): unbounded LQDs.
      'lqdBeta','lqdBetaLen',
    ].forEach(n => { UF[n] = gl.getUniformLocation(fractalProg, 'u_' + n); });

    // ---- Uniform locations (sphere program) ----
    const US = {};
    ['mvp','fractal','fractalCenter','maskHalfExtent','farFieldColor','rimDarken','eyeDir']
      .forEach(n => { US[n] = gl.getUniformLocation(sphereProg, 'u_' + n); });

    // ---- Uniform locations (overlay program) ----
    const UO = {};
    ['mvp','color'].forEach(n => { UO[n] = gl.getUniformLocation(overlayProg, 'u_' + n); });

    // ---- Full-screen triangle VBO (for fractal pass) ----
    const fracVBO = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, fracVBO);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1,-1, 3,-1, -1, 3]), gl.STATIC_DRAW);
    const fracAPos = gl.getAttribLocation(fractalProg, 'a_pos');

    // ---- Sphere mesh VBO + IBO ----
    const SC = (typeof SphereCommon !== 'undefined') ? SphereCommon
      : (typeof window !== 'undefined' && window.SphereCommon ? window.SphereCommon : null);
    if (!SC) { console.error('sphere-webgl: SphereCommon not loaded.'); return null; }

    const mesh = SC.buildSphereMesh(96, 48);
    const sphVBO = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, sphVBO);
    gl.bufferData(gl.ARRAY_BUFFER, mesh.positions, gl.STATIC_DRAW);
    const sphIBO = gl.createBuffer();
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, sphIBO);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, mesh.indices, gl.STATIC_DRAW);
    const sphAPos = gl.getAttribLocation(sphereProg, 'a_pos');
    const ovAPos  = gl.getAttribLocation(overlayProg, 'a_pos');

    // ---- Overlay VBO (dynamic — rebuilt in setPhi) ----
    const ovVBO = gl.createBuffer();
    let ovLineStart  = 0;  // boundary start (in vertex index)
    let ovLineCount  = 0;  // boundary segment vertices (2 per segment)
    let ovMarkStart  = 0;  // pole + north-pole marker start
    let ovMarkCount  = 0;  // marker segment vertices

    // ---- Phi state (matches phiState in schwarz-webgl.js) ----
    const phiState = {
      unbounded:    false,
      familyId:     0,
      w0:           new Float32Array(2),
      c:            0,
      gamma:        new Float32Array(2),
      z0:           new Float32Array(2),
      absZ0:        0,
      rInfConj:     new Float32Array(2),
      polyA:        new Float32Array(MAX_LAURENT * 2),
      polyALen:     0,
      branchZ:      new Float32Array(MAX_BRANCHES * 2),
      branchA:      new Float32Array(MAX_BRANCHES * MAX_K * 2),
      branchACount: new Int32Array(MAX_BRANCHES),
      nBranches:    0,
      lqdBeta:      new Float32Array(MAX_BETA * 2),
      lqdBetaLen:   0,
      mask:         null,
      maskCenter:   new Float32Array(2),
      maskHalfExtent: new Float32Array(2),
      escapeR:      1e10,
      // Fractal view extent: covers the mask region (stored scalar).
      fractalR:     1,
      fractalCx:    0,
      fractalCy:    0,
    };

    // ---- Render params ----
    const rp = {
      maxIter:   64,
      scaleMode: 'smooth',
      modK:      8,
      texSize:   1024,
    };

    // ---- Display params ----
    const dp = {
      rimDarken:      0.3,
      showBoundary:   true,
      showPoles:      true,
      showNorthPole:  true,
      boundaryColor:  [1.0, 1.0, 1.0, 1.0],
      poleColor:      [1.0, 0.9, 0.3, 1.0],
    };

    // Far-field color: same "escape" gray as the Schwarz renderer.
    const FAR_FIELD_COLOR = [80/255, 80/255, 90/255, 1.0];

    // ---- Colormap + FBO textures ----
    let colormapTex  = H.buildColormapTexture(gl, H.pickColormap('magma'));
    let activeColormap = 'magma';

    let fractalFBO = null;     // { fbo, tex, texSize }
    let fractalDirty = true;   // re-render fractal on next render() call

    // 2×2 default gray texture used when no φ has been captured yet. The sphere
    // shader still does its inverse-stereographic + texture sample, but the
    // sample returns a uniform mid-gray, giving a flat-shaded "placeholder
    // sphere" that the user can orbit/zoom to confirm 3D interaction.
    const defaultFractalTex = (function () {
      const tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      const px = new Uint8Array([
        180,180,184,255,  180,180,184,255,
        180,180,184,255,  180,180,184,255,
      ]);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 2, 2, 0, gl.RGBA, gl.UNSIGNED_BYTE, px);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      return tex;
    })();
    let hasPhi = false;   // becomes true after a successful setPhi()

    // =========================================================================
    // setPhi — pack phi data + build mask + boundary overlay geometry
    // =========================================================================
    function setPhi(phi, opts) {
      opts = opts || {};
      const polyPts = opts.boundaryPts || [];

      // γ-branch merge for unboundedLQD_singular (HANDOFF #24 / #26 / #28):
      // append {z: phi.z0, A: phi.lqdGamma} to the uploaded branches so the
      // shader's branch loop picks up the synthetic-branch contribution.
      // Mirrors withSyntheticBranch in schwarz-common.js and
      // _phiWithSyntheticBranch in solver-uqd-lqd-singular.js.
      const gammaArr = phi.lqdGamma || [];
      const effBranches =
        (phi.family === 'unboundedLQD_singular' && gammaArr.length > 0 && phi.z0)
          ? (phi.branches || []).concat([{ z: phi.z0, A: gammaArr }])
          : (phi.branches || []);

      // Validate capacity (same as schwarz-webgl.js).
      const nb    = effBranches.length;
      const polyA = phi.polyA || phi.F || [];
      const beta  = phi.lqdBeta || [];
      if (nb > MAX_BRANCHES || polyA.length > MAX_LAURENT) return false;
      if (beta.length > MAX_BETA) return false;
      for (let j = 0; j < nb; j++) {
        if (effBranches[j].A.length > MAX_K) return false;
      }

      // Pack phi state.
      phiState.unbounded = !!phi.unbounded;
      switch (phi.family) {
        case 'boundedLQD':            phiState.familyId = 2; break;
        case 'boundedLQD_singular':   phiState.familyId = 3; break;
        case 'unboundedLQD':          phiState.familyId = 4; break;
        case 'unboundedLQD_singular': phiState.familyId = 5; break;
        default: phiState.familyId = phi.unbounded ? 1 : 0;
      }
      phiState.w0[0]    = phi.w0    ? phi.w0.re    : 0;
      phiState.w0[1]    = phi.w0    ? phi.w0.im    : 0;
      phiState.gamma[0] = phi.gamma ? phi.gamma.re : 0;
      phiState.gamma[1] = phi.gamma ? phi.gamma.im : 0;
      phiState.z0[0]    = phi.z0    ? phi.z0.re    : 0;
      phiState.z0[1]    = phi.z0    ? phi.z0.im    : 0;
      phiState.absZ0    = phi.z0    ? Math.hypot(phi.z0.re, phi.z0.im) : 0;
      phiState.c        = phi.c != null ? phi.c : 0;
      phiState.polyALen = polyA.length;
      phiState.polyA.fill(0);
      for (let l = 0; l < polyA.length; l++) {
        phiState.polyA[2*l]   = polyA[l].re;
        phiState.polyA[2*l+1] = polyA[l].im;
      }
      // Polynomial-h β-correction (HANDOFF #22 / #26 / #28): copy β into the
      // GPU-side Float32Array. The shader sees B(1/z) = Σ_l β_l/z^l in φ's
      // exp argument for families 4 (unboundedLQD) and 5 (unboundedLQD_singular).
      phiState.lqdBeta.fill(0);
      phiState.lqdBetaLen = beta.length;
      for (let l = 0; l < beta.length; l++) {
        phiState.lqdBeta[2*l]   = beta[l].re;
        phiState.lqdBeta[2*l+1] = beta[l].im;
      }
      // rInfConj must reflect the MERGED branches so the shader's
      // `sum - cconj(u_rInfConj)` term sees the correct ∞-baseline including
      // the γ-branch's contribution at infinity.
      if (phiState.familyId === 4 || phiState.familyId === 5) {
        const rInf = H.jsRHashAtInfinity({ branches: effBranches });
        phiState.rInfConj[0] =  rInf.re;
        phiState.rInfConj[1] = -rInf.im;
      } else {
        phiState.rInfConj[0] = phiState.rInfConj[1] = 0;
      }
      phiState.nBranches = nb;
      phiState.branchZ.fill(0); phiState.branchA.fill(0); phiState.branchACount.fill(0);
      for (let j = 0; j < nb; j++) {
        const br = effBranches[j];
        phiState.branchZ[2*j]   = br.z.re;
        phiState.branchZ[2*j+1] = br.z.im;
        phiState.branchACount[j] = br.A.length;
        for (let k = 0; k < br.A.length; k++) {
          phiState.branchA[2*(j*MAX_K+k)]   = br.A[k].re;
          phiState.branchA[2*(j*MAX_K+k)+1] = br.A[k].im;
        }
      }

      // Build polygon mask texture.
      const padFactor = phiState.unbounded ? 5.0 : 2.4;
      if (phiState.mask) gl.deleteTexture(phiState.mask);
      if (polyPts.length) {
        const m = H.buildMaskTexture(gl, polyPts, padFactor);
        phiState.mask             = m.tex;
        phiState.maskCenter[0]    = m.maskCenter.re;
        phiState.maskCenter[1]    = m.maskCenter.im;
        phiState.maskHalfExtent[0] = m.maskHalfExtent.x;
        phiState.maskHalfExtent[1] = m.maskHalfExtent.y;
      } else {
        phiState.mask = null;
        phiState.maskCenter[0] = phiState.maskCenter[1] = 0;
        phiState.maskHalfExtent[0] = phiState.maskHalfExtent[1] = 1;
      }

      phiState.escapeR = opts.escapeR ||
        (phiState.unbounded ? phiState.maskHalfExtent[0] * 6.0 : 1e10);

      // The fractal texture covers [cx−R, cx+R]² in w-space.
      // R is the mask half-extent (already includes the pad factor).
      phiState.fractalR  = phiState.maskHalfExtent[0];
      phiState.fractalCx = phiState.maskCenter[0];
      phiState.fractalCy = phiState.maskCenter[1];

      // Build overlay geometry (boundary + poles + north-pole marker).
      _buildOverlayGeometry(phi, polyPts);

      fractalDirty = true;
      hasPhi = true;
      return true;
    }

    // =========================================================================
    // setRenderParams
    // =========================================================================
    function setRenderParams(params) {
      params = params || {};
      let dirty = false;
      if (params.maxIter   != null && params.maxIter   !== rp.maxIter)   { rp.maxIter   = params.maxIter;   dirty = true; }
      if (params.scaleMode != null && params.scaleMode !== rp.scaleMode) { rp.scaleMode = params.scaleMode; dirty = true; }
      if (params.modK      != null && params.modK      !== rp.modK)      { rp.modK      = params.modK;      dirty = true; }
      if (params.colormap  != null && params.colormap  !== activeColormap) {
        activeColormap = params.colormap;
        gl.deleteTexture(colormapTex);
        colormapTex = H.buildColormapTexture(gl, H.pickColormap(activeColormap));
        dirty = true;
      }
      if (params.texSize != null && params.texSize !== rp.texSize) {
        rp.texSize = params.texSize;
        dirty = true;
      }
      if (dirty) fractalDirty = true;
    }

    function setDisplayParams(params) {
      params = params || {};
      if (params.rimDarken     != null) dp.rimDarken     = params.rimDarken;
      if (params.showBoundary  != null) dp.showBoundary  = params.showBoundary;
      if (params.showPoles     != null) dp.showPoles     = params.showPoles;
      if (params.showNorthPole != null) dp.showNorthPole = params.showNorthPole;
      if (params.boundaryColor != null) dp.boundaryColor = params.boundaryColor;
      if (params.poleColor     != null) dp.poleColor     = params.poleColor;
    }

    // =========================================================================
    // render  — main entry point called by sphere-ui.js each frame
    // =========================================================================
    // camera = { azimuth, elevation, distance }  (azimuth/elevation in radians)
    // size   = { W, H }  in physical pixels (canvas.width / canvas.height)
    function render(camera, size) {
      const W = size ? size.W : canvas.clientWidth;
      const H = size ? size.H : canvas.clientHeight;
      if (canvas.width !== W || canvas.height !== H) {
        canvas.width = W; canvas.height = H;
      }

      // --- Compute eye position and MVP matrices ----------------------------
      const cosEl = Math.cos(camera.elevation);
      const sinEl = Math.sin(camera.elevation);
      const cosAz = Math.cos(camera.azimuth);
      const sinAz = Math.sin(camera.azimuth);
      const dist  = camera.distance;
      const eye   = [dist*cosEl*cosAz, dist*cosEl*sinAz, dist*sinEl];

      const viewMat = SC.mat4lookAt(eye, [0,0,0], [0,0,1]);
      const projMat = SC.mat4perspective(Math.PI / 3, W / Math.max(H, 1), 0.05, 20);
      const mvp     = SC.mat4multiply(projMat, viewMat);

      // Normalized eye direction (for rim shading in sphere fragment shader).
      const eyeLen  = Math.sqrt(eye[0]*eye[0] + eye[1]*eye[1] + eye[2]*eye[2]);
      const eyeDir  = [eye[0]/eyeLen, eye[1]/eyeLen, eye[2]/eyeLen];

      // --- Pass 1: Fractal to FBO (cached) ----------------------------------
      // Skipped entirely when no φ has been captured — the sphere shader will
      // fall back to the 2×2 default gray texture so the user still sees a
      // flat-shaded placeholder sphere they can orbit.
      if (hasPhi && fractalDirty) {
        _renderFractalToFBO(rp.texSize);
        fractalDirty = false;
      }

      // --- Pass 2 + 3: Sphere + overlay to canvas ---------------------------
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, W, H);
      gl.clearColor(0.08, 0.08, 0.12, 1.0);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

      // Pass 2: textured sphere.
      gl.enable(gl.DEPTH_TEST);
      gl.depthFunc(gl.LESS);
      gl.enable(gl.CULL_FACE);
      gl.cullFace(gl.BACK);

      gl.useProgram(sphereProg);
      gl.uniformMatrix4fv(US.mvp, false, mvp);
      gl.uniform3fv(US.eyeDir, eyeDir);
      gl.uniform2fv(US.fractalCenter, [phiState.fractalCx, phiState.fractalCy]);
      gl.uniform1f(US.maskHalfExtent, phiState.fractalR);
      gl.uniform4fv(US.farFieldColor,  FAR_FIELD_COLOR);
      gl.uniform1f(US.rimDarken, dp.rimDarken);

      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D,
        (hasPhi && fractalFBO) ? fractalFBO.tex : defaultFractalTex);
      gl.uniform1i(US.fractal, 0);

      gl.bindBuffer(gl.ARRAY_BUFFER, sphVBO);
      gl.enableVertexAttribArray(sphAPos);
      gl.vertexAttribPointer(sphAPos, 3, gl.FLOAT, false, 0, 0);
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, sphIBO);
      gl.drawElements(gl.TRIANGLES, mesh.nTris * 3, gl.UNSIGNED_SHORT, 0);

      // Pass 3: overlay geometry.
      const anyOverlay = (dp.showBoundary && ovLineCount > 0) ||
                         ((dp.showPoles || dp.showNorthPole) && ovMarkCount > 0);
      if (anyOverlay) {
        gl.disable(gl.CULL_FACE);
        gl.depthFunc(gl.LEQUAL);

        gl.useProgram(overlayProg);
        gl.uniformMatrix4fv(UO.mvp, false, mvp);

        gl.bindBuffer(gl.ARRAY_BUFFER, ovVBO);
        gl.enableVertexAttribArray(ovAPos);
        gl.vertexAttribPointer(ovAPos, 3, gl.FLOAT, false, 0, 0);

        // Boundary polyline (white or user-chosen color).
        if (dp.showBoundary && ovLineCount > 0) {
          gl.uniform4fv(UO.color, dp.boundaryColor);
          gl.drawArrays(gl.LINES, ovLineStart, ovLineCount);
        }
        // Pole markers (yellow or user-chosen color).
        if ((dp.showPoles || dp.showNorthPole) && ovMarkCount > 0) {
          gl.uniform4fv(UO.color, dp.poleColor);
          gl.drawArrays(gl.LINES, ovMarkStart, ovMarkCount);
        }
      }
    }

    // =========================================================================
    // Pass 1: render fractal to FBO
    // =========================================================================
    function _renderFractalToFBO(texSize) {
      // Create or resize FBO.
      if (!fractalFBO || fractalFBO.texSize !== texSize) {
        if (fractalFBO) {
          gl.deleteFramebuffer(fractalFBO.fbo);
          gl.deleteTexture(fractalFBO.tex);
        }
        fractalFBO = _createFBO(texSize);
        if (!fractalFBO) { console.error('sphere-webgl: FBO creation failed'); return; }
      }

      gl.bindFramebuffer(gl.FRAMEBUFFER, fractalFBO.fbo);
      gl.viewport(0, 0, texSize, texSize);
      gl.useProgram(fractalProg);

      // View: a square of [cx−R, cx+R]² covering the boundary region.
      const cx = phiState.fractalCx, cy = phiState.fractalCy;
      const R  = phiState.fractalR;
      const pxPerUnit = texSize / (2 * Math.max(R, 1e-9));

      gl.uniform2f(UF.viewCenter, cx, cy);
      gl.uniform1f(UF.pxPerUnit, pxPerUnit);
      gl.uniform2f(UF.canvasSize, texSize, texSize);

      // Phi state uniforms.
      gl.uniform1i(UF.unbounded, phiState.unbounded ? 1 : 0);
      gl.uniform1i(UF.family,    phiState.familyId);
      gl.uniform2fv(UF.w0,        phiState.w0);
      gl.uniform2fv(UF.gamma,     phiState.gamma);
      gl.uniform2fv(UF.z0,        phiState.z0);
      gl.uniform1f(UF.absZ0,      phiState.absZ0);
      gl.uniform2fv(UF.rInfConj,  phiState.rInfConj);
      gl.uniform1f(UF.c,          phiState.c);
      gl.uniform2fv(UF.polyA,     phiState.polyA);
      gl.uniform1i(UF.polyALen,   phiState.polyALen);
      gl.uniform2fv(UF.lqdBeta,   phiState.lqdBeta);
      gl.uniform1i(UF.lqdBetaLen, phiState.lqdBetaLen);
      gl.uniform2fv(UF.branchZ,   phiState.branchZ);
      gl.uniform2fv(UF.branchA,   phiState.branchA);
      gl.uniform1iv(UF.branchACount, phiState.branchACount);
      gl.uniform1i(UF.nBranches,  phiState.nBranches);
      gl.uniform1i(UF.maxIter,    rp.maxIter);
      gl.uniform1f(UF.escapeR,    phiState.escapeR);
      gl.uniform1i(UF.scaleMode,  H.SCALE_MODE_ID[rp.scaleMode] | 0);
      gl.uniform1i(UF.modK,       Math.max(2, rp.modK | 0));

      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, phiState.mask);
      gl.uniform1i(UF.mask, 0);
      gl.uniform2fv(UF.maskCenter,     phiState.maskCenter);
      gl.uniform2fv(UF.maskHalfExtent, phiState.maskHalfExtent);

      gl.activeTexture(gl.TEXTURE1);
      gl.bindTexture(gl.TEXTURE_2D, colormapTex);
      gl.uniform1i(UF.colormap, 1);

      // Full-screen triangle.
      gl.bindBuffer(gl.ARRAY_BUFFER, fracVBO);
      gl.enableVertexAttribArray(fracAPos);
      gl.vertexAttribPointer(fracAPos, 2, gl.FLOAT, false, 0, 0);
      gl.disable(gl.DEPTH_TEST);
      gl.disable(gl.CULL_FACE);
      gl.drawArrays(gl.TRIANGLES, 0, 3);

      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    }

    function _createFBO(texSize) {
      const tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, texSize, texSize, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      const fbo = gl.createFramebuffer();
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
      const ok = gl.checkFramebufferStatus(gl.FRAMEBUFFER) === gl.FRAMEBUFFER_COMPLETE;
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      if (!ok) { gl.deleteTexture(tex); gl.deleteFramebuffer(fbo); return null; }
      return { fbo, tex, texSize };
    }

    // =========================================================================
    // Build overlay geometry: boundary polyline + pole/north-pole markers
    // =========================================================================
    function _buildOverlayGeometry(phi, boundaryPts) {
      // --- Boundary line segments ---
      const boundVerts = [];
      if (boundaryPts && boundaryPts.length >= 2) {
        const pts = boundaryPts;
        const N   = pts.length;
        // Line strip: consecutive pairs.
        for (let i = 0; i < N; i++) {
          const pA = pts[i];
          const pB = pts[(i + 1) % N];
          const sA = SC.projectToSphere({ re: pA.re, im: pA.im });
          const sB = SC.projectToSphere({ re: pB.re, im: pB.im });
          boundVerts.push(sA.x, sA.y, sA.z, sB.x, sB.y, sB.z);
        }
      }

      // --- Pole markers (star shape ✸ on the sphere surface) ---
      const markVerts = [];
      const STAR_SIZE  = 0.03;   // arm length (world units)
      const STAR_N_POL = 3;      // North-pole star
      const STAR_N_POL_SIZE = 0.05;

      // Finite poles from phi.branches (the a_j = φ(z_j) values).
      // Project z_j → a_j using evalPhi if available; otherwise use the
      // pole positions from the Schwarz handle (hData poles).
      // Simplest: draw stars at the quadrature nodes (a_j values from phi.branches)
      // projected to sphere. These are in w-space; project w → sphere.
      // Note: phi.branches store the z_j (preimage in 𝔻), not the a_j (pole in w-space).
      // The a_j lives in hData. We receive phi but not hData here.
      // Fallback: project the boundary centroid — no, that's wrong.
      // Instead, use the quadrature data's pole positions stored in the
      // phi object (if available). Classical QDs store poles in the inverse-Faber
      // form; for the overlay we just mark the center of the boundary polygon
      // as an approximation. The proper approach is to accept poles as a param.
      // See the sphereState.polesSnapshot stored in sphere-ui.js.
      // For now we'll accept poles via a separate setPolePts() call from sphere-ui.js.

      // North-pole marker (always, if poles are shown).
      const npStar = _buildStar([0, 0, 1], STAR_N_POL_SIZE, STAR_N_POL);
      markVerts.push(...npStar);

      // Build VBO data: boundary first, markers second.
      ovLineStart = 0;
      ovLineCount = boundVerts.length / 3;  // each entry is one coord; vertices = length/3
      ovMarkStart = ovLineCount;
      ovMarkCount = markVerts.length / 3;

      const all = new Float32Array(boundVerts.length + markVerts.length);
      all.set(boundVerts, 0);
      all.set(markVerts, boundVerts.length);

      gl.bindBuffer(gl.ARRAY_BUFFER, ovVBO);
      gl.bufferData(gl.ARRAY_BUFFER, all, gl.DYNAMIC_DRAW);
    }

    // setPolePts: called from sphere-ui.js with the actual quadrature pole
    // positions in w-space so we can draw markers at the correct locations.
    function setPolePts(poles) {
      if (!poles || !poles.length) return;
      const STAR_SIZE  = 0.03;
      const STAR_N_ARM = 3;
      const STAR_N_POL_SIZE = 0.05;

      // Re-fetch existing boundary geometry from VBO isn't practical.
      // We'll store an "extra mark verts" array and re-upload the whole thing.
      // This requires knowing the current boundary verts too.
      // Simpler: call _buildOverlayGeometry_withPoles — but we'd need to store
      // boundaryPts. Add a stateful cache.
      // For now, build pole mark geometry and append to existing ovVBO.
      // Since this is called once per capture, just rebuild entirely.
      // Store the pole positions for the next rebuild.
      _storedPolePts = poles;
      _rebuildMarkVerts();
    }

    let _storedPolePts = [];
    let _storedBoundaryPts = [];

    function _rebuildMarkVerts() {
      const STAR_SIZE      = 0.03;
      const STAR_N_ARM     = 3;
      const STAR_POLE_SIZE = 0.05;
      const markVerts = [];

      // Finite poles.
      for (const p of _storedPolePts) {
        const sp = SC.projectToSphere({ re: p.re, im: p.im });
        const star = _buildStar([sp.x, sp.y, sp.z], STAR_SIZE, STAR_N_ARM);
        markVerts.push(...star);
      }
      // North-pole ✸.
      const npStar = _buildStar([0, 0, 1], STAR_POLE_SIZE, STAR_N_ARM);
      markVerts.push(...npStar);

      // Re-build boundary section.
      const boundVerts = [];
      const pts = _storedBoundaryPts;
      const N   = pts.length;
      for (let i = 0; i < N; i++) {
        const pA = pts[i], pB = pts[(i + 1) % N];
        const sA = SC.projectToSphere({ re: pA.re, im: pA.im });
        const sB = SC.projectToSphere({ re: pB.re, im: pB.im });
        boundVerts.push(sA.x, sA.y, sA.z, sB.x, sB.y, sB.z);
      }

      ovLineStart = 0;
      ovLineCount = boundVerts.length / 3;
      ovMarkStart = ovLineCount;
      ovMarkCount = markVerts.length / 3;

      const all = new Float32Array(boundVerts.length + markVerts.length);
      all.set(boundVerts, 0);
      all.set(markVerts, boundVerts.length);
      gl.bindBuffer(gl.ARRAY_BUFFER, ovVBO);
      gl.bufferData(gl.ARRAY_BUFFER, all, gl.DYNAMIC_DRAW);
    }

    // Build a ✸ star pattern around a point on the sphere.
    // center: [x, y, z] unit vector; size: arm length; nArms: number of arms.
    // Returns flat [x,y,z, x,y,z, ...] line-pair vertices (2 verts per arm).
    function _buildStar(center, size, nArms) {
      const cx = center[0], cy = center[1], cz = center[2];
      // Tangent basis at center using Gram-Schmidt.
      let t1x, t1y, t1z;
      if (Math.abs(cy) > 0.9) { t1x=1; t1y=0; t1z=0; } // default
      else                     { t1x=0; t1y=1; t1z=0; }
      const d1 = t1x*cx + t1y*cy + t1z*cz;
      t1x -= d1*cx; t1y -= d1*cy; t1z -= d1*cz;
      const l1 = Math.sqrt(t1x*t1x + t1y*t1y + t1z*t1z);
      t1x /= l1; t1y /= l1; t1z /= l1;
      // t2 = cross(center, t1)
      const t2x = cy*t1z - cz*t1y;
      const t2y = cz*t1x - cx*t1z;
      const t2z = cx*t1y - cy*t1x;

      const verts = [];
      for (let i = 0; i < nArms; i++) {
        const angle = Math.PI * i / nArms;
        const cos   = Math.cos(angle), sin = Math.sin(angle);
        const dx = cos*t1x + sin*t2x;
        const dy = cos*t1y + sin*t2y;
        const dz = cos*t1z + sin*t2z;
        verts.push(
          cx - size*dx, cy - size*dy, cz - size*dz,
          cx + size*dx, cy + size*dy, cz + size*dz,
        );
      }
      return verts;
    }

    // =========================================================================
    // destroy
    // =========================================================================
    function destroy() {
      if (phiState.mask)        gl.deleteTexture(phiState.mask);
      if (colormapTex)          gl.deleteTexture(colormapTex);
      if (defaultFractalTex)    gl.deleteTexture(defaultFractalTex);
      if (fractalFBO) {
        gl.deleteFramebuffer(fractalFBO.fbo);
        gl.deleteTexture(fractalFBO.tex);
      }
      gl.deleteBuffer(fracVBO); gl.deleteBuffer(sphVBO);
      gl.deleteBuffer(sphIBO);  gl.deleteBuffer(ovVBO);
      gl.deleteProgram(fractalProg);
      gl.deleteProgram(sphereProg);
      gl.deleteProgram(overlayProg);
      [fracVS, fracFS, sphVS, sphFS, ovVS, ovFS].forEach(s => { if (s) gl.deleteShader(s); });
    }

    // =========================================================================
    // Public helper: store boundary pts reference so _rebuildMarkVerts works.
    // Called internally from setPhi; also called by sphere-ui when adding poles.
    // =========================================================================
    const rendererPublic = {
      available: true,
      setPhi(phi, opts) {
        _storedBoundaryPts = (opts && opts.boundaryPts) ? opts.boundaryPts : [];
        return setPhi(phi, opts);
      },
      setPolePts,
      setRenderParams,
      setDisplayParams,
      render,
      destroy,
      markFractalDirty() { fractalDirty = true; },
    };
    return rendererPublic;
  }

  Sphere.createRenderer = createSphereRenderer;

})(typeof globalThis !== 'undefined' ? globalThis : this);
