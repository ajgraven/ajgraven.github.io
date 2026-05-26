// =============================================================================
// sphere-common.js  — Pure math kernel for the Riemann-sphere visualization.
//
// NO DOM. NO WebGL. Safe to load in node-test.js.
//
// Exports (window.SphereCommon in browser; module.exports.SphereCommon in Node):
//
//   Stereographic projection (north-pole convention: ∞ ↔ north pole (0,0,1)):
//     projectToSphere(w)          {re,im} → {x,y,z}
//     unprojectFromSphere(p,eps)  {x,y,z} → {re,im} | null  (null near north pole)
//
//   UV-sphere mesh:
//     buildSphereMesh(nLon, nLat) → { positions:Float32Array(nVerts×3),
//                                     indices:Uint16Array(nTris×3),
//                                     nVerts, nTris, nLon, nLat }
//       Orientation: j=0 ↔ (0,0,+1) (north pole = ∞);
//                    j=nLat ↔ (0,0,−1) (south pole = origin in C).
//       Fits in Uint16Array for nLon ≤ 255, nLat ≤ 254.
//
//   4×4 matrix math (column-major, WebGL convention):
//     mat4identity()
//     mat4multiply(a, b)
//     mat4lookAt(eye, target, up)      — right-handed view matrix
//     mat4perspective(fovY, aspect, near, far)
//     mat4invertRigid(m)               — fast inverse for rigid-body transforms
//
// =============================================================================

(function (global) {
  'use strict';

  // ===========================================================================
  // Stereographic projection   (north-pole projection)
  // ===========================================================================
  // Convention: project from the north pole so that ∞ maps to (0, 0, +1).
  //
  // Forward (C → S²):
  //   w = u + iv,  r² = u² + v²
  //   x = 2u/(1 + r²),   y = 2v/(1 + r²),   z = (r²−1)/(r²+1)
  //
  //   ∞    → (0, 0, +1)   (north pole)
  //   0    → (0, 0, −1)   (south pole)
  //   |w|=1 → z=0          (equator)
  //
  // Inverse (S² \ {north pole} → C):
  //   u = x/(1−z),   v = y/(1−z)
  //   Singular when z = 1 (north pole). Returns null for |1−z| < eps.

  function projectToSphere(w) {
    const u = w.re, v = w.im;
    const r2 = u * u + v * v;
    const d  = 1.0 + r2;
    return { x: 2 * u / d, y: 2 * v / d, z: (r2 - 1) / (r2 + 1) };
  }

  function unprojectFromSphere(p, eps) {
    if (eps == null) eps = 1e-9;
    const denom = 1.0 - p.z;
    if (Math.abs(denom) < eps) return null;
    if (p.z > 0) {
      // When z is close to 1 (large |w|), computing x/(1-z) suffers catastrophic
      // cancellation.  Use the algebraically equivalent form
      //   u = x*(1+z) / (x²+y²)  (since x²+y² = 1−z² = (1−z)(1+z))
      // which avoids subtracting two nearly-equal numbers.
      const r2sq = p.x * p.x + p.y * p.y;   // = 1 − z²  on the unit sphere
      if (r2sq === 0) return { re: 0, im: 0 };
      const fac = (1.0 + p.z) / r2sq;
      return { re: p.x * fac, im: p.y * fac };
    }
    return { re: p.x / denom, im: p.y / denom };
  }

  // ===========================================================================
  // UV-sphere mesh
  // ===========================================================================
  // Returns vertices as flat Float32Array of (x,y,z) triples and a Uint16Array
  // index buffer for triangle rendering.  All vertex positions lie on the unit
  // sphere (|v| = 1 exactly by construction via sin/cos).
  //
  // The seam at longitude=0 and longitude=2π has duplicate vertices so UV
  // parameterization is consistent (u=0 and u=1 produce the same 3D position
  // but are stored separately).

  function buildSphereMesh(nLon, nLat) {
    nLon = nLon || 96;
    nLat = nLat || 48;

    const nVerts  = (nLon + 1) * (nLat + 1);
    const nTris   = nLon * nLat * 2;
    const positions = new Float32Array(nVerts * 3);
    const indices   = new Uint16Array(nTris * 3);

    // Vertex positions
    let vi = 0;
    for (let j = 0; j <= nLat; j++) {
      const theta = Math.PI * j / nLat;   // 0 (north) … π (south)
      const sinT  = Math.sin(theta);
      const cosT  = Math.cos(theta);
      for (let i = 0; i <= nLon; i++) {
        const phi  = 2 * Math.PI * i / nLon;
        positions[vi++] = sinT * Math.cos(phi);   // x
        positions[vi++] = sinT * Math.sin(phi);   // y
        positions[vi++] = cosT;                    // z
      }
    }

    // Triangle indices (CCW from the outside)
    let ii = 0;
    for (let j = 0; j < nLat; j++) {
      for (let i = 0; i < nLon; i++) {
        const a = j * (nLon + 1) + i;
        const b = a + 1;
        const c = a + (nLon + 1);
        const d = c + 1;
        indices[ii++] = a; indices[ii++] = c; indices[ii++] = b;
        indices[ii++] = b; indices[ii++] = c; indices[ii++] = d;
      }
    }

    return { positions, indices, nVerts, nTris, nLon, nLat };
  }

  // ===========================================================================
  // 4×4 matrix math  (column-major, gl.uniformMatrix4fv-compatible)
  // ===========================================================================
  // Storage: m[col * 4 + row].  Equivalently, the first 4 floats are column 0,
  // the next 4 are column 1, etc.

  function mat4identity() {
    return new Float64Array([1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]);
  }

  // C = A * B  (both column-major)
  function mat4multiply(a, b) {
    const out = new Float64Array(16);
    for (let col = 0; col < 4; col++) {
      for (let row = 0; row < 4; row++) {
        let s = 0;
        for (let k = 0; k < 4; k++) s += a[k * 4 + row] * b[col * 4 + k];
        out[col * 4 + row] = s;
      }
    }
    return out;
  }

  // Right-handed view matrix: eye, target, up are [x,y,z] arrays.
  // forward = normalize(target - eye)
  // right   = normalize(cross(forward, up))
  // vup     = cross(right, forward)   (re-orthogonalised)
  //
  // Returns a Float64Array (full double precision).  gl.uniformMatrix4fv
  // accepts any TypedArray and converts to float32 internally, so this is
  // safe to pass directly to WebGL.
  function mat4lookAt(eye, target, up) {
    const f = _normalize(_sub(target, eye));
    const s = _normalize(_cross(f, up));
    const u = _cross(s, f);
    const m = new Float64Array(16);
    m[ 0] =  s[0]; m[ 4] =  s[1]; m[ 8] =  s[2]; m[12] = -_dot(s, eye);
    m[ 1] =  u[0]; m[ 5] =  u[1]; m[ 9] =  u[2]; m[13] = -_dot(u, eye);
    m[ 2] = -f[0]; m[ 6] = -f[1]; m[10] = -f[2]; m[14] =  _dot(f, eye);
    m[ 3] =     0; m[ 7] =     0; m[11] =     0; m[15] =  1;
    return m;
  }

  // Standard symmetric perspective frustum.
  // fovY in radians; aspect = width / height.
  // Returns a Float64Array; see mat4lookAt note above.
  function mat4perspective(fovY, aspect, near, far) {
    const f        = 1.0 / Math.tan(fovY * 0.5);
    const rangeInv = 1.0 / (near - far);
    const m = new Float64Array(16);
    m[ 0] = f / aspect;
    m[ 5] = f;
    m[10] = (near + far) * rangeInv;
    m[11] = -1;
    m[14] = 2 * near * far * rangeInv;
    return m;
  }

  // Fast inverse for a rigid-body (rotation + translation) matrix.
  // M = [R | Rt; 0 | 1]  →  M⁻¹ = [R^T | −R^T·Rt; 0 | 1]
  // Works for any matrix returned by mat4lookAt.
  function mat4invertRigid(m) {
    const out = new Float64Array(16);
    // Transpose the 3×3 rotation block.
    out[ 0] = m[0]; out[ 4] = m[1]; out[ 8] = m[2];
    out[ 1] = m[4]; out[ 5] = m[5]; out[ 9] = m[6];
    out[ 2] = m[8]; out[ 6] = m[9]; out[10] = m[10];
    // Translation: −R^T · t.
    const tx = m[12], ty = m[13], tz = m[14];
    out[12] = -(out[0]*tx + out[4]*ty + out[ 8]*tz);
    out[13] = -(out[1]*tx + out[5]*ty + out[ 9]*tz);
    out[14] = -(out[2]*tx + out[6]*ty + out[10]*tz);
    out[15] = 1;
    return out;
  }

  // ---- Internal 3-vector helpers (not exported) ----------------------------
  function _sub(a, b) { return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]; }
  function _dot(a, b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
  function _cross(a, b) {
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
  }
  function _normalize(a) {
    const len = Math.sqrt(_dot(a, a));
    if (len < 1e-300) return [0, 0, 1];
    return [a[0]/len, a[1]/len, a[2]/len];
  }

  // ===========================================================================
  // Export
  // ===========================================================================
  const SphereCommon = {
    projectToSphere,
    unprojectFromSphere,
    buildSphereMesh,
    mat4identity,
    mat4multiply,
    mat4lookAt,
    mat4perspective,
    mat4invertRigid,
  };

  if (typeof window !== 'undefined') {
    window.SphereCommon = SphereCommon;
  } else if (typeof module !== 'undefined') {
    module.exports = { SphereCommon };
  }

})(typeof globalThis !== 'undefined' ? globalThis : this);
