// Preset dictionary for the Logarithmic Potential Explorer.
// Each entry has:
//   pPos    – 3 canvas [x,y] positions for P1, P2, P3 (positive charges / zeros)
//   nPos    – 3 canvas [x,y] positions for N1, N2, N3 (negative charges / poles)
//   pActive – 3 booleans, whether each positive charge is active
//   nActive – 3 booleans, whether each negative charge is active
//   zoom    – initial zoom level
//   center  – initial plot-space center [cx, cy]
//
// Canvas coords for default zoom=0.5, center=[0,0]:
//   canvas(x,y) = (plot_x * 0.5 + 1, plot_y * 0.5 + 1)
var pt_presets = {
  "single_pos": {
    pPos:    [[1.0, 1.0], [-0.5, 1.0], [-0.5, 0.5]],
    nPos:    [[2.5, 0.5], [2.5, 1.0], [2.5, 1.5]],
    pActive: [true, false, false],
    nActive: [false, false, false],
    zoom:    0.5,
    center:  [0, 0],
  },
  "dipole": {
    pPos:    [[0.75, 1.0], [-0.5, 1.0], [-0.5, 0.5]],
    nPos:    [[1.25, 1.0], [2.5, 1.0], [2.5, 1.5]],
    pActive: [true, false, false],
    nActive: [true, false, false],
    zoom:    0.5,
    center:  [0, 0],
  },
  "two_pos": {
    pPos:    [[0.7, 1.0], [1.3, 1.0], [-0.5, 0.5]],
    nPos:    [[2.5, 0.5], [2.5, 1.0], [2.5, 1.5]],
    pActive: [true, true, false],
    nActive: [false, false, false],
    zoom:    0.5,
    center:  [0, 0],
  },
  "quadrupole": {
    pPos:    [[0.75, 1.25], [1.25, 0.75], [-0.5, 0.5]],
    nPos:    [[1.25, 1.25], [0.75, 0.75], [2.5, 1.5]],
    pActive: [true, true, false],
    nActive: [true, true, false],
    zoom:    0.5,
    center:  [0, 0],
  },
  "full_hex": {
    pPos:    [[1.0, 1.35], [0.65, 0.825], [1.35, 0.825]],
    nPos:    [[0.65, 1.175], [1.35, 1.175], [1.0, 0.65]],
    pActive: [true, true, true],
    nActive: [true, true, true],
    zoom:    0.5,
    center:  [0, 0],
  },
};


function addarrays(a, b)      { return a.map((e, i) => e + b[i]); }
function subtractarrays(a, b) { return a.map((e, i) => e - b[i]); }
function arraymult(arr, m)    { return arr.map(e => e * m); }


var iniscript = function(preset, res) {
  return `
  use("CindyGL");

  zoom   = ${preset.zoom};
  center = [${preset.center[0]}, ${preset.center[1]}];

  p1_on = ${preset.pActive[0]};
  p2_on = ${preset.pActive[1]};
  p3_on = ${preset.pActive[2]};
  n1_on = ${preset.nActive[0]};
  n2_on = ${preset.nActive[1]};
  n3_on = ${preset.nActive[2]};

  // Convert HSV (h,s,v each in [0,1]) to an RGB triple
  hsvToRGB(h, s, v) := (
    regional(j, fr, p, q, t);
    h = (h - floor(h)) * 6;
    j = floor(h);
    fr = h - j;
    p = 1.0 - s;
    q = 1.0 - s * fr;
    t = 1.0 - s * (1.0 - fr);
    if(j == 0, [1, t, p],
    if(j == 1, [q, 1, p],
    if(j == 2, [p, 1, t],
    if(j == 3, [p, q, 1],
    if(j == 4, [t, p, 1],
    [1, p, q]))))) * v
  );

  // Standard domain coloring:
  //   hue   = arg(w) / (2 pi)
  //   value = modulated logarithmically by |w|, creating rings around zeros/poles
  domainColor(w) := (
    regional(n, logw, zfract, grey1, grey2);
    n = 12;
    logw = log(w) / (2 * pi);
    zfract = n * logw - floor(n * logw);
    grey1 = im(zfract);
    grey2 = re(zfract);
    hsvToRGB(im(logw), 1.0, 0.5 + 0.5 * re(sqrt(grey1 * grey2)))
  );

  createimage("potential", ${res}, ${res});
  `;
};


// Move script: executed on every frame.
// Steps before colorplot run on CPU and set GLSL uniforms.
// colorplot runs on GPU using those uniforms.
var movescript = [
  // Convert P1,P2,P3 canvas coords -> plot coords
  'z1 = (re(complex(P1.xy))-1)/zoom + center_1 + i*((im(complex(P1.xy))-1)/zoom + center_2);',
  'z2 = (re(complex(P2.xy))-1)/zoom + center_1 + i*((im(complex(P2.xy))-1)/zoom + center_2);',
  'z3 = (re(complex(P3.xy))-1)/zoom + center_1 + i*((im(complex(P3.xy))-1)/zoom + center_2);',
  // Convert N1,N2,N3 canvas coords -> plot coords
  'n1 = (re(complex(N1.xy))-1)/zoom + center_1 + i*((im(complex(N1.xy))-1)/zoom + center_2);',
  'n2 = (re(complex(N2.xy))-1)/zoom + center_1 + i*((im(complex(N2.xy))-1)/zoom + center_2);',
  'n3 = (re(complex(N3.xy))-1)/zoom + center_1 + i*((im(complex(N3.xy))-1)/zoom + center_2);',
  // GPU domain-color pass
  'colorplot(',
  '  [center_1 - 1/zoom, center_2 - 1/zoom],',
  '  [center_1 + 1/zoom, center_2 - 1/zoom],',
  '  "potential",',
  '  domainColor(',
  '    if(p1_on, complex(#) - z1, 1) *',
  '    if(p2_on, complex(#) - z2, 1) *',
  '    if(p3_on, complex(#) - z3, 1) /',
  '    if(n1_on, complex(#) - n1, 1) /',
  '    if(n2_on, complex(#) - n2, 1) /',
  '    if(n3_on, complex(#) - n3, 1)',
  '  )',
  ');',
  'drawimage([0,0],[2,0],"potential");',
  // Draw labels over draggable charge points
  'if(p1_on, drawtext(P1+(0.04,0.04),"+",color->(0.9,0,0),size->14));',
  'if(p2_on, drawtext(P2+(0.04,0.04),"+",color->(0.9,0,0),size->14));',
  'if(p3_on, drawtext(P3+(0.04,0.04),"+",color->(0.9,0,0),size->14));',
  'if(n1_on, drawtext(N1+(0.04,0.04),"-",color->(0,0,0.9),size->14));',
  'if(n2_on, drawtext(N2+(0.04,0.04),"-",color->(0,0,0.9),size->14));',
  'if(n3_on, drawtext(N3+(0.04,0.04),"-",color->(0,0,0.9),size->14));',
].join('\n');


class PotentialPlot {
  constructor(varName, preset, canvasName, canvasID, canvasWidth = 600, canvasHeight = 600, res = 600) {
    this._varName      = varName;
    this._canvasName   = canvasName;
    this._canvasID     = canvasID;
    this._canvasWidth  = canvasWidth;
    this._canvasHeight = canvasHeight;
    this._res          = res;
    this._center       = preset.center.slice();
    this._zoom         = preset.zoom;
    this._pActive      = preset.pActive.slice();
    this._nActive      = preset.nActive.slice();
    this._mousepos     = [0, 0];
    this._anyChargeAtMouse = false;

    var vn = varName;

    var keydownscript = 'javascript("' + vn + '.keypress(\'"+"\\" + key()+"\'.charCodeAt(0))");';

    // Only pan when no charge point is under the mouse
    var mousedragscript =
      'javascript("if (!(' + vn + '.anyChargeAtMouse)) {' +
        vn + '.shift(arraymult(' + vn + '.mouseshift,1/' + vn + '.zoom));' +
      '}");';

    var mousedownscript = 'javascript("' + vn + '.mouseshift === null;");';

    this._cindy = CindyJS({
      canvasname: canvasName,
      scripts: {
        init:      iniscript(preset, res),
        move:      movescript,
        keydown:   keydownscript,
        mousedrag: mousedragscript,
        mousedown: mousedownscript,
      },
      geometry: [
        { name: "P1", type: "Free", pos: preset.pPos[0], color: [0.9, 0, 0], pinned: false, size: 5 },
        { name: "P2", type: "Free", pos: preset.pPos[1], color: [0.9, 0, 0], pinned: false, size: 5 },
        { name: "P3", type: "Free", pos: preset.pPos[2], color: [0.9, 0, 0], pinned: false, size: 5 },
        { name: "N1", type: "Free", pos: preset.nPos[0], color: [0, 0, 0.9], pinned: false, size: 5 },
        { name: "N2", type: "Free", pos: preset.nPos[1], color: [0, 0, 0.9], pinned: false, size: 5 },
        { name: "N3", type: "Free", pos: preset.nPos[2], color: [0, 0, 0.9], pinned: false, size: 5 },
      ],
      ports: [{
        id:        canvasID,
        width:     canvasWidth,
        height:    canvasHeight,
        transform: [{ visibleRect: [0, 2, 2, 0] }],
      }],
    });
  }

  evokeCS(code) { this._cindy.evokeCS(code); }

  // Returns true if any charge point is currently under the mouse cursor.
  // Result is cached via a side-channel write from CindyScript.
  get anyChargeAtMouse() {
    this.evokeCS(
      'javascript("' + this._varName + '._anyChargeAtMouse = " + ' +
      '(contains(elementsatmouse(),P1)|contains(elementsatmouse(),P2)|contains(elementsatmouse(),P3)|' +
      'contains(elementsatmouse(),N1)|contains(elementsatmouse(),N2)|contains(elementsatmouse(),N3)) + ";");'
    );
    return this._anyChargeAtMouse;
  }

  get mousepos() {
    this.evokeCS(`javascript("${this._varName}._mousepos = "+mouse());`);
    return this._mousepos;
  }

  get mouseshift() {
    this._mouseshift = subtractarrays(this._mousepos, this.mousepos);
    return this._mouseshift;
  }

  zoomIn(ratio) { this.zoom = this._zoom * ratio; }

  shift(vec) { this.center = addarrays(this._center, vec); }

  keypress(key) {
    switch (key) {
      case 187: this.zoomIn(2);                        break; // +
      case 189: this.zoomIn(0.5);                      break; // -
      case 38:  this.shift([0,  1/(this._zoom * 4)]);  break; // up
      case 40:  this.shift([0, -1/(this._zoom * 4)]);  break; // down
      case 39:  this.shift([ 1/(this._zoom * 4), 0]);  break; // right
      case 37:  this.shift([-1/(this._zoom * 4), 0]);  break; // left
    }
  }

  // Toggle active state for one charge. isPos=true for positive (P), false for negative (N).
  // idx is 0-based. Returns the new boolean state.
  toggleCharge(idx, isPos) {
    if (isPos) {
      this._pActive[idx] = !this._pActive[idx];
      var flag = 'p' + (idx + 1) + '_on';
      this.evokeCS(flag + ' = ' + this._pActive[idx] + ';');
      return this._pActive[idx];
    } else {
      this._nActive[idx] = !this._nActive[idx];
      var flag = 'n' + (idx + 1) + '_on';
      this.evokeCS(flag + ' = ' + this._nActive[idx] + ';');
      return this._nActive[idx];
    }
  }

  applyPreset(preset) {
    this._pActive = preset.pActive.slice();
    this._nActive = preset.nActive.slice();

    // Move geometry points
    this.evokeCS('P1.xy=[' + preset.pPos[0][0] + ',' + preset.pPos[0][1] + '];');
    this.evokeCS('P2.xy=[' + preset.pPos[1][0] + ',' + preset.pPos[1][1] + '];');
    this.evokeCS('P3.xy=[' + preset.pPos[2][0] + ',' + preset.pPos[2][1] + '];');
    this.evokeCS('N1.xy=[' + preset.nPos[0][0] + ',' + preset.nPos[0][1] + '];');
    this.evokeCS('N2.xy=[' + preset.nPos[1][0] + ',' + preset.nPos[1][1] + '];');
    this.evokeCS('N3.xy=[' + preset.nPos[2][0] + ',' + preset.nPos[2][1] + '];');

    // Update active flags
    this.evokeCS(
      'p1_on = ' + preset.pActive[0] + '; ' +
      'p2_on = ' + preset.pActive[1] + '; ' +
      'p3_on = ' + preset.pActive[2] + '; ' +
      'n1_on = ' + preset.nActive[0] + '; ' +
      'n2_on = ' + preset.nActive[1] + '; ' +
      'n3_on = ' + preset.nActive[2] + ';'
    );

    // Update view
    this.zoom   = preset.zoom;
    this.center = preset.center.slice();
  }

  exportImage(imageName) { this._cindy.exportPNG(imageName); }

  set zoom(zval) {
    this._zoom = zval;
    this.evokeCS(`zoom = ${zval};`);
  }

  set center(cval) {
    this._center = cval;
    this.evokeCS(`center_1 = ${cval[0]}; center_2 = ${cval[1]};`);
  }

  set res(rval) {
    this._res = parseInt(rval);
    this.evokeCS(`createimage("potential", ${this._res}, ${this._res});`);
  }

  get zoom()    { return this._zoom;    }
  get center()  { return this._center;  }
  get res()     { return this._res;     }
  get pActive() { return this._pActive; }
  get nActive() { return this._nActive; }
}


// ── UI helpers ───────────────────────────────────────────────────────────────

var getCenterInput = () => document.getElementById("inpcenter").value.split(",").map(parseFloat);
var getZoomInput   = () => parseFloat(document.getElementById("inpzoom").value);
var getResInput    = () => parseInt(document.getElementById("inpres").value);

var setCenterInput = (v) => {
  document.getElementById("inpcenter").value =
    v.map(x => parseFloat(x.toPrecision(6))).join(",");
};
var setZoomInput = (v) => { document.getElementById("inpzoom").value = v; };

var apply_changes = function() {
  plot.center = getCenterInput();
  plot.zoom   = getZoomInput();
  plot.res    = getResInput();
};

var apply_preset = function(presetKey) {
  var p = pt_presets[presetKey];
  setCenterInput(p.center);
  setZoomInput(p.zoom);
  plot.applyPreset(p);

  // Update positive charge button styles
  for (var i = 0; i < 3; i++) {
    var btn = document.getElementById('btn-p' + (i + 1));
    if (btn) btn.style.backgroundColor = p.pActive[i] ? '#cc3333' : '#888888';
  }
  // Update negative charge button styles
  for (var i = 0; i < 3; i++) {
    var btn = document.getElementById('btn-n' + (i + 1));
    if (btn) btn.style.backgroundColor = p.nActive[i] ? '#3333cc' : '#888888';
  }
};

// chargeId: 'p1', 'p2', 'p3', 'n1', 'n2', 'n3'
var toggle = function(chargeId) {
  var isPos = chargeId[0] === 'p';
  var idx   = parseInt(chargeId[1]) - 1;
  var active = plot.toggleCharge(idx, isPos);
  var btn = document.getElementById('btn-' + chargeId);
  if (btn) {
    btn.style.backgroundColor = active
      ? (isPos ? '#cc3333' : '#3333cc')
      : '#888888';
  }
};
