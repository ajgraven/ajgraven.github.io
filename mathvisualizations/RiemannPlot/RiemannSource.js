// Utility functions shared across Riemann visualization tools.
// The main conformal map logic lives inline in riemann_plot.html.

function addarrays(a, b)      { return a.map((e, i) => e + b[i]); }
function subtractarrays(a, b) { return a.map((e, i) => e - b[i]); }
function arraymult(arr, m)    { return arr.map(e => e * m); }

// Preset dictionary: each entry has f (CindyScript expression), center, zoom
var presets = {
  "z^2":         { f: "z^2",          center: [0, 0], zoom: 0.5  },
  "z^3":         { f: "z^3",          center: [0, 0], zoom: 0.5  },
  "1/z":         { f: "1/z",          center: [0, 0], zoom: 0.5  },
  "1/(z^2+1)":   { f: "1/(z^2+1)",   center: [0, 0], zoom: 0.5  },
  "sin(z)":      { f: "sin(z)",       center: [0, 0], zoom: 0.35 },
  "cos(z)":      { f: "cos(z)",       center: [0, 0], zoom: 0.35 },
  "tan(z)":      { f: "tan(z)",       center: [0, 0], zoom: 0.35 },
  "exp(z)":      { f: "exp(z)",       center: [0, 0], zoom: 0.15 },
  "log(z)":      { f: "log(z)",       center: [0, 0], zoom: 0.35 },
  "sqrt(z)":     { f: "sqrt(z)",      center: [0, 0], zoom: 0.5  },
  "z+1/z":       { f: "z+1/z",        center: [0, 0], zoom: 0.35 },
  "(z-1)/(z+1)": { f: "(z-1)/(z+1)", center: [0, 0], zoom: 0.5  },
};



var iniscript = function(preset, res) {
  return `
  use("CindyGL");

  zoom   = ${preset.zoom};
  center = [${preset.center[0]}, ${preset.center[1]}];

  f(z) := (${preset.f});

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
  //   hue     = arg(w) / (2 pi)
  //   value   = modulated logarithmically by |w|, creating rings around zeros/poles
  domainColor(w) := (
    regional(n, logw, zfract, grey1, grey2);
    n = 12;
    logw = log(w) / (2 * pi);
    zfract = n * logw - floor(n * logw);
    grey1 = im(zfract);
    grey2 = re(zfract);
    hsvToRGB(im(logw), 1.0, 0.5 + 0.5 * re(sqrt(grey1 * grey2)))
  );

  createimage("domain", ${res}, ${res});
  `;
};


class DomainPlot {
  constructor(varName, preset, canvasName, canvasID, canvasWidth = 600, canvasHeight = 600, res = 600) {
    this._varName      = varName;
    this._canvasName   = canvasName;
    this._canvasID     = canvasID;
    this._canvasWidth  = canvasWidth;
    this._canvasHeight = canvasHeight;
    this._res          = res;
    this._f            = preset.f;
    this._center       = preset.center.slice();
    this._zoom         = preset.zoom;
    this._mousepos     = [0, 0];

    var movescript = [
      'colorplot(',
      '  [center_1 - 1/zoom, center_2 - 1/zoom],',
      '  [center_1 + 1/zoom, center_2 - 1/zoom],',
      '  "domain",',
      '  domainColor(f(complex(#)))',
      ');',
      'drawimage([0,0],[2,0],"domain");',
    ].join('\n');

    var keydownscript   = 'javascript("' + varName + '.keypress(\'"+"\\" + key()+"\'.charCodeAt(0))");';
    var mousedragscript = 'javascript("' + varName + '.shift(arraymult(' + varName + '.mouseshift,1/' + varName + '.zoom));");';
    var mousedownscript = 'javascript("' + varName + '.mouseshift === null;");';

    this._cindy = CindyJS({
      canvasname: canvasName,
      scripts: {
        init:      iniscript(preset, res),
        move:      movescript,
        keydown:   keydownscript,
        mousedrag: mousedragscript,
        mousedown: mousedownscript,
      },
      geometry: [],
      ports: [{
        id:        canvasID,
        width:     canvasWidth,
        height:    canvasHeight,
        transform: [{ visibleRect: [0, 2, 2, 0] }],
      }],
    });
  }

  evokeCS(code) { this._cindy.evokeCS(code); }

  ApplyPreset(preset) {
    this.f      = preset.f;
    this.center = preset.center.slice();
    this.zoom   = preset.zoom;
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

  exportImage(imageName) { this._cindy.exportPNG(imageName); }

  get mousepos() {
    this.evokeCS(`javascript("${this._varName}._mousepos = "+mouse());`);
    return this._mousepos;
  }

  get mouseshift() {
    this._mouseshift = subtractarrays(this._mousepos, this.mousepos);
    return this._mouseshift;
  }

  set f(fval) {
    this._f = fval;
    this.evokeCS(`f(z) := (${fval});`);
  }

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
    this.evokeCS(`createimage("domain", ${this._res}, ${this._res});`);
  }

  get zoom()   { return this._zoom;   }
  get center() { return this._center; }
  get f()      { return this._f;      }
  get res()    { return this._res;    }
}


// ── UI helpers ───────────────────────────────────────────────────────────────

var getFInput      = () => document.getElementById("inpf").value;
var getCenterInput = () => document.getElementById("inpcenter").value.split(",").map(parseFloat);
var getZoomInput   = () => parseFloat(document.getElementById("inpzoom").value);
var getResInput    = () => parseInt(document.getElementById("inpres").value);

var setFInput = (v) => { document.getElementById("inpf").value = v; };
var setCenterInput = (v) => {
  document.getElementById("inpcenter").value =
    v.map(x => parseFloat(x.toPrecision(6))).join(",");
};
var setZoomInput = (v) => { document.getElementById("inpzoom").value = v; };

var apply_preset = function(presetKey) {
  var p = presets[presetKey];
  setFInput(p.f);
  setCenterInput(p.center);
  setZoomInput(p.zoom);
  plot.ApplyPreset(p);
};

var apply_changes = function() {
  plot.f      = getFInput();
  plot.center = getCenterInput();
  plot.zoom   = getZoomInput();
  plot.res    = getResInput();
};
