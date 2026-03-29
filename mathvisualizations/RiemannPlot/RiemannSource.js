var plot_preset_dict = {
  "disk":       {f:"z", zoom: 1, center: [0,0]}
};


var iniscript = function(preset,res) {
  return `
  use("CindyGL");

  // Fixed constants
  e = 2.71828182845904523536028747; // Euler's number

  reim(z) := [re(z),im(z)];

  lwZeroApprox(z) := ( // Approximation of the principal branch of the Lambert W function near the origin
  ezsqrt = sqrt(1 + e*z);
      (12*ezsqrt*(45*sqrt(2) + 32*ezsqrt))/(sqrt(e)*(623 + 83*e*z + 372*sqrt(2)*ezsqrt))-1;
  );

  lwInftyApprox(z) := log(z)-log(log(z))+log(log(z))/log(z); // Approximation of the principal branch of the Lambert W function near infinity

  lambertw(z) :=( // Custom implementation of the Lambert W function
      if(abs(z)<1.7,
          w=lwZeroApprox(z),
          w=lwInftyApprox(z)
      );
      repeat(5,
          w=(w^2+z/exp(w))/(w+1);
      );
      w
  );

  arg(z) := arctan2(reim(z));


  // initial values of parameters
  zoom = ${preset.zoom};
  center = [${preset.center}];

  f(z) := (
      ${preset.f};
  );

  // Coordinate transformations
  PltToCanvX(x) := (x-center_1)*zoom+1; // plot coordinate to canvas coordinate
  PltToCanvY(y) := (y-center_2)*zoom+1; // plot coordinate to canvas coordinate
  PltToCanvXY(XY) := [PltToCanvX(XY_1),PltToCanvY(XY_2)];
  PltToCanvZ(z) := PltToCanvX(re(z))+i*PltToCanvY(im(z)); // plot coordinate to canvas coordinate

  CanvToPltX(x) := (x-1)/zoom+center_1; // canvas coordinate to plot coordinate
  CanvToPltY(y) := (y-1)/zoom+center_2; // canvas coordinate to plot coordinate
  CanvToPltXY(XY) := [CanvToPltX(XY_1),CanvToPltY(XY_2)];
  CanvToPltZ(z) := CanvToPltX(re(z))+i*CanvToPltY(im(z)); // canvas coordinate to plot coordinate

  //Colors
  Z0.color  = (1,1,1);
  colorFcn(u) := (
      if(u==n,(0,0,0),
          u = u/n;
          u = (3*u/(2*u+1));
          (4*u,1.3*u,(1-u)^2*.7);
      );
  );
  `
}

function addarrays(a,b){ // add a pair of arrays elementwise
    return a.map((e,i) => e + b[i]);
}

function subtractarrays(a,b){ // subtract a pair of arrays elementwise
    return a.map((e,i) => e - b[i]);
}

function arraymult(arr,m){ // multiply an array by a scalar
    return arr.map((e,i) => e*m);
}


class RiemannPlot {
  constructor(varName, paramDict, canvasName, canvasID, callbacks = {}, canvasWidth = 500, canvasHeight = 500, res = 500) {
    this._varName = varName; // object instance variable name (string)
    this._canvasName = canvasName; // name of javascript canvas (string)
    this._canvasID = canvasID; // canvas ID (string)
    this._canvasWidth = canvasWidth; // canvas width in pixels (int)
    this._canvasHeight = canvasHeight; // canvas height in pixels (int)
    this._res = res; // resolution (res x res) (int)
    this._callbacks = callbacks;
    this._f = paramDict.f; // mapping function (string) f(z)="this._f" 
    this._center = paramDict.center; // center of plot (array)
    this._zoom = paramDict.zoom; // default zoom level
    this._mousepos = [0,0]; // current position of mouse (in canvas coordinates)
    this._mouseshift = [0,0]; //shift vector for change in mouse position (in canvas coordinates)
    this._isPtSelected = false;
    this._movescript = ''; // list of lines of JS code (strings) to be executed when any element of the canvas moves
    this._keydownscript = ''; // list of lines of JS code (strings) to be executed when a key is pressed
    this._mousemovescript = ''; // list of lines of JS code (strings) to be executed when the mouse is moved
    this._mouseclickscript = ''; // list of lines of JS code (strings) to be executed when the mouse is clicked
    this._mousedownscript = ''; // list of lines of JS code (strings) to be executed when the mouse is clicked down
    this._mousedragscript = ''; // list of lines of JS code (strings) to be executed repeatedly when the mouse is being dragged
    this._movescript = 'colorplot([center_1-1/zoom,center_2-1/zoom],[center_1+1/zoom,center_2-1/zoom],"julia",colorFcn(dynIter(complex(#), c)));' +
        'drawimage([0,0],[2,0], "julia");' +
        'connect(apply(dynKIter(CanvToPltZ(complex(Z0.xy)),c,nplot-1),reim(PltToCanvZ(#))),color->[1,1,1],size->1.8);' +
        'drawtext(Z0+(.025,.025), "z0="+CanvToPltZ(complex(Z0.xy)), color->[1,1,1],size->15);';
    this._keydownscript   = 'javascript("' + this._varName + '.keypress(\'"+"\\" + key()+"\'.charCodeAt(0))");';
    this._mousedragscript = `javascript("
      if (!(${this._varName}.z0AtMouse)) {
        ${this._varName}.shift(arraymult(${this._varName}.mouseshift,1/${this._varName}.zoom));
      }");
      `;
    this._mousedownscript = `javascript("${this._varName}.mouseshift === null");`; // fixes timing issue with click and drag feature
    if ("move" in this._callbacks) {
      this._movescript += GenCindyJSCode(this._callbacks.move);
    }
    if ("keydown" in this._callbacks) {
      this._keydownscript += GenCindyJSCode(this._callbacks.keydown);
    }
    if ("mousedrag" in this._callbacks) {
      this._mousedragscript += GenCindyJSCode(this._callbacks.mousedrag);
    }
    if ("mousedown" in this._callbacks) {
      this._mousedownscript += GenCindyJSCode(this._callbacks.mousedown);
    }
    if ("mousemove" in this._callbacks) {
      this._mousemovescript += GenCindyJSCode(this._callbacks.mousemove);
    }
    if ("mouseclick" in this._callbacks) {
      this._mouseclickscript += GenCindyJSCode(this._callbacks.mouseclick);
    }
    this._cindy = CindyJS({
      canvasname: this._canvasName,
      scripts: {
        init:       iniscript(paramDict,this._res),
        move:       this._movescript,
        keydown:    this._keydownscript,
        mousedrag:  this._mousedragscript,
        mousedown:  this._mousedownscript,
        mousemove:  this._mousemovescript,
        mouseclick: this._mouseclickscript
      },
      //geometry: [{name: "Z0", kind: "P", type: "Free", pos: this.PlotToCanv(this._z0), size: 3 }],
      ports: [{
        id: this._canvasID,
        width: this._canvasWidth,
        height: this._canvasHeight,
        transform: [{ visibleRect: [0,2,2,0] }],},]
    });
  }

  evokeCS(cscode) {
    this._cindy.evokeCS(cscode);
  }

  ApplyPreset(preset) {
    this.center = preset.center;
    this.zoom = preset.zoom;
    this.f = preset.f;
  }

  CanvToPlot(z) {
    return [(z[0]-1)/this._zoom+this._center[0], (z[1]-1)/this._zoom+this._center[1]];
  }

  PlotToCanv(z) {
    return [(z[0]-this._center[0])*this._zoom+1, (z[1]-this._center[1])*this._zoom+1];
  }

  zoomIn(ratio) {
    this.zoom = this.zoom*ratio;
  }

  shift(vec) {
    this.center = addarrays(this._center,vec);
  }


  keypress(key) {
    switch(key) {
      case 187:  this.zoomIn(2); break;
      case 189:  this.zoomIn(1/2); break;
      case 38:   this.shift([0,1/(this._zoom*4)]); break;
      case 40:   this.shift([0,-1/(this._zoom*4)]); break;
      case 39:   this.shift([1/(this._zoom*4),0]); break;
      case 37:   this.shift([-1/(this._zoom*4),0]); break;
    }
  }

  exportImage(imageName) { // Export image as PNG
    this.cindy.exportPNG(imageName);
  }
    

  set isPtSelected(isSelected) {
    this._isPtSelected = isSelected;
  }


  set f(fval) {
    this._f = fval;
    this.evokeCS(`f(z,c) := (${this._f});`);
  }

  set zoom(zoomval) {
    this._zoom = zoomval;
    this.evokeCS(`zoom=${zoomval};`);
  }

  set center(centerval) {
    this._center = centerval;
    var str = `center_1=${this._center[0]};
               center_2=${this._center[1]};`;
    this.evokeCS(str);
  }

  set res(resVal) {
    this._res = resVal;
    this.evokeCS(`createimage("julia", ${this._res}, ${this._res})`);
  }

  get mousepos() {
    this.evokeCS(`javascript("${this._varName}._mousepos = "+mouse());`);
    return this._mousepos;
  }

  get mouseshift() {
    this._mouseshift = subtractarrays(this._mousepos,this.mousepos);
    return this._mouseshift;
  }

  get zoom() {
    return this._zoom;
  }

  get center() {
    return this._center;
  }

  get f() {
    return this._f;
  }

  get cindy() {
    return this._cindy;
  }

  get range() {
    return [this._center[0]-1/this._zoom,this._center[0]+1/this._zoom,
            this._center[1]-1/this._zoom,this._center[1]+1/this._zoom];
  }

  get res() {
    return this._res;
  }

}


/////////////////////

// Generate string that calls the lines of JS code
// takes in an array of strings where each is a line
// of JS code (w/o semicolons)
var GenCindyJSCode = function(arr) {
  return  'javascript("' + arr.join("; ") +';");';
}




////////////////////

var complex = function(zArr) {
  if (zArr[1] >= 0) {
    return zArr[0] + '+i*' + zArr[1];
  } else {
    return zArr[0] + '-i*' + (-zArr[1]);
  }
}

var re = function(z) {
  out = z.split(/\-|\+/);
  if (z[0] == '-') {
    out = (-1)*parseFloat(out[1]);
  } else {
    out = parseFloat(out[0]);
  }
  return out;
}

var im = function(z) {
  if (z.includes('-i*') || z.includes('+i*')) {
    imag = z.split(/\-i\*|\+i\*/)[1];
    out = parseFloat(z[z.length-imag.length-3]+'1')*parseFloat(imag);
  } else {
    imag = z.split(/\-|\+/).at(-1);
    out = parseFloat(z[z.length-imag.length-1]+'1')*parseFloat(imag.slice(0,imag.length-1));
  }
  return out;
}

var reim = function(z) {
  return [re(z),im(z)];
}












// Javascript-HTML interfacing functions.
var getFInput = function() {
  return document.getElementById("inpf").value;
}

var getCenterInput = function() {
  return document.getElementById("inpcenter").value.split(",").map(parseFloat);
}

var getZoomInput = function() {
  return parseFloat(document.getElementById("inpzoom").value);
}

//////////
/*
var setNInput = function(nval) {
  document.getElementById("inpmn").value = nval;
  document.getElementById("inpjn").value = nval;
}

var setFInput = function(fval) {
  document.getElementById("inpf").value = fval;
}

var setDEscInput = function(escval) {
  document.getElementById("inpje").value = escval;
}

var setPEscInput = function(escval) {
  document.getElementById("inpme").value = escval;
}

var setPCenterInput = function(centerval) {
  centerval = centerval.map(vals => vals.toPrecision(6)).map(parseFloat);
  document.getElementById("inpparamcenter").value = centerval;
}

var setPZoomInput = function(zoomval) {
  document.getElementById("inpparamzoom").value = zoomval;
}

var setDCenterInput = function(centerval) {
  centerval = centerval.map(vals => vals.toPrecision(6)).map(parseFloat);
  document.getElementById("inpdyncenter").value = centerval;
}

var setDZoomInput = function(zoomval) {
  document.getElementById("inpdynzoom").value = zoomval;
}


var getPresetDicts = function() { 
  return [{f:getFInput(), c:getCInput(),                      n:getPNInput(), nplot:getPNPlotInput(), escape:getPEscInput(), zoom: getPZoomInput(), center: getPCenterInput()},
          {f:getFInput(), c:julia_fract.c, z0:julia_fract.z0, n:getDNInput(), nplot:getDNPlotInput(), escape:getDEscInput(), zoom: getDZoomInput(), center: getDCenterInput()}];
}

var setInputs = function(preset_val) {
  setPCenterInput(param_preset_dict[preset_val].center);
  setPZoomInput(param_preset_dict[preset_val].zoom);
  setDCenterInput(dyn_preset_dict[preset_val].center);
  setDZoomInput(dyn_preset_dict[preset_val].zoom);
  setCInput(param_preset_dict[preset_val].c);
  setNInput(param_preset_dict[preset_val].n);
  setFInput(param_preset_dict[preset_val].f);
  setPEscInput(param_preset_dict[preset_val].escape);
  setDEscInput(dyn_preset_dict[preset_val].escape);
}

var apply_preset = function(preset_val) {
  setInputs(preset_val);
  julia_fract.ApplyPreset(    dyn_preset_dict[preset_val]);
  parameter_fract.ApplyPreset(param_preset_dict[preset_val]);
}

var apply_changes = function() {
  julia_fract.ApplyPreset(    getPresetDicts()[1]);
  parameter_fract.ApplyPreset(getPresetDicts()[0]);
  parameter_fract.res = document.getElementById('inpdres').value;
  julia_fract.res = document.getElementById('inpjres').value;
}
*/



// Utility Javascript functions

// Save canvas to image file
// Possible filetypes include: png,jpeg,svg,pdf
var saveCanvasAs = function(canvas,filename,filetype) { 
  var downloadLink = document.createElement('a');
  downloadLink.download = filename + "." + filetype;
  downloadLink.href = canvas.toDataURL("image/" + filetype);
  downloadLink.click();
  downloadLink.remove();
  return true;
}