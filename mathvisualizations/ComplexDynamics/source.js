var param_preset_dict = {
  "mandelbrot":     {f:"z^2+c", c:'-.7-.4*i', n:'100', escape:"abs(z)>2", zoom: .75, center: [-.75,0]},
  "tricorn":        {f:"conjugate(z^2)+c",  c:'.1091+i*.502',    n:'50', escape:"abs(z)>2", zoom: .55, center: [-.25,0]},
  "burning ship":   {f:"(abs(re(z))+i*abs(im(z)))^2-c", c:'1.6185+i*0.0471', n:'100', escape:"abs(z)>2", center:[1.7,.05], zoom: 10},
  "butterfly":      {f:"conjugate(z^2)+c*re(1/z)", c:'.252+i*0', n:'150', escape:"abs(z)>4", center:[0,0], zoom:.7},
  "exponential map":{f:"c*e^(z-1)", c:'1+0*i', n:'25', escape:'re(log(c)+z)>3000', center:[3.5,0], zoom: .08},
  "exp Schwarz":    {f:"c0 = c^2/z;c1 = lambertw(-c0);conjugate(c0/exp(c1+c^2/c1));", c:'1+0*i', n:'50', escape:"abs(lambertw(-c^2/f(z,c)))>abs(c)",center: [1,0], zoom: .5}
};
var dyn_preset_dict = {
  "mandelbrot":     {f:"z^2+c", c:'.2541-.0333*i', n:'100', escape:"abs(z)>2", zoom: .65, center: [0,0]},
  "tricorn":        {f:"conjugate(z^2)+c",  c:'.2541-i*0.2302',    n:'50', escape:"abs(z)>2", zoom: .65, center: [0,0]},
  "burning ship":   {f:"(abs(re(z))+i*abs(im(z)))^2-c", c:'-.8217+i*0.1233', n:'100', escape:"abs(z)>2", center:[0,0], zoom: .5},
  "butterfly":      {f:"conjugate(z^2)+c*re(1/z)", c:'-.4547-i*.7733', n:'150', escape:"abs(z)>4", center:[0,0], zoom:.9},
  "exponential map":{f:"c*e^(z-1)", c:'1.418-i*.119', n:'25', escape:'re(log(c)+z)>3000', zoom:.13, center:[8,0]},
  "exp Schwarz":    {f:"c0 = c^2/z;c1 = lambertw(-c0);conjugate(c0/exp(c1+c^2/c1));", c:'2.92-.48*i', n:'50', escape:"if(re(z)<-5,true,if(abs(lambertw(-c^2/f(z,c)))>abs(c),true,false));", center:[6,0], zoom:.15}
};

var iniscript = function(preset) {
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

  // initial values of parameters
  n = ${preset.n}; // number of iterates to use when generating image
  c = ${preset.c}; // initial value of c
  nplot = 6; // number of iterates to plot
  zoom = ${preset.zoom};
  center = [${preset.center}];

  f(z, c) := (
      ${preset.f};
  );

  escape(z, c) := (
      ${preset.escape};
  );


  // Coordinate transformations
  PltToCanvX(x) := (x-center_1)*zoom+1; // plot coordinate to canvas coordinate
  PltToCanvY(y) := (y-center_2)*zoom+1; // plot coordinate to canvas coordinate
  PltToCanvZ(z) := PltToCanvX(re(z))+i*PltToCanvY(im(z)); // plot coordinate to canvas coordinate

  CanvToPltX(x) := (x-1)/zoom+center_1; // canvas coordinate to plot coordinate
  CanvToPltY(y) := (y-1)/zoom+center_2; // canvas coordinate to plot coordinate
  CanvToPltZ(z) := CanvToPltX(re(z))+i*CanvToPltY(im(z)); // canvas coordinate to plot coordinate


  // Iteration functions
  preIter(z, c) := ( //returns the number of iterates to escape for dynamical plane, or n (max) if escape fails
      kmax=0;
      repeat(n,k,
          if(not(escape(z,c)),
              z = f(z,c);
              kmax=k;
          );
      );
  );

  preKIter(z,c,k) := ( // returns first k iterates of dynamical plane starting at z
      zs = [z];
      repeat(k,l,
          if(not(escape(z,c)),
              z = f(z,c);
              zs=append(zs,z);
          );
      );
      append(zs,f(z,c));
  );

  paramIter(z)  := (preIter(z,z)); // iterator for parameter space
  dynIter(z,c)  := (preIter(z,c)); // iterator for dynamical plane

  paramKIter(z,k) := (preKIter(z,z,k)); // generates first k iterates for parameter space
  dynKIter(z,c,k) := (preKIter(z,c,k)); // generates first k iterates for dynamical plane

  //Colors
  Z0.color  = (1,1,1);
  colorFcn(u) := (
      if(u==n,(0,0,0),
          u = u/n;
          u = (3*u/(2*u+1));
          (4*u,1.3*u,(1-u)^2*.7);
      );
  );

  //Generate image
  createimage("julia", 500, 500);
  `
}


class FractalPlot {
  constructor(varName, paramDict, canvasName, canvasID, callbacks = {}, fractType = "dyn", canvasWidth = 500, canvasHeight = 500) {
    this._varName = varName;
    this._canvasName = canvasName;
    this._canvasID = canvasID;
    this._fractType = fractType;
    this._canvasWidth = canvasWidth;
    this._canvasHeight = canvasHeight;
    this._callbacks = callbacks;
    this._c = paramDict.c;
    this._f = paramDict.f;
    this._n = paramDict.n;
    this._esc = paramDict.escape;
    this._center = paramDict.center;
    this._zoom = paramDict.zoom;
    this._z0 = reim(paramDict.c);
    this._movescript = null;
    switch(this._fractType) {
      case "dyn":
        this._movescript = 'colorplot([center_1-1/zoom,center_2-1/zoom],[center_1+1/zoom,center_2-1/zoom],"julia",colorFcn(dynIter(complex(#), c)));' +
        'drawimage([0,0],[2,0], "julia");' +
        'connect(apply(dynKIter(CanvToPltZ(complex(Z0.xy)),c,nplot),reim(PltToCanvZ(#))),color->[1,1,1],size->1.8);' +
        'drawtext(Z0+(.015,.015), "z="+CanvToPltZ(complex(Z0.xy)), color->[1,1,1],size->15);';
        break;
      case "param":
        this._movescript = 'colorplot([center_1-1/zoom,center_2-1/zoom],[center_1+1/zoom,center_2-1/zoom],"julia",colorFcn(paramIter(complex(#))));' +
        'drawimage([0,0],[2,0], "julia");' +
        'connect(apply(paramKIter(CanvToPltZ(complex(Z0.xy)),nplot),reim(PltToCanvZ(#))),color->[1,1,1],size->1.8);' +
        'drawtext(Z0+(.015,.015), "z="+CanvToPltZ(complex(Z0.xy)), color->[1,1,1],size->15);';
        break;
    }
    if ("move" in this._callbacks) {
      this._movescript += 'javascript("' + this._callbacks.move.join("; ") +';");';
    }
    this._keydownscript = 'javascript("' + this._varName + '.keypress(\'"+"\\" + key()+"\'.charCodeAt(0))");';
    if ("keydown" in this._callbacks) {
      this._keydownscript += 'javascript("' + this._callbacks.keydown.join("; ") +';");';
    }
    this._cindy = CindyJS({
      canvasname: this._canvasName,
      scripts: {
        init: iniscript(paramDict),
        move: this._movescript,
        keydown: this._keydownscript
      },
      geometry: [{name: "Z0", kind: "P", type: "Free", pos: this.PlotToCanv(this._z0), size: 3 }],
      ports: [{
        id: this._canvasID,
        width: this._canvasWidth,
        height: this._canvasHeight,
        transform: [{ visibleRect: [0,2,2,0] }],},]
    });
  }

  ApplyPreset(preset) {
    this.center = preset.center;
    this.zoom = preset.zoom;
    this.c = preset.c;
    this.f = preset.f;
    this.n = preset.n;
    this.esc = preset.escape;
    this.z0 = reim(this.c);
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
    this.center = [this._center[0] + vec[0],this._center[1] + vec[1]];
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


  set z0(z0Val) {
    this._z0 = z0Val;
    this._cindy.evokeCS('Z0.xy=[' + this.PlotToCanv(this._z0) + '];');
  }

  set c(cval) {
    this._c = cval;
    this._cindy.evokeCS('c=' + this._c + ';');
  }

  set f(fval) {
    this._f = fval;
    this._cindy.evokeCS('f(z,c) := (' + this._f + ');');
  }

  set esc(escval) {
    this._esc = escval;
    this._cindy.evokeCS('escape(z,c) := (' + this._esc + ');');
  }

  set n(nval) {
    this._n = nval;
    this._cindy.evokeCS('n=' + this._n + ';');
  }

  set zoom(zoomval) {
    this._zoom = zoomval;
    this._cindy.evokeCS('zoom=' + zoomval + ';');
  }

  set center(centerval) {
    this._center = centerval;
    var str = 'center_1='+this._center[0]+';'+
              'center_2='+this._center[1]+';';
    this._cindy.evokeCS(str);
  }


  get z0() {
    this._cindy.evokeCS('javascript("' + this._varName + '._z0 = ' + this._varName + '.CanvToPlot("+Z0.xy+")");');
    return this._z0;
  }

  get zoom() {
    return this._zoom;
  }

  get center() {
    return this._center;
  }

  get c() {
    return this._c;
  }

  get f() {
    return this._f;
  }

  get esc() {
    return this._esc;
  }

  get n() {
    return this._n;
  }

  get cindy() {
    return this._cindy;
  }

  get range() {
    return [this._center[0]-1/this._zoom,this._center[0]+1/this._zoom,
            this._center[1]-1/this._zoom,this._center[1]+1/this._zoom];
  }

}

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

var getCInput = function() {
  return document.getElementById("inpc").value;
}

var getNInput = function() {
  return parseInt(document.getElementById("inpn").value);
}

var getFInput = function() {
  return document.getElementById("inpf").value;
}

var getDEscInput = function() {
  return document.getElementById("inpje").value;
}

var getPEscInput = function() {
  return document.getElementById("inpme").value;
}

var getPCenterInput = function() {
  return document.getElementById("inpparamcenter").value.split(",").map(parseFloat);
}

var getPZoomInput = function() {
  return parseFloat(document.getElementById("inpparamzoom").value);
}

var getDCenterInput = function() {
  return document.getElementById("inpdyncenter").value.split(",").map(parseFloat);
}

var getDZoomInput = function() {
  return parseFloat(document.getElementById("inpdynzoom").value);
}

//////////

var setCInput = function(cval) {
  document.getElementById("inpc").value = cval;
}

var setNInput = function(nval) {
  document.getElementById("inpn").value = nval;
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


var apply_changes = function() {
  julia_fract.ApplyPreset(    getPresetDicts()[1]);
  parameter_fract.ApplyPreset(getPresetDicts()[0]);
}

var getPresetDicts = function() {
  return [{f:getFInput(), c:getCInput(), n:getNInput(), escape:getPEscInput(), zoom: getPZoomInput(), center: getPCenterInput()},
          {f:getFInput(), c:getCInput(), n:getNInput(), escape:getDEscInput(), zoom: getDZoomInput(), center: getDCenterInput()}];
}


var apply_preset = function(preset_val) {
  setPCenterInput(param_preset_dict[preset_val].center);
  setPZoomInput(param_preset_dict[preset_val].zoom);
  setDCenterInput(dyn_preset_dict[preset_val].center);
  setDZoomInput(dyn_preset_dict[preset_val].zoom);
  setCInput(param_preset_dict[preset_val].c);
  setNInput(param_preset_dict[preset_val].n);
  setFInput(param_preset_dict[preset_val].f);
  setPEscInput(param_preset_dict[preset_val].escape);
  setDEscInput(dyn_preset_dict[preset_val].escape);
  julia_fract.ApplyPreset(    dyn_preset_dict[preset_val]);
  parameter_fract.ApplyPreset(param_preset_dict[preset_val]);
}

