var preset_dict = {
  "mandelbrot":     {f:"z^2+c",             c:'-.7-.4*i', n:'50', escape:"abs(z)>2", paramescape:"abs(z)>2",
                     plimits:[-2.1,.8,-1.45,1.45],  jlimits:[-2,2,-2,2]},
  "tricorn":        {f:"conjugate(z^2)+c",  c:'-0.488 + i*0.128',    n:'50', escape:"abs(z)>2", paramescape:"abs(z)>2",
                     plimits:[-2,2,-2,2],           jlimits:[-2,2,-2,2]},
  "burning ship":   {f:"(abs(re(z))+i*abs(im(z)))^2-c", c:'1.6185 + i*0.0471', n:'50', escape:"abs(z)>2", paramescape:"abs(z)>2", plimits:[1.6,1.8,-.05,.15], jlimits:[-2,2,-2,2]},
  "exponential map":{f:"c*e^(z-1)", c:'1', n:'25', escape:'re(log(c)+z)>3000', paramescape:'re(log(c)+z)>3000',
                     plimits:[-1,5,-3,3],           jlimits:[-1,11,-6,6]},
  "exp Schwarz":    {f:"c0 = c^2/z;c1 = lambertw(-c0);conjugate(c0/exp(c1+c^2/c1));",             c:'1+0*i',    n:'50', escape:"if(re(z)<-5,true,if(abs(lambertw(-c^2/f(z,c)))>abs(c),true,false));",
                     paramescape:"abs(lambertw(-c^2/f(z,c)))>abs(c)",
                     plimits:[-1,5,-3,3],           jlimits:[-1,11,-6,6]}
};




var gslp = [
    { name: "C", kind: "P", type: "Free", pos: [2/3,1], size: 3 },
];

var gslp2 = [
    { name: "I", kind: "P", type: "Free", pos: [2/3,1], size: 3 },
];


var mcdy = CindyJS({
    canvasname: "MCSCanvas",
    scripts: "mcs*",
    initscript:"csinit",
    geometry: gslp,
    ports: [
        {
            id: "MCSCanvas",
            width: 500,
            height: 500,
            transform: [{ visibleRect: [0,0,2,2] }],
        },
    ],
});


var jcdy = CindyJS({
    canvasname: "JCSCanvas",
    scripts: "jcs*",
    initscript:"csinit",
    geometry: gslp2,
    ports: [
        {
            id: "JCSCanvas",
            width: 500,
            height: 500,
            transform: [{ visibleRect: [0,0,2,2] }],
        },
    ],
});

var update_cindyvars=function(f,c,n,mesc,jesc){
  mcdy.evokeCS('moveto(C,reim(mPltToCanvZ(' + c + ')));');
  mcdy.evokeCS('f(z,c) := (' + f + ');');
  jcdy.evokeCS('f(z,c) := (' + f + ');');
  mcdy.evokeCS('mandelescape(z,c) := (' + mesc + ');');
  jcdy.evokeCS('juliaescape(z,c) := (' + jesc + ');');
  mcdy.evokeCS('n=' + n + ';');
  jcdy.evokeCS('n=' + n + ';');
}

var apply_changes=function(){
  c = document.getElementById('inpc').value;
  n = document.getElementById('inpn').value;
  f = document.getElementById('inpf').value;
  mesc = document.getElementById('inpme').value;
  jesc = document.getElementById('inpje').value;
  update_cindyvars(f,c,n,mesc,jesc);
}

var apply_preset=function(){
  preset_choice = document.getElementById('fractal_presets').value;
  params = preset_dict[preset_choice];
  console.log('setJuliaWindow([' + params.jlimits + ']);');
  jcdy.evokeCS('setJuliaWindow([' + params.jlimits + ']);');
  mcdy.evokeCS('setMandelWindow([' + params.plimits + ']);');
  document.getElementById('inpc').value = params.c;
  document.getElementById('inpn').value = params.n;
  document.getElementById('inpf').value = params.f;
  document.getElementById('inpme').value = params.paramescape;
  document.getElementById('inpje').value = params.escape;
  apply_changes();
}

var updatec=function(cval){
    jcdy.evokeCS('c=complex([' + cval[0] + ',' + cval[1] + ']);');
}



var setf=function(e, b){
  var chCode=e.which ? e.which:e.keyCode;
  //if(chCode==13){
    mcdy.evokeCS('f(z,c) := (' + b.value + ');');
  //}
}

var setparameterescape=function(e, b){
  var chCode=e.which ? e.which:e.keyCode;
  //if(chCode==13){
    mcdy.evokeCS('mandelescape(z,c) := (' + b.value + ');');
  //}
}

var setdynplaneescape=function(e, b){
  var chCode=e.which ? e.which:e.keyCode;
  //if(chCode==13){
    jcdy.evokeCS('juliaescape(z,c) := (' + b.value + ');');
  //}
}

var setc=function(b){
    mcdy.evokeCS('moveto(C,reim(mPltToCanvZ(' + b.value + ')));');
}

var setn=function(e, b){
  var chCode=e.which ? e.which:e.keyCode;
  //if(chCode==13){
    mcdy.evokeCS('n=' + b.value + ';');
    jcdy.evokeCS('n=' + b.value + ';');
  //}
}

// Parameter space coordinate transformations

var paramleft=function(){
  //if(chCode==13){
    mcdy.evokeCS('mShift(-1/8,"horizontal");');
  //}
}

var paramright=function(){
  //if(chCode==13){
    mcdy.evokeCS('mShift(1/8,"horizontal");');
  //}
}
var paramup=function(){
  //if(chCode==13){
    mcdy.evokeCS('mShift(1/8,"vertical");');
  //}
}

var paramdown=function(){
  //if(chCode==13){
    mcdy.evokeCS('mShift(-1/8,"vertical");');
  //}
}

var paramzoomin=function(){
  //if(chCode==13){
    mcdy.evokeCS("mZoom(2/3);");
  //}
}

var paramzoomout=function(){
  //if(chCode==13){
    mcdy.evokeCS("mZoom(3/2);");
  //}
}

// Dynamical plane coordinate transformations

var jleft=function(){
  //if(chCode==13){
    jcdy.evokeCS('jShift(-1/8,"horizontal");');
  //}
}

var jright=function(){
  //if(chCode==13){
    jcdy.evokeCS('jShift(1/8,"horizontal");');
  //}
}
var jup=function(){
  //if(chCode==13){
    jcdy.evokeCS('jShift(1/8,"vertical");');
  //}
}

var jdown=function(){
  //if(chCode==13){
    jcdy.evokeCS('jShift(-1/8,"vertical");');
  //}
}

var jzoomin=function(){
  //if(chCode==13){
    jcdy.evokeCS("jZoom(2/3);");
  //}
}

var jzoomout=function(){
  //if(chCode==13){
    jcdy.evokeCS("jZoom(3/2);");
  //}
}