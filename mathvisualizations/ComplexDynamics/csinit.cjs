use("CindyGL");

// Fixed constants
e = 2.71828182845904523536028747; // Euler's number

// Initial values of constants
n = 50; // number of iterates to use when generating image
c = -.7-.4*i; // initial value of c
nplot = 6; // number of iterates to plot





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

f(z, c) := (
    z^2+c;
);

color(u) := (
    //hue(re(log(1+(e-1)*(u/n))));
    u = n-u;
    (re(log(u+1)/log(n+1)),(u/n)^2,.3)
);



// Parameter space plot coordinate transformations
mandelminx = -2.1;
mandelminy = -1.45;
mandelmaxx = .8;
mandelmaxy = 1.45;
mandelcanvminx = 0;
mandelcanvminy = 0;
mandelcanvmaxx = 2;
mandelcanvmaxy = 2;


mPltToCanvX(x) := (x-mandelminx)*((mandelcanvmaxx-mandelcanvminx)/(mandelmaxx-mandelminx))+mandelcanvminx; // plot coordinate to canvas coordinate
mPltToCanvY(y) := (y-mandelminy)*((mandelcanvmaxy-mandelcanvminy)/(mandelmaxy-mandelminy))+mandelcanvminy; // plot coordinate to canvas coordinate
mPltToCanvZ(z) := mPltToCanvX(re(z))+i*mPltToCanvY(im(z)); // plot coordinate to canvas coordinate

mCanvToPltX(x) := (x-mandelcanvminx)*((mandelmaxx-mandelminx)/(mandelcanvmaxx-mandelcanvminx))+mandelminx; // canvas coordinate to plot coordinate
mCanvToPltY(y) := (y-mandelcanvminy)*((mandelmaxy-mandelminy)/(mandelcanvmaxy-mandelcanvminy))+mandelminy; // canvas coordinate to plot coordinate
mCanvToPltZ(z) := mCanvToPltX(re(z))+i*mCanvToPltY(im(z)); // canvas coordinate to plot coordinate

setMandelWindow(arr) := (
    mandelminx = arr_1;
    mandelmaxx = arr_2;
    mandelminy = arr_3;
    mandelmaxy = arr_4;
);



// Dynamical plane plot coordinate transformations
jminx = -2;
jminy = -2;
jmaxx = 2;
jmaxy = 2;
jcanvminx = 0;
jcanvminy = 0;
jcanvmaxx = 2;
jcanvmaxy = 2;


jPltToCanvX(x) := (x-jminx)*((jcanvmaxx-jcanvminx)/(jmaxx-jminx))+jcanvminx; // plot coordinate to canvas coordinate
jPltToCanvY(y) := (y-jminy)*((jcanvmaxy-jcanvminy)/(jmaxy-jminy))+jcanvminy; // plot coordinate to canvas coordinate
jPltToCanvZ(z) := jPltToCanvX(re(z))+i*jPltToCanvY(im(z)); // plot coordinate to canvas coordinate

jCanvToPltX(x) := (x-jcanvminx)*((jmaxx-jminx)/(jcanvmaxx-jcanvminx))+jminx; // canvas coordinate to plot coordinate
jCanvToPltY(y) := (y-jcanvminy)*((jmaxy-jminy)/(jcanvmaxy-jcanvminy))+jminy; // canvas coordinate to plot coordinate
jCanvToPltZ(z) := jCanvToPltX(re(z))+i*jCanvToPltY(im(z)); // canvas coordinate to plot coordinate

setJuliaWindow(arr) := (
    jminx = arr_1;
    jmaxx = arr_2;
    jminy = arr_3;
    jmaxy = arr_4;
);


juliaescape(z, c) := (
    abs(z)>2;
);

mandelescape(z, c) := (
    abs(z)>2;
);

juliaiter(z, c) := ( //returns the number of iterates to escape for dynamical plane, or n (max) if escape fails
    kmax=0;
    repeat(n,k,
        if(not(juliaescape(z,c)),
            z = f(z,c);
            kmax=k;
        );
    );
);


mandeliter(z, c) := ( //returns the number of iterates to escape for dynamical plane, or n (max) if escape fails
    kmax=0;
    repeat(n,k,
        if(not(mandelescape(z,c)),
            z = f(z,c);
            kmax=k;
        );
    );
);

mandelniter(c,k) := ( // returns first k iterates of dynamical plane starting at z
    zs = [c];
    z = c;
    repeat(k,l,
        if(not(mandelescape(z,c)),
            z = f(z,c);
            zs=append(zs,z);
        );
    );
    append(zs,f(z,c));
);

julianiter(z,c,k) := ( // returns first k iterates of dynamical plane starting at z
    zs = [z];
    repeat(k,l,
        if(not(juliaescape(z,c)),
            z = f(z,c);
            zs=append(zs,z);
        );
    );
    append(zs,f(z,c));
);



//Generate images
if(C==C,
    C.xy = reim(mPltToCanvZ(c));
    // parameter space script
    mShift(ratio,dir) := ( //Translates parameter space window by ratio of window width/height
        Cpos = C.xy;
        if(dir=="horizontal",
            shift = (mandelmaxx-mandelminx)*ratio;
            mandelminxtemp = mandelminx + shift;
            mandelmaxx = mandelmaxx + shift;
            mandelminx = mandelminxtemp;
            C.x = Cpos_1,
            shift = (mandelmaxy-mandelminy)*ratio;
            mandelminytemp = mandelminy + shift;
            mandelmaxy = mandelmaxy + shift;
            mandelminy = mandelminytemp;
            C.y = Cpos_2;
        );
    );
    mZoom(ratio) := ( //Zoom parameter space window by ratio of window width/height
        Cpos = C.xy;
        cx = (mandelmaxx+mandelminx)/2;
        cy = (mandelmaxy+mandelminy)/2;
        mxminc = (mandelminx-cx)*ratio;
        mxmaxc = (mandelmaxx-cx)*ratio;
        myminc = (mandelminy-cy)*ratio;
        mymaxc = (mandelmaxy-cy)*ratio;
        mandelminx = mxminc + cx;
        mandelmaxx = mxmaxc + cx;
        mandelminy = myminc + cy;
        mandelmaxy = mymaxc + cy;
        C.xy = Cpos;
    );
    c = mCanvToPltZ(complex(C));
    C.color  = (0,0,1);
    createimage("mandel", 800, 800);, 
    // dynamical plane script
    jShift(ratio,dir) := ( //Translate parameter space window by ratio of window width/height
        Ipos = I.xy;
        if(dir=="horizontal",
            shift = (jmaxx-jminx)*ratio;
            jminxtemp = jminx + shift;
            jmaxx = jmaxx + shift;
            jminx = jminxtemp;
            I.x = Ipos_1,
            shift = (jmaxy-jminy)*ratio;
            jminytemp = jminy + shift;
            jmaxy = jmaxy + shift;
            jminy = jminytemp;
            I.y = Ipos_2;
        );
    );
    jZoom(ratio) := ( //Zooms parameter space window by ratio of window width/height
        Ipos = I.xy;
        cx = (jmaxx+jminx)/2;
        cy = (jmaxy+jminy)/2;
        jxminc = (jminx-cx)*ratio;
        jxmaxc = (jmaxx-cx)*ratio;
        jyminc = (jminy-cy)*ratio;
        jymaxc = (jmaxy-cy)*ratio;
        jminx = jxminc + cx;
        jmaxx = jxmaxc + cx;
        jminy = jyminc + cy;
        jmaxy = jymaxc + cy;
        I.xy = Ipos;
    );
    I.color  = (0,0,1);
    createimage("julia", 800, 800) 
);