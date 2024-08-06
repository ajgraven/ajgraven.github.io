use("CindyGL");

// Fixed constants
e = 2.71828182845904523536028747; // Euler's number

// Initial values of constants
n = 100; // number of iterates to use when generating image
c = -.7-.4*i; // initial value of c
nplot = 8; // number of iterates to plot





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




// Dynamical plane plot coordinate transformations
jminx = -2;
jminy = -2;
jmaxx = 2;
jmaxy = 2;
jcanvminx = 0;
jcanvminy = 0;
jcanvmaxx = 2;
jcanvmaxy = 2;


jPltToCanvX(x) := (x-jminx)*(2/(jmaxx-jminx)); // plot coordinate to canvas coordinate
jPltToCanvY(y) := (y-jminy)*(2/(jmaxy-jminy)); // plot coordinate to canvas coordinate
jPltToCanvZ(z) := jPltToCanvX(re(z))+i*jPltToCanvY(im(z)); // plot coordinate to canvas coordinate

jCanvToPltX(x) := x*(jmaxx-jminx)/2+jminx; // canvas coordinate to plot coordinate
jCanvToPltY(y) := y*(jmaxy-jminy)/2+jminy; // canvas coordinate to plot coordinate
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

juliaiter(z, c) := ( //returns the number of iterates to escape for dynamical plane, or n (max) if escape fails
    kmax=0;
    repeat(n,k,
        if(not(juliaescape(z,c)),
            z = f(z,c);
            kmax=k;
        );
    );
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



jShift(ratio,dir) := ( //Translate parameter space window by ratio of window width/height
    Z0pos = Z0.xy;
    if(dir=="horizontal",
        shift = (jmaxx-jminx)*ratio;
        jminxtemp = jminx + shift;
        jmaxx = jmaxx + shift;
        jminx = jminxtemp;
        Z0.x = Z0pos_1,
        shift = (jmaxy-jminy)*ratio;
        jminytemp = jminy + shift;
        jmaxy = jmaxy + shift;
        jminy = jminytemp;
        Z0.y = Z0pos_2;
    );
);

jZoom(ratio) := ( //Zooms parameter space window by ratio of window width/height
    Z0pos = Z0.xy;
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
    Z0.xy = Z0pos;
);


//Generate images
Z0.color  = (1,1,1);
createimage("julia", 800, 800);