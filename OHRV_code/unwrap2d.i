/*
  UNWRAP2D.I
  
  Functions for phase unwrapping 2D data 


  P.M. Celliers
  LLNL
  March 2010
  
    
*/

require, "unwrap2D_baldi.i";
require, "unwrap2D_GERI.i";
//require, "fft_r2r.i";
require, "image.i";

func unwrap2d(img, tsize, &weights, weight=, seed=, ltest=, split=,  plot=,
              w=, timing=, skip_inner=, bmethod=, method=, mask=, NIT=, NIT_lp=)
/* DOCUMENT unwrap2d, tsize, weight=, seed=, ltest=, split=,  plot=,
              w=, timing=, skip_inner=, bmethod=, method=, mask=

   ARGUMENTS:
      tsize - tile size in pixels

   KEYWORDS:
      method= 1 -> GERI reliability method
              2 -> Baldi region growing method
              3 -> DCT algorithm 1
              4 -> DCT algorithm 2
              5 -> DCT algorithm 3
              6 -> Lp norm algorithm
     
   SEE ALSO:
 */
{
  if (is_void(method)) method= 1;
  //print, "unwrap2d:: method = ", method;
  
  if (method == 1 || method == "reliability") {
    if(is_void(mask)) {
      return img_copy(img, data= unwrap2D_GERI(img_data(img)));
    } else {
      return img_copy(img, data= unwrap2D_GERI(img_data(img), mask));
    }
    
  } else if (method == 2 || method == "Baldi") {
    return unwrap2d_Baldi(img, tsize, weight= weight, seed= seed, ltest= ltest,
                          split= split, plot= plot, w= w, timing= timing, skip_inner= skip_inner,
                          method= bmethod);
    
  } else if (method == 3) {
    return img_copy(img, data= _unwrap2d_dct0(img_data(img)));
    
  } else if (method == 4) {
    return img_copy(img, data= _unwrap2d_dct1(img_data(img), img_data(weight), NIT= NIT));
    
  } else if (method == 5) {
    return img_copy(img, data= _unwrap2d_dct2(img_data(img), img_data(weight), NIT= NIT));
    
  } else if (method == 6) {
    return img_copy(img, data= _unwrap2d_LPnorm(img_data(img), weights, mask= mask, NIT= NIT, NIT_lp= NIT_lp));
    
  } else error, "Invalid method";
}

func _unwrap2d_dct0(phi)
/* DOCUMENT _unwrap2d_dct0(phi)

   Apply the discrete cosine transformation method to determine an
   unwrapped 2D phase from a wrapped phase.  This algorithm is not too
   useful for noisy or inconsistent data.

   Implements algorithm 1 in the article "Robust two-dimensional
   weighted and unweighted phase unwrapping that uses fast transforms
   and iterative methods", D.C. Ghiglia and L.A. Romero,
   J.Opt.Soc.Am. vol 11, p 107 (1994)
   
   SEE ALSO:
     _unwrap2d_dct1, _unwrap2d_dct2, fftw_r2r
 */
{
  extern phiw, phiwdx, phiwdy, rho;
  
  dims= dimsof(phi);
  n= dims(2); m= dims(3);
  
  //Some checks on the input data
  if (dims(1) != 2) error, "Can only unwrap 2D arrays";
  if (phi(*)(min) < -pi || phi(*)(max) > pi) error, "Data not in wrapped range [-pi, pi]";
  if (n < 3 || m < 3) error, "Need at least 3x3 array";

  //Working storage
  //deltax= deltay= rho= array(0.0, dimsof(phi));
  phiw= array(0.0, n+2, m+2);
  igrid= indgen(0:n-1)(,-:1:m);
  jgrid= indgen(0:m-1)(-:1:n,);

  //Populate phi working array
  phiw(2:n+1, 2:m+1)= phi;
  phiw(1,)= phiw(2,);
  phiw(n+2,)= phiw(n+1,);
  phiw(,1)= phiw(,2);
  phiw(,m+2)= phiw(,m+1);
   
  //Compute the delta differences
  phiwdx= _wrap(phiw(dif,));
  phiwdy= _wrap(phiw(,dif));
  rho= phiwdx(2:n+1,2:m+1) - phiwdx(1:n,2:m+1) +
         phiwdy(2:n+1,2:m+1) - phiwdy(2:n+1,1:m);

  //Compute rho
  //deltax(1:n-1,)= _wrap(phi(dif,));
  //deltay(,1:m-1)= _wrap(phi(,dif));
  //rho(1,1)= deltax(1,1) + deltay(1,1);        //first element
  //rho(1,2:m)= deltax(1,2:m) + deltay(1,dif);  //first x row
  //rho(2:n,1)= deltax(dif,1) + deltay(2:n,1);  //first y row
  //rho(2:n,2:m)= deltax(dif,2:m)+deltay(2:n,dif);  //remaining elements

  //Solve for phi (equations 13 & 17)
  //RHO= dctw(rho,1)/n/m;
  RHO= fftw_r2r(rho,"REDFT10")/n/m/2;
  den= cos(pi*igrid/n) + cos(pi*jgrid/m) - 2.0;
  den(1,1)= 1.0;
  PHI= RHO/den/4;  //equation 13
  //return dctw(PHI,-1);
  return fftw_r2r(PHI,"REDFT01");
}

func _unwrap2d_dct1(phi, weight, NIT=)
/* DOCUMENT _unwrap2d_dct1(phi, weight, NIT=)

   Apply the discrete cosine transformation method to determine an
   unwrapped 2D phase from a wrapped phase and an array of weights
   associated with the wrapped phase data.

   Implements algorithm 2 in the article "Robust two-dimensional
   weighted and unweighted phase unwrapping that uses fast transforms
   and iterative methods", D.C. Ghiglia and L.A. Romero,
   J.Opt.Soc.Am. vol 11, p 107 (1994)
   
   SEE ALSO:
     _unwrap2d_dct0, unwrap2d_dct1, dct, dctw
 */
{
  dims= dimsof(phi);
  n= dims(2); m= dims(3);
  if(is_void(NIT)) NIT= 10;
  if (is_void(weight)) weight= array(1.0, dimsof(phi));
  
  //Some checks on the input data
  if (rankof(phi) != rankof(weight)) error, "Unmatched array ranks";
  if (dims(1) != 2) error, "Can only unwrap 2D arrays";
  if (phi(*)(min) < -pi || phi(*)(max) > pi) error, "Data not in wrapped range [-pi, pi]";
  if (n < 3 || m < 3) error, "Need at least 3x3 array";
  if (!allof(dimsof(phi) == dimsof(weight))) error, "Unmatched array dimensions";

  //Scale weights to [0,1]
  weight/= weight(*)(max);
  clip, weight, 0.0, 1.0;

  //Working storage
  cij= rho= deltax= deltay= array(0.0, dimsof(phi));
  //deltax= deltay= array(0.0, n+1, m+1);
  w= phiw= array(0.0, n+2, m+2);
  igrid= indgen(0:n-1)(,-:1:m);
  jgrid= indgen(0:m-1)(-:1:n,);

  //Weighting factors
  w(2:n+1,2:m+1)= weight;
  w(1,)= w(2,); w(n+2,)= w(n+1,);
  w(,1)= w(,2); w(,m+2)= w(,m+1);
  w2= w*w;
  wx1= [w2(3:n+2,2:m+1), w2(2:n+1,2:m+1)](..,min);
  wx2= [w2(2:n+1,2:m+1), w2(1:n,2:m+1)](..,min);
  wy1= [w2(2:n+1,3:m+2), w2(2:n+1,2:m+1)](..,min);
  wy2= [w2(2:n+1,2:m+1), w2(2:n+1,1:m)](..,min);

  //Populate phi working array for this iteration
  phiw(2:n+1, 2:m+1)= phi;
  phiw(1,)= phiw(2,);
  phiw(,1)= phiw(,2);
  phiw(n+2,)= phiw(n+1,);
  phiw(,m+2)= phiw(,m+1);
    
  //Compute the delta differences
  deltax= _wrap(phiw(dif,))(2:,2:m+1);
  deltay= _wrap(phiw(,dif))(2:n+1,2:);
  deltaxm1= roll(deltax, [1,0]); deltaxm1(1,)= 0.0;
  deltaym1= roll(deltay, [0,1]); deltaym1(,1)= 0.0;

  //Compute cij with weighting factors (equation 34)
  cij= wx1*deltax - wx2*deltaxm1 + wy1*deltay - wy2*deltaym1;

  //Initial guess on phi
  phi0= array(0.0, dimsof(phi));
  k= 0;
  
  do {

    //Initialize phi
    phiw(2:n+1, 2:m+1)= phi0;
    phiw(1,)= phiw(2,);
    phiw(,1)= phiw(,2);
    phiw(n+2,)= phiw(n+1,);
    phiw(,m+2)= phiw(,m+1);
    
    //Compute phi difference terms
    dphix= phiw(dif,);
    dphiy= phiw(,dif);
    dphix1= dphix(2:n+1,2:m+1);
    dphix2= dphix(1:n,2:m+1);
    dphiy1= dphiy(2:n+1,2:m+1);
    dphiy2= dphiy(2:n+1,1:m);

    //Compute rho iterate (equation 35)
    rho= cij - ((wx1 - 1)*dphix1 - (wx2-1)*dphix2 + (wy1 - 1)*dphiy1 - (wy2-1)*dphiy2);

    //Solve for phi (equation 17)
    RHO= fftw_r2r(rho,"REDFT10")/n/m/2;
    den= cos(pi*igrid/n) + cos(pi*jgrid/m) - 2.0;
    den(1,1)= 1.0;
    phi0= fftw_r2r(RHO/den/4, "REDFT01");
    
  } while (k++ < NIT);
  
  return phi0;
}

func _unwrap2d_dct2(phi, weight, NIT=)
/* DOCUMENT _unwrap2d_dct2(phi, weight, NIT=)

   Apply the discrete cosine transformation method to determine an
   unwrapped 2D phase from a wrapped phase (phi) and an array of
   weights (weight) associated with the wrapped phase data.

   Implements algorithm 3 in the article "Robust two-dimensional
   weighted and unweighted phase unwrapping that uses fast transforms
   and iterative methods", D.C. Ghiglia and L.A. Romero,
   J.Opt.Soc.Am. vol 11, p 107 (1994).

   This algorithm is iterative and employs a preconditioned conjugate
   gradient method to accelerate the convergence. Nevertheless it may
   take several dozen iterations.
   
   SEE ALSO:
     _unwrap2d_dct0, unwrap2d_dct1, dct, dctw
 */
{
  extern w2, w, wght, rk, phik;
  
  dims= dimsof(phi);
  n= dims(2); m= dims(3);
  if(is_void(NIT)) NIT= 10;
  if (is_void(weight)) weight= array(1.0, dimsof(phi));
  
  //Some checks on the input data
  if (rankof(phi) != rankof(weight)) error, "Unmatched array ranks";
  if (dims(1) != 2) error, "Can only unwrap 2D arrays";
  if (phi(*)(min) < -pi || phi(*)(max) > pi) error, "Data not in wrapped range [-pi, pi]";
  if (n < 3 || m < 3) error, "Need at least 3x3 array";
  if (!allof(dimsof(phi) == dimsof(weight))) error, "Unmatched array dimensions";

  //Scale weights to [0,1]
  wght= weight/(median(weight(*)) + 2*weight(*)(rms));
  clip, wght, 0.0, 1.0;

  //Working storage
  cij= rho= deltax= deltay= array(0.0, dimsof(phi));
  //deltax= deltay= array(0.0, n+1, m+1);
  w= phiw= array(0.0, n+2, m+2);
  igrid= indgen(0:n-1)(,-:1:m);
  jgrid= indgen(0:m-1)(-:1:n,);

  //Compute weighting factors
  w(2:n+1,2:m+1)= wght;
  w(1,)= w(2,); w(n+2,)= w(n+1,);
  w(,1)= w(,2); w(,m+2)= w(,m+1);
  w2= w*w;
  wx1= [w2(3:n+2,2:m+1), w2(2:n+1,2:m+1)](..,min);
  wx2= [w2(2:n+1,2:m+1), w2(1:n,2:m+1)](..,min);
  wy1= [w2(2:n+1,3:m+2), w2(2:n+1,2:m+1)](..,min);
  wy2= [w2(2:n+1,2:m+1), w2(2:n+1,1:m)](..,min);
  w2= [wx1, wx2, wy1, wy2];

  //Populate phi working array for this iteration
  phiw(2:n+1, 2:m+1)= phi;
  phiw(1,)= phiw(2,);
  phiw(,1)= phiw(,2);
  phiw(n+2,)= phiw(n+1,);
  phiw(,m+2)= phiw(,m+1);
    
  //Compute the delta differences
  deltax= _wrap(phiw(dif,))(2:,2:m+1);
  deltay= _wrap(phiw(,dif))(2:n+1,2:);
  deltaxm1= roll(deltax, [1,0]); deltaxm1(1,)= 0.0;
  deltaym1= roll(deltay, [0,1]); deltaym1(,1)= 0.0;

  //Algorithm 3, step 1
  //Compute cij with weighting factors (equation 34)
  rk= cij= wx1*deltax - wx2*deltaxm1 + wy1*deltay - wy2*deltaym1;

  //Initial guess on phi
  phik= array(0.0, dimsof(phi));
  k= 0;
  
  do {

    //Algoithm 3, step 2: Solve for zk (equation 17)
    RHO= fftw_r2r(rk,"REDFT10")/n/m/2;
    den= cos(pi*igrid/n) + cos(pi*jgrid/m) - 2.0;
    den(1,1)= 1.0;
    zk= fftw_r2r(RHO/den/4, "REDFT01");

    k++;  //step 3
    
    if (k == 1) {    //step 4
      rkm1= rk;
      zkm1= zk;
      pk= zkm1;
      betak= 0.0;
    } else {         //step 5
      zkm2= zkm1;
      rkm2= rkm1;
      rkm1= rk;
      zkm1= zk;
      pkm1= pk;
      betak= rkm1(*)(+)*zkm1(*)(+)/(rkm2(*)(+)*zkm2(*)(+));
      pk= zkm1 + betak*pkm1
    }

    //step 6
    phikm1= phik;
    Qpk= _dct2_Qop(pk,w2);  //Apply the weighted Laplacian operator to pk
    alphak= rkm1(*)(+)*zkm1(*)(+)/(pk(*)(+)*Qpk(*)(+));
    phik= phikm1 + alphak*pk;
    rk= rkm1 - alphak*Qpk;

    //write, format=" k = %02d, rk(rms) = %e, rk(*)(max) = %e, rk(*)(avg) = %e, alphak= %e, betak= %e\n",
    //  k, rk(*)(rms), rk(*)(max), rk(*)(avg), alphak, betak;

    //window, 60; fma; pli, rk;
    
  } while (k < NIT);
  
  return phik;
}

func _dct2_Qop(p, w2)
/* DOCUMENT _dct2_Qop(p, w2)

   Applies the weighted Laplacian operator to a 2D phase data set
     
   SEE ALSO:
     _unwrap2d_dct2
 */
{
  extern n, m, phiw;

  //Initialize phi
  phiw(2:n+1, 2:m+1)= p;
  phiw(1,)= phiw(2,);
  phiw(,1)= phiw(,2);
  phiw(n+2,)= phiw(n+1,);
  phiw(,m+2)= phiw(,m+1);
    
  //Compute phi difference terms
  dphix= phiw(dif,);
  dphiy= phiw(,dif);
  dphix1= dphix(2:n+1,2:m+1);
  dphix2= dphix(1:n,2:m+1);
  dphiy1= dphiy(2:n+1,2:m+1);
  dphiy2= dphiy(2:n+1,1:m);

  //Apply weighted Laplacian operator
  return w2(..,1)*dphix1 - w2(..,2)*dphix2 + w2(..,3)*dphiy1 - w2(..,4)*dphiy2;
}

func _unwrap2d_dct2_A(cij, w2, NIT=)
/* DOCUMENT _unwrap2d_dct2(cij, UV, NIT=)

   Apply the discrete cosine transformation method to determine an
   unwrapped 2D phase from a wrapped phase (phi) and an array of
   weights (weight) associated with the wrapped phase data.

   Implements algorithm 3 in the article "Robust two-dimensional
   weighted and unweighted phase unwrapping that uses fast transforms
   and iterative methods", D.C. Ghiglia and L.A. Romero,
   J.Opt.Soc.Am. vol 11, p 107 (1994).

   This algorithm is iterative and employs a preconditioned conjugate
   gradient method to accelerate the convergence. Nevertheless it may
   take several dozen iterations.
   
   SEE ALSO:
     _unwrap2d_dct0, unwrap2d_dct1, dct, dctw
 */
{
  extern rk, phik;
  
  dims= dimsof(cij);
  n= dims(2); m= dims(3);
  if(is_void(NIT)) NIT= 25;
  if(is_void(NIT_lp)) NIT_lp= 0;
  write, format="-----------> NIT=%02d\r", NIT;

  //Some checks on the input data
  if (dims(1) != 2) error, "Can only unwrap 2D arrays";
  if (n < 3 || m < 3) error, "Need at least 3x3 array";
  if (!allof(dimsof(cij)(2:3) == dimsof(w2)(2:3))) error, "Unmatched array dimensions";

  //Working storage
  phiw= array(0.0, n+2, m+2);
  igrid= indgen(0:n-1)(,-:1:m);
  jgrid= indgen(0:m-1)(-:1:n,);

  //Algorithm 3, step 1
  //Compute cij with weighting factors (equation 34)
  rk= cij;

  //Initial guess on phi
  phik= array(0.0, dimsof(cij));
  k= 0;
  
  do {

    //Algoithm 3, step 2: Solve for zk (equation 17)
    RHO= fftw_r2r(rk,"REDFT10")/n/m/2;
    den= cos(pi*igrid/n) + cos(pi*jgrid/m) - 2.0;
    den(1,1)= 1.0;
    zk= fftw_r2r(RHO/den/4, "REDFT01");

    k++;  //step 3
    
    if (k == 1) {    //step 4
      rkm1= rk;
      zkm1= zk;
      pk= zkm1;
      betak= 0.0;
    } else {         //step 5
      zkm2= zkm1;
      rkm2= rkm1;
      rkm1= rk;
      zkm1= zk;
      pkm1= pk;
      betak= rkm1(*)(+)*zkm1(*)(+)/(rkm2(*)(+)*zkm2(*)(+));
      pk= zkm1 + betak*pkm1
    }

    //step 6
    phikm1= phik;
    Qpk= _dct2_Qop(pk,w2);  //Apply the weighted Laplacian operator to pk
    alphak= rkm1(*)(+)*zkm1(*)(+)/(pk(*)(+)*Qpk(*)(+));
    phik= phikm1 + alphak*pk;
    rk= rkm1 - alphak*Qpk;

    write, format="l= %02d/%02d, k=%02d/%02d, rk(rms)=%e, rk(*)(max)=%e, rk(*)(avg)=%e, alphak=%e, betak=%e\r",
      l, NIT_lp, k, NIT, rk(*)(rms), rk(*)(max), rk(*)(avg), alphak, betak;

    //window, 60; fma; pli, rk;
    
  } while (k < NIT && rk(*)(max) > 1e-3);
  
  return phik;
}

func _unwrap2d_LPnorm(psi, &UV, phi0=, NIT_lp=, NIT=, p=, mask=)
/* DOCUMENT _unwrap2d_LPnorm(psi, &weights, NIT=, p=)

   Apply the minimum Lp-norm phase unwrapping algorithm of Ghiglia &
   Romero. 

   Implements the algorithm described in the article "Minimum Lp-norm
   two-dimensional phase unwrapping", D.C. Ghiglia and L.A. Romero,
   J.Opt.Soc.Am. vol 13, p 1999 (1996).

   This algorithm solves the weighted least-squares unwrapping

   KEYWORDS:
     phi0= <doesn't seem to do anything>
     NIT_lp= number of of lp-norm outer iterations
     NIT= number of dct2 inner iterations
     p=
     mask=
   
   SEE ALSO:
     _unwrap2d_dct2_A, _LP_weights
 */
{
  local f, g;
  extern R;
  
  dims= dimsof(psi);
  n= dims(2); m= dims(3);
  if (is_void(NIT)) NIT= 25;
  if (is_void(NIT_lp)) NIT_lp= 20;
  if (is_void(p)) p= 0.0;
  
  //Some checks on the input data
  if (dims(1) != 2) error, "Can only unwrap 2D arrays";
  if (psi(*)(min) < -pi || psi(*)(max) > pi)
    error, "Data not in wrapped range [-pi, pi]";
  if (n < 3 || m < 3) error, "Need at least 3x3 array";
  if (!is_void(phi0)) {
    if (rankof(phi0) != 2 || !allof(dimsof(psi) == dimsof(phi0)))
      error, "phi0::Incommensurate dimensions";
  }
  if (!is_void(mask)) {
    if (rankof(mask) != 2 || !allof(dimsof(psi) == dimsof(mask)))
      error, "mask::Incommensurate dimensions";
    write, format="%s\n", "Applying a mask.";
  }

  //Storage for weighting array
  w2= array(0.0, [3,n,m,2]);

  //Phase differences from the data set
  f= _wrap(psi(dif,));
  g= _wrap(psi(,dif));
  f1= f2= g1= g2= array(0.0, dimsof(psi));
  f1(1:n-1,)= f;  f2(2:n,)= f;
  g1(,1:m-1)= g;  g2(,2:m)= g;
  
  //Step 2: Initial guess on phi
  phil= array(0.0, dimsof(psi));
  l= 0;
  
  do {
    
    //Step 3
    R= _wrap(psi - phil);
    //Step 4
    if (_test_resid(R)) {
      //Step 5
      Ruw= _unwrap2d_dct0(R);
      phil+= Ruw;
      UV= _LP_weights(phil);
      write, format=" %s\n", "Converged.";
      return phil;
    }

    //Step 6
    UV= _LP_weights(phil);
    
    //Step 7
    w2= array(0.0, [3,n,m,4]);
    w2(..,1)= UV(..,1);
    w2(2:n,1:m,2)= UV(1:n-1,1:m,1);
    w2(..,3)= UV(..,2);
    w2(1:n,2:m,4)= UV(1:n,1:m-1,2);  
    cij= f1*w2(..,1) - f2*w2(..,2) + g1*w2(..,3) - g2*w2(..,4);
    // cij= array(0.0,n,m);
    // cij(1:n-1,)+=f*UV(1:n-1,,1);
    // cij(2:n,)-=f*UV(2:n,,1);
    // cij(,1:m-1)+=g*UV(,2:m-1,2);
    // cij(,2:m,)-=g*UV(,2:m,2);

    //Step 8
    phil= _unwrap2d_dct2_A(cij, w2, NIT= NIT);
    write, format=" l = %02d\r", l;
    window, 60; fma; pli, phil;
    window, 61; fma; pli, cij;
    
  //Steps 9, 10
  } while (l++ <= NIT_lp);
  write, format=" %s\n", "Done.";

  R= _wrap(psi - phil);
  return phil;
}

func _LP_weights(phi, p, eps0)
/* DOCUMENT _LP_weights, phi, p, eps0

   Returns the weights as prescribed by the Ghiglia-Romero LP-norm
   algorithm.
     
   SEE ALSO:
     _unwrap2D_LPnorm
 */
{
  //extern f, g, n, m, u, v;

  if (is_void(eps0)) eps0= 0.01;
  U= V= array(0.0, n, m);

  //Equations (41) & (42) in JOSA v13p1999 (1996)
  u= abs(phi(dif,) - f);
  v= abs(phi(,dif) - g);

  if (is_void(p) || p == 0.0) {
    U(1:n-1,)= eps0/(u*u + eps0);
    V(,1:m-1)= eps0/(v*v + eps0);
  } else {
    U(1:n-1,)= eps0/(u^(2.-p) + eps0);
    V(,1:m-1)= eps0/(v^(2.-p) + eps0);    
  }
  if (!is_void(mask)) { U*=mask; V*= mask; }
  return [U, V];
}

func _test_resid(phi, &resid)
/* DOCUMENT r= _test_resid(phi, resid)

   Test a phase map for inconsistencies (phase residuals along closed
   paths).  Returns 0 if no inconsistencies exist, 1 if there is at
   least one inconsistency.
     
   SEE ALSO:
     _unwrap2d_LPnorm
 */
{
  local f, g;
  dims= dimsof(phi);
  n= dims(2); m= dims(3);

  //Compute the delta differences
  f= _wrap(phi(dif,));  //= fij
  g= _wrap(phi(,dif));  //= gij

  resid= g(1:n-1,) + f(,2:m) - g(2:n,) - f(,1:m-1);
  return noneof(abs(resid) > 1e-14);
}

func _wrap(x)
/* DOCUMENT phiw= _wrap(phi)

   For an input phase function returns a wrapped phase bounded within
   the interval [-pi, pi].  This is the inverse of the unwrap function
   
   SEE ALSO:
     unwrap
 */
{
  if (numberof(x) > 1) {
    lz= x < 0;
    wl= where(lz);
    wp= where(!lz);
    pm= x;
    if(is_array(wl)) pm(wl)= (x(wl) - pi)%(2*pi) + pi;
    if(is_array(wp)) pm(wp)= (x(wp) + pi)%(2*pi) - pi;
    return pm;
  } else {
    if (x < 0) return (x - pi)%(2*pi) + pi;
    else return (x + pi)%(2*pi) - pi;
  }
}

func _unwrap_test_data(n, smax, kbkg, noise, .., type=)
/* DOCUMENT       _unwrap_test_data, n, smax, kbkg, noise, type=
            -or-  _unwrap_test_data, n, smax, kbkg, noise, sg, type=

   Generates a synthetic phase map for testing unwrapping algorithms,
   consisting of a supergaussian added onto a background phase ramp.

   Parameters are:
     n - numberof pixels per side
     smax - maximum extent of scale
     kbkg - [kx, ky]
     noise - noise level
     sg - [ampl, r0, sg_power]
     
   KEYWORDS:
     type= "sg" [default], "ramp", "ramp+shear"
     
   SEE ALSO:
 */
{

  if (is_void(noise)) noise= 0.0;
  if (is_void(type)) type= "sg";
  p= IMG_DAT();
  p.nx= p.ny= n;
  p.xscale= p.yscale= &span(-smax, smax, n);
  x= (*p.xscale)(,-:1:p.ny);
  y= (*p.yscale)(-:1:p.nx,);
  p.data= &(kbkg(1)*x + kbkg(2)*y);
  p.x_label="Position";
  p.y_label="Position";
  p.z_label="Phase";
  p.x_unit="pixel";
  p.y_unit="pixel";
  p.z_unit="rad";
  
  if (type == "ramp") {
    *p.data+= noise*random_n(dimsof(x));
    *p.data= _wrap(*p.data);
    return p;
  } else if (type =="ramp+shear") {
    xl= where(x < 0.0);
    (*p.data)(xl)= -kbkg(1)*x(xl) -kbkg(2)*y(xl);
    *p.data+= noise*random_n(dimsof(x));
    *p.data= _wrap(*p.data);
    dx= (*p.xscale)(dif)(avg);
    w= array(1.0, dimsof(*p.data));
    w(where(x >= 0.0 & x < 1.5*dx))= 0.0;
    wght= img_copy(p, data= w);
    //return save(phi=p, weight= wght);
    return p;
  } else if (type == "sg") {
    sg= next_arg();
    x2= x*x;
    y2= y*y;
    r2= x2 + y2;
    p.data= &(sg(1)*exp(-(r2/sg(2)/sg(2))^(sg(3)/2.0)));
    *p.data+= kbkg(1)*x + kbkg(2)*y;
    *p.data+= noise*random_n(dimsof(x));
    *p.data= _wrap(*p.data);
    return p;
  }
}


func _unwrap2d_test(..)
/* DOCUMENT _unwrap2d_test

   Perform a suite of unwrapping tests using various unwrapping
   methods.
     
   SEE ALSO:
     unwrap2d
 */
{
  p00= _unwrap_test_data(256, 15.0, 2*pi*[1./5.5, 1./15.5], 0.0, [5.5, 12.0, 9]);
  p05= _unwrap_test_data(256, 15.0, 2*pi*[1./5.5, 1./15.5], 0.5, [5.5, 12.0, 9]);
  p10= _unwrap_test_data(256, 15.0, 2*pi*[1./5.5, 1./15.5], 1.0, [5.5, 12.0, 9]);
  p15= _unwrap_test_data(256, 15.0, 2*pi*[1./5.5, 1./15.5], 1.5, [5.5, 12.0, 9]);
  p20= _unwrap_test_data(256, 15.0, 2*pi*[1./5.5, 1./15.5], 2.0, [5.5, 12.0, 9]);

  mask= array(1.0f, dimsof(img_data(p05)));
  xg= img_grid(p05, 1);
  yg= img_grid(p05,2);
  mask(where(sqrt(xg*xg + yg*yg) < 5.0))= 0.0;

  //Unwrapping with no mask
  sh, 10, unwrap2d(p05, method= 1), notitle= 1;
  pltitle, "GERI unwrapper, noise=0.5 rad";
  sh, 11, unwrap2d(p10, method= 1), notitle= 1;
  pltitle, "GERI unwrapper, noise=1.0 rad";
  sh, 12, unwrap2d(p15, method= 1), notitle= 1;
  pltitle, "GERI unwrapper, noise=1.5 rad";

  sh, 20, unwrap2d(p05, method= 6), notitle= 1;
  pltitle, "Minimum LP-norm unwrapper, noise=0.5 rad";
  sh, 21, unwrap2d(p10, method= 6), notitle= 1;
  pltitle, "Minimum LP-norm unwrapper, noise=1.0 rad";
  sh, 22, unwrap2d(p15, method= 6), notitle= 1;
  pltitle, "Minimum LP-norm unwrapper, noise=1.5 rad";

  //Unwrapping with central mask
  sh, 13, unwrap2d(p05, method= 1, mask= mask), notitle= 1;
  pltitle, "GERI unwrapper, noise=0.5 rad + mask";
  sh, 14, unwrap2d(p10, method= 1, mask= mask), notitle= 1;
  pltitle, "GERI unwrapper, noise=1.0 rad + mask";
  sh, 15, unwrap2d(p15, method= 1, mask= mask), notitle= 1;
  pltitle, "GERI unwrapper, noise=1.5 rad + mask";

  sh, 23, unwrap2d(p05, method= 6, mask= mask), notitle= 1;
  pltitle, "Minimum LP-norm unwrapper, noise=0.5 rad + mask";
  sh, 24, unwrap2d(p10, method= 6, mask= mask), notitle= 1;
  pltitle, "Minimum LP-norm unwrapper, noise=1.0 rad + mask";
  sh, 25, unwrap2d(p15, method= 6, mask= mask), notitle= 1;
  pltitle, "Minimum LP-norm unwrapper, noise=1.5 rad + mask";

  //ramp+shear
  q00= _unwrap_test_data(256, 15.0, 2*pi*[1./5.5, 1./15.5], 0.0, type="ramp+shear");
  q05= _unwrap_test_data(256, 15.0, 2*pi*[1./5.5, 1./15.5], 0.5, type="ramp+shear");
  q10= _unwrap_test_data(256, 15.0, 2*pi*[1./5.5, 1./15.5], 1.0, type="ramp+shear");
  q15= _unwrap_test_data(256, 15.0, 2*pi*[1./5.5, 1./15.5], 1.5, type="ramp+shear");
  q20= _unwrap_test_data(256, 15.0, 2*pi*[1./5.5, 1./15.5], 2.0, type="ramp+shear");

  mask= array(1.0f, dimsof(img_data(q05)));
  xg= img_grid(q05, 1);
  yg= img_grid(q05,2);
  mask(where(sqrt(xg*xg + yg*yg) < 5.0))= 0.0;

  //Unwrapping with no mask
  sh, 30, unwrap2d(q05, method= 1), notitle= 1;
  pltitle, "Ramp/shear, GERI unwrapper, noise=0.5 rad";
  sh, 31, unwrap2d(q10, method= 1), notitle= 1;
  pltitle, "Ramp/shear, GERI unwrapper, noise=1.0 rad";
  sh, 32, unwrap2d(q15, method= 1), notitle= 1;
  pltitle, "Ramp/shear, GERI unwrapper, noise=1.5 rad";

  sh, 40, unwrap2d(q05, method= 6), notitle= 1;
  pltitle, "Ramp/shear, Minimum LP-norm unwrapper, noise=0.5 rad";
  sh, 41, unwrap2d(q10, method= 6), notitle= 1;
  pltitle, "Ramp/shear, Minimum LP-norm unwrapper, noise=1.0 rad";
  sh, 42, unwrap2d(q15, method= 6), notitle= 1;
  pltitle, "Ramp/shear, Minimum LP-norm unwrapper, noise=1.5 rad";

  //Unwrapping with central mask
  sh, 33, unwrap2d(q05, method= 1, mask= mask), notitle= 1;
  pltitle, "Ramp/shear, GERI unwrapper, noise=0.5 rad + mask";
  sh, 34, unwrap2d(q10, method= 1, mask= mask), notitle= 1;
  pltitle, "Ramp/shear, GERI unwrapper, noise=1.0 rad + mask";
  sh, 35, unwrap2d(q15, method= 1, mask= mask), notitle= 1;
  pltitle, "Ramp/shear, GERI unwrapper, noise=1.5 rad + mask";

  sh, 43, unwrap2d(q05, method= 6, mask= mask), notitle= 1;
  pltitle, "Ramp/shear, Minimum LP-norm unwrapper, noise=0.5 rad + mask";
  sh, 44, unwrap2d(q10, method= 6, mask= mask), notitle= 1;
  pltitle, "Ramp/shear, Minimum LP-norm unwrapper, noise=1.0 rad + mask";
  sh, 45, unwrap2d(q15, method= 6, mask= mask), notitle= 1;
  pltitle, "Ramp/shear, Minimum LP-norm unwrapper, noise=1.5 rad + mask";
}
