/*
  DVVMISC.I
  
  Definitions for miscellaneous functions related to dVV/OHRV analysis
  tasks

  P.M. Celliers
  LLNL
  December, 2007

*/

/*

  dvv_Normalize - normalize the magnitude, preserve phase
  dvv_fitPhase - fit phase to a 2D polynomial
  dvv_fitSurf - retuns a fitted surface in an image
  dvv_AzAvg - azimuthal average over a 2D data set
  dvv_FillMask - 
  dvv_Apodize - apply a supergaussian apodization to an image
  dvv_AngledLineout
  dvv_DiffMTF
  dvv_GetPhaseResults
  dvv_PhaseResidues
  dvv_PhaseWeight
  dvv_RippleLinout
  dvv_mkBlankFrameSet
  dvv_TestWarp
  dvv_ModelShock

 */

func dvv_Normalize(img)
/* DOCUMENT dvv_Normalize, img

   Normalizes a complex valued image by setting the magnitude
   of the image array to 1.0, while preserving the phase.
   Zero-valued elements are not changed (undefined phase).

   KEYWORDS:
     
   SEE ALSO:
 */
{
  ni= img_copy(img);
  w0= where(*ni.data == 0);
  ni1= (*ni.data);
  ni1a= abs(*ni.data);
  if(is_array(w0)) ni1a(w0)= 1e-6;
  ni1/= ni1a;
  ni1(w0)= 0+0i;
  ni.data= &ni1
  return ni;
}

// func dvv_unwrap(img) - use img_wrap instead
// /* DOCUMENT dvv_unwrap, img

//    Rough cut phase unwrapping - normalizes all phases in the
//    range from -pi to +pi
     
//    SEE ALSO:
//  */
// {
//   uw= img_copy(img);
//   do {
//     wl= where(*uw.data < -pi);
//     wg= where(*uw.data >= pi);
//     if (numberof(wl) > 0) (*uw.data)(wl)+= 2.0*pi;
//     if (numberof(wg) > 0) (*uw.data)(wg)-= 2.0*pi;
//   } while (numberof(wl) + numberof(wg) > 0);
//   return uw;
// }

//func dvv_bkgPhase(img, reg=, mask=, box=, degree=, weight=, rej_mask=, sel_mask=)
func dvv_fitPhase(img, reg=, mask=, box=, degree=, weight=, rej_mask=, sel_mask=)
/* DOCUMENT bkgPhi= dvv_bkgPhase(img, reg=, box=, degree=, weight=)
     
   SEE ALSO:
 */
{
  //if(is_void(reg)) reg= REGION(x1= -200, x2= 200, y1= -200, y2= 200);
  if(is_void(box)) box= 210;
  if(!is_void(reg)) {
    X= img_extract(img, reg); 
    if(!is_void(weight)) Xw= img_extract(weight, reg);
  } else {
    X= img_copy(img);
    if(!is_void(weight)) Xw= img_copy(weight);
  }
  
  if(!is_void(rej_mask)) {
    rmsk= load_filter(rej_mask);
    //print, "rmsk = ", rmsk;
    mrej= img_GenMaskSet(X, rmsk);
  }
  
  if(!is_void(sel_mask)) {
    smsk= load_filter(sel_mask);
    //print, "smsk = ", smsk;
    msel= img_GenMaskSet(X, smsk);
  }
  
  if(is_void(msel) && !is_void(mrej)) m= mrej;
  else if (!is_void(msel) && is_void(mrej)) m= msel;
  else if (!is_void(msel) && !is_void(mrej)) m= img_CombineMasks(msel, mrej);
  if(!is_void(m)) img_ApplyMask, X, m, maskarr;
  
  //info, maskarr; info, m;
  //img_ShowMask, X, m;
  Xuw= unwrap2d(X, box, weight= Xw, mask= maskarr);
  //img_display, Xuw, w= 20;
  xg= (*X.xscale)(,-:1:X.ny);
  yg= (*X.yscale)(-:1:X.nx,);
  dd= *Xuw.data;
  //window, 21; pli, dd;
  ddm= array(1n, dimsof(dd));
  if(!is_void(m)) ddm(m)= 0n;
  //window, 22; pli, dd*ddm;
  af= fitsurf(*Xuw.data, yg, xg, degree= degree, mask= where(ddm));
  return dvv_mkSurf(img, af, type= "polysurf");
}

func dvv_fitSurf(img, reg=, degree=)
/* DOCUMENT s= dvv_fitSurf(img, &af, reg=, degree=)

   Returns an IMG_DAT data structure containing a 
     
   SEE ALSO:
 */
{
  if(is_void(degree)) degree= 2;
  if(is_void(reg)) reg= REGION(x1= -350, x2= 350, y1= -350, y2= 350);
  X= img_extract(img, reg);
  xg= (*X.xscale)(,-:1:X.ny);
  yg= (*X.yscale)(-:1:X.nx,);
  af= fitsurf(*X.data, yg, xg, degree= degree);
  return dvv_mkSurf(img, af, type= "polysurf");
}

func dvv_AzAvg(ft, mask=)
/* DOCUMENT dvv_AzAvg, ft, mask=

    Computes an azimuthal average over a 2D data set;
    assumes that the data set is centered on the origin.
    The resulting vector is
        z(r) = integral_(-pi)^(+pi) z(r,phi) dphi
    is returned as a N x 2 element vector: [r, z]

   KEYWORDS:
    mask=
     
   SEE ALSO:
 */
{
  mg= zg= sg= rg= (*ft.xscale)(where(*ft.xscale >= 0.0));
  nr= numberof(rg);
  xg= (*ft.xscale)(,-:1:ft.ny);
  yg= (*ft.yscale)(-:1:ft.nx,);
  // Masked (i.e. filtered) spectrum
  if(!is_void(mask)) {
    msk= array(1.0, dimsof(*ft.data));
    msk(mask)= 0.0;
  }
  for(i= 1; i<= nr; i++) {
    //if(i%100 == 0) print, "i= ", i;
    yv= span(-pi, pi, int(2*pi*i)+1);
    yc= rg(i)*sin(yv(1:-1));
    xc= rg(i)*cos(yv(1:-1));
    // For masked spectra computed a weighted average
    if(!is_void(mask)) {
      mgd= interp2(yc, xc, msk, yg, xg);
      zgd= interp2(yc, xc, (*ft.data), yg, xg);
      zg(i)= (mgd * zgd)(avg);
      mg(i)= mgd(avg);
      wm= where(mgd == 1.0);
      if(numberof(wm) > 0) sg(i)= zgd(wm)(rms);
      else sg(i)= 0.0;
    } else {
      zgd= interp2(yc, xc, (*ft.data), yg, xg);
      zg(i)= zgd(avg);
      sg(i)= zgd(rms);
    }
  }
  if(!is_void(mask)) return [rg, zg/mg, sg];
  else return [rg, zg, sg];
}

func dvv_FillMask(img, az, mask)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  ft= dvv_fft2D(ft);
  xsc= *ft.xscale;
  ysc= *ft.yscale;
  ix= (mask - 1)%ft.nx + 1;
  iy= (mask - 1)/ft.ny + 1;
  N= numberof(mask);

  xx= xsc(ix);
  yy= ysc(iy);
  rr= sqrt(xx*xx + yy*yy);
  I= 0.0+1i;

  dt= *ft.data;
  fill= interp(az(,2), az(,1), rr(1:N/2));
  phs= exp(I*2*pi*random(N/2) - pi);
  fill= fill*phs;
  dt(ix(1:N/2), iy(1:N/2))= fill;
  fill.im= -fill.im;
  dt(ix(N/2+1:0)(0:1:-1), iy(N/2+1:0)(0:1:-1))= fill;
  
  ft.data= &dt;
  
  return ftn;
}

func dvv_Apodize(img, x0, y0, R, order=)
/* DOCUMENT dvv_Apodize, img, x0, y0, R, order=

   Apodizes a region of an image data set by multiplying
   the data by a supergaussian function (default is 6th order).
   Returns the apodized image.

   KEYWORDS:
     order=  defines the supergaussian order (default = 6)
     
   SEE ALSO:
   
 */
{
  ni= img_copy(img);
  if(is_void(order)) order= 6;
  xg= (*img.xscale)(,-:1:img.ny);
  yg= (*img.yscale)(-:1:img.nx,);
  r2= (xg - x0)*(xg - x0) + (yg - y0)*(yg - y0);
  sg= exp(-(r2/R/R)^(order/2));
  ni.data= &((*img.data)*sg);
  return ni;
}

func dvv_AngledLineout(img, pt, theta, lwidth=, w=, wlo=, marks=, width=, color=)
/* DOCUMENT dvv_AngledLineout, img, pt, theta, lwidth=, w=, wlo=, marks=, width=, color=





   SEE ALSO:
 */
{
  //print, ">>>>>>> dvv_AngledLineout ...";
  if(is_void(pt)) pt= [0.0, 0.0];
  if(numberof(theta) == 2) theta= 180.0*atan(theta(2) - pt(2), theta(1) - pt(1))/pi;
  nx= img.nx;
  dx= (*img.xscale)(dif)(avg);
  if(is_void(lwidth)) lwidth= dx;
  n= 1 + int(lwidth/dx);
  theta= theta/180.0*pi;
  //print, "n= ", n, "dx= ", dx, "lwidth= ", lwidth;
  a= span(-lwidth/2.0, lwidth/2.0, n);
  //print, a;
  if(is_void(w)) w= 0;
  s= *img.xscale;
  xv= s*cos(theta);
  yv= s*sin(theta);
  //yv= sqrt(s*s - xv*xv)*sign(xv);
  //window, 0; plg, yv, xv, color= "red";
  xv= xv(,-:1:n) + (pt(1) - (a*sin(theta))(-:1:nx,));
  yv= yv(,-:1:n) + (pt(2) + (a*cos(theta))(-:1:nx,));
  //print, transpose([xv,yv]);
  //ldtab= tl2cub(,*img.xscale, *img.yscale, *img.data);
  //ldt= tl2cub(ldtab, xv, yv);
  ldt= interp2d(yv, xv, *img.data, *img.yscale, *img.xscale, fill= 0.0);
  window, w, legends= 0;
  plg, ldt(,avg), s(,1), marks= marks, width= width, color= color;
  if(!is_void(wlo)) {
    window, wlo; plsys, 1;
    plg, yv(,1), xv(,1), color= "red", marks= 0;
    plg, yv(,0), xv(,0), color= "blue", marks= 0;
    //plg, [yv(1,1), yv(1,0)], [xv(1,1), xv(1,0)], color= "red", marks= 0;
    //plg,[yv(0,1), yv(0,0)], [xv(0,1), xv(0,0)], color= "red", marks= 0;
  }
  if(!am_subroutine()) return [s(,1), ldt(,avg)];
}

func dvv_DiffMTF(f, f0)
/* DOCUMENT dvv_DiffMTF, f, f0

   Computes the MTF function for a well-corrected, diffraction-limited
   lens as a function of spatial frequency f.  The cutoff frequency
   f0 is related to the F# of the lens by

   f0 = 1/(2*lambda*F#)

   e.g. F/3 at 0.5 microns gives f0 = 333 [1/mm units]; the maximum
   frequency that can be detected with such a lens is 2*f0 (667 [1/mm]).
     
   SEE ALSO:
 */
{
  y= f;
  v= f < 2*f0 & f >= 0.0;
  wl= where(v);
  wg= where(!v);
  f2= f(wl)/2./f0;
  y(wl)= (2/pi)*(acos(f2) - f2*sqrt(1 - f2*f2));
  if(numberof(wg) > 0) y(wg)= 1e-6;
  return y;
}

func dvv_GetPhaseResults(phase_set)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern dtp, dtpUW, ampl, vel;
  extern rfp, rfpUW;
  extern dtV, rfV, dtVaz, rfVaz;
  extern MLO, HLO, VLO;

  p= dvv_PhaseSet();

  p.MLO= &MLO;
  p.HLO= &HLO;
  p.VLO= &VLO;
  p.dtVaz= &dtVaz;
  p.rfVaz= &rfVaz;

  p.dtp= dtp;
  p.dtpUW= dtpUW;
  p.rfp= rfp;
  p.rfpUW= rfpUW;
  p.ampl= ampl;
  p.vel= vel;
  p.dtV= dtV;
  
  return p;
}

func dvv_PhaseResidues(phi)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  r= img_copy(phi);
  r.nx-= 1;
  r.ny-= 1;
  r.xscale= &((*r.xscale)(zcen));
  r.yscale= &((*r.yscale)(zcen));
  
  vd= (*phi.data)(,dif);
  hd= (*phi.data)(dif,);
  wvg= where(vd >  pi);
  wvl= where(vd < -pi);
  whg= where(hd >  pi);
  whl= where(hd < -pi);
  vd1= array(0.0, dimsof(vd));
  hd1= array(0.0, dimsof(hd));
  
  vd1(wvg)= 1.;
  vd1(wvl)= -1.;
  hd1(whg)= 1.;
  hd1(whl)= -1.;

  r.data= &(vd1(dif,) - hd1(,dif));
  return r;
}

func dvv_PhaseWeight(dt, k)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _CSPHI1, _SNPHI1;

  if(is_void(k)) k= 5;
  lga= dvv_Ampl(dt);
  wl= where(*lga.data <= 0.0);
  if(numberof(wl) > 0) (*lga.data)(wl)= 1e-6;
  dx= (*lga.xscale)(dif)(avg);
  lga.data= &log(*lga.data);
  bx= dvv_mkBox(lga, (k+0.01)*dx);
  cm= dvv_Correlate(img_mul(lga, _CSPHI1), bx);
  sm= dvv_Correlate(img_mul(lga, _SNPHI1), bx);
  ss= img_add(img_mul(cm, cm), img_mul(sm, sm));
  sm.data= &(sqrt(*ss.data)/k/k);
  return sm;
}

func dvv_RippleLineout(pUW, mode, box, radius=, w=, wlo=, norm=, shot=, rge=,
                       lwidth=, marks=, width=)
/* DOCUMENT dvv_RippleLineout, pUW, mode, box, w=, wlo=, radius=, norm=, shot=,
                rge=, lwidth=, marks=, width=
     
   SEE ALSO:
 */
{
  if(is_void(norm)) norm= 1;
  if(is_void(shot)) shot= "";
  if(is_void(radius)) radius= 0.005;
  if(is_void(w)) w= 0;
  if(is_void(wlo)) w= 1;
  
  V= dvv_Vel(img_extract(pUW, box), scale= 1000.0);
  Vamp= dvv_AMP(V, norm= norm, wndw= "hanning", nodc= 1);
  if(norm == 1) Vamp.z_unit= "m/s";
  else Vamp.z_unit= "m/s-!mm"
  Vpsd= dvv_PSD(V, norm= norm, wndw= "hanning", nodc= 1);
  mode1= dvv_GenFilterMask(V, _(mode, radius), type= "mode_select");
  mode_amplitude= dvv_PSD_rms(Vpsd, norm= norm, mask= mode1);
  write, format="Mode amplitude: %12.5e rms, %12.5e %s peak-to-peak\n",
    mode_amplitude, sqrt(2)*2*mode_amplitude, V.z_unit;
  if(is_void(rge)) dvv_DisplayImg, Vamp, w= wlo;
  else dvv_DisplayImg, Vamp, rge(1), rge(2), w= wlo;
  plsys, 1;
  dvv_plCircle, mode, radius, color= "red", width= 2, nolimits= 1;
  dvv_plCircle, -mode, radius, color= "red", width= 2, nolimits= 1;
  lo1= dvv_AngledLineout(Vamp, ,mode, wlo= wlo, w= w, lwidth= lwidth, marks= marks, width= width);
  window, w;
  logxy, 1, 0;
  xytitles, "Spatial frequency ("+Vamp.x_unit+")",
    "Velocity spectral density ("+Vamp.z_unit+")", [-0.01,0.0];
  if(is_void(rge)) limits, 0.05, 0.5;
  else limits, 0.01, 0.5, rge(1), rge(2);
  f= open(shot+"_angle_lineout.txt", "w");
  write, f, format="# Angled lineout through a pre-imposed mode for shot %s\n", shot;
  write, f, format="#  Mode is [%6.4f, %6.4f], selection radius, %6.4f\n", mode(1), mode(2), radius;
  write, f, format="#  Mode wavelentgh is %5.2f\n", 1.0/sqrt((mode^2)(sum));
  write, f, format="#  Sampling box is (x1= %6.1f, x2= %6.1f, y1= %6.1f, y2= %6.1f)\n",
            box.x1,box.x2,box.y1,box.y2;
  write, f, format="#  Mode amplitude: %12.5e rms, %12.5e %s peak-to-peak\n",
            mode_amplitude, sqrt(2)*2*mode_amplitude, V.z_unit;
  write, f, format="# Spatial frequency (%s)  Mode amplitude (%s)\n", Vamp.x_unit, Vamp.z_unit;
  write, f, lo1(,1), lo1(,2);
  close, f;
}

func dvv_mkBlankFrameSet(&fs0, rf=, fm=, warp_delta=, phi_delta=, dwmethod=, data=, n=, anoise=)
/* DOCUMENT dvv_mkBlankFrameSet, rf=, fm=, warp_delta=, phi_delta=, dwmethod=, data=, n=, anoise=

   SEE ALSO:
     dvv_GenerateFrameSet
*/
{
  extern _phi_data;

  if (is_void(anoise)) anoise= 0.0;
  if (is_void(n)) n= 2100;
  if (is_void(fm)) fm= [0.025, -0.03];
  if (is_void(phi_delta)) delta= pi/2; else delta= phi_delta;
  if (is_void(warp_delta)) warp_delta= 0.0;
  if (is_void(order)) order= 3;
  if (is_void(dwmethod)) dwmethod= "polywarp";
  fs= dvv_FrameSet();

  mapping= dvv_SetMapping();

  gc= array(POINT, 5);
  gc.x= [0.0,-300,-300,300,300];
  gc.y= [0.0,-300,300,-300,300];

  pts= ptsf= array(POINT, 81);
  pts.x= reform((span(-4, 4, 9)*100)(,-:1:9), [1,81]);
  pts.y= reform((span(-4, 4, 9)*100)(-:1:9,), [1,81]);
  mapping.warp_refpts= &pts;

  mapping.ch1S= dvv_FrameMapping(gridList= &gc, angle= 0.0, sclx= 1.0, scly= 1.0);
  mapping.ch2S= dvv_FrameMapping(gridList= &gc, angle= 0.0, sclx= 1.0, scly= 1.0);
  mapping.ch1P= dvv_FrameMapping(gridList= &gc, angle= 0.0, sclx= 1.0, scly= 1.0);
  mapping.ch2P= dvv_FrameMapping(gridList= &gc, angle= 0.0, sclx= 1.0, scly= 1.0);

  ptsf= pts;
  ptsf.x+= random_n(dimsof(pts))*warp_delta;
  ptsf.y+= random_n(dimsof(pts))*warp_delta;
  polywarp, ptsf.x, ptsf.y, pts.x, pts.y, order, kx1s, ky1s;
  mapping.ch1S.warp_pts= &ptsf;
  mapping.ch1S.kx= &kx1s;
  mapping.ch1S.ky= &ky1s;

  ptsf= pts;
  ptsf.x+= random_n(dimsof(pts))*warp_delta;
  ptsf.y+= random_n(dimsof(pts))*warp_delta;
  polywarp, ptsf.x, ptsf.y, pts.x, pts.y, order, kx2s, ky2s;
  mapping.ch2S.warp_pts= &ptsf;
  mapping.ch2S.kx= &kx2s;
  mapping.ch2S.ky= &ky2s;

  ptsf= pts;
  ptsf.x+= random_n(dimsof(pts))*warp_delta;
  ptsf.y+= random_n(dimsof(pts))*warp_delta;
  polywarp, ptsf.x, ptsf.y, pts.x, pts.y, order, kx1p, ky1p;
  mapping.ch1P.warp_pts= &ptsf;
  mapping.ch1P.kx= &kx1p;
  mapping.ch1P.ky= &ky1p;

  ptsf= pts;
  ptsf.x+= random_n(dimsof(pts))*warp_delta;
  ptsf.y+= random_n(dimsof(pts))*warp_delta;
  polywarp, ptsf.x, ptsf.y, pts.x, pts.y, order, kx2p, ky2p;
  mapping.ch2P.warp_pts= &ptsf;
  mapping.ch2P.kx= &kx2p;
  mapping.ch2P.ky= &ky2p;

  fs.map= mapping;

  xscl= yscl= span(-445.0, 445.0, n);
  tmp= img_new(xscl, yscl, array(0.0, n, n), xlabel= "Position",
	       ylabel= "Position", zlabel= "Signal", xunit= "!mm",yunit= "!mm",
               zunit= "counts", title= "Test"); 

  xg= img_grid(tmp, 1);
  yg= img_grid(tmp, 2);

  p= 2*pi*(fm(1)*xg + fm(2)*yg);
  if (data) {
    phi= img_data(dvv_ModelShock(tmp, 1.83, 21.0));
  } else phi= 0.0;

  if (rf) rfdat= img_data(dvv_frameSum(rf, average= 1)); else rfdat= 1.0;
  
  fs.ch1S= img_copy(tmp, data= rfdat*(1 + 0.5*cos(p + phi)) + anoise*(random(dimsof(p)) - 0.5));
  fs.ch2S= img_copy(tmp, data= rfdat*(1 + 0.5*cos(p + phi + pi)) + anoise*(random(dimsof(p)) - 0.5));
  fs.ch1P= img_copy(tmp, data= rfdat*(1 + 0.5*cos(p + phi + delta)) + anoise*(random(dimsof(p)) - 0.5));
  fs.ch2P= img_copy(tmp, data= rfdat*(1 + 0.5*cos(p + phi + delta + pi)) + anoise*(random(dimsof(p)) - 0.5));
  //fs.ch1S= img_copy(tmp, data= rfdat*(1 + 0.5*sin(p + phi)));
  //fs.ch2S= img_copy(tmp, data= rfdat*(1 - 0.5*sin(p + phi)));
  //fs.ch1P= img_copy(tmp, data= rfdat*(1 - 0.5*sin(p + phi + delta)));
  //fs.ch2P= img_copy(tmp, data= rfdat*(1 + 0.5*sin(p + phi + delta)));
  _phi_data= p+phi;
  print, "dvv_mkBlankFrameSet - fscopy";  
  fs0= dvv_FSCopy(fs);

  fs.ch1S.shotid+= " - ch1S";
  fs.ch2S.shotid+= " - ch2S";
  fs.ch1P.shotid+= " - ch1P";
  fs.ch2P.shotid+= " - ch2P";
  
  if (warp_delta) dvv_ApplyWarp, fs, dwmethod= dwmethod;
  
  return fs;
}

func dvv_TestWarp(.., rf=, dwmethod=, phi_delta=, warp_delta=, fm=, rge=, prlo=, title=, domaps=)
/* DOCUMENT dvv_TestWarp, rf=, dwmethod=, phi_delta=, warp_delta=, fm=, rge=, prlo=, title=, domaps=
     
   SEE ALSO:
 */
{
  //extern shot;
  extern _center, _box;
  extern _DVV_AMPL_MAP, _DVV_PHASE_MAP;
  extern pV, pVaz, tspUW, ts, ts0, tspbkg;

  if (is_void(_center)) _center= [0.0,0.0];
  if (is_void(_box) || _box < 50) _box= 100;
  if (is_void(rge)) rge= [0,200];
  if (is_void(title)) title="";
  if (is_void(phi_delta)) phi_delta= pi/2;
  
  shot= strtr(title,'_','-');

  ts= dvv_mkBlankFrameSet(ts0, rf= rf, dwmethod= dwmethod, phi_delta= phi_delta, warp_delta= warp_delta, fm= fm);

  if (domaps) {
    dvv_SetPA_Maps, ts, phiexpect= phi_delta, method= 2;
  } else {
    _DVV_AMPL_MAP= img_copy(_DVV_AMPL_MAP, data= 1.0);
    _DVV_PHASE_MAP= img_copy(_DVV_PHASE_MAP, data= phi_delta);
  }
  
  tsp= unwrap2d(dvv_Phase(ts, [pi/2, 1.0, 0.0, 0.0], method= 1, invert= 1));
  tspbkg= dvv_fitPhase(tsp, reg= 400*[-1,1,-1,1], degree= 3);
  tspUW= unwrap2d(img_sub(tsp, tspbkg));
  //dvv_DisplaySet, ts, w= next_window(); palette, "yarg.gp";
  dvv_ShowSpectV, tspUW, pV, pVaz, w= next_window(), bxctr= _center, az0= az0, prlo= prlo,
      box= REGION(x1=-_box/2, x2= _box/2, y1=-_box/2, y2= _box/2), rge= rge, norm= 2,
      nodc= nodc, opfilt= opfilt;
  palette, "yarg.gp";
  wpdf, current_window(), "TestWarp-2D-"+title;

  window, next_window();
  plg, pVaz(,2), pVaz(,1), color= "red", width= 3, marks= 0;
  xytitles, "Spatial Frequency (!mm^-1^)", "Velocity density (m/s-!mm)";
  logxy, 1, 1;
  gridxy, 1, 1;
  limits, 0.005, 0.5, 10, 500;
  pltitle, shot;
  wpdf, current_window(), "TestWarp-1D-"+title;
}

func dvv_ModelShock(img, vpf, vel, R0, xy00=)
{
  if (is_void(R0)) R0= 350.0;
  if (is_void(xy00)) {
    x0= y0= 0.0;
  } else {
    x0= xy00(1); y0= xy00(2);
  }
  xg= img_grid(img, 1);
  yg= img_grid(img, 2);
  rg= sqrt(xg*xg + yg*yg);
  phi= 2*pi*vel*(1 - ((xg-x0)^4 + (yg-y0)^4 + 2*(xg-x0)^2*(yg-y0)^2)/R0^4)/vpf;
  wg= where(rg > 350);
  if (is_array(wg)) phi(wg)= 0.0;
  return img_copy(img, data= phi);
}

func dvv_plQuadComponents(A, phi, sd, pd, x, y, invert=)
{
  Av= img_zvalue(A, y, x);
  pv= img_zvalue(phi, y, x);
  if (invert) pv*= -1;
  ar= img_zvalue(_DVV_AMPL_MAP, y, x);
  dd= img_zvalue(_DVV_PHASE_MAP, y, x);
  sdv= img_zvalue(sd, y, x);
  pdv= img_zvalue(pd, y, x);
  sd1= ar*sdv;
  pd1= (pdv - ar*sdv*cos(dd))/sin(dd);
  plg, [Av*sin(pv), 0.0], [Av*cos(pv), 0.0], width= 3, marks= 0, type= 1, color= "black";
  plmk, Av*sin(pv), Av*cos(pv), msize= 0.5, width= 10, marker= 5, color= "magenta";
  plg, [pdv, pdv], [-50000, 50000], width= 3, type= 2, marks= 0, color= "blue";
  plg, [-50000, 50000], [sdv, sdv], width= 3, type= 2, marks= 0, color= "red";
  plg, [pd1, pd1], [-50000, 50000], width= 3, type= 4, marks= 0, color= "blue";
  plg, [-50000, 50000], [sd1, sd1], width= 3, type= 4, marks= 0, color= "red";
  //limits, -1.1*Av, 1.1*Av, -1.1*Av, 1.1*Av;
  limits, sd1-0.1*Av, sd1+0.1*Av, pd1-0.1*Av, pd1+0.1*Av;
  plt, swrite(format="a= %4.2f", sqrt(sd1*sd1 + pd1*pd1)/Av), 0.44, 0.89, tosys= 0;
}

func dvv_QuadComponents(fr)
{
  sdiff= img_sub(fr.ch2S, fr.ch1S);
  pdiff= img_sub(fr.ch2P, fr.ch1P);
  csd= img_operate(_DVV_PHASE_MAP, cos);
  snd= img_operate(_DVV_PHASE_MAP, sin);
  sd1= img_mul(sdiff, _DVV_AMPL_MAP);
  pd1= img_div(img_add(pdiff, img_mul(sdiff, _DVV_AMPL_MAP, csd)), snd);
  ampl= img_operate(img_add(img_mul(sd1, sd1), img_mul(pd1, pd1)), sqrt);
  ampl0= img_operate(img_add(img_mul(sdiff, sdiff), img_mul(pdiff, pdiff)), sqrt);
  return save(csphi= sd1, snphi= pd1, ampl= ampl, ampl0= ampl0);
}
