/*
  DVVSPECTRUM.I
  
  Definitions for routines to perform spectral operations on dVV/OHRV
  data

  P.M. Celliers
  LLNL
  December, 2007

*/

/*
  Provides the following functions:

    dvv_PSD        - compute 2D power spectrum
    dvv_PSD_rms    - extract RMS measure from 2D power spectrum
    dvv_PSD_1D     - collapse 2D power spectrum to 1D representation by summing
    dvv_PSD_1D_rms -
    
    dvv_AMP        - computes the SQRT of the 2D power spectrum
    dvv_AMP_rms    - extract RMS measure from 2D rms amplitude spectrum
    dvv_AMP_1D     - collapse 2D amplitude spectrum to 1D representation by summing
    dvv_AMP_1D_rms -

    dvv_rms - 

    dvv_fft2D      - compute complex-valued 2D fft (forwards & backwards)
    dvv_Filter     - 
    dvv_Mode       -
    dvv_GenFilterMask -
    dvv_GenFilterMaskSet -
    dvv_ShowFilterMask -
    dvv_CombineMasks -
    dvv_ORMasks    -
    mode_vec
    dvv_mkModes
    dvv_fromCX
    dvv_mkNoiseModel
    dvv_OHRVNoiseModel
    dvv_White

*/

func dvv_PSD(img, p=, wndw=, nodc=, norm=)
/* DOCUMENT dvv_PSD, img, p=, wndw=, nodc=, norm=

   Computes the fft power spectrum for the data set contained in img.
   The power spectrum is normalized such that an appropriate sum of
   the non-zero frequency modes matches the rms of the input data set.

   See section 13.4 in "Numerical Recipes in C"

   There are two normalizations available for the Fourier
   coefficients:

     norm= 1(default): "Amplitude^2 per Fourier Mode": with this
       normalization the RMS = SQRT(Sum over all PSD coefficients,
       excluding DC term). This means that when the sampling interval
       is changed, so will the PSD coefficients.

     norm=2: "Amplitude^2 per unit frequency^2": with this normalization
       the "box size" of the sample is factored into the coefficients.
       This will produce, on average, the same spectral density
       independent of the sampling details. To get the RMS one must
       average over the spectrum; in this process the frequency binning
       units have to be divided out (i.e. the box size parameters
       get divided out):

       RMS=SQRT(Average_over_all_PSD_coefficients_except_DC_Term * Lx*Ly) 

   The function dvv_PSD_rms calculates the correct RMS for each type of
   normalization.

   KEYWORDS:
     wndw= "hanning", "hamming", "bartlett", "triangular", "edge7%"
     p=      fftw plan
     nodc= [] (void), no adjustment of background or DC offset
           0 - sets the dc average to zero, useful for data sets
              with very large DC offsets
           1,2,3 - subtracts n-order polynomial fit to the data before
              taking the transform (removes the lowest order modes).
     norm=  Normalization:
              1) PSD_k,l = |C_k,l|^2 / (nx^2*ny^2)
              2) PSD_k,l |C_k,l|^2 * Lx*Ly / (nx*ny)
   
   SEE ALSO:
     dvv_PSD_rms

*/
{
  extern _NO_FFTW;
  local ft, dat;
  ft= img_copy(img);
  wfac= 1.0;

  dat= *img.data;
  if (!is_void(nodc)) {
    if (nodc == 0) {
      dat-= dat(*)(avg);  /* Subtract DC offset */
    } else if (nodc > 0 && nodc < 4) {
      fit= img_fitsurf(img, degree= nodc);
      dat-= img_data(fit);  /* Subtract low order polynomial fit (< 4th order) */
    } else error, "nodc= keyword should be void, or 0, 1, 2 or 3"; 
  }
  
  nxy= double(ft.nx)*double(ft.ny);

  //Factor in the window function and compute its normalization factor
  if(!is_void(wndw)) {
    xg= span(-pi, pi, img.nx)(,-:1:img.ny);
    yg= span(-pi, pi, img.ny)(-:1:img.nx,);
    if(wndw == "hanning") {
      wndw= 0.25*(1.0 + cos(xg))*(1.0 + cos(yg));
      wfac= (wndw*wndw)(*)(sum);
      dat*= wndw;
    } else if (wndw == "hamming") {
      wndw= (0.54 + 0.46*cos(xg))*(0.54 + 0.46*cos(yg));
      wfac= (wndw*wndw)(*)(sum);
      dat*= wndw;
    } else if (wndw == "triangular") {
      wndw= (1.0 - abs(xg)/pi)*(1.0 - abs(yg)/pi);
      wfac= (wndw*wndw)(*)(sum);
      dat*= wndw;
    } else if (wndw == "edge7%") {
      ww= (abs(xg) > 0.93*pi) | (abs(yg) > 0.93*pi);
      wndw= array(1.0, dimsof(xg));
      wndw(where(ww))= 0.0;
      wfac= (wndw*wndw)(*)(sum);
      dat*= wndw;
    } else error, "Unknown window function";
  } else {
    wfac= nxy;
  }
  //print, "wfac=", wfac;

  ft.data= &dat;
  //dvv_DisplayImg, ft, w= 0;

  if(_NO_FFTW) {
    ft.data= &roll(abs(fft(dat, +1)));
  } else {  
    //if(is_void(p)) p= fftw_plan(dimsof(dat), +1);
    //ft.data= &roll(abs(fftw(dat, p)));
    //if(is_void(p)) p= fftw_plan(dimsof(dat), +1);
    ft.data= &roll(abs(fftw(dat, +1)));
  }
  dx= (*img.xscale)(dif)(avg);
  dy= (*img.yscale)(dif)(avg);
  ft.xscale= &roll(fft_indgen(ft.nx)/double(ft.nx)/dx);
  ft.yscale= &roll(fft_indgen(ft.ny)/double(ft.ny)/dy);
  if (!(ft.nx%2)) (*ft.xscale)(1)= - abs((*ft.xscale)(1));
  if (!(ft.ny%2)) (*ft.yscale)(1)= - abs((*ft.yscale)(1));
  ft.x_unit= ft.x_unit+"^-1^";
  ft.y_unit= ft.y_unit+"^-1^";
  ft.x_label= "frequency";
  ft.y_label= "frequency";
  
  //Normalize the fft according to the mesh size,
  //window function and normalization type
  if(is_void(norm) || norm == 1) {
    //Normalization type 1, power per mode
    *ft.data*= *ft.data;
    *ft.data/= wfac*nxy;
    //ft.z_unit= ft.z_unit+"^2^";
  } else if (norm == 2) {
    //Normalization type 2, power per unit frequency
    //requires box size Lx, Ly
    Lx= (*img.xscale)(0) - (*img.xscale)(1);
    Ly= (*img.yscale)(0) - (*img.yscale)(1);
    //ft.z_unit= img.z_unit+"^2^ "+"-"+img.x_unit+"^2^";
    *ft.data*= *ft.data*Lx*Ly/wfac/nxy;
  }
  return ft;
}

func dvv_rms(img, fl, fu, mask=)
/* DOCUMENT  rms_sum= dvv_rms(img)
        -or- rms_sum= dvv_rms(img, mask= mask)
        -or- rms_sum= dvv_rms(img, fl, fu, mask= mask)
        -or- rms_sum= dvv_rms(img, [fl, fu], mask= mask)

   Evaluates an rms sum over a 2D image data set. In the first calling
   version rms is calculated directly from the image.  A filter mask
   may be applied to filter out parts of the spectrum from being
   included in the sum additional modes. In addition, band pass
   filtering limits given by fl and fu can be specified. If both the
   mask and band pass limits are specified, the mask= is OR-ed with the
   bandpass filter defined by the limits, f[l, fu].

   The input data set must be a 2D IMG_DAT() data structure containing
   a real-space image.

   KEYWORDS: mask= use a filter mask defined by dvv_GenFilterMask
     before evaluating the rms
     
   SEE ALSO:
     dvv_GenFilterMask, dvv_Filter
 */
{
  //write, format="dvv__rms: %s\n", " ";
  if(!is_void(fl)) {
    if(numberof(fl) == 2 && is_void(fu)) {
      fu= fl(2); fl= fl(1);
    } else if (is_void(fu)) error, "Incorrectly specified frequency limits.";
    mask1= dvv_GenFilterMask(img, [fl, fu], type= "band_pass");
    if(is_void(mask)) mask= mask1;
    else mask= dvv_ORmasks(mask, mask1);
  }
  
  if(is_void(mask)) return (*img.data)(*)(rms);
  else {
    flt= dvv_Filter(img, mask);
    return (*flt.data)(*)(rms);
  }
}

func dvv_PSD_rms(psd, fl, fu, norm=, mask=)
/* DOCUMENT  rms_sum= dvv_PSD_rms(psd, norm= norm)
        -or- rms_sum= dvv_PSD_rms(psd, fl, fu, norm= norm, mask= mask)
        -or- rms_sum= dvv_PSD_rms(psd, [fl, fu], norm= norm, mask= mask)
        -or- rms_sum= dvv_PSD_rms(psd, norm= norm, mask= mask)

   Evaluates an rms sum over a 2D power spectrum. The rms sum is
   calculated by summing over all of the modes within the frequency
   limits given by fl and fu. In addition a filter mask can be applied
   to filter out other parts of the spectrum from being included in
   the sum additional modes.  If the mask is specified, it is OR-ed
   with the bandpass filter defined by the limits, [fl, fu]

   The input power spectrum must be a 2D IMG_DAT() structure as
   returned by dvv_PSD.

   The result depends on the normalization used in computing the PSD.
   The same normalization (=1 or =2) must be specified as was used
   when computing the psd data.

   KEYWORDS:  
     norm= Normalization:
           1) |C_k|^2 / (nx^2*ny^2)
           2) |C_k|^2 *Lx*Ly / (nx*ny)^2
           
     mask= define this with a filter mask which is applied to the data
         before the summation; for extracting PSD over a subset
         of the frequency spectrum
            
   SEE ALSO:
     dvv_PSD, dvv_PSD_1D, dvv_GenFilterMask
 */
{
  //write, format="dvv_PSD_rms: norm= %d\n", norm;
  if(!is_void(fl)) {
    if(numberof(fl) == 2 && is_void(fu)) {
      fu= fl(2); fl= fl(1);
    } else if (is_void(fu)) error, "Incorrectly specified frequency limits.";
    //write, format="dvv_PSD_rms: fl= %e, fu= %e\n", fl, fu;
    xg= (*psd.xscale)(,-:1:psd.ny);
    yg= (*psd.yscale)(-:1:psd.nx,);
    r= sqrt(xg*xg + yg*yg);
    mask1= where(r < fl | r > fu);
    if(is_void(mask)) mask= mask1;
    else mask= dvv_ORmasks(mask, mask1);
    //info, mask;
  }

  if(is_void(norm) || norm == 1) {
    //Normalization type == 1, sum over all modes except DC
    dat= *psd.data;
    dat(psd.nx/2+1,psd.ny/2+1)= 0.0;
    if(!is_void(mask) && numberof(mask) > 0) dat(mask)= 0.0;
    return sqrt(dat(*)(sum));   
  } else if(norm == 2) {
    //Normalization type == 2, average over the PSD
    //and divide by the box area
    Lx= 1./(*psd.xscale)(dif)(avg)*(psd.nx-1)/psd.nx;
    Ly= 1./(*psd.yscale)(dif)(avg)*(psd.ny-1)/psd.ny;
    dat= *psd.data;
    dat(psd.nx/2+1,psd.ny/2+1)= 0.0;
    if(!is_void(mask) && numberof(mask) > 0) dat(mask)= 0.0;
    return sqrt(dat(*)(sum)/(Lx*Ly));
  }   
}

func dvv_PSD_1D(psd, qavg=)
/* DOCUMENT psd_1d= dvv_PSD_1D(psd)

   Collapses a 2D PSD data set into an equivalent 1D PSD, by assigning
   to mode N the quadrature sum (i.e. total power) of all mode falling
   within k0*n and k0*(n+1) where k0 is the fundamental mode.  The
   quadrature sum of the resulting spectrum should give the same total
   power as the sum over the 2D spectrum.

   The result is set of three vectors in a 3 x N array; the first
   vector is the frequency of each bin (frequency scale), the second
   is the quadrature sum, and the last is the number of modes assigned
   to the bin.


   KEYWORDS:
     qavg=,  when non-zero the returned result is the quadrature average
   
   SEE ALSO:
     dvv_PSD_1D_rms, dvv_PSD, dvv_PSD_rms
 */
{

  kx0= (*psd.xscale)(dif)(avg);
  ky0= (*psd.yscale)(dif)(avg);

  //Unequal frequency binning in the x- and y- dimensions
  //The WYKO data sets have this problem. Select the modes
  //according to the correct frequencies
  if(abs((kx0 - ky0)/kx0) > 1.0e-12) {
    //print, "Unequal kx, ky grid ... ";
    kk= sqrt(kx0*ky0);
    kxg= (*psd.xscale)(,-:1:psd.ny);
    kyg= (*psd.yscale)(-:1:psd.nx,);
    //klim= min((*psd.xscale)(0), (*psd.yscale)(0));
    //nk= int(klim/kk);
    //ksc= kk*indgen(1:nk);
    //ii= int(sqrt(kxg*kxg + kyg*kyg)/kk);
    ii= round(sqrt(kxg*kxg + kyg*kyg)/kk);
    //wp= where(ii >= 0);
    wp= where(ii > 0);
    NModes= histogram(ii(wp));
    if(qavg) psd1d= histogram(ii(wp), (*psd.data)(wp))/NModes;
    else psd1d= histogram(ii(wp), (*psd.data)(wp));
    //ksc= (indgen(1:numberof(psd1d))-0.5)*kk;
    ksc= indgen(1:numberof(psd1d))*kk;
    
  //Equal frequency binning in the x- and y- dimensions
  } else {
    ix= (*psd.xscale)(,-:1:psd.ny)/kx0;
    iy= (*psd.yscale)(-:1:psd.nx,)/ky0;
    //ii= int(sqrt(ix*ix + iy*iy));
    ii= round(sqrt(ix*ix + iy*iy));
    //wp= where(ii >= 0);
    wp= where(ii > 0);
    NModes= histogram(ii(wp));
    if (qavg) psd1d= histogram(ii(wp), (*psd.data)(wp))/NModes;
    else psd1d= histogram(ii(wp), (*psd.data)(wp));
    //ksc= (indgen(1:numberof(psd1d))-0.5)*kx0;
    ksc= indgen(1:numberof(psd1d))*kx0;
  }
  return [ksc, psd1d, NModes];
}

func dvv_PSD_1D_rms(psd_1d, fl, fu, qavg=, norm=)
/* DOCUMENT   rms_sum= dvv_PSD_1D_rms(psd_1d, fl, fu, norm= norm, qavg= qavg)
         -or- rms_sum= dvv_PSD_1D_rms(psd_1d, frange, norm= norm, qavg= qavg)
         -or- rms_sum= dvv_PSD_1D_rms(psd_1d, norm= norm, qavg= qavg)

   Evaluates an rms sum over a collapsed 1D power spectrum as
   generated by dvv_PSD_1D. The rms sum is calculated by summing over
   all of the modes within the frequency limits given by fl and fu.
   
    Case 1 (norm= 1, default):
       rms_density = sqrt(QUADRATURE_MODE_SUM)
    Case 2 (norm= 2):
       rms_density = sqrt(QUADRATURE_MODE_SUM)/L

   Either form of spectrum normalization is acceptable to use.  The
   second form expresses the spectral density in physical units,
   i.e. rms amplitude per unit frequency
   
   KEYWORDS:
   
   norm= Normalization:
           1) |C_k|^2 / (nx^2*ny^2)
           2) |C_k|^2 *Lx*Ly / (nx*ny)^2
               
   SEE ALSO:
     dvv_PSD_1D, dvv_PSD, dvv_PSD_rms
 */
{
  ds= dimsof(psd_1d);
  if (ds(1) != 2) error, "Invalid data structure for the 1D power spectrum";

  if(!is_void(fl)) {
    if(numberof(fl) == 2 && is_void(fu)) {
      fu= fl(2); fl= fl(1);
    } else if (is_void(fu)) error, "Incorrectly specified frequency limits.";
  } else {
    fl= psd_1d(1,1);
    fu= psd_1d(0,1);
  }
  
  ii= digitize([fl,fu],psd_1d(,1))-1;
  N= numberof(psd_1d(,1));
  L= 1./psd_1d(dif,1)(avg)*(N-1)/N;
  //L= 1./psd_1d(dif,1)(avg);

  //1D spectrum is a quadrature sum
  if(is_void(qavg) || !qavg) prms= sqrt(psd_1d(,2)(ii(1):ii(2))(sum));
  //1D spectrum is azimuthally averaged
  else {
    if(ds(3) < 3) nvec= 2*pi*indgen(1:N); else nvec= psd_1d(,3);
    prms= sqrt((nvec*psd_1d(,2))(ii(1):ii(2))(sum));
  }
  if(is_void(norm) || norm == 1) return prms;
  else if (norm == 2) return prms/L;
  else error, "Invalid value for keyword: norm";
}

func dvv_AMP(img, p=, wndw=, nodc=, norm=)
/* DOCUMENT dvv_AMP, img, p=, wndw=, nodc=, norm=

   Computes the fft amplitude spectrum for the data set contained in
   img.  The amplitude spectrum is the square root of the power
   spectrum.  It normalized such that an appropriate sum of the
   non-zero frequency modes matches the rms of the input data set.

   See section 13.4 in "Numerical Recipes in C"

   There are two normalizations available for the Fourier
   coefficients:

     norm= 1(default): "Amplitude^2 per Fourier Mode": with this
       normalization the RMS = SQRT(Sum over all PSD coefficients,
       excluding DC term). This means that when the sampling interval
       is changed, so will the PSD coefficients.

     norm=2: "Amplitude^2 per unit frequency^2": with this normalization
       the "box size" of the sample is factored into the coefficients.
       This will produce, on average, the same spectral density
       independent of the sampling details. To get the RMS one must
       average over the spectrum; in this process the frequency binning
       units have to be divided out (i.e. the box size parameters
       get divided out):

       RMS=SQRT(Average_over_all_PSD_coefficients_except_DC_Term * Lx*Ly) 

   The function dvv_PSD_rms calculates the correct RMS for each type of
   normalization.

   KEYWORDS:
     wndw= "hanning", "hamming", "bartlett", "triangular", "edge7%"
     p=      fftw plan
     nodc= [] (void), no adjustment of background or DC offset
           0 - sets the dc average to zero, useful for data sets
              with very large DC offsets
           1,2,3 - subtracts n-order polynomial fit to the data before
              taking the transform (removes the lowest order modes).
     norm=   Normalization:
               1) PSD_k,l = |C_k,l|^2 / (nx^2*ny^2)
               2) PSD_k,l |C_k,l|^2 * Lx*Ly / (nx*ny)
   
   SEE ALSO:
     dvv_PSD, dvv_PSD_rms
 */
{
  psd= dvv_PSD(img, wndw= wndw, nodc= nodc, norm= norm);
  psd.data= &sqrt(*psd.data);
  return psd;
}
             
func dvv_AMP_rms(amp, fl, fu, norm=, mask=)
/* DOCUMENT  rms_sum= dvv_AMP_rms(psd, fl, fu, norm= norm, mask= mask)
        -or- rms_sum= dvv_AMP_rms(psd, [fl, fu], norm= norm, mask= mask)
        -or- rms_sum= dvv_AMP_rms(psd, norm= norm, mask= mask)

   Retuns the equivalent RMS sum for a given amplitude spectrum.  The
   input power spectrum can be a 2D IMG_DAT() structure as returned by
   dvv_AMP.  The rms sum is calculated by summing over all of the
   modes within the frequency limits given by fl and fu. In addition a
   filter mask can be applied to filter out other parts of the
   spectrum from being included in the sum additional modes.  If the
   mask is specified, it is OR-ed with the bandpass filter defined by
   the limits, f[l, fu]

   The result depends on the normalization used in computing the PSD.
   The same normalization (=1 or =2) must be specified as was used
   when computing the psd data.


   KEYWORDS:
   
   norm= Normalization:
           1) |C_k|^2 / (nx^2*ny^2)
           2) |C_k|^2 *Lx*Ly / (nx*ny)^2
           
   mask= define this with a filter mask which is applied to the data
         before the summation; for extracting PSD over a subset
         of the frequency spectrum

   SEE ALSO:
     dvv_PSD_rms, dvv_AMP, dvv_AMP_1D, dvv_GenFilterMask
 */
{
  psd= img_copy(amp);
  *psd.data*= (*psd.data);
  return dvv_PSD_rms(psd, fl, fu, norm= norm, mask= mask);
}

func dvv_AMP_1D(amp, qavg=)
/* DOCUMENT amp_1d= dvv_AMP_1D(amp, qavg=)

   Collapses a 2D amplitude spectrum data set into an equivalent 1D
   amplitude spectrum, by assigning to mode N average of the amplitude
   of all modes falling within k0*n and k0*(n+1) where k0 is the
   fundamental mode.  A sum over the square of resulting spectrum
   should give the same total power as the sum over the 2D spectrum.

   The result is set of three vectors in a 3 x N array; the first
   vector is the frequency of each bin (frequency scale), the second
   is the square root of the quadrature sum, and the last is the
   number of modes assigned to the bin.

   This function makes use of dvv_PSD_1D to do most of the
   computation.


   KEYWORDS:
     qavg=, when non-zero the returned result is the square
     root of the average power of all the modes falling within each
     bin
   
   SEE ALSO:
     dvv_PSD_1D_rms, dvv_PSD, dvv_PSD_rms
 */
{
  //Reconstruct PSD from amplitude
  psd= img_copy(amp);
  (*psd.data)*= (*psd.data);
  amp_1d= dvv_PSD_1D(psd, qavg= qavg);
  //if(qavg) amp_1d(,2)= sqrt(amp_1d(,2));
  //else amp_1d(,2)= sqrt(amp_1d(,2));
  amp_1d(,2)= sqrt(amp_1d(,2));
  return amp_1d;

  // {
  //   psd= img_copy(amp);
  //   //if(is_void(alt) || !alt) *psd.data*= (*psd.data);
  //   kx0= (*psd.xscale)(dif)(avg);
  //   ky0= (*psd.yscale)(dif)(avg);

  //   //Unequal frequency binning in the x- and y- dimensions
  //   //The WYKO data sets have this problem. Select the modes
  //   //according to the correct frequencies
  //   if(abs((kx0 - ky0)/kx0) > 1.0e-12) {
  //     //print, "Unequal kx, ky grid ... ";
  //     kk= sqrt(kx0*ky0);
  //     kxg= (*psd.xscale)(,-:1:psd.ny);
  //     kyg= (*psd.yscale)(-:1:psd.nx,);
  //     //klim= min((*psd.xscale)(0), (*psd.yscale)(0));
  //     //nk= int(klim/kk);
  //     //ksc= kk*indgen(1:nk);
  //     ii= int(sqrt(kxg*kxg + kyg*kyg)/kk);
  //     wp= where(ii > 0);
  //     NModes= histogram(ii(wp));
  // //     if(is_void(alt) || !alt) {
  // //       if(qavg) psd1d= sqrt(histogram(ii(wp), (*psd.data)(wp)))/NModes;
  // //       else psd1d= sqrt(histogram(ii(wp), (*psd.data)(wp)));
  // //     } else {
  //       if(qavg) psd1d= histogram(ii(wp), (*psd.data)(wp))/NModes;
  //       else psd1d= histogram(ii(wp), (*psd.data)(wp));
  // //    }
  //     ksc= indgen(1:numberof(psd1d))*kk;
    
  //     //Equal frequency binning in the x- and y- dimensions
  //   } else {
  //     ix= (*psd.xscale)(,-:1:psd.ny)/kx0;
  //     iy= (*psd.yscale)(-:1:psd.nx,)/ky0;
  //     ii= int(sqrt(ix*ix + iy*iy));
  //     wp= where(ii > 0);
  //     NModes= histogram(ii(wp));
  // //     if(is_void(alt) || !alt) {
  // //       if (qavg) psd1d= sqrt(histogram(ii(wp), (*psd.data)(wp)))/NModes;
  // //       else psd1d= sqrt(histogram(ii(wp), (*psd.data)(wp)));
  // //     } else {
  //       if (qavg) psd1d= histogram(ii(wp), (*psd.data)(wp))/NModes;
  //       else psd1d= histogram(ii(wp), (*psd.data)(wp));
  // //    }
  //     ksc= indgen(1:numberof(psd1d))*kx0;
  //   }
  //   return [ksc, psd1d, NModes];
  // }

}

func dvv_AMP_1D_rms(amp_1d, fl, fu, qavg=, norm=)
/* DOCUMENT   rms_sum= dvv_AMP_1D_rms(amp_1d, fl, fu, norm= norm, qavg= qavg)
         -or- rms_sum= dvv_AMP_1D_rms(amp_1d, frange, norm= norm, qavg= qavg)
         -or- rms_sum= dvv_AMP_1D_rms(amp_1d, norm= norm, qavg= qavg)

   Evaluates an rms sum over a collapsed 1D power spectrum as
   generated by dvv_AMP_1D. The rms sum is calculated by summing over
   all of the modes within the frequency limits given by fl and fu.
   
    Case 1 (norm= 1, default):
       rms_density = sqrt(QUADRATURE_MODE_SUM)
    Case 2 (norm= 2):
       rms_density = sqrt(QUADRATURE_MODE_SUM)/L

   Either form of spectrum normalization is acceptable to use.  The
   second form expresses the spectral density in physical units,
   i.e. rms amplitude per unit frequency
   
   KEYWORDS:
   
   norm= Normalization:
           1) |C_k|^2 / (nx^2*ny^2)
           2) |C_k|^2 *Lx*Ly / (nx*ny)^2
               
   SEE ALSO:
     dvv_AMP_1D, dvv_AMP, dvv_AMP_rms
 */
/*{
 *  ds= dimsof(amp_1d);
 * if (ds(1) != 2) error, "Invalid data structure for the 1D amplitude spectrum";
 * psd_1d= amp_1d;
 * psd_1d(,2)*= psd_1d(,2);
 * return dvv_PSD_1D_rms(psd_1d, fl, fu, qavg= qavg, norm= norm);
 *}
*/
{
  ds= dimsof(amp_1d);
  if (ds(1) != 2) error, "Invalid data structure for the 1D amplitude spectrum";

  if(!is_void(fl)) {
    if(numberof(fl) == 2 && is_void(fu)) {
      fu= fl(2); fl= fl(1);
    } else if (is_void(fu)) error, "Incorrectly specified frequency limits.";
  } else {
    fl= amp_1d(1,1);
    fu= amp_1d(0,1);
  }
  
  ii= digitize([fl,fu],amp_1d(,1))-1;
  N= numberof(amp_1d(,1));
  L= 1./amp_1d(dif,1)(avg)*(N-1)/N;
  //L= 1./amp_1d(dif,1)(avg);

  //1D spectrum is a quadrature sum for each mode
  if(ds(3) < 3) nvec= 2*pi*indgen(1:N); else nvec= amp_1d(,3);
  //if(is_void(qavg) || !qavg) prms= sqrt((amp_1d(,2)^2/nvec)(ii(1):ii(2))(sum));
  if(is_void(qavg) || !qavg) prms= sqrt((amp_1d(,2)^2)(ii(1):ii(2))(sum));
  //1D spectrum is azimuthally averaged
  else {
    prms= sqrt((nvec*(amp_1d(,2))^2)(ii(1):ii(2))(sum));
    //prms= sqrt(2*pi*integ(amp_1d(,2)^2, amp_1d(,1), [fl, fu])(dif));
  }
  if(is_void(norm) || norm == 1) return prms;
  else if (norm == 2) return prms/L;
  else error, "Invalid value for keyword: norm";
}

func dvv_fft2D(img, p=, wndw=, nodc=, dir=)
/* DOCUMENT dvv_fft2D, img, p=, wndw=, nodc=, dir=

   Computes the complex fft spectrum for the data set
   contained in img.

   See section 13.4 in "Numerical Recipes in C"

   KEYWORDS:
     wndw= "hanning", "hamming", "bartlett", "triangular", "edge7%"
     p=      fftw plan
     nodc=   1, sets the dc average to zero
     dir=    +1 for forward transform
             -1 for reverse transform
   
   SEE ALSO:
     dvv_PSD
 */
{
  extern _NO_FFTW;

  ft= img_copy(img);
  if(is_void(dir)) dir= 1;
  //if(dir > 0) print, "Forward transform ...";
  //if(dir < 0) print, "Backward transform ...";
  
  dat= *img.data;
  if(nodc) dat-= dat(*)(avg);  /* Subtract DC offset */
  nxy= double(ft.nx)*double(ft.ny);

  if(!is_void(wndw)) {
    xg= span(-pi, pi, img.nx)(,-:1:img.ny);
    yg= span(-pi, pi, img.ny)(-:1:img.nx,);
    if(wndw == "hanning") {
      wndw= 0.25*(1.0 + cos(xg))*(1.0 + cos(yg));
      wfac= (wndw*wndw)(sum,sum);
    } else if (wndw == "hamming") {
      wndw= (0.54 + 0.46*cos(xg))*(0.54 + 0.46*cos(yg));
      wfac= (wndw*wndw)(sum,sum);
    } else if (wndw == "triangular" || wndw == "bartlett") {
      wndw= (1.0 - abs(xg)/pi)*(1.0 - abs(yg)/pi);
      wfac= (wndw*wndw)(sum,sum);
    } else if (wndw == "edge7%") {
      ww= (abs(xg) > 0.93*pi) | (abs(yg) > 0.93*pi);
      wndw= array(1.0, dimsof(xg));
      wndw(where(ww))= 0.0;
      wfac= (wndw*wndw)(sum,sum);
    } else error, "Unknown window function";
  } else {
    wfac= nxy;
  }
  if(dir == 1 && !is_void(wndw)) dat*= wndw;

  //if(is_void(p) && !_NO_FFTW) p= fftw_plan(dimsof(dat), dir);

  //The roll function doesn't seem to work perfectly for odd-numbered
  //sample sizes, so we need to fix it in the case of backward
  //transforms.  The fix is applied to the data, but not to the scales(?)
  if(_NO_FFTW) {
    if(dir > 0) {
      ft.data= &roll(fft(dat,+1));
    } else if(dir < 0) {
      ft.data= &fft(roll(dat,[ft.nx/2+ft.nx%2,ft.ny/2+ft.ny%2]), -1);
    } else error, "Invalid fft direction parameter.";
  } else {
    if(dir > 0) {
      ft.data= &roll(fftw(dat, p));
    } else if(dir < 0) {
      ft.data= &fftw(roll(dat,[ft.nx/2+ft.nx%2,ft.ny/2+ft.ny%2]), p);
    } else error, "Invalid fft direction parameter.";
  }
  dx= (*img.xscale)(dif)(avg);
  dy= (*img.yscale)(dif)(avg);
  if(dir < 0) {
    //ft.xscale= &roll(fft_indgen(ft.nx)/double(ft.nx)/dx, ft.nx/2+ft.nx%2);
    //ft.yscale= &roll(fft_indgen(ft.ny)/double(ft.ny)/dy, ft.ny/2+ft.ny%2);
    ft.xscale= &roll(fft_indgen(ft.nx)/double(ft.nx)/dx, ft.nx/2);
    ft.yscale= &roll(fft_indgen(ft.ny)/double(ft.ny)/dy, ft.ny/2);
    if (!(ft.nx%2)) (*ft.xscale)(1)= - abs((*ft.xscale)(1));
    if (!(ft.ny%2)) (*ft.yscale)(1)= - abs((*ft.yscale)(1));
    if(!is_monotonic(*ft.xscale)) error, "Internal roll error on xscale";
    if(!is_monotonic(*ft.yscale)) error, "Internal roll error on yscale";
  } else {
    ft.xscale= &roll(fft_indgen(ft.nx)/double(ft.nx)/dx);
    ft.yscale= &roll(fft_indgen(ft.ny)/double(ft.ny)/dy);
    if (!(ft.nx%2)) (*ft.xscale)(1)= - abs((*ft.xscale)(1));
    if (!(ft.ny%2)) (*ft.yscale)(1)= - abs((*ft.yscale)(1));
    if(!is_monotonic(*ft.xscale)) error, "Internal roll error on xscale";
    if(!is_monotonic(*ft.yscale)) error, "Internal roll error on yscale";
  }
  ft.x_unit= ft.x_unit+"^-1^";
  ft.y_unit= ft.y_unit+"^-1^";
  ft.x_label= "frequency";
  ft.y_label= "frequency";

  //print, "wfac, nxy= ", wfac, nxy;
  if(dir > 0) *ft.data/= sqrt(nxy);
  else if(dir < 0) *ft.data/= sqrt(nxy);
  if(!is_void(wndw) && dir < 0) {
    wz= where(wndw==0.0);
    if(numberof(wz)) wndw(wz)= 1e10;
    ft.data= &(*ft.data/wndw);
  }
  if(dir < 0) ft.data= &((*ft.data).re);
  return ft;
}

func dvv_Mode(psd, kx, ky, n=)
/* DOCUMENT dvv_Mode, psd, fx, fy

   Returns the psd components for the mode with frequency (fx, fy).
   Since the mode at (fx,fy) has psd components at (fx,fy) and
   (-fx,-fy), both terms are returned.

   KEYWORDS:
     n=  returns the psd components are frequencies up to to n modes away
     
   SEE ALSO:
 */
{
  kx0= (*psd.xscale)(dif)(avg);
  ky0= (*psd.yscale)(dif)(avg);
  ix= digitize([-kx,kx], *psd.xscale);
  iy= digitize([-ky,ky], *psd.yscale);
  kxg= (*psd.xscale)(,-:1:psd.ny);
  kyg= (*psd.yscale)(-:1:psd.nx,);
  kg= sqrt(kxg*kxg + kyg*kyg);
  //print, "kmag= ",kg(ix(1)-1, iy(1)-1), kg(ix(2), iy(2));
  //print, "kmag= ",kg(ix(1)-1, iy(1)-1), kg(ix(2)-1, iy(2)-1);
  //if(is_void(n)) return [(*psd.data)(ix(1)-1, iy(1)-1),
  //        (*psd.data)(ix(2), iy(2))];
  if(is_void(n)) return [(*psd.data)(ix(1)-1, iy(1)-1),
          (*psd.data)(ix(2)-1, iy(2)-1)];
  else return [(*psd.data)(ix(1)-1-n:ix(1)-1+n, iy(1)-1-n:iy(1)-1+n),
          (*psd.data)(ix(2)-1-n:ix(2)-1+n, iy(2)-1-n:iy(2)-1+n)];
}

func dvv_Filter(img, mask, az, p=, wndw=, rnfill=)
  /* DOCUMENT      dvv_Filter, img, mask, p=, wndw=, rnfill=
              -or- dvv_Filter, img, mask, az, p=, wndw=, rnfill=
              -or- dvv_Filter, img, fmodel, p=, wndw=

  1st form: Filters a data set by zeroing elements of the DFT of the
    image array.  The elements to be zeroed are given on the array
    mask, which contains a a list of indices indicating which elements
    of the DFT of the function need to be zeroed
     
  2nd form: Performs the same function as the first form with the
    additional step of filling in the masked elements with values
    corresponding to the azimuthal average of the remaining (unmasked)
    elements.  The magnitudes of these elements are set to the
    magnitude of the azimuthal average, and the phase is
    random. Because this is a DFT of a real-valued function the mask
    must obey point-reflection symmetry (i.e. for each [kx,ky]
    element, there is a [-kx,-ky] element in the mask list. These
    elements are then filled with such that the symmetric elements are
    complex conjugates of each other, i.e.  [kx, ky] <--
    Azavg*exp(I*R), [-kx,-ky] <-- ComplexConjugate([kx,ky]).

  3rd form: Applies a filter function defined in the complex plane to
    DFT of the data

  Filter masks are generated by dvv_GenFilterMask.

     
  SEE ALSO:
    dvv_GenFilterMask, dvv_CombineMasks, dvv_ShowFilterMask
  */
{
  extern _NO_FFTW;
  ft= img_copy(img);
  
  nxy= double(img.nx)*double(img.ny);
  wfac= [];

  if(!is_void(wndw)) {
    xg= span(-pi, pi, img.nx)(,-:1:img.ny);
    yg= span(-pi, pi, img.ny)(-:1:img.nx,);

    if(wndw == "hanning") {
      win= 0.25*(1.0 + cos(xg))*(1.0 + cos(yg));
      wfac= (win*win)(sum,sum);
      dat= (*img.data) * win;
    } else if (wndw == "hamming") {
      win= (0.54 + 0.46*cos(xg))*(0.54 + 0.46*cos(yg));
      wfac= (win*win)(sum,sum);
      dat= (*img.data) * win;
    } else if (wndw == "triangular") {
      win= (1.0 - abs(xg)/pi)*(1.0 - abs(yg)/pi);
      wfac= (win*win)(sum,sum);
      dat= (*img.data) * win;
    } else if (wndw == "edge7%") {
      ww= (abs(xg) > 0.93*pi) | (abs(yg) > 0.93*pi);
      win= array(1.0, dimsof(xg));
      win(where(ww))= 0.0;
      wfac= (win*win)(sum,sum);
      dat= (*img.data) * win;
    } else error, "Unknown window function";
  } else {
    dat= *img.data;
  }

  //print, "nxy= ", nxy, "wfac= ", wfac;
  //filtered= img_copy(img);
  //print, "Sum: ", (*img.data)(*)(sum);
  if(_NO_FFTW) {
    //print, "Using Yorick fft (+1) ... ";
    ft.data= &(roll(fft(dat, +1)));
  } else {
    //pf= fftw_plan(dimsof(*img.data), +1);
    //ft.data= &(roll(fftw(dat, pf)));
    ft.data= &(roll(fftw(dat, +1)));
  }
  if(typeof(mask) == "long") (*ft.data)(mask)= 0.0 + 0.0i;
  else if (structof(mask) == IMG_DAT) {
    if(mask.nx != ft.nx || mask.ny != img.ny) error, "invalid mask structure";
    else (*ft.data)*= *mask.data;
  } else if (typeof(mask) == "string") {
    mask= dvv_GenFilterMask(img, mask);
    (*ft.data)(mask)= 0.0 + 0.0i;
  } else error, "invalid mask structure";
  
  //print, (*ft.data)(ft.nx/2-10:ft.nx/2+10, ft.ny/2);
  if (!is_void(az)) {
    //print,"Applying fill to filter mask ... ";
    dx= (*img.xscale)(dif)(avg);
    dy= (*img.yscale)(dif)(avg);
    xsc= roll(fft_indgen(ft.nx)/double(ft.nx)/dx);
    ysc= roll(fft_indgen(ft.ny)/double(ft.ny)/dy);
    if (!(ft.nx%2)) xsc(1)*= -1.0;
    if (!(ft.ny%2)) ysc(1)*= -1.0;
    ix= (mask - 1)/ft.nx + 1;
    iy= (mask - 1)%ft.ny + 1;
    
    //For even numbered data sets, ignore the first index
    //in the rolled transform data (corresponds to Nyquist
    //frequency), remove the F_Nyquist vakues from the mask
    if(!(ft.nx%2) && !(ft.ny%2)) wx= where(ix != 1 & iy != 1);
    else if (!(ft.nx%2)) wx= where(ix != 1);
    else if (!(ft.ny%2)) wx= where(iy != 1);
    //Shift the indices to +ve and -ve
    //Sorted indices are now in ixs, iys arrays
    if(numberof(wx) > 0) {
      //Sort the mask values according to x-, and y-indices
      print, "Removing ix= 1, iy= 1 ...";
      is=msort(ix, iy);
      ms= mask(is(wx));
      ixs= ix(is(wx));
      iys= iy(is(wx));
    } else {
      is=msort(ix, iy);
      ms= mask(is);
      ixs= ix(is);
      iys= iy(is);
    }
    
    N= numberof(ixs);
    ixs1= ixs - (ft.nx/2 + 1);
    iys1= iys - (ft.ny/2 + 1);
    //print, "allof x", allof(ixs1(1:N/2) == -ixs1(N/2+1:0)(0:1:-1));
    //info, where(ixs1(1:N/2) == -ixs1(N/2+1:0)(0:1:-1));
    //print, "allof y", allof(iys1(1:N/2) == -iys1(N/2+1:0)(0:1:-1));
    info, where(iys1(1:N/2) == -iys1(N/2+1:0)(0:1:-1));
    xx= xsc(ixs);
    yy= ysc(iys);
    rr= sqrt(xx*xx + yy*yy);
    I= 0.0+1.0i;
    fill= interp(sqrt(az(,2)), az(,1), rr(1:N/2));
    if(!is_void(rnfill)) fill+= fill*random_n(dimsof(fill));
    //if(!is_void(wndw)) fill*= win(ms(1:N/2));
    phs= exp(I*pi*(2*random(N/2) - 1.0));
    window, 0;
    plmk, yy(1:N/2), xx(1:N/2), msize= 0.3;
    window, 1;
    plmk, fill, rr(1:N/2), msize= 0.2;
    fill= fill*phs;
    //Fill in the x > 0 half-plane
    info, fill;
    print, "ixs= ", ixs(1:20)-ft.nx/2-1, ixs(-19:0)-ft.nx/2-1;
    print, "iys= ", iys(1:20)-ft.ny/2-1, iys(-19:0)-ft.ny/2-1;
    window, 4;
    zim= abs(*ft.data);
    wzm= where(zim <= 0.0);
    if(numberof(wzm) > 0) zim(wzm)= 1e-20;
    pli, bytscl(log10(zim), cmin= -6, cmax= -3);
    print, (*ft.data)(ms(1:N/2))(-99:0);
    (*ft.data)(ms(1:N/2))= fill;
    //print, (*ft.data)(ms(1:N/2))(-99:0);
    //print, rr(1:N/2)(-99:0);
    //print, phs(1:N/2)(-99:0);
    //print, abs(phs(1:N/2))(-99:0);
    //Fill in the x < 0 half-plane with the complex conjugate
    fill.im= -fill.im;
    (*ft.data)(ms(N/2+1:0)(0:1:-1))= fill;
    
    ixx= (ms - 1)/ft.nx + 1;
    iyy= (ms - 1)%ft.ny + 1;
    for(i= 1; i<= 10; i++) {
      print, ixx(i), xsc(ixx(i)), iyy(i), ysc(iyy(i)), (*ft.data)(ms(i));
      print, ixx(-i+1), xsc(ixx(-i+1)), iyy(-i+1), ysc(iyy(-i+1)), (*ft.data)(ms(-i+1));
    }
  }
  //zim= abs(*ft.data);
  //wzm= where(zim <= 0.0);
  //if(numberof(wzm) > 0) zim(wzm)= 1e-20;
  //ft.data= &zim;
  //dvv_DisplayImg, ft, w= 2, lgscl= 1, 1e-5, 1e-2;
  //pli, bytscl(log10(zim), cmin= -6, cmax= -3);
  ft.data= &roll(*ft.data, [ft.nx/2+ft.nx%2,ft.ny/2+ft.ny%2]);
  if(_NO_FFTW) {
    //print, "Using Yorick fft (-1) ... ";
    ft.data= &(fft(*ft.data, -1)/nxy);
  } else {
    //pr= fftw_plan(dimsof(*ft.data), -1);
  //  filtered.data= &((fftw(*ft.data, pr).re)/(ft.nx*ft.ny));
    //ft.data= &(fftw(*ft.data, pr)/nxy);
    ft.data= &(fftw(*ft.data, -1)/nxy);
  }
  
  if(!is_void(wfac)) {
    wz= where(win==0.0);
    if(numberof(wz)) win(wz)= 1e10;
    ft.data= &((*ft.data)/win);
  }
  //print, "real:",((*ft.data).re)(sum,sum);
  //print, "imag:",((*ft.data).im)(sum,sum);
  ft.data= &((*ft.data).re);
  return ft;
}

func dvv_GenFilterMask(img, p, type=)
/* DOCUMENT      dvv_GenFilterMask, img, p, type=
            -or- dvv_GenFilterMask, img, maskspec

   Generates a binary filter mask of various types based on the
   specified mask specification. The mask is an array of indices into
   the frequency components of the 2D fft of the data that are set to
   0 during filtering operations.

   The mask specification consists of an identifying string and a list
   of parameters.  In the first calling form the type is given by the
   type=<type> keyword and the parameter list is given in the second
   argument. In the second calling form the type and parameters are
   input in the second argument as (1) a single string: "<type> <p1>
   <p2> <p3> ... "; with tokens delimited by spaces, and/or commas
   and/or tabs; or, (2) as a string array of the individual tokens.

   For example:
   
        msk= dvv_GenFilterMask(img, [0.023, 0.19, 0.05], type="mode_reject")
   -or- msk= dvv_GenFilterMask(img, "mode_reject 0.023 0.19 0.05")
   -or- msk= dvv_GenFilterMask(img, ["mode_reject", "0.023", "0.19", "0.05"])
      
   All of these calling forms generate a mode reject filter for the
   mode centered at kx=0.023 ky=0.19 and band reject radius of 0.05
   centered on this mode.

   Available filter types are: [<type>, <param1>, <param2>, ... ]
   
     mode_reject, kx, ky, band_radius
     mode_select, kx, ky, band_radius
     low_pass, kcut
     high_pass, kcut
     band_pass, klow, khigh
     high_pass_1D, kx_coeff, ky_coeff, kcut
     low_pass_1D, kx_coeff, ky_coeff, kcut
     high_pass_slot, kx_coeff, ky_coeff, kcut, hp_cutoff
     offset_high_pass_slot, kx_coeff, ky_coeff, width, hp_cutoff, offset
     
   SEE ALSO:
     dvv_CombineMasks, dvv_ShowFilterMask, dvv_Filter
*/
{ 
  if(typeof(p) == "string") {
    if(numberof(p) <= 1) {
      mm= parse_line(p, ", \t",10);
      type= mm(1);
      p= strread(mm(2:0), "%e");
    } else {
      type= p(1);
      p= strread(p(2:0), "%e");
    }
  }

  if (is_void(type)) error, "No filter type specified";
  if (is_void(p)) error, "No filter parameters specified";
  
  dx= (*img.xscale)(dif)(avg);
  dy= (*img.yscale)(dif)(avg);
  xsc= roll(fft_indgen(img.nx)/double(img.nx)/dx, img.nx/2);
  ysc= roll(fft_indgen(img.ny)/double(img.ny)/dy, img.ny/2);
  if (!(img.nx%2)) xsc(1)*= -1.0;
  if (!(img.ny%2)) ysc(1)*= -1.0;
  if(!is_monotonic(xsc)) error, "Internal roll error on xscale";
  if(!is_monotonic(ysc)) error, "Internal roll error on yscale";
  xg= xsc(,-:1:img.ny);
  yg= ysc(-:1:img.nx,);
  
  if (type == "mode_reject") {
    return where((((xg - p(1))*(xg - p(1)) + (yg - p(2))*(yg - p(2)) < p(3)*p(3)) |
                  ((xg + p(1))*(xg + p(1)) + (yg + p(2))*(yg + p(2)) < p(3)*p(3))));
  } else if (type == "mode_select") {
    return where((((xg - p(1))*(xg - p(1)) + (yg - p(2))*(yg - p(2)) > p(3)*p(3)) &
                  ((xg + p(1))*(xg + p(1)) + (yg + p(2))*(yg + p(2)) > p(3)*p(3))));
  } else if (type == "mode_pass") {
    w1= ((xg - p(1))*(xg - p(1)) + (yg - p(2))*(yg - p(2)) > p(3)*p(3));
    w2= ((xg + p(1))*(xg + p(1)) + (yg + p(2))*(yg + p(2)) > p(3)*p(3));
    //print, "Sum w1= ", w1(sum, sum), "Sum w2=", w2(sum, sum), "Sum w1&w2= ", (w1 & w2)(sum, sum);
    return where(w1 & w2);
  } else if (type == "low_pass") {
    // first parameter is the cutoff frequency
    return where(xg*xg + yg*yg > p(1)*p(1));
  } else if (type == "high_pass") {
    // first parameter is the cutoff frequency
    return where(xg*xg + yg*yg < p(1)*p(1));
  } else if (type == "band_pass") {
    // first parameter is the cutoff frequency
    return where(xg*xg + yg*yg > p(2)*p(2) | xg*xg + yg*yg < p(1)*p(1));
  } else if (type == "high_pass_SQ") {
    // first parameter is the cutoff frequency
    return where(abs(xg) < p(1) & abs(yg) < p(1));
  } else if (type == "high_pass_1D") {
    // p(1) is x-coeff, p(2) is y-coeff, p(3) is cutoff
    return where((p(1)*xg + p(2)*yg < p(3)) & (- p(1)*xg - p(2)*yg < p(3)));
  } else if (type == "low_pass_1D") {
    // p(1) is x-coeff, p(2) is y-coeff, p(3) is cutoff
    return where((p(1)*xg + p(2)*yg > p(3)) | (- p(1)*xg - p(2)*yg > p(3)));
  } else if (type == "high_pass_slot") {
    print, "xg min:max: ", min(xg), max(xg);
    print, "yg min:max: ", min(yg), max(yg);
    pr= sqrt((p(1:2)^2)(sum));
    xc= p(2)/pr; yc= -p(1)/pr;
    return where(((xc*xg + yc*yg < p(3)) &
                  (- xc*xg - yc*yg < p(3))) &
                 (xg*xg + yg*yg > p(4)*p(4)));
  } else if (type == "offset_high_pass_slot") {
    print, "xg min:max: ", min(xg), max(xg);
    print, "yg min:max: ", min(yg), max(yg);
    pr= sqrt((p(1:2)^2)(sum));
    xc= p(2)/pr; yc= -p(1)/pr;
    xoff= 1.0/xc;
    yoff= 1.0/yc;
    offN= sqrt(xoff*xoff + yoff*yoff);
    xoff= xoff*p(5)/offN;
    yoff= yoff*p(5)/offN;
    wminus= where( ((xc*(xg-xoff) + yc*(yg-yoff) < p(3)) &
                   (- xc*(xg-xoff) - yc*(yg-yoff) < p(3))) &
                  ((xg-xoff)*(xg-xoff) + (yg-yoff)*(yg-yoff) > p(4)*p(4)));
    wplus= where( ((xc*(xg+xoff) + yc*(yg+yoff) < p(3)) &
                   (- xc*(xg+xoff) - yc*(yg+yoff) < p(3))) &
                  ((xg+xoff)*(xg+xoff) + (yg+yoff)*(yg+yoff) > p(4)*p(4)));
    return _(wminus, wplus);
  } else {
    error, swrite(format="%s is an unrecognized filter mask type", type);
  }
}

func dvv_GenFilterMaskSet(img, p, type=)
/* DOCUMENT      dvv_GenFilterMaskSet, img, p, type=
            -or- dvv_GenFilterMaskSet, img, maskspec

   Generates compound filter masks based on dvv_GenFilterMask.
   Multiple filters may be specified and combined at once in a single
   call to dvv_GenFilterMaskSet to take the place of multiple calls to
   dvv_GenFilterMask and dvv_CombineMasks.

   (1) For a specification consisting of a list of N filter masks all
   of the same type one can use:
   
       m= dvv_GenFilterMaskSet(img, p, type=<type>). For example,

       p= [[0.139, 0.0359, 0.0085],[0.1005, 0.0831, 0.0085],
           [0.1030, 0.0998, 0.0085]];
       m= dvv_GenFilterMaskSet(img, p, type= "mode_reject");

   The parameter list p must be an [M x N] array where M is the
   numberof parameters in the filter definition, N is the number of
   filters, and the key word type must be defined.

   (2) For a specification consisting of a list of N miscellaneous
   filters one must use:
   
      m= dvv_GenFilterMaskSet(img, [maskspec1, maskspec2, maskspec3 ...]
      where maskspec<n> is a string, or a string array as defined in
      the 2nd and 3rd calling forms of dvv_GenFilerMask. For example,

       farr= ["mode_reject 0.139 0.0359 0.0085",
              "mode_reject 0.1005 0.0831 0.0085",
              "mode_reject 0.1030 0.0998 0.0085",
              "band_pass 0.02 0.5"];
       m= dvv_GenFilterMaskSet(img, farr);

   where maskspec<n> is a string, or a string array as defined in the
   2nd and 3rd calling forms of dvv_GenFilerMask.

   SEE ALSO:
     dvv_GenFilterMask, dvv_ShowFilterMask, dvv_Filter
 */
{
  if(typeof(p) == "string") {
    if(dimsof(p)(1) == 2) {
      NF= dimsof(p)(0);
      marr= array(pointer, NF);
      for(i= 1; i<= NF; i++) marr(i)= &dvv_GenFilterMask(img, p(,i));
    } else if (dimsof(p)(1) == 1) {
      NF= numberof(p);
      marr= array(pointer, NF);
      for(i= 1; i<= NF; i++) marr(i)= &dvv_GenFilterMask(img, p(i));
    } else {
      if(dimsof(p)(1) != 2) error, "Need a 2D array for the filter parameters";
      NF= dimsof(p)(0);
      marr= array(pointer, NF);
      for(i= 1; i<= NF; i++) marr(i)= &dvv_GenFilterMask(img, p(,i), type= type);
    }
  }
  m= *marr(1);
  for(i= 2; i<= NF; i++) m= dvv_CombineMasks(m, *marr(i));
  return m;
}

func dvv_ShowFilterMask(img, m, w=)
/* DOCUMENT dvv_ShowFilterMask, img, m, w=
     
   SEE ALSO:
     dvv_GenFilterMask
 */
{
  dx= (*img.xscale)(dif)(avg);
  dy= (*img.yscale)(dif)(avg);
  xsc= roll(fft_indgen(img.nx)/double(img.nx)/dx, img.nx/2);
  ysc= roll(fft_indgen(img.ny)/double(img.ny)/dy, img.ny/2);
  if (xsc(1) > 0.0) xsc(1)*= -1.0;
  if (ysc(1) > 0.0) ysc(1)*= -1.0;
  if(!is_monotonic(xsc)) error, "Internal roll error on xscale";
  if(!is_monotonic(ysc)) error, "Internal roll error on yscale";
  xg= xsc(,-:1:img.ny);
  yg= ysc(-:1:img.nx,);
  zg= array(1.0, dimsof(xg));
  zg(m)= 0.0;
  m= img_copy(img);
  m.xscale= &xsc;
  m.yscale= &ysc;
  m.data= &zg;
  m.x_unit= img.x_unit+"^-1^";
  m.y_unit= img.y_unit+"^-1^";
  dvv_DisplayImg, m, w= w, 0, 1;
  if(!am_subroutine()) return m;
}

func dvv_CombineMasks(m1, m2, ..)
/* DOCUMENT dvv_CombineMasks, m1, m2, ..
     
   SEE ALSO:
     dvv_GenFilterMask
 */
{
  if(is_void(m2) && typeof(m1) == "pointer") {
    if (numberof(m1) > 1) {
      mm= [];
      for (i= 1; i<= numberof(m1); i++) mm= _(mm, *m1(i));
    } else {
      mm= *m1(1);
    }
  } else {
    //Combine the mask elements into one big array;
    mm= _(m1,m2);
    while(more_args()) mm= _(mm, next_arg());
  }
  //info, mm;
  mm= mm(sort(mm));
  
  //Sort and compute pairwise differences
  //print, mm(dif)(1:10);
  d= _(1,mm(dif));
  //info, d;

  //Removes repeated elements
  return mm(where(d));
}

func dvv_ORMasks(m1, ..)
/* DOCUMENT dvv_ORMasks, m1, m2, ..
     
   SEE ALSO:
     dvv_GenFilterMask
 */
{
  //Combine the mask elements into one big array;
  ma= m1;
  do {
    mb= next_arg();
    mm= _(ma, mb);
    //Sort and compute pairwise differences
    mm= mm(sort(mm));
    d= _(1,mm(dif));

    //Extract repeated indices
    ma= mm(where(!d));
  } while (more_args());
  return ma;
}

func mode_vec(N, min_wav, max_wav, min_ampl, max_ampl)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  am= lx= ly= [];
  do {
    lx= _(lx, max_wav*(random(2*N) - 0.5));
    ly= _(ly, max_wav*(random(2*N) - 0.5));
    wm= where(abs(lx) > min_wav & abs(ly) > min_wav);
  }
  while (numberof(wm) < N);
  if(is_void(min_ampl) || is_void(max_ampl)) am= array(1.0, N);
  else do {
    am= _(am, max_ampl*(random(2*N) - 0.5));
    wa= where(abs(am) > min_ampl);
  }
  while (numberof(wa) < N);
  return [am(wa)(1:N), lx(wm)(1:N), ly(wm)(1:N)];
}

func dvv_mkModes(img, mv, coherent=)
/* DOCUMENT m= dvv_mkModes(img, mv)

   Generates a synthetic data set consisting of an array of
   pure sinusoidal modes, as specified by the mv argument.

   ARGUMENTS:
     mv - mode vector, an nx3 vector of the form:
             [[a1, a2, ...],[lx1, lx2, ...],[ly1, ly2, ...]]
          where the a's are the amplitudes, lx's & ly's are the
          x- and y-components of the wavelength. For example
          the mode vector:
            mv = [[0.01],[15.0],[6.0]]
          produces the mode:
            z = 0.01*cos(2*pi*(x/15 + y/6.0));
   SEE ALSO:
 */
{
  mm= img_copy(img);
  xg= (*mm.xscale)(,-:1:mm.ny);
  yg= (*mm.yscale)(-:1:mm.nx,);
  zg= array(0.0, dimsof(xg));
  for(i= 1; i<= dimsof(mv)(2); i++) {
    if(!coherent) dp= random()*2*pi; else dp= 0.0;
    if (abs(mv(i,2)) < 1e-4) {
      zg+= mv(i,1)*cos(2*pi*yg/mv(i,3) + dp);
    } else if (abs(mv(i,3)) < 1e-4) {
      zg+= mv(i,1)*cos(2*pi*xg/mv(i,2) + dp);
    } else {
      zg+= mv(i,1)*cos(2*pi*(xg/mv(i,2) + yg/mv(i,3)) + dp);
    }
  }
  mm.data= &zg;
  return mm;
}

func dvv_fromCX(img, type=, nodc=, norm=)
/* DOCUMENT dvv_fromCX, img, type=, nodc= , norm=

   Available types
     ampl - amplitude
     phase
     phaseV
     spectP
     spectA
     spectX
     spectV
     shifted

   KEYWORDS:
     type=
     nodc=
     
   SEE ALSO:
 */
{
  extern _FRINGE_MODE;
  extern _WAV, _TAU, _ETALON_DELTA, _VPF;
    
  f= img_copy(img);
  if(type == "ampl") {
    f.data= &abs(*f.data);
    return f;
  } else if (type == "phase") {
    f.data= &atan((*f.data).im, (*f.data).re);
    f.z_unit = "rad";
    f.z_label = "Phase";
    return f;
  } else if (type == "phaseV") {
    if(typeof(*f.data) == "complex") f.data= &(atan((*f.data).im, (*f.data).re)*_VPF/2/pi*1000.0);
    else if(typeof(*f.data) == "double") f.data= &((*f.data)*_VPF/2/pi*1000.0);
    f.z_unit = "m/s";
    f.z_label = "Velocity";
    return f;
  } else if (type == "spectP") {
    return dvv_PSD(f, wndw= "hanning", nodc= nodc, norm= norm);
  } else if (type == "spectA") {
    return dvv_AmplSpectr(dvv_PSD(f, wndw= "hanning", nodc= nodc, norm= norm), norm= norm);
  } else if (type == "spectX") {
    a= dvv_AmplSpectr(dvv_PSD(f, wndw= "hanning", nodc= nodc, norm= norm), norm= norm);
    a.z_unit = "!mm";
    if(norm == 2) a.z_unit = "!mm-!mm";
    a.z_label = "Displacement";
    a.data= &((*a.data)*_WAV/2/pi);
    return a;
  } else if (type == "spectV") {
    a= dvv_AmplSpectr(dvv_PSD(f, wndw= "hanning", nodc= nodc, norm= norm), norm= norm);
    a.z_unit = "m/s";
    if(norm == 2) a.z_unit = "m/s-!mm";
    a.z_label = "Velocity per mode";
    a.data= &((*a.data)*_VPF/2/pi*1000.0);
    return a;
  } else if (type == "shifted") {
    print, "dvv_fromCX, shifted", _FRINGE_MODE;
    xg= (*f.xscale)(,-:1:f.ny);
    yg= (*f.yscale)(-:1:f.nx,);
    md= exp(-(0.+1.i)*2*pi*(xg*_FRINGE_MODE(1) + yg*_FRINGE_MODE(2)));
    f.data= &((*f.data)*md);
    return f;
  } else error, "need to specify operation type.";
}

func dvv_mkNoiseModel(img, rmsval, mag, model)
/* DOCUMENT dvv_mkNoiseModel, img, rmsval, model

   Applies a passband filter to a white spectrum to create a spectrum
   matching a given noise model


   The distribution is scaled so that its integral (summed over
   all pixel elements) is unity.

   Arguments:
     img
     model
     params
     
   SEE ALSO:
     dvv_mkBox, dvv_mkCross
 */
{
  local cx, cf, psd0, rms1, nm;
  
  cx= img_copy(img);
  cx.data= &random_n(dimsof(*cx.data));
  psd0= dvv_PSD(cx);
  nm= model(psd0, mag);
  cf= dvv_Filter(cx, nm);
  rms1= dvv_rms(cf);
  cf.data= &(*cf.data*rmsval/rms1);
  return cf;
}

func dvv_OHRVNoiseModel(psd, mag, rmsval=)
/* DOCUMENT     dvv_OHRVNoiseModel(psd, mag, rmsval=);
           -or- dvv_OHRVNoiseModel(f, mag, rmsval=);

   Returns the frequency-dependent rms noise amplitude with magnitude
   matching the noise floor observed in the OHRV interferometer
   system.  In the first calling form the first argument represents a
   2D power spectral density (PSD) structure that is used as a
   template to fill in the rms noise amplitude for all 2D frequency
   components.

   In the second calling form the spectral frequencies are given
   by the first argument.

   An optional argument mag is used to scale the frequencies, such
   that  f -> f/mag
   
   SEE ALSO:
 */
{
  aa= [2.40051,4.36256,2.37368,30.1376];

  if(is_void(mag)) mag= 1.0;
  if (structof(psd) == IMG_DAT) {
    nm= img_copy(psd);
    xg= (*nm.xscale)(,-:1:nm.ny);
    yg= (*nm.yscale)(-:1:nm.nx,);
    ff= sqrt(xg*xg + yg*yg)/abs(mag);
    nm.data= &(aa(1)*exp(-abs(aa(2))*ff) + aa(3)*exp(-abs(aa(4))*ff));
    //Zero out the Nyquist-frequency terms for even-numbered
    //meshes
    //    if(nm.nx%2 == 0) (*nm.data)(1,)= 0+0i;
    //    if(nm.ny%2 == 0) (*nm.data)(,1)= 0+0i;
    if(!is_void(rmsval)) {
      rmsv= dvv_AMP_rms(nm, norm= 1);
      nm.data= &(*nm.data*rmsval/rmsv);
    }
    return nm;
  } else {
    ff= psd/abs(mag);
    nm= aa(1)*exp(-abs(aa(2))*ff) + aa(3)*exp(-abs(aa(4))*ff);
    if(!is_void(rmsval)) {
      rmsv= dvv_AMP_1D_rms([ff, nm], norm= 1, qavg= 1);
      nm*= rmsval/rmsv;
    }
    return nm;
  }
}

/* White spectrum */
func dvv_White(N, L, stdev)
/* DOCUMENT      dvv_White, N, L, stdev
            -or- dvv_White, [Nx, Ny], L, stdev
            -or- dvv_White, [Nx, Ny], [Lx, Ly], stdev
     
   SEE ALSO:
 */
{
  p= IMG_DAT();
  
  if(numberof(N) == 1) {Nx= N; Ny= N;}
  else {Nx= N(1); Ny= N(2);}

  if(numberof(L) == 1) {Lx= L; Ly= L;}
  else {Lx= L(1); Ly= L(2);}

  p.nx= Nx;
  p.ny= Ny;
  p.xscale= &span(-Lx/2., Lx/2., Nx);
  p.yscale= &span(-Ly/2., Ly/2., Ny);
  x= (*p.xscale)(,-:1:p.ny);
  y= (*p.yscale)(-:1:p.nx,);
  p.data= &(stdev*random_n([2,p.nx,p.ny]));
  p.x_label="Position";
  p.y_label="Position";
  p.z_label="Phase";
  p.x_unit="pixel";
  p.y_unit="pixel";
  p.z_unit="rad";
  return p;
}
