/*
  DVVPROCESS.I
  
  Definitions for processing dVV/OHRV data sets

  $Id$

  P.M. Celliers
  LLNL
  December, 2007

*/

/*

  dvv_Quadrature   - produce a quadrature pair from a frameset
  dvv_QuadFit      - use fitting methods to extract the quadrature pair
  dvv_Phase        - compute phase from quadrature pair
  dvv_GhostCorrect - attempt to determine a correction for ghost reflections
  dvv_Ampl         - return fringing amplitude from a frame set
  dvv_Contrast     - return fringing contrast from a frame set
  dvv_ContrastScan
  dvv_PhaseScan
  dvv_Intensity    - return sum of the 4-frame intensity
  dvv_SetPA_Maps
  dvv_MapPhaseDiff
  dvv_GetPhaseDiff
  dvv_mkMap
  dvv_mkSurf
  dvv_FitPhase
  dvv_fitProfile
  dvv_Vel
  dvv_ShiftShock
  dvv_ImgPhase
  dvv_ImgMag
  dvv_FieldShift
  dvv_optContrast
  dvv_optPhase
  dvv_optLissajous
  dvv_plPhasePoints
  dvv_plQuadrature
  dvv_jumpList
  dvv_adjustPhase
  dvv_delPhi
  dvv_PhaseAmplScan
  dvv_OffsetScan
  
 */

func dvv_Quadrature(frdat, params, &weight, highpass=, nocorr=, prefilter=,
                    filtered=, wndw=, reg=, drop1s=, method=, median_filter=)
/* DOCUMENT dvv_Quadrature(frdat, params, weight, highpass=, nocorr=, prefilter=,
                      filtered=, wndw=, reg=, drop1s=, method=, median_filter=)

  Processes a fringing frameset to produce a quadrature pair corrected
  for the amplitude and phase maps determined by dvv_SetPA_Maps

  ARGUMENTS:
    frdat - fringe frame set, registered, dewarped & flat-fielded
    params - Lissajous parameters (or null), as determined
             by dvv_PhaseScan or dvv_ContrastScan
    weight - [output] weight function for unwrapping

  KEYWORDS:
    highpass=
    nocorr=
    filtered=  string containing s-channel filter def (see dvv_GenFilterMaskSet)
    wndw=
    reg=
    drop1S=
    method=  0 or undefined - compute phase from standard sdiff/pdiff pair
             1 from 6 sdiff/pdiff pairs compute phase by linear fit to
               (dx, dy) data and compute angle from the slope.
             2 from 6 sdiff/pdiff pairs normalized using stot, ptot and (stot+ptot)
               compute phase using arctan and take average angle
    median_filter= apply a median_filter
           
  SEE ALSO:
    dvv_ContrastScan, dvv_SetPA_Maps
*/
{
  //  stot= img_add(frdat.ch1S, frdat.ch2S, average= 1);
  //  ptot= img_add(frdat.ch1P, frdat.ch2P, average= 1);
  extern _DVV_PHASE_MAP, _DVV_AMPL_MAP, _DVV_XOFF_MAP, _DVV_YOFF_MAP;
  extern _CSPHI, _SNPHI, _CSPHI0, _SNPHI0, _CSPHI1, _SNPHI1;
  extern _DVV_CH1S_FILTER, _DVV_CH2S_FILTER;
  extern _DVV_CH1P_FILTER, _DVV_CH2P_FILTER;
  extern _DVV_CMMN_FILTER, _DVV_CMMN_SUM_FILTER;
  extern _DVV_CHP_FILTER, _DVV_CHS_FILTER;
  extern _DVV_HF_FILTER, _DVV_SELECT_MODE;
  extern _ms, _mp;
  extern _DC_LIMIT;
  extern sdiff, pdiff;

  extern delta;

  if(!is_void(params)) {
    if(numberof(params) != 4) error, "Invalid parameter set.";
    delta= params(1);
    amp_ratio= params(2);
    xoff= params(3);
    yoff= params(4);
  }

  nx= frdat.ch1S.nx;
  ny= frdat.ch1S.ny;
  xscl= *frdat.ch1S.xscale;
  yscl= *frdat.ch1S.yscale;
  
  if(is_void(delta) && !is_void(_DVV_PHASE_MAP)) {
    delta= img_copy(_DVV_PHASE_MAP);
    //if (frdat.binning > 1 && frdat.ch1S.nx != _DVV_PHASE_MAP.nx)
    //  delta= img_rebin(delta, frdat.binning, average= 1);
    if (nx != _DVV_PHASE_MAP.nx) delta= img_resample(delta, xscl, yscl);
  }
  
  if(is_void(amp_ratio) && !is_void(_DVV_AMPL_MAP)) {
    amp_ratio= img_copy(_DVV_AMPL_MAP);
    //if (frdat.binning > 1 && frdat.ch1S.nx != _DVV_AMPL_MAP.nx)
    //  amp_ratio= img_rebin(amp_ratio, frdat.binning, average= 1);
    if (nx != _DVV_AMPL_MAP.nx) amp_ratio= img_resample(amp_ratio, xscl, yscl);
  }

  if(is_void(xoff) && !is_void(_DVV_XOFF_MAP)) {
    xoff= img_copy(_DVV_XOFF_MAP);
    //if (frdat.binning > 1 && frdat.ch1S.nx != _DVV_XOFF_MAP.nx)
    //  xoff= img_rebin(xoff, frdat.binning, average= 1);
    if (nx != _DVV_XOFF_MAP.nx) xoff= img_resample(xoff, xscl, yscl);

  }

  if(is_void(yoff) && !is_void(_DVV_YOFF_MAP)) {
    yoff= img_copy(_DVV_YOFF_MAP);
    //if (frdat.binning > 1 && frdat.ch1S.nx != _DVV_YOFF_MAP.nx)
    //  yoff= img_rebin(yoff, frdat.binning, average= 1);
    if (nx != _DVV_YOFF_MAP.nx) yoff= img_resample(yoff, xscl, yscl);
  }

  if(is_void(reg)) {
    ch1S= img_copy(frdat.ch1S);
    ch2S= img_copy(frdat.ch2S);
    ch1P= img_copy(frdat.ch1P);
    ch2P= img_copy(frdat.ch2P);
  } else {
    ch1S= img_extract(frdat.ch1S, reg);
    ch2S= img_extract(frdat.ch2S, reg);
    ch1P= img_extract(frdat.ch1P, reg);
    ch2P= img_extract(frdat.ch2P, reg);
  }

  //Apply a prefiltering operation to each frame prior to phase processing
  if (prefilter) {
    pf= parse_line(prefilter);
    print, "pf = ", pf;
    if (pf(1) == "smooth") {
      fwhm= tonum(pf(2));
      psf= dvv_mkGaussian(ch1S, fwhm);
      ch1S= dvv_Smooth(ch1S, psf);
      ch2S= dvv_Smooth(ch2S, psf);
      ch1P= dvv_Smooth(ch1P, psf);
      ch2P= dvv_Smooth(ch2P, psf);
    } else {
      prefilt= dvv_GenFilterMask(ch1S, prefilter);
      ch1S= dvv_Filter(ch1S, prefilt, wndw= wndw);
      ch2S= dvv_Filter(ch2S, prefilt, wndw= wndw);
      ch1P= dvv_Filter(ch1P, prefilt, wndw= wndw);
      ch2P= dvv_Filter(ch2P, prefilt, wndw= wndw);
    }
  }

  stot= img_add(ch1S, ch2S, average= 1);
  ptot= img_add(ch1P, ch2P, average= 1);
  tot= img_add(stot, ptot, average= 1);

  if (is_void(method) || method == 0) {

    // Channel 1S is corrupted, use only the other 3 channels
    if (drop1s) {
      stot= img_copy(ptot);
      ch1S= img_sub(stot, ch2S);
    }
  
    _DVV_HF_FILTER= _(_DVV_CHS_FILTER, _DVV_CHP_FILTER);
  
    if(filtered && (!is_void(_DVV_CMMN_FILTER) || !is_void(_DVV_HF_FILTER))) {
      dvv_SetupCommonSumFilter;
      _csum= dvv_GenFilterMaskSet(ch1S, _(_DVV_CMMN_SUM_FILTER, _DVV_HF_FILTER));
      stot= dvv_Filter(stot, _csum, wndw= wndw);
      ptot= dvv_Filter(ptot, _csum, wndw= wndw);
    }

    //s-pol is cos(phi)
    //sdiff= img_sub(ch2S, ch1S);
    sdiff= img_div(img_sub(ch2S, ch1S), stot);
    //Apply filters to s-channel data
    _CSPHI0= img_copy(sdiff);
    if (filtered) {
      if (is_void(_DVV_CH1S_FILTER) && is_void(_DVV_CH2S_FILTER) &&
          is_void(_DVV_CMMN_FILTER) && is_void(_DVV_HF_FILTER) && is_void(_DVV_SELECT_MODE))
        print, "WARNING: Ch1S, Ch2s Filters are not specified, s-output NOT filtered ...";
      else {
        _ms= dvv_GenFilterMaskSet(sdiff, _(_DVV_CH1S_FILTER, _DVV_CH2S_FILTER,
                                           _DVV_HF_FILTER, _DVV_CMMN_FILTER, _DVV_SELECT_MODE));
        //    m= dvv_GenFilterMaskSet(sdiff, s_filter);
        sdiff= dvv_Filter(sdiff, _ms, wndw= wndw);
        //dvv_DisplayImg, csphi, w= 1, -1, 1;
      }
    }

    //p-pol is cos(phi + delta)
    //pdiff= img_sub(ch2P, ch1P);
    pdiff= img_div(img_sub(ch2P, ch1P), ptot);
    _SNPHI0= img_copy(pdiff);
    //Apply filters to p-channel data
    if (filtered) {
      if (is_void(_DVV_CH1P_FILTER) && is_void(_DVV_CH2P_FILTER) &&
          is_void(_DVV_CMMN_FILTER) && is_void(_DVV_HF_FILTER) && is_void(_DVV_SELECT_MODE))
        print, "WARNING: Ch1p, Ch2p Filters are not specified, p-output NOT filtered ...";
      else {
        _mp= dvv_GenFilterMaskSet(pdiff, _(_DVV_CH1P_FILTER, _DVV_CH2P_FILTER,
                                           _DVV_HF_FILTER, _DVV_CMMN_FILTER, _DVV_SELECT_MODE));
        pdiff= dvv_Filter(pdiff, _mp, wndw= wndw);
      }
    }

    //Apply a high pass filter to subtract phase offsets
    if (!is_void(highpass)) {
      print, "Applying high pass filter ..., highpass, dc_limit =", highpass, _DC_LIMIT;
      if(is_void(_DC_LIMIT)) _DC_LIMIT= highpass;
      mask= dvv_GenFilterMask(csphi, [_DC_LIMIT], type= "high_pass_SQ");
      sdiff= dvv_Filter(sdiff, mask);
      pdiff= dvv_Filter(pdiff, mask);
    }
    
  } else if (method == 1 || method == 2) {
    //Apply a prefiltering operation to each frame prior to phase processing
    print, "dvv_Quadrature - fscopy";
    frcpy= dvv_FSCopy(frdat);
    if (prefilter) {
      pf= parse_line(prefilter);
      print, "pf = ", pf;
      if (pf(1) == "smooth") {
        fwhm= tonum(pf(2));
        psf= dvv_mkGaussian(ch1S, fwhm);
        ch1S= dvv_Smooth(ch1S, psf);
        ch2S= dvv_Smooth(ch2S, psf);
        ch1P= dvv_Smooth(ch1P, psf);
        ch2P= dvv_Smooth(ch2P, psf);
      } else {
        prefilt= dvv_GenFilterMask(ch1S, prefilter);
        ch1S= dvv_Filter(ch1S, prefilt, wndw= wndw);
        ch2S= dvv_Filter(ch2S, prefilt, wndw= wndw);
        ch1P= dvv_Filter(ch1P, prefilt, wndw= wndw);
        ch2P= dvv_Filter(ch2P, prefilt, wndw= wndw);
      }
      frcpy.ch1S= ch1S;
      frcpy.ch2S= ch2S;
      frcpy.ch1P= ch1P;
      frcpy.ch2P= ch2P;
    }
    spdiff= dvv_QuadFit(frcpy, reg= reg, method= method, median_filter= median_filter);
    restore, spdiff, pdiff, sdiff;
  } else {
    error, "Invalid method.";
  }
    
  _CSPHI= img_copy(sdiff);
  _SNPHI= img_copy(pdiff);

  // if(!is_void(amp_ratio)) {
  //   if(struct_type(amp_ratio) == "IMG_DAT") {
  //     if(is_void(reg)) ar= img_data(amp_ratio);
  //     else ar= img_data(img_extract(amp_ratio, reg));
  //     sdiff.data= &((*sdiff.data)*ar);
  //   } else if(numberof(amp_ratio) > 1) {
  //     xg= (*sdiff.xscale)(,-:1:sdiff.ny);
  //     yg= (*sdiff.yscale)(-:1:sdiff.nx,);
  //     ar= polysurf(amp_ratio, yg, xg);
  //     *sdiff.data*= ar;
  //   } else {
  //     ar= amp_ratio;
  //     *sdiff.data*= ar;
  //   }
  // }

  xg= img_grid(sdiff, 1);
  yg= img_grid(sdiff, 2);
  
  if(struct_type(delta) == "IMG_DAT") {
    if (is_void(reg)) ddlta= img_data(delta);
    else ddlta= img_data(img_extract(delta, reg));
  } else if (numberof(delta) > 1) {
    //xg= (*sdiff.xscale)(,-:1:sdiff.ny);
    //yg= (*sdiff.yscale)(-:1:sdiff.nx,);
    ddlta= polysurf(delta, yg, xg);
  } else ddlta= delta;

  if (!is_void(xoff)) {
    if (struct_type(xoff) == "IMG_DAT") {
      //if (is_void(reg)) *sdiff.data-= img_data(xoff)*img_data(stot);
      if (is_void(reg)) *sdiff.data-= img_data(xoff);
      //else *sdiff.data-= img_data(img_extract(xoff,reg))*img_data(img_extract(stot,reg));
      else *sdiff.data-= img_data(img_extract(xoff,reg));
      //} else if (numberof(xoff) == 1 && xoff != 0.0) *sdiff.data-= xoff*(*stot.data);
    } else if (numberof(xoff) == 1 && xoff != 0.0) *sdiff.data-= xoff;
  }
  
  if (!is_void(yoff)) {
    if (struct_type(yoff) == "IMG_DAT") {
      //if (is_void(reg)) *pdiff.data-= img_data(yoff)*img_data(ptot);
      if (is_void(reg)) *pdiff.data-= img_data(yoff);
      //else *pdiff.data-= img_data(img_extract(yoff,reg))*img_data(img_extract(ptot, reg));
      else *pdiff.data-= img_data(img_extract(yoff,reg));
      //} else if (numberof(yoff) == 1 && yoff != 0.0) *pdiff.data-= yoff*(*ptot.data);
    } else if (numberof(yoff) == 1 && yoff != 0.0) *pdiff.data-= yoff;
  }

  if(!is_void(amp_ratio)) {
    if(struct_type(amp_ratio) == "IMG_DAT") {
      if(is_void(reg)) ar= img_data(amp_ratio);
      else ar= img_data(img_extract(amp_ratio, reg));
      sdiff.data= &((*sdiff.data)*ar);
    } else if(numberof(amp_ratio) > 1) {
      //xg= (*sdiff.xscale)(,-:1:sdiff.ny);
      //yg= (*sdiff.yscale)(-:1:sdiff.nx,);
      ar= polysurf(amp_ratio, yg, xg);
      *sdiff.data*= ar;
    } else {
      ar= amp_ratio;
      *sdiff.data*= ar;
    }
  }
  
  //Construct sin(phi) from the data
  csphi= img_copy(sdiff);
  if (nocorr)
    snphi= img_mul(pdiff, -1.0);
  else
    snphi= img_div(img_sub(pdiff, img_mul(sdiff, cos(ddlta))), sin(ddlta));
  //snphi= img_div(img_sub(img_mul(sdiff, cos(ddlta)), pdiff), sin(ddlta));

  //Subtract estimated ghost signals here
  //   if(!is_void(ghost_params)) {  
  //     //Get amplitude [ gh(1) ] and phase [ gh(2) ] of the ghost image
  //     gh= dvv_Ghost(ghost_params, reg= reg);
  //     //Adjust quadrature components by subtracting the ghost terms
  //     csphi= img_sub(csphi,img_mul(img_operate(gh(2), cos), gh(1)));
  //     snphi= img_sub(snphi,img_mul(img_operate(gh(2), sin), gh(1)));
  //   }

  _CSPHI1= img_copy(csphi);
  _SNPHI1= img_copy(snphi);

  return save(csphi= _CSPHI1, snphi= _SNPHI1, stot= stot, ptot= ptot, tot= tot);  
}

func dvv_QuadFit(fs, reg=, method=, median_filter=)
/* DOCUMENT dvv_QuadFit, fs, reg=, method=, median_filter=

   Computes the quadrature phase from an average determined from a
   set of estimates of the quadrature pair values from a given pixel.

   Both methods in this algorithm begin with a set of data values
   determined for each pixel in the frameset: xa= N(ch2P)-1, xb=
   N(ch1P)-1, ya= N(ch1P)-1, yb= N(ch2P)-1.  Here the signals are
   normalized N(frame) by dividing by the average total intensity in
   the frame.  The average total intensity is estimated 3 ways: By
   summing the 2 s-channel frames, summing the 2 p-channel frames and
   summing all 4 frames.  Thus there are 3 estimates of the normalized
   intensity for each point on each frame.  From this set of values
   there are six "quadrature pair estimates", i.e. [xa(1), ya(1)],
   [xa(2), ya(2)], ...  [xb(3), yb(3)].
   
   KEYWORDS:
     reg=    Perform fit only over the data within reg
     method= 1 - (default) use a least quares slope fit to all QP estimates
             2 - use angles from individual measurements
     median_filter=
             apply a median filter to the intensity-normalized data; this
             will reduce the effect of phase dropouts from low intensity
             regions of the data
     
   SEE ALSO:
 */
{
  local tot, stot, ptot;

  if (!is_void(reg)) {
    print, "dvv_QuadFit - fscopy";
    fs= dvv_FSCopy(fs, reg= reg);
  }
  tot= dvv_frameSum(fs, stot, ptot, average= 1);
  ddims= img_dims(tot);
  if(is_void(method)) method= 2;

  if(method == 1) {
    //Method 1 computes the phase from a linear fit to determine the slope
    //of a data set given by [(xa, ya), (0,0), (xb, yb)]
    if(median_filter) {
      svec= [img_data(img_median(img_div(fs.ch2S, stot))),
             img_data(img_median(img_div(fs.ch2S, ptot))),
             img_data(img_median(img_div(fs.ch2S, tot))),
             img_data(img_median(img_div(fs.ch1S, stot))),
             img_data(img_median(img_div(fs.ch1S, ptot))),
             img_data(img_median(img_div(fs.ch1S, tot))),
             array(1.,ddims), array(1.,ddims), array(1.,ddims)]-1;
    
      pvec= [img_data(img_median(img_div(fs.ch1P, stot))),
             img_data(img_median(img_div(fs.ch1P, ptot))),
             img_data(img_median(img_div(fs.ch1P, tot))),
             img_data(img_median(img_div(fs.ch2P, stot))),
             img_data(img_median(img_div(fs.ch2P, ptot))),
             img_data(img_median(img_div(fs.ch2P, tot))),
             array(1.,ddims), array(1.,ddims), array(1.,ddims)]-1;
    } else {
      svec= [img_data(img_div(fs.ch2S, stot)),
             img_data(img_div(fs.ch2S, ptot)),
             img_data(img_div(fs.ch2S, tot)),
             img_data(img_div(fs.ch1S, stot)),
             img_data(img_div(fs.ch1S, ptot)),
             img_data(img_div(fs.ch1S, tot)),
             array(1.,ddims), array(1.,ddims), array(1.,ddims)]-1;
    
      pvec= [img_data(img_div(fs.ch1P, stot)),
             img_data(img_div(fs.ch1P, ptot)),
             img_data(img_div(fs.ch1P, tot)),
             img_data(img_div(fs.ch2P, stot)),
             img_data(img_div(fs.ch2P, ptot)),
             img_data(img_div(fs.ch2P, tot)),
             array(1.,ddims), array(1.,ddims), array(1.,ddims)]-1;
    }

    //Construct a linear least squares estimate of the quadrature pair
    //slope (following Bevington, section 6-3, equations 6-9, p. 104).
    //The fit is applied to y versus x and x versus y, and the results
    //from both fits are averaged.
    //print, "dimsof svec, pvec:: ", dimsof(svec);
    N=dimsof(svec)(0);
    //print, "dvv_QuadFit::N = ", N;
    DDx=N*(svec*svec)(..,sum) - (svec(..,sum))^2;
    Mx= (N*(svec*pvec)(..,sum) - svec(..,sum)*pvec(..,sum))/DDx;
    
    DDy=N*(pvec*pvec)(..,sum) - (pvec(..,sum))^2;
    My= (N*(svec*pvec)(..,sum) - svec(..,sum)*pvec(..,sum))/DDy;

    // //Remove zero slope cases
    // cMx0= Mx == 0.0; wMx0= where(cMx0); wMx1= where(!cMx0);
    // cMy0= My == 0.0; wMy0= where(cMy0); wMy1= where(!cMy0);
    // if (is_array(wMy0))
    //   Mxavg= merge(100.0*sign(Mx(wMy0)), 0.5*(Mx(wMy1) + 1/My(wMy1)), cMy0);
    // else
    //   Mxavg= 0.5*(Mx + 1./My);
    // if (is_array(wMy0))
    //   Myavg= merge(100.0*sign(My(wMx0)), 0.5*(My(wMx1) + 1/Mx(wMx1)), cMx0);
    // else
    //   Myavg= 0.5*(My + 1./Mx);

    //Remove zero slope cases
    cMx0= Mx == 0.0; wMx0= where(cMx0); wMx1= where(!cMx0);
    cMy0= My == 0.0; wMy0= where(cMy0); wMy1= where(!cMy0);
    if (is_array(wMx0)) Mx(wMx0)= 1e-6*sign(My(wMx0));
    if (is_array(wMy0)) My(wMy0)= 1e-6*sign(Mx(wMy0));
    
    //Bias near-zero slope cases, if any
    cMx0= abs(Mx) < 0.05 & svec(..,rms) > pvec(..,rms); wMx0= where(cMx0); wMx1= where(!cMx0);
    cMy0= abs(My) < 0.05 & pvec(..,rms) > svec(..,rms); wMy0= where(cMy0); wMy1= where(!cMy0);
    if (is_array(wMx0)) Myavg= merge(1./Mx(wMx0), My(wMx1), cMx0);
    if (is_array(wMy0)) Mxavg= merge(1./My(wMy0), Mx(wMy1), cMy0);
    
    //Compute average slope base on Mx and My results
    Mxx = sign(Mxavg)*sqrt((Mxavg*Mxavg + 1./(Myavg*Myavg))/2.0);
    Myy = sign(Myavg)*sqrt((Myavg*Myavg + 1./(Mxavg*Mxavg))/2.0);
    
    //Use sign from X(Y) fit for |slope| > 1.0, from Y(X) fit for |slope| < 1.0
    cS1 = Mxx > 1.0; wCs1 = where(cS1); wCs0= where(!cS1);
    MM= merge(1./Myy(wCs1),Mxx(wCs0),cS1);

    /*  this code works but produces outliers ....
    sgnxx= sign(svec(..,1:3)(..,avg) - svec(..,4:6)(..,avg));
    sgnyy= sign(pvec(..,1:3)(..,avg) - pvec(..,4:6)(..,avg));

    sdiff= sgnxx*1.0/sqrt(1.+MM*MM);
    pdiff= MM*sdiff;
    */

    //This code produces outliers -- need to investigate ...
    sgnxx= sign(svec(..,1:3)(..,avg) - svec(..,4:6)(..,avg));
    sdiffx= sgnxx(wCs0)/sqrt(1.+Mxx(wCs0)*Mxx(wCs0));
    pdiffx= Mxx(wCs0)*sdiffx;

    sgnyy= sign(pvec(..,1:3)(..,avg) - pvec(..,4:6)(..,avg));
    pdiffy= sgnyy(wCs1)/sqrt(1.+Myy(wCs1)*Myy(wCs1));
    sdiffy= Myy(wCs1)*pdiffy;

    sdiff= merge(sdiffy, sdiffx, cS1);
    pdiff= merge(pdiffy, pdiffx, cS1);
    
    return save(sdiff= img_copy(tot, data= sdiff*img_data(tot)),
                pdiff= img_copy(tot, data= -pdiff*img_data(tot)),
                angle= img_copy(tot, data= atan(pdiff, sdiff)),
                // cS1= img_copy(tot, data= cS1),
                // DDx= img_copy(tot, data= DDx),
                // DDy= img_copy(tot, data= DDy),
                // Mx= img_copy(tot, data= Mx),
                // My= img_copy(tot, data= My),
                // Mxx= img_copy(tot, data= Mxx),
                // Myy= img_copy(tot, data= Myy),
                // sgnxx= img_copy(tot, data= sgnxx),
                // sgnyy= img_copy(tot, data= sgnyy),
                MM= img_copy(tot, data= MM));
    
  } else if (method == 2) {
    //Method 2 computes the phase from the atan function applied to each of the
    //six quadrature pair estimates.
    if (median_filter) {
      svec= [img_data(img_median(img_div(fs.ch2S, stot))),
             img_data(img_median(img_div(fs.ch2S, ptot))),
             img_data(img_median(img_div(fs.ch2S, tot))),
             img_data(img_median(img_div(fs.ch1S, stot))),
             img_data(img_median(img_div(fs.ch1S, ptot))),
             img_data(img_median(img_div(fs.ch1S, tot)))]-1;
    
      pvec= [img_data(img_median(img_div(fs.ch1P, stot))),
             img_data(img_median(img_div(fs.ch1P, ptot))),
             img_data(img_median(img_div(fs.ch1P, tot))),
             img_data(img_median(img_div(fs.ch2P, stot))),
             img_data(img_median(img_div(fs.ch2P, ptot))),
             img_data(img_median(img_div(fs.ch2P, tot)))]-1;
    } else {
      svec= [img_data(img_div(fs.ch2S, stot)),
             img_data(img_div(fs.ch2S, ptot)),
             img_data(img_div(fs.ch2S, tot)),
             img_data(img_div(fs.ch1S, stot)),
             img_data(img_div(fs.ch1S, ptot)),
             img_data(img_div(fs.ch1S, tot))]-1;
    
      pvec= [img_data(img_div(fs.ch1P, stot)),
             img_data(img_div(fs.ch1P, ptot)),
             img_data(img_div(fs.ch1P, tot)),
             img_data(img_div(fs.ch2P, stot)),
             img_data(img_div(fs.ch2P, ptot)),
             img_data(img_div(fs.ch2P, tot))]-1;
    }

    angles= atan(pvec,svec);
    angles= _wrap(angles + ([2,2,2,1,1,1](-:1:ddims(2),-:1:ddims(3),)*pi)) + 2*pi;
    //da= median(angles,3)(..,-:1:6) - angles;
    //ww= where(abs(da) > pi);
    //if(is_array(ww)) angles(ww)-= 2*pi*sign(da(ww));
    aa= angles(..,avg)-(2*pi);
    return save(sdiff= img_copy(tot, data= cos(aa)*img_data(tot)),
                pdiff= img_copy(tot, data= -sin(aa)*img_data(tot)),
                angle= img_copy(tot, data= aa),
                MM= img_copy(tot, data= tan(aa))); 
  }
}

func dvv_Phase(frdat, params, &weight, unwrap=, cplex=, highpass=, prefilter=,
               shot=, filtered=, wndw=, reg=, invert=, drop1s=,
               method=, median_filter=)
/* DOCUMENT dvv_Phase(frdat, weight, unwrap=, cplex=, highpass=, prefilter=,
                shot=, filtered=, wndw=, reg=, invert=, drop1s=, method=,
                median_filter=

  Extracts and returns the wrapped phase from a 2D visar, 4-frame data
  set using the lissajous parameters params.  If the parameter weight
  is defined also computes a weighting function to be used for phase
  unwrapping.

  ARGUMENTS:
    frdat - fringe frame set, registered, dewarped & flat-fielded
    params - Lissajous parameters (or null), as determined
             by dvv_PhaseScan or dvv_ContrastScan
    weight - [output] weight function for unwrapping

  KEYWORDS:
    unwrap=
    cplex=
    highpass=
    prefilter=
    shot=
    filtered=  string containing s-channel filter def (see dvv_GenFilterMaskSet)
    wndw=
    reg=
    invert=    invert the phase
    // ghost_params= --removed--
    method=
    median_filter=
           
  SEE ALSO:
    dvv_ContrastScan, dvv_SetPA_Maps
*/
{
  //  stot= img_add(frdat.ch1S, frdat.ch2S, average= 1);
  //  ptot= img_add(frdat.ch1P, frdat.ch2P, average= 1);
  extern _DVV_PHASE_MAP, _DVV_AMPL_MAP, _DVV_XOFF_MAP, _DVV_YOFF_MAP;
  extern _CSPHI, _SNPHI, _CSPHI0, _SNPHI0, _CSPHI1, _SNPHI1;
  extern _DVV_CH1S_FILTER, _DVV_CH2S_FILTER;
  extern _DVV_CH1P_FILTER, _DVV_CH2P_FILTER;
  extern _DVV_CMMN_FILTER, _DVV_CMMN_SUM_FILTER;
  extern _DVV_CHP_FILTER, _DVV_CHS_FILTER;
  extern _DVV_HF_FILTER, _DVV_SELECT_MODE;
  extern _ms, _mp;
  extern _DC_LIMIT;

  extern delta;

  //message, "dvv_Phase::median_filter = ", median_filter;
  //message, "dvv_Phase::phase_method = ", method;
  
  qph= dvv_Quadrature(frdat, params, &weight, highpass= highpass, prefilter= prefilter,
                      filtered= filtered, wndw= wndw, reg= reg, drop1s= drop1s,
                      method= method, median_filter= median_filter);
  restore, qph, csphi, snphi, stot, ptot;
  
  if (invert) snphi= img_mul(snphi, -1.0);

  //Construct phase from arctangent
  phs= img_copy(csphi);
  if(is_void(cplex)) {
    phs.data= &atan(*snphi.data, *csphi.data);
    if(!is_void(shot)) phs.shotid=shot+"-phase";
    phs.z_unit= "rad";

    //Unwrap the phase
    if(!is_void(unwrap)) {
      phs= img_unwrap(phs, or= "h", dpass= 1, force= 1, midp= 1);
      *phs.data-= (*phs.data)(phs.nx/2, phs.ny/2);
    }
  } else {
    phs.data= &(img_data(csphi) + (0.0 + 1.0i) * img_data(snphi));
  }
  if(!is_void(shot)) phs.shotid=shot+"-phase";
  phs.z_unit="rad";
  phs.z_label="Phase ";

  //Construct weight function from log amplitude
  if(!is_void(weight)) {
    weight= img_add(stot, ptot);
    w= where(*weight.data <= 0.0);
    if(numberof(w) > 0) (*weight.data)(w)= 1e-6;
    weight.data= &log(*weight.data);
  }
  return phs;
}

func dvv_GhostCorrect(frdat, rfdat, phase, params, invert=, gphi=, glevel=, poff=, method=, median_filter=)
/* DOCUMENT frdat, rfdat, phase, params, gphi=, glevel=, poff=, method=, median_filter=

   Generates corrected data sets for 

   KEYWORDS:
     gphi=
     glevel=
     poff=    phase offset
     method=
     median_filter=
     
   SEE ALSO:
 */
{
  extern _DVV_FIELD;
  field= (*frdat.ch1S.xscale)(0)/1.1125;
  //field= _DVV_FIELD/1.1125;
  
  if (is_void(poff)) poff= 0.0;
  if (is_void(method)) method= 1;
  if (!is_void(median_filter)) medf= median_filter; else medf = 0;
  if (is_void(invert)) sgn= 1.0; else sgn= -1.0;
  //if (is_void(params)) params= _DVV_LISS_PARAMS;
  
  rr= img_zclip(dvv_Contrast(frdat, params), [0., 1.]);
  if (medf > 0) {print, "median filter :", medf; rr= img_median(rr, medf);}
  cc= img_zclip(dvv_Contrast(rfdat, params), [0., 1.]);
  refcontr= img_data(img_extract(cc, REGION(x1= -field/4, x2= field/4, y1= -field, y2= -field/2)))(avg, avg);
  cc= img_copy(cc, data= array(refcontr, img_dims(cc)));
  ii= dvv_Intensity(frdat);

  cphi= img_operate(img_add(phase,poff), cos); //cosine of data phase
  sphi= img_operate(img_add(phase,poff), sin); //sine of data phase
  c2phi= img_operate(img_mul(img_add(phase, poff), 2), cos);  //cos(2 x  data phase)

  if (method == 1) {
    num1a= img_add(img_mul(cc, cc), img_mul(rr, rr, -1));  //numerator of A expression
    num1f= img_add(img_mul(cphi, cc, rr, 2), img_mul(cc, cc, -1),
                   img_mul(rr, rr, c2phi, -1));            //Denominator of both expressions
    den1= img_add(img_mul(cc, cc), img_mul(rr, rr), img_mul(cc, rr, cphi, -2));

    aa= img_div(num1a, den1);          //A expression
    cphis= img_div(num1f, den1);       //cosine of the phase
    acphis= img_operate(img_zclip(cphis, [-1.0, 1.0]), acos); //Arccosine of cos(phase)
    sphisd= img_data(phase);            //Begin construct sin(phase) from unwrapped phase
    pn= _wrap(sphisd) < 0;               //We need to wrap the phase  
    wn= where(pn);  wp= where(!pn);
    if(is_array(wp)) sphisd(wp)= sgn*sin(img_data(acphis)(wp));  //Use +ve result if wrapped(phi) > 0
    if(is_array(wn)) sphisd(wn)= -sgn*sin(img_data(acphis)(wn)); //Use -ve result if wrapped(phi) < 0
    sphis= img_copy(cphis, data= sphisd);
    xx= img_zclip(img_div(rr, cc), [0.0, 1.0]);
    
  } else if (method == 2) {
    xx= img_zclip(img_div(rr, cc), [0.0, 1.0]);
    xxsq= img_mul(xx, xx);
    cphi2x= img_mul(xx, cphi, 2);
    num3a= img_mul(img_add(xxsq, -1), -1);
    num3f= img_add(cphi2x, img_mul(xxsq, c2phi, -1), -1);
    den3= img_sub(img_add(xxsq, 1), cphi2x);
    
    aa=img_div(num3a, den3);          //A expression
    cphis= img_div(num3f, den3);       //cosine of the phase
    acphis= img_operate(img_zclip(cphis, [-1.0, 1.0]), acos); //Arccosine of cos(phase)
    sphisd= img_data(phase);            //Construct sin(phase)
    pn= _wrap(sphisd) < 0;
    wn= where(pn); wp= where(!pn);
    if(is_array(wp)) sphisd(wp)= sgn*sin(img_data(acphis)(wp)); //Use +ve result if wrapped(phi) > 0
    if(is_array(wn)) sphisd(wn)= -sgn*sin(img_data(acphis)(wn));//Use  -ve result if wrapped(phi) < 0
    sphis= img_copy(cphis, data= sphisd);

  } else if (method == 3) {
    xx= img_zclip(img_div(rr, cc), [0.0, 1.0]);
    xxsq= img_mul(xx, xx);
    cphi2x= img_mul(xx, cphi, 2);
    num3a= img_mul(img_add(xxsq, -1), -1);
    num3f= img_add(cphi2x, img_mul(xxsq, c2phi, -1), -1);
    den3= img_sub(img_add(xxsq, 1), cphi2x);
    
    aa=img_div(num3a, den3);          //A expression
    cphis= img_div(num3f, den3);       //cosine of the phase
    acphis= img_operate(img_zclip(cphis, [-1.0, 1.0]), acos); //Arccosine of cos(phase)
    sphisd= img_data(phase);            //Construct sin(phase)
    pn= wrap(sphisd) < 0;
    wn= where(pn); wp= where(!pn);
    if(is_array(wp)) sphisd(wp)= -sgn*sin(img_data(acphis)(wp)); //Use -ve result if wrapped(phi) > 0
    if(is_array(wn)) sphisd(wn)= sgn*sin(img_data(acphis)(wn)); //Use +ve result if wrapped(phi) < 0
    sphis= img_copy(cphis, data= sphisd);
  }

  //dvv_DisplayImg, acphis, w= 60;
  //dvv_DisplayImg, cphis, w= 61;
  //dvv_DisplayImg, sphis, w= 62;
  phis= img_copy(cphi, data=atan(img_data(sphis), img_data(cphis)));
  if(is_void(gphi)) print, ">>> dvv_GhostCorrect: gphi is VOID";

  //Mask out areas where the phase is very close to the null phase
  //within the threshold given by gphi.
  if (!is_void(gphi) && gphi > 0.0) {
    write, format=">>> dvv_GhostCorrect: ****** Applying null phase correction: gphi=%12.5e\n", gphi;
    if(is_void(glevel)) glevel= 0.0;
    write, format=">>> dvv_GhostCorrect: ****** glevel=%12.5e\n", glevel;
    uwd= img_data(img_wrap(phase));
    phisd= img_data(img_wrap(phis));
    //fx= where(-gphi < phisd & phisd < gphi);
    fx= where(-gphi < uwd & uwd < gphi);
    phisd(fx)= uwd(fx);
    phis= img_copy(phis, data= phisd);
    ad= img_data(aa);
    if (glevel != 0.0) {
      gg= img_copy(phis, data= glevel);
      afix= img_data(img_div(gg, img_sub(ii, gg)));
      ad(fx)= afix(fx);
    } else {
      ad(fx)= 0.0;
    }
    aa= img_copy(aa, data= ad);
  } else {
    print, ">>> dvv_GhostCorrect: ****** No null phase correction *******>>>>>>>";
  }

  isig= img_div(ii, img_add(aa, 1));
  ighost= img_div(img_mul(ii, aa), img_add(aa, 1));
  isig.shotid+="signal";
  ighost.shotid+="ghost";
  return [phis, isig, ighost, aa, xx];
}

func dvv_Ampl(frdat, params, shot=, reg=, highpass=, prefilter=,
               filtered=, wndw=, invert=, drop1s=,
               method=, median_filter=)
/* DOCUMENT ampl= dvv_Ampl(frdat, shot=, reg=)
      -or-  ampl= dvv_Ampl(frdat, params, shot=, reg=, highpass=,
               prefilter=, filtered=, wndw=, invert=, drop1s=,
               method=, median_filter=)

  Returns an IMG_DAT data structure containing the amplitude of the
  fringing component of the signal

  frdat - fringe frame set, registered, dewarped & flat-fielded
  params - 

  KEYWORDS:
    shot=
    reg=
           
  SEE ALSO:
    dvv_Phase, dvv_SetPA_Maps, dvv_Contrast
*/
{
  extern _DVV_PHASE_MAP, _DVV_AMPL_MAP;
  extern _CSPHI, _SNPHI, _CSPHI1, _SNPHI1;

  qph= dvv_Quadrature(frdat, params, &weight, highpass= highpass, prefilter= prefilter,
                      filtered= filtered, wndw= wndw, reg= reg, drop1s= drop1s,
                      method= method, median_filter= median_filter);
  restore, qph, csphi, snphi, tot;
  
  //Construct amplitude from the data
  tt= img_add(img_mul(csphi, csphi), img_mul(snphi, snphi));
  tt.data= &(sqrt(img_data(tt))*img_data(tot));
  if(!is_void(shot)) tt.shotid=shot+"-amplitude";
  tt.z_unit= "counts";
  tt.z_label= "Amplitude ";
  return tt;
}

func dvv_Contrast(frdat, params, nocorr=,  shot=, reg=, highpass=,
                  prefilter=, filtered=, wndw=, invert=, drop1s=,
                  method=, median_filter=)
/* DOCUMENT c= dvv_Contrast(params, nocorr=,  shot=, reg=, highpass=,
                   prefilter=, filtered=, wndw=, invert=, drop1s=,
                   method=, median_filter=)

  Returns an IMG_DAT data structure containing the contrast function
  of the interferogram set.

  ARGUMENTS:
    frdat - fringe data in a frame set
    params - lissajous parameters for phase extraction

  KEYWORDS:
    nocorr=
    reg=
    method=
           
  SEE ALSO:
    dvv_Phase, dvv_Ampl, dvv_Intensity
*/
{
  extern _DVV_PHASE_MAP, _DVV_AMPL_MAP;
  extern _CSPHI, _SNPHI, _CSPHI0, _SNPHI0, _CSPHI1, _SNPHI1;
  extern _DVV_CH1S_FILTER, _DVV_CH2S_FILTER;
  extern _DVV_CH1P_FILTER, _DVV_CH2P_FILTER;
  extern _DVV_CMMN_FILTER, _DVV_CMMN_SUM_FILTER;
  extern _DVV_CHP_FILTER, _DVV_CHS_FILTER;
  extern _DVV_HF_FILTER;
  extern _ms, _mp;
  extern _DC_LIMIT;

  qph= dvv_Quadrature(frdat, params, &weight, highpass= highpass, prefilter= prefilter,
                      nocorr= nocorr, filtered= filtered, wndw= wndw, reg= reg,
                      drop1s= drop1s, method= method, median_filter= median_filter);

  restore, qph, csphi, snphi, stot, ptot;

  //return img_div(img_operate(img_add(img_mul(csphi,csphi), img_mul(snphi,snphi)), sqrt),
  //               img_add(stot, ptot, average= 1));
  return img_add(img_mul(csphi,csphi), img_mul(snphi,snphi));
}

func dvv_ContrastScan(frdat, params, N=, w=, delta=, ar=, xoff=, yoff=, title=, reg=)
/* DOCUMENT new_params= dvv_ContrastScan, frdat, init_params,
   N=, w= , delta=, ar=, xoff= , yoff= ,title=, reg=

   Performs a parameter scan of the contrast function while varying
   one of four parameters. Four parameters are normally examined, and
   each is scanned while holding the other 3 fixed.  Each scan is
   fitted with a quadratic polynomial and the extremum of the
   polynomial is determined and returned as an updated value for the
   given parameter. The scan results are also plotted in a 4-frame
   plot.

   The quadrature components of the lissajous pattern are given by the
   cosine and delta-corrected sine terms:
   
   x-component = AR * (sdiff + XOFF)
   y-component = ((sdiff + XOFF) * cos(DELTA) - (pdiff + YOFF)) / sin(DELTA)

   Small adjustments are then made to the lissajous components
   according to the parameters:
     (1) delta - phase offset between s-and p-channels (approx -pi/2 = -1.5),
         cos(DELTA) ~ 0 and sin(DELTA) ~ 1
     (2) ar - normalized amplitude ratio between the p-fringes and
         s-fringes (nominal value 1.0)
     (3) xoff - s-channel offset or XOFF of the lissajous pattern
     (4) yoff - p-channel offset or YOFF of the lissajous pattern

   The merit functions examined during the scan are the fourier
   amplitudes of the contrast function at the fringe mode and the
   second harmonic of the fringe-mode.  A well-adjusted parameter set
   should minimize the spectral content at these two modes.  The phase
   and amplitude terms are adjusted to minimize the 2nd harmonic; the
   offsets are adjusted to minimize the 1st harmonic.
     
   SEE ALSO:
     dvv_Contrast
 */
{
  if(is_void(N)) N= 5;
  if(is_void(w)) w= 0;
  if(!is_void(params)) {
    if(numberof(params) != 4) error, "Invalid parameter set.";
    delta= params(1);
    ar= params(2);
    xoff= params(3);
    yoff= params(4);
    mult= 1.0;
  } else {
    if(is_void(delta)) delta= -pi/2;
    if(is_void(ar)) ar= 1.0;
    if(is_void(xoff)) xoff= 0.0;
    if(is_void(yoff)) yoff= 0.0;
    params= [delta, ar, xoff, yoff];
    mult= 5.0;
  }

  tot= dvv_frameSum(frdat, stot, ptot, average= 1);
  ssig= img_zclip(img_div(img_sub(frdat.ch2S, frdat.ch1S), tot), [-1.5, 1.5]);
  fm= dvv_getFringeMode(ssig);
  
  //mode1= swrite(format="mode_select %7.5f %7.5f %7.5f", _FRINGE_MODE(1), _FRINGE_MODE(2), 0.004);
  //mode2= swrite(format="mode_select %7.5f %7.5f %7.5f", 2*_FRINGE_MODE(1), 2*_FRINGE_MODE(2), 0.005);
  mode1= swrite(format="mode_select %7.5f %7.5f %7.5f", fm(1), fm(2), 0.004);
  mode2= swrite(format="mode_select %7.5f %7.5f %7.5f", 2*fm(1), 2*fm(2), 0.005);
  if(is_void(reg)) {
    m1= dvv_GenFilterMask(frdat.ch1S, mode1);
    m2= dvv_GenFilterMask(frdat.ch1S, mode2);
  } else {
    m1= dvv_GenFilterMask(img_extract(frdat.ch1S, reg), mode1);
    m2= dvv_GenFilterMask(img_extract(frdat.ch1S, reg), mode2);
  }

  window, w, style= "square2x2vg.gs", legends= 0; fma;
  if(!is_void(delta)) {
    d1= d2= dd= delta+mult*span(-0.045, 0.045, N) ;
    q= params;
    for(i= 1; i<= N; i++) {
      q(1)= dd(i);
      c= dvv_Contrast(frdat, q, reg= reg);
      //d1(i)= dvv_PSD_rms(dvv_PSD(c, wndw= "hanning"), mask= m1);
      d2(i)= dvv_PSD_rms(dvv_PSD(c, wndw= "hanning"), mask= m2);
    }
    plsys, 1; limits, -1.55, -1.45;
    plmk, d2, dd, msize= 0.5;
    cc= poly1_fit(d2, dd, 2);
    xx= delta + span(-0.1, 0.1, 1000);
    plg, poly1(xx, cc), xx, color= "red", marks= 0;
    delta0= interp(xx, poly1_deriv(xx, cc), 0.0);
    //delta0= max(xx(1), min(delta0, xx(0)));
    plmk, poly1(delta0, cc), delta0;
    //print, "delta = ", delta0;
    plt, swrite(format="Phase offset (rad): %6.4f", delta0),
      0.14, 0.9, justify="LH", height=12;
    plsys, 1; limits, delta0 - 0.04, delta0 + 0.04;
  }
    
  if(!is_void(ar)) {
    d1= d2= dd= ar+mult*span(-0.06, 0.06, N) ;
    print, "ar = ", ar, "dd= ", dd;
    q= params; q(1)= delta0;
    for(i= 1; i<= N; i++) {
      q(2)= dd(i);
      c= dvv_Contrast(frdat, q, reg= reg);
      //d1(i)= dvv_PSD_rms(dvv_PSD(c, wndw= "hanning"), mask= m1);
      d2(i)= dvv_PSD_rms(dvv_PSD(c, wndw= "hanning"), mask= m2);
    }
    plsys, 2; limits, 0.98, 1.08;
    plmk, d2, dd, msize= 0.5;
    cc= poly1_fit(d2, dd, 2);
    xx= ar + span(-0.1, 0.1, 1000);
    plg, poly1(xx, cc), xx, color= "red", marks= 0;
    ar0= interp(xx, poly1_deriv(xx, cc), 0.0);
    plmk, poly1(ar0, cc), ar0;
    //print, "ar = ", ar0;
    plsys, 0;
    plt, swrite(format="Amplitude ratio: %6.4f", ar0),
      0.49, 0.9, justify="LH", height=12;
    plsys, 2; limits, ar0 - 0.05, ar0 + 0.05;
  }

  if(!is_void(xoff)) {
    d1= d2= dd= xoff+mult*span(-0.03, 0.03, N) ;
    q= params; q(1)= delta0; q(2)= ar0;
    for(i= 1; i<= N; i++) {
      q(3)= dd(i);
      c= dvv_Contrast(frdat, q, reg= reg);
      d1(i)= dvv_PSD_rms(dvv_PSD(c, wndw= "hanning"), mask= m1);
      //d2(i)= dvv_PSD_rms(dvv_PSD(c, wndw= "hanning"), mask= m2);
    }
    plsys, 3; limits, -0.05, 0.05;
    plmk, d1, dd, msize= 0.5;
    cc= poly1_fit(d1, dd, 2);
    xx= xoff + span(-0.03, 0.03, 10*N);
    plg, poly1(xx, cc), xx, color= "red", marks= 0;
    xoff0= interp(xx, poly1_deriv(xx, cc), 0.0);
    plmk, poly1(xoff0, cc), xoff0;
    //print, "xoff = ", xoff0;
    plsys, 0;
    plt, swrite(format="X-offset: %7.5f", xoff0),
      0.14, 0.52, justify="LH", height=12;
    plsys, 3; limits, xoff0 - 0.05, xoff0 + 0.05;
  }

  if(!is_void(yoff)) {
    d1= d2= dd= yoff+mult*span(-0.03, 0.03, N) ;
    q= params; q(1)= delta0; q(2)= ar0; q(3)= xoff0;
    for(i= 1; i<= N; i++) {
      q(4)= dd(i);
      c= dvv_Contrast(frdat, q, reg= reg);
      d1(i)= dvv_PSD_rms(dvv_PSD(c, wndw= "hanning"), mask= m1);
      //d2(i)= dvv_PSD_rms(dvv_PSD(c, wndw= "hanning"), mask= m2);
    }
    plsys, 4; limits, -0.05, 0.05;
    plmk, d1, dd, msize= 0.5;
    cc= poly1_fit(d1, dd, 2);
    xx= yoff + span(-0.03, 0.03, 10*N);
    plg, poly1(xx, cc), xx, color= "red", marks= 0;
    yoff0= interp(xx, poly1_deriv(xx, cc), 0.0);
    plmk, poly1(yoff0, cc), yoff0;
    //print, "yoff = ", yoff0;
    plsys, 0;
    plt, swrite(format="Y-offset: %7.5f", yoff0),
      0.49, 0.52, justify="LH", height=12;
    plsys, 4; limits, yoff0 - 0.05, yoff0 + 0.05;
  }

  if(!is_void(title)) {plsys, 0; plt, "Contrast: "+title, 0.1, 0.96, justify="LH", height= 18;}
  if(!am_subroutine()) return [delta0, ar0, xoff0, yoff0]; 
}

func dvv_PhaseScan(frdat, params, N=, w=, delta=, ar=, xoff=, yoff=, title=, reg=)
/* DOCUMENT new_params= dvv_PhaseScan, frdat, params,
                        N=, w= , delta=, ar=, xoff= , yoff= ,title=, reg=

   Performs a parameter scan of the phase function amplitude at the
   first and second harmonics of the fringe function while varying one
   of four parameters. Four parameters are normally examined, and each
   is scanned while holding the other 3 fixed.  Each scan is fitted
   with a quadratic polynomial and the extremum of the polynomial is
   determined and returned as an updated value for the given
   parameter. The scan results are also plotted in a 4-frame plot.

   The quadrature components of the lissajous pattern are given by the
   cosine and delta-corrected sine terms:
   
   x-component = AR * sdiff + XOFF
   y-component = (sdiff * cos(DELTA) - pdiff) / sin(DELTA) + YOFF

   Small adjustments are then made to the lissajous components
   according to the parameters:
     (1) delta - phase offset between s-and p-channels (approx -pi/2 = -1.5)
     (2) ar - normalized amplitude ratio between the p-fringes and
         s-fringes (nominal value 1.0)
     (3) xoff - s-channel offset or XOFF of the lissajous pattern
     (4) yoff - p-channel offset or YOFF of the lissajous pattern

   The merit functions examined during the scan are the fourier
   amplitudes of the unwrapped phase at the second harmonic of the
   fringe-mode.  A well-adjusted parameter set should minimize the
   spectral content at these two modes.
     
   SEE ALSO:
     dvv_ContrastScan, dvv_Contrast
 */
{
  if(is_void(N)) N= 5;
  if(is_void(w)) w= 0;
  if(!is_void(params)) {
    if(numberof(params) != 4) error, "Invalid parameter set.";
    delta= params(1);
    ar= params(2);
    xoff= params(3);
    yoff= params(4);
    mult= 3.0;
  } else {
    if(is_void(delta)) delta= -pi/2;
    if(is_void(ar)) ar= 1.0;
    if(is_void(xoff)) xoff= 0.0;
    if(is_void(yoff)) yoff= 0.0;
    params= [delta, ar, xoff, yoff];
    mult= 5.0;
  }
  
  tot= dvv_frameSum(frdat, stot, ptot, average= 1);
  ssig= img_zclip(img_div(img_sub(frdat.ch2S, frdat.ch1S), tot), [-1.5, 1.5]);
  fm= dvv_getFringeMode(ssig);
  
  mode1= swrite(format="mode_select %7.5f %7.5f %7.5f", fm(1), fm(2), 0.004);
  mode2= swrite(format="mode_select %7.5f %7.5f %7.5f", 2*fm(1), 2*fm(2), 0.005);
  //mode1= swrite(format="mode_select %7.5f %7.5f %7.5f", _FRINGE_MODE(1), _FRINGE_MODE(2), 0.004);
  //mode2= swrite(format="mode_select %7.5f %7.5f %7.5f", 2*_FRINGE_MODE(1), 2*_FRINGE_MODE(2), 0.005);
  if(is_void(reg)) {
    m1= dvv_GenFilterMask(frdat.ch1S, mode1);
    m2= dvv_GenFilterMask(frdat.ch1S, mode2);
  } else {
    m1= dvv_GenFilterMask(img_extract(frdat.ch1S, reg), mode1);
    m2= dvv_GenFilterMask(img_extract(frdat.ch1S, reg), mode2);
  }
  
  window, w, style= "square2x2vg.gs", legends= 0; fma;
  if(!is_void(delta)) {
    d1= d2= dd= delta+mult*span(-0.025, 0.025, N) ;
    q= params;
    for(i= 1; i<= N; i++) {
      q(1)= dd(i);
      rfp= dvv_Phase(frdat, q, reg= reg);
      //rfbkg= dvv_bkgPhase(rfp, degree=2, box= 100, reg= reg);
      rfbkg= img_extract(dvv_bkgPhase(frdat), reg);
      //rfpUW= dvv_unwrap(img_sub(rfp, rfbkg));
      rfpUW= unwrap2d(img_sub(rfp, rfbkg));
      d1(i)= dvv_PSD_rms(dvv_PSD(rfpUW, wndw= "hanning"), mask= m1);
      d2(i)= dvv_PSD_rms(dvv_PSD(rfpUW, wndw= "hanning"), mask= m2);
    }
    window, w; plsys, 1; limits, -1.6, -1.4;
    plmk, d1, dd, msize= 0.5, width= 3;
    plmk, d2, dd, msize= 0.5, color= "red", width= 11;
    cc= poly1_fit(d2, dd, 2);
    xx= delta + span(-0.25, 0.25, 100*N);
    plg, poly1(xx, cc), xx, color= "red", marks= 0;
    delta0= interp(xx, poly1_deriv(xx, cc), 0.0);
    //delta0= max(xx(1), min(delta0, xx(0)));
    plmk, poly1(delta0, cc), delta0;
    print, "delta = ", delta0;
    plt, swrite(format="Phase offset (rad): %6.4f", delta0),
      0.14, 0.9, justify="LH", height=12;
    plsys, 1; limits, delta0 - 0.08, delta0 + 0.08;
  }
    
  if(!is_void(ar)) {
    d1= d2= dd= ar+mult*span(-0.025, 0.025, N) ;
    q= params; q(1)= delta0;
    for(i= 1; i<= N; i++) {
      q(2)= dd(i);
      rfp= dvv_Phase(frdat, q, reg= reg);
      //rfbkg= dvv_bkgPhase(rfp, degree=2, box= 100, reg= reg);
      rfbkg= img_extract(dvv_bkgPhase(frdat), reg);
      //rfpUW= dvv_unwrap(img_sub(rfp, rfbkg));
      rfpUW= unwrap2d(img_sub(rfp, rfbkg));
      //d1(i)= dvv_PSD_rms(dvv_PSD(rfpUW, wndw= "hanning"), mask= m1);
      d2(i)= dvv_PSD_rms(dvv_PSD(rfpUW, wndw= "hanning"), mask= m2);
    }
    window, w; plsys, 2; limits, 0.98, 1.08;
    plmk, d2, dd, msize= 0.5;
    cc= poly1_fit(d2, dd, 2);
    xx= ar + span(-0.25, 0.25, 100*N);
    plg, poly1(xx, cc), xx, color= "red", marks= 0;
    ar0= interp(xx, poly1_deriv(xx, cc), 0.0);
    plmk, poly1(ar0, cc), ar0;
    print, "ar = ", ar0;
    plsys, 0;
    plt, swrite(format="Amplitude ratio: %6.4f", ar0),
      0.49, 0.9, justify="LH", height=12;
    plsys, 2; limits, ar0 - 0.05, ar0 + 0.05;
  }

  if(!is_void(xoff)) {
    d1= d2= dd= xoff+mult*span(-0.06, 0.06, N) ;
    q= params; q(1)= delta0; q(2)= ar0;
    for(i= 1; i<= N; i++) {
      q(3)= dd(i);
      rfp= dvv_Phase(frdat, q, reg= reg);
      //rfbkg= dvv_bkgPhase(rfp, degree=2, box= 100, reg= reg);
      rfbkg= img_extract(dvv_bkgPhase(frdat), reg);
      //rfpUW= dvv_unwrap(img_sub(rfp, rfbkg));
      rfpUW= unwrap2d(img_sub(rfp, rfbkg));
      d1(i)= dvv_PSD_rms(dvv_PSD(rfpUW, wndw= "hanning"), mask= m1);
      //d2(i)= dvv_PSD_rms(dvv_PSD(rfpUW, wndw= "hanning"), mask= m2);
    }
    window, w; plsys, 3; limits, -0.05, 0.05;
    plmk, d1, dd, msize= 0.5;
    cc= poly1_fit(d1, dd, 2);
    xx= xoff + span(-0.23, 0.23, 100*N);
    plg, poly1(xx, cc), xx, color= "red", marks= 0;
    xoff0= interp(xx, poly1_deriv(xx, cc), 0.0);
    plmk, poly1(xoff0, cc), xoff0;
    print, "xoff = ", xoff0;
    plsys, 0;
    plt, swrite(format="X-offset: %7.5f", xoff0),
      0.14, 0.52, justify="LH", height=12;
    plsys, 3; limits, xoff0 - 0.05, xoff0 + 0.05;
  }

  if(!is_void(yoff)) {
    d1= d2= dd= yoff+mult*span(-0.06, 0.06, N) ;
    q= params; q(1)= delta0; q(2)= ar0; q(3)= xoff0;
    for(i= 1; i<= N; i++) {
      q(4)= dd(i);
      rfp= dvv_Phase(frdat, q, reg= reg);
      //rfbkg= dvv_bkgPhase(rfp, degree=2, box= 100, reg= reg);
      rfbkg= img_extract(dvv_bkgPhase(frdat), reg);
      //rfpUW= dvv_unwrap(img_sub(rfp, rfbkg));
      rfpUW= unwrap2d(img_sub(rfp, rfbkg));
      d1(i)= dvv_PSD_rms(dvv_PSD(rfpUW, wndw= "hanning"), mask= m1);
      //d2(i)= dvv_PSD_rms(dvv_PSD(rfpUW, wndw= "hanning"), mask= m2);
    }
    window, w; plsys, 4; limits, -0.05, 0.05;
    plmk, d1, dd, msize= 0.5;
    cc= poly1_fit(d1, dd, 2);
    xx= yoff + span(-0.23, 0.23, 100*N);
    plg, poly1(xx, cc), xx, color= "red", marks= 0;
    yoff0= interp(xx, poly1_deriv(xx, cc), 0.0);
    plmk, poly1(yoff0, cc), yoff0;
    print, "yoff = ", yoff0;
    plsys, 0;
    plt, swrite(format="Y-offset: %7.5f", yoff0),
      0.49, 0.52, justify="LH", height=12;
    plsys, 4; limits, yoff0 - 0.05, yoff0 + 0.05;
  }
  if(!is_void(title)) {plsys, 0; plt, "Phase: "+title, 0.1, 0.96, justify="LH", height= 18;}
  if(!am_subroutine()) return [delta0, ar0, xoff0, yoff0]; 
}

func dvv_Intensity(frdat, shot=, smooth=, filter=)
/* DOCUMENT dvv_Field, frdat, shot=, smooth=, filter=
     
   SEE ALSO:
 */
{
  extern _DVV_PHASE_MAP, _DVV_AMPL_MAP;
  extern _CSPHI, _SNPHI, _CSPHI1, _SNPHI1;
  extern _FRINGE_MODE, _LAMBDA;

  //s- and p-channel intensities
  stot= img_add(frdat.ch1S, frdat.ch2S, average= 1);
  ptot= img_add(frdat.ch1P, frdat.ch2P, average= 1);

  //Construct amplitude from the data
  tt= img_add(stot, ptot, average= 1);
  ttavg= (*tt.data)(*)(avg);
  *tt.data-=ttavg;
  if(!is_void(smooth)) tt= dvv_Smooth(tt, dvv_mkGaussian(tt, smooth));
  if(!is_void(filter)) tt= dvv_Filter(tt, dvv_GenFilterMask(tt, filter), wndw= "hanning");
  *tt.data+=ttavg;
  if(!is_void(shot)) tt.shotid=shot+"-intensity";
  tt.z_unit= "a.u.";
  tt.z_label= "Intensity ";
  return tt;
}

func dvv_SetPA_Maps(rf, phiexpect=, ampexpect=, dbg=, type=, method=, nrg=, degree= )
/* DOCUMENT dvv_SetPA_Maps, rf, phiexpect=, ampexpect=, dbg=, type=, method=, nrg=, degree=

   Defines a map of phase differences and amplitude ratio from a
   reference data set.  These are needed to correct for slight
   deviations of the quadrature separation between the 4 frames, where
   the S- versus P- phase difference is slightly different from pi/2,
   the amplitude ratios are slightly different from 1, and there might
   be slight offsets from zero in the S- and P-channels. The maps are
   generated over a series of tiles, typically around 100x100 pixels
   covering the field of view, and then interpolated over the entire
   image.These maps are used as correction terms in the phase
   extraction algorithm.

   The result of this calculation is to define a set of global symbols
   which contain the corrections in the form of IMG_DAT() data sets.
     
     _DVV_AMPL_MAP, _DVV_PHASE_MAP, _DVV_XOFF_MAP, _DVV_YOFF_MAP
     _DVV_AMPL_DATA, _DVV_PHASE_DATA, _DVV_XOFF_DATA, _DVV_YOFF_DATA

   KEYWORDS:
     phiexpect=  expected value of the phase offset, phi
     ampexpect=  expected value of the amplitude ratio, a
     dbg= show debugging output
     type= specify interpolation type for generating the map
     method=  1 - find amp/phase from the fringe mode peak in the FFT space
                  (dvv_MapPhaseDiff) -- default
              2 - find amp/phase from sunusoidal fits to the quadrature pair
                  (dvv_MapPhaseDiff)
              3 - find amp/phase by mulit-parameter minimization of the 1st
                  and 2nd harmonics of the mode peaks in the 2D spectrum of
                  the unwrapped phase (dvv_optPhase)
              4 - find amp/phase by mulit-parameter minimization of the 1st
                  and 2nd harmonics of the mode peaks in the 2D spectrum of
                  the contrast function (dvv_optContrast)
              5 - find amp/phase by mulit-parameter minimization of the 1st
                  and 2nd harmonics of the mode peaks in the 2D spectrum of
                  the contrast function (dvv_optLissajous)
              6 - finds phase map by using FTM on S- and P- separately
     
   SEE ALSO:
     dvv_MapPhaseDiff, dvv_optPhase, dvv_optContrast, dvv_optLissajous
 */
{
  extern _DVV_MODE, _DVV_AMPL_MAP, _DVV_PHASE_MAP, _DVV_XOFF_MAP, _DVV_YOFF_MAP;
  extern _DVV_AMPL_DATA, _DVV_PHASE_DATA, _DVV_XOFF_DATA, _DVV_YOFF_DATA;
  extern _DVV_MASK_DATA;
  extern _FRINGE_MODE;
  extern _pdiffs, _aratio, _xoffs, _yoffs, _rmask, _rmask1;

  if(is_void(phiexpect)) phiexpect= 1.55;
  if(is_void(ampexpect)) ampexpect= 1.05;
  if(is_void(method)) method= 1;

  if (method < 6) {
    if(rf.binning > 1)
      _DVV_MODE= dvv_GenFilterMask(img_rebin(_DVV_FF.ch1S, rf.binning),
                                   _(_FRINGE_MODE, 0.05), type= "mode_pass");
    else
      _DVV_MODE= dvv_GenFilterMask(_DVV_FF.ch1S, _(_FRINGE_MODE, 0.05), type= "mode_pass");

    pmap= dvv_MapPhaseDiff(rf, _DVV_MODE, pdiffs, aratio, xoffs, yoffs, rmask, nrg= nrg,
                           midp= 0, dbg= dbg, method= method, phiexpect= phiexpect);

    window, 60;
    pli, rmask;
    _pdiffs= pdiffs;
    _aratio= aratio;
    _rmask= rmask;
    _xoffs= xoffs;
    _yoffs= yoffs;

    //Further processing of masked off points
    if(!is_void(rmask)) {
      //Find median values of unmasked data points
      pmed= median(pdiffs(..,1)(where(!rmask)));
      amed= median(aratio(..,1)(where(!rmask)));
      xmed= median(xoffs(..,1)(where(!rmask)));
      ymed= median(yoffs(..,1)(where(!rmask)));

      //Apply a further masking step to eliminate points that are too
      //far away from the medians
      rmask1= rmask |
        abs(pdiffs(..,1) - pmed) > 0.15 |
        abs(aratio(..,1) - amed) > 0.15 |
        abs(xoffs(..,1) - xmed) > 0.15 |
        abs(yoffs(..,1) - ymed) > 0.15;
      _rmask1= rmask1;
      print, "Total number of mask points:", rmask1(*)(sum);
    }
  
    asurf= dvv_mkMap(rf.ch1S, aratio, type= type, outlier= 0.10, pdexpect= ampexpect, mask= rmask1);
    _DVV_AMPL_DATA= img_copy(_mkMapData);
  
    psurf= dvv_mkMap(rf.ch1S, pdiffs, type= type, outlier= 0.15, pdexpect= phiexpect, mask= rmask1);
    _DVV_PHASE_DATA= img_copy(_mkMapData);

    xsurf= dvv_mkMap(rf.ch1S, xoffs, type= type, outlier= 0.05, pdexpect= 0.0, mask= rmask1);
    _DVV_XOFF_DATA= img_copy(_mkMapData);

    ysurf= dvv_mkMap(rf.ch1S, yoffs, type= type, outlier= 0.05, pdexpect= 0.0, mask= rmask1);
    _DVV_YOFF_DATA= img_copy(_mkMapData);
    _DVV_MASK_DATA= img_copy(_mkMapData, data= rmask1);

    _DVV_AMPL_MAP= img_copy(asurf);
    _DVV_PHASE_MAP= img_copy(psurf);
    _DVV_XOFF_MAP= img_copy(xsurf);
    _DVV_YOFF_MAP= img_copy(ysurf);
  
  } else if (method == 6) {

    print, "dvv_SetPA_Maps: Executing method 6";
    FOV= 300.0;  //Fit out to +/- 300 micron field
    if (is_void(degree)) degree= 1;   //Use first order polynomial surface fit
    fwhm= 100;   //Convolution box size

    //Phase
    tot= dvv_frameSum(rf, stot, ptot, average= 1);
    sfr= img_zclip(img_div(img_sub(rf.ch2S, rf.ch1S), stot), [-1,1]*1.5);
    pfr= img_zclip(img_div(img_sub(rf.ch2P, rf.ch1P), ptot), [-1,1]*1.5);
    mode_filter= swrite(format="mode_select %f %f %f", _FRINGE_MODE(1), _FRINGE_MODE(2), 0.01);
    sfrF= img_zclip(dvv_Filter(sfr, mode_filter, wndw= "hanning"), [-1,1]*1.5);
    pfrF= img_zclip(dvv_Filter(pfr, mode_filter, wndw= "hanning"), [-1,1]*1.5);
    box= img_floatData(img_mk2Dtophat(tot, fwhm));

    //Phase difference
    sph= img_fringePhase(sfrF, or= "h", wndw= "hanning");
    pph= img_fringePhase(pfrF, or= "h", wndw= "hanning");
    phi= unwrap2d(img_wrap(img_sub(sph, pph)));
    cc= dvv_Smooth(img_zclip(img_add(img_mul(sfrF, sfrF), img_mul(pfrF, pfrF)), [0, 1]), 25.0);
    cc= img_mask(cc, [-1,1,-1,1]*FOV);
    img_fitsurf, phi, af, weight= img_mul(tot,cc), reg= [-1,1,-1,1]*FOV, degree= degree, bin= 2;
    _DVV_PHASE_DATA= img_copy(phi);
    _DVV_PHASE_MAP= img_mksurf(phi, af, type= "polysurf");

    //Offsets
    _DVV_XOFF_DATA= img_convolve(sfr, box);
    _DVV_YOFF_DATA= img_convolve(pfr, box);
    img_fitsurf, _DVV_XOFF_DATA, ac, weight= stot, reg= [-1,1,-1,1]*FOV, degree= degree, bin= 2;
    img_fitsurf, _DVV_YOFF_DATA, pc, weight= ptot, reg= [-1,1,-1,1]*FOV, degree= degree, bin= 2;
    _DVV_XOFF_MAP= img_mksurf(sfr, ac, type= "polysurf");
    _DVV_YOFF_MAP= img_mksurf(pfr, pc, type= "polysurf");

    //Amplitude ratio
    sa= img_copy(sfr, data= abs(img_data(img_sub(sfrF, _DVV_XOFF_MAP))));
    pa= img_copy(pfr, data= abs(img_data(img_sub(pfrF, _DVV_YOFF_MAP))));
    sav= img_convolve(sa, box);
    pav= img_convolve(pa, box);
    _DVV_AMPL_DATA= img_div(pav, sav);
    img_fitsurf, _DVV_AMPL_DATA, aa, weight= tot, reg= [-1,1,-1,1]*FOV, degree= degree, bin= 2;
    _DVV_AMPL_MAP= img_mksurf(sa, aa, type= "polysurf");
    
  } else {
    error, "Unknown method";
  }
}

func dvv_MapPhaseDiff(frdat, mask, &pdiffs, &aratio, &xoffs, &yoffs, &rmask, dbg=, nrg=, cutoff=, ampl=, midp=, method=, phiexpect=)
/* DOCUMENT dvv_MapPhaseDiff, frdat, mask, pdiffs, aratio,  &rmask, dbg=, nrg=, cutoff=, ampl=, midp=, method=, phiexpect=

   Maps phase and amplitude differences measured over a reference data set of
   OHRV data.

   KEYWORDS:

   dbg= turn on debugging output
     nrg= numberof of tiles to use over the 800 micron FOV
     cutoff= <not used?>
     ampl= <not used?>
     midp= <not used, obsolete?>
     method= 1, 2, 3, 4 or 5 (see dvv_SetPA_Maps)
     phiexpect= expected phi difference (radians)
     
   SEE ALSO:
     dvv_SetPA_Maps, dvv_GetPhaseDiff, dvv_Phase
 */
{
  extern _iS, _iP, _fw;
  extern _DVV_FIELD;
  field= (*frdat.ch1S.xscale)(0)/1.1125;
  //field= _DVV_FIELD/1.1125;
  
  if(is_void(_fw)) fw= 425.0; else fw= _fw;
  if(is_void(phiexpect)) phiexpect= 1.5;
  
  //stot= img_add(frdat.ch1S, frdat.ch2S, average= 1);
  //ptot= img_add(frdat.ch1P, frdat.ch2P, average= 1);
  //stot= img_add(frdat.ch1S, frdat.ch2S);
  //ptot= img_add(frdat.ch1P, frdat.ch2P);
  //tot= img_add(stot, ptot, average= 1);

  tot= dvv_frameSum(frdat, stot, ptot, average= 1);
  ssig= img_zclip(img_div(img_sub(frdat.ch2S, frdat.ch1S), tot), [-1.5, 1.5]);
  psig= img_zclip(img_div(img_sub(frdat.ch2P, frdat.ch1P), tot), [-1.5, 1.5]);
  if(!is_void(mask)) {
    _iS= dvv_Filter(ssig, mask);
    _iP= dvv_Filter(psig, mask);
  } else {
    _iS= img_copy(ssig);
    _iP= img_copy(psig);
  }
  
  pmap= img_copy(_iS);
  *pmap.data= 0.0;
  if(is_void(nrg)) nrg= 9;
  xrg= span(-field, field, nrg);
  yrg= span(-field, field, nrg);
  ar= pg= array(1.0, nrg-1, nrg-1);
  xoff= yoff= array(0.0, nrg-1, nrg-1);
  rmask= array(0, nrg-1, nrg-1);
  ooc= save();
  np= 10;
  yarr= array(0.0, np);

  print, "dvv_MapPhaseDiff::method = ", method;
  fm= dvv_getFringeMode(ssig);
  mode1= swrite(format="mode_select %7.5f %7.5f %7.5f", fm(1), fm(2), 0.007);
  mode2= swrite(format="mode_select %7.5f %7.5f %7.5f", 2*fm(1), 2*fm(2), 0.007);
  params0= [1.55, 1.0, 0.0, 0.0];
  
  for(i= 1; i< nrg; i++) {
    if (!is_void(dbg)) print, "i= ", i, "nx= ", nx, "ny= ", ny;
    for(j= 1; j< nrg; j++) {
      reg= REGION(x1= xrg(i), x2= xrg(i+1), y1= yrg(j), y2= yrg(j+1));
      //if (dbg) print, "i= ", i, "j= ", j, "region= ", reg;
      write, format= "i= %d, j= %d, region= [[%f,%f], [%f,%f]]\r", i, j, reg.x1, reg.x2, reg.y1, reg.y2;
      if (method < 3) {
        if (dbg) {print, "dvv_MapPhaseDiff: Executing method 1 or 2"; reg; method;}
        d= dvv_GetPhaseDiff(_iS, _iP, reg, zpg, arx, err,
                            dbg= dbg, cutoff=  cutoff, ampl= ampl, midp= midp, method= method);
        if (err) {
          print, "i= ", i, "j= ", j, "region= ", reg, "TRYING with unfiltered inputs ... ";
          d= dvv_GetPhaseDiff(ssig, psig, reg, zpg, arx, err,
                              dbg= dbg, cutoff=  cutoff, ampl= ampl, midp= midp, method= method);
          if (err) rmask(i,j)= 1;
        }
        pmap= img_insert(pmap, reg, d);
        pg(i,j)= zpg;
        ar(i,j)= arx;
        
      } else if (method == 3) {

        if (dbg) {print, "dvv_MapPhaseDiff: Executing method 3"; reg;}
        m1= dvv_GenFilterMask(img_extract(frdat.ch1S, reg), mode1);
        m2= dvv_GenFilterMask(img_extract(frdat.ch1S, reg), mode2);
        rbkg= img_extract(dvv_bkgPhase(frdat), reg);
        save, ooc, frdat= frdat, fm= fm, reg= reg, m1= m1, m2= m2, rbkg= rbkg, np= np;
        params= params0;
        fit= lmfit(dvv_optContrast, ooc, params, yarr, tol= 1e-5, gain= 2.0);
        //if (fit.niter >= 100) params= [1.5, 1.0, 0.0, 0.0];
        dev= (params - params0)(rms);
        if (fit.niter >= 100 || dev > 0.15) {
          params= params0;
          print, "dvv_MapPhaseDiff: Method 3, region=", reg;
          rmask(i,j)= 1;
        }
        pg(i,j)= params(1);
        ar(i,j)= params(2);
        xoff(i,j)= params(3);
        yoff(i,j)= params(4);
        
      } else if (method == 4) {

        if (dbg) {print, "dvv_MapPhaseDiff: Executing method 4"; reg;}
        m1= dvv_GenFilterMask(img_extract(frdat.ch1S, reg), mode1);
        m2= dvv_GenFilterMask(img_extract(frdat.ch1S, reg), mode2);
        rbkg= img_extract(dvv_bkgPhase(frdat), reg);
        save, ooc, frdat= frdat, fm= fm, reg= reg, m1= m1, m2= m2, rbkg= rbkg, np= np;
        params= params0;
        //fit= lmfit(dvv_optPhase, ooc, params, yarr, fit= [1,2], tol= 1e-6, gain= 2.);
        //fit= lmfit(dvv_optPhase, ooc, params, yarr, fit= [3,4], tol= 1e-6, gain= 2.);
        fit= lmfit(dvv_optPhase, ooc, params, yarr, tol= 1e-6, gain= 2.);
        dev= (params - params0)(rms);
        if (fit.niter >= 100 || dev > 0.15) {
          print, "dvv_MapPhaseDiff: Method 4, region=", reg, "dev = ", dev, "params= ", params;
          params= params0;
          rmask(i,j)= 1;
        }
        pg(i,j)= params(1);
        ar(i,j)= params(2);
        xoff(i,j)= params(3);
        yoff(i,j)= params(4);
        
      } else if (method == 5) {

        if (dbg) {print, "dvv_MapPhaseDiff: Executing method 5"; reg;}
        m1= dvv_GenFilterMask(img_extract(frdat.ch1S, reg), mode1);
        m2= dvv_GenFilterMask(img_extract(frdat.ch1S, reg), mode2);
        rbkg= img_extract(dvv_bkgPhase(frdat), reg);
        save, ooc, frdat= frdat, fm= fm, reg= reg, m1= m1, m2= m2, rbkg= rbkg, np= np;
        params= params0;
        fit= lmfit(dvv_optLissajous, ooc, params, yarr, tol= 1e-5, gain= 2.0);
        dev= (params - params0)(rms);
        if (fit.niter >= 100 || dev > 0.15) {
          params= params0;
          print, "dvv_MapPhaseDiff: Method 5, region=", reg;
          rmask(i,j)= 1;
        }
        //if (fit.niter >= 100) params= [1.5, 1.0, 0.0, 0.0]; 
        pg(i,j)= params(1);
        ar(i,j)= params(2);
        xoff(i,j)= params(3);
        yoff(i,j)= params(4);
      } else error, "Invalid method";
    }
  }
  xg= xrg(zcen)(,-:1:nrg-1);
  yg= yrg(zcen)(-:1:nrg-1,);
  pdiffs= [_wrap(pg), yg, xg];
  aratio= [ar, yg, xg];
  xoffs= [xoff, yg, xg];
  yoffs= [yoff, yg, xg];
  return pmap;
}

func dvv_GetPhaseDiff(iS, iP, reg, &zd0, &ar0, &err, or=, dbg=, cutoff=, ampl=, dpass=, midp=, method=)
/* DOCUMENT dvv_GetPhaseDiff, iS, iP, reg, zd0, ar0, err, or=, dbg=, cutoff=, ampl=, dpass=, midp=, method=

   Implements methods 1 and 2 of dvv_MapPhaseDiff.  Method 1 extracts
   the phases directly from the fft of the fringe mode.  Method 2
   extracts the phase by fitting sinusoids to the data and computing
   the phase difference between the quadrature pairs averaged over the
   tile segment.


  INPUTS:
    iS
    iP
    reg

  OUTPUTS:
    zd0
    ar0
    err

  KEYWORDS:
    method= 1 or 2 (default)
    or=   <fringe orientation ["h" or "v"] for 1D phase unwrapping, obsolete>
    dbg=  set to enable debugging output
    cutoff= <not used?>
    ampl= <not used?>
    dpass= <for 1D phase unwrapping, obsolete>
    midp= <for 1D phase unwrapping, obsolete>
     
  SEE ALSO:
*/
{
  extern _na, _nb;
  extern _xgrad, _ygrad;

  if (is_void(method)) method= 2;
  
  if (method == 1) {
    
    iSx= img_extract(iS, reg);
    iPx= img_extract(iP, reg);
    fms= dvv_getFringeMode(iSx);
    fmp= dvv_getFringeMode(iPx);
    sphs= img_cplxPhase(img_fft2D(iSx));
    pphs= img_cplxPhase(img_fft2D(iPx));
    samp= img_cplxAmplitude(img_fft2D(iSx));
    pamp= img_cplxAmplitude(img_fft2D(iPx));
    sp= dvv_Mode(sphs, fms(1), fms(2));
    pp= dvv_Mode(pphs, fmp(1), fmp(2));
    sa= dvv_Mode(samp, fms(1), fms(2));
    pa= dvv_Mode(pamp, fmp(1), fmp(2));
    ar0= (pa/sa)(1);
    zd0= (sp - pp)(1);
    err= 0;
    return array(zd0, dimsof(img_data(iSx)));

  } else if (method == 2) {
    
    if(is_void(or)) or= "v";
    if(is_void(cutoff)) cutoff= 25.5;
    if(is_void(_nb)) _nb= 1;
    if(is_void(_na)) _na= 2;
    if(is_void(ampl)) ampl= 1.0;
    if(is_void(dpass)) dpass= 1;
    if(!is_void(midp) && midp == 0) xmidp= []; else xmidp= 1;

    err= 0;
    if(!is_void(dbg)) {
      winkill, 0;
      winkill, 1;
      winkill, 2;
      winkill, 3;
      winkill, 4;
      winkill, 5;
      winkill, 6;
      winkill, 8;
      winkill, 9;
      winkill, 10;
      winkill, 11;
    }

    nb2= _nb*_nb;
    na2= _na*_na;

    iSx= img_extract(iS, reg);
    iPx= img_extract(iP, reg);
    x0= (*iSx.xscale)(avg);
    y0= (*iSx.yscale)(avg);
  
    xg= (*iSx.xscale)(,-:1:iSx.ny) - x0;
    yg= (*iSx.yscale)(-:1:iSx.nx,) - y0;
  
    if(!is_void(dbg)) print, "iSx.nx= ", iSx.nx, "iSx.ny= ", iSx.ny;

    if (or == "h") or1= "v"; else or1= "h";
    iSx1= img_resample(iSx, fft_good(iSx.nx), fft_good(iSx.ny));
    iPx1= img_resample(iPx, fft_good(iPx.nx), fft_good(iPx.ny));
    xg1= (*iSx1.xscale)(,-:1:iSx1.ny) - x0;
    yg1= (*iSx1.yscale)(-:1:iSx1.nx,) - y0;
    filt= FILT(type= "phase", f1= 1./iSx.nx, f2= 1.0);
    if(!is_void(dbg)) filt;
    //  iSxp= img_unwrap(img_phase(iSx1, filt, or= or),
    //                   or= or1, dpass= dpass, force= 1, midp= xmidp);
    //  iPxp= img_unwrap(img_phase(iPx1, filt, or= or),
    //                   or= or1, dpass= dpass, force= 1, midp= xmidp);
    iSxp= unwrap2d(img_phase(iSx1, filt, or= or));
    iPxp= unwrap2d(img_phase(iPx1, filt, or= or));

    if(!is_void(dbg)) {
      window, 8; fma; sh, 8, img_phase(iSx1, filt, or= or);
      window, 9; fma; sh, 9, img_phase(iPx1, filt, or= or);
      window, 10; fma; sh, 10, iSxp;
      window, 11; fma; sh, 11, iPxp;
      pause, 1000;
    }

    // Generate the initial guess from the central 10% of the unwrapped
    // data
    inx1= iSx1.nx;
    iny1= iSx1.ny;
    xReg= iREGION(ix1= inx1/2 - inx1/10, ix2= inx1/2 + inx1/10,
                  iy1= iny1/2 - iny1/10, iy2= iny1/2 + iny1/10);

    iSxpx= img_extract(iSxp, xReg);
    iPxpx= img_extract(iPxp, xReg);
  
    if(!is_void(dbg)) {
      window, 12; fma; sh, 12, iSxpx;
      window, 13; fma; sh, 13, iPxpx;
    }
    x00= (*iSxpx.xscale)(avg);
    y00= (*iSxpx.yscale)(avg);
    xg1x= (*iSxpx.xscale)(,-:1:iSxpx.ny) - x00;
    yg1x= (*iSxpx.yscale)(-:1:iSxpx.nx,) - y00;

    afs= fitsurf(*iSxpx.data, yg1x, xg1x, degree= _na-1);
    afp= fitsurf(*iPxpx.data, yg1x, xg1x, degree= _na-1);
    afs= _(array(0.0, nb2), reform(afs, [1,na2]));
    afp= _(array(0.0, nb2), reform(afp, [1,na2]));
    afs(1)= afp(1)= ampl;
    if(!is_void(dbg)) {
      print, "afs init=", afs;
      print, "afp init=", afp;
    }
    afpi= afp;
    afsi= afs;

    //Now use the initial fit as input to a direct non-linear fit
    //on the S- image
    r= lmfit(dvv_FitPhase, [xg, yg], afs, *iSx.data, deriv= 1);
    if(r.niter >= 100) {
      afs= afsi;
      err= 1;
      print, "Reverting to initial fit ...";
    }
    if(!is_void(dbg)) afs;
    //Now use the initial fit as input to a direct non-linear fit
    //on the S- image
    r= lmfit(dvv_FitPhase, [xg, yg], afp, *iPx.data, deriv= 1);
    if(r.niter >= 100) {
      afp= afpi;
      err= 1;
      print, "Reverting to initial fit ...";
    }
    if(!is_void(dbg)) afp;

    //Construct the fitted functions
    zps= polysurf(reform(afs(nb2+1:0), [2,_na,_na]), yg, xg);
    zpp= polysurf(reform(afp(nb2+1:0), [2,_na,_na]), yg, xg);
    as= polysurf(reform(afs(1:nb2), [2,_nb,_nb]), yg, xg);
    ap= polysurf(reform(afp(1:nb2), [2,_nb,_nb]), yg, xg);

    if(!is_void(dbg)) {
      iSf= img_copy(iSx);
      iPf= img_copy(iPx);
      // iSf.data= &(afs(1)*cos(zps));
      // iPf.data= &(afp(1)*cos(zpp));
      iSf.data= &(afs(1)*sin(zps));
      iPf.data= &(afp(1)*sin(zpp));
      sh, 0, img_sub(iSx, iSf);
      sh, 1, img_sub(iPx, iPf);
      sh, 2, iSx; sh, 3, iSf;
      sh, 4, iPx; sh, 5, iPf;
    }

    //Phase difference over the area
    zdiff= zpp-zps;
    //wp= where(zdiff > pi);
    //wn= where(zdiff < pi);
    //if(!is_null(wp)) zdiff(wp)= zdiff(wp) - int((zdiff(wp) + pi)/(2*pi));
    //if(!is_null(wn)) zdiff(wn)= zdiff(wn) + int((zdiff(wn) - pi)/(2*pi));
    zd0= zdiff(avg,avg);
    if(zd0 > pi) zd0= zd0 - 2*pi;
    ar0= (ap/as)(avg,avg);
    //print, "zd0= ", zd0, "ar0= ", ar0;
    if(!is_void(dbg)) {
      print, "zd0= ", zd0, "ar0= ", ar0;
      window, 6;
      pli, zdiff;
      print, "zdiff=", zdiff(avg,avg);
    }
    return zdiff; 
  }
}

func dvv_mkMap(img, pdiffs, type=, outlier=, pdexpect=, mask=)
/* DOCUMENT dvv_mkMap, img, pdiffs, type=, outlier=, pdexpect=

   Generate a 2D 
     
   SEE ALSO:
 */
{
  extern _mkMapData;
    
  if (is_void(type)) type= "polysurf";

  if (!is_void(mask)) {
    wm= where(mask);
    wg= where(!mask);
    af= fitsurf(pdiffs(..,1)(wg), pdiffs(..,2)(wg), pdiffs(..,3)(wg), degree= 1);
    zz= pdiffs(..,1);
    if (numberof(wm)) {
      zz(wm)= polysurf(af, pdiffs(..,2)(wm), pdiffs(..,3)(wm));
      pdiffs(..,1)= zz;
    }
    
  } else if(!is_void(outlier)) {
    pdmed= median(pdiffs(..,1)(*)); 
    if (!is_void(pdexpect)) {
      if (pdexpect != 0 && abs((abs(pdmed)-abs(pdexpect))/pdexpect) > 0.2) {
        write(format="WARNING -- OUT OF RANGE: Observed data median=%f, expected=%f\n", pdmed, pdexpect);
        pdmed= pdexpect;
      } else if (pdexpect == 0 && abs(pdmed) > 0.1) {
        write(format="WARNING -- OUT OF RANGE: Observed data median=%f, expected=%f\n", pdmed, pdexpect);
        pdmed= pdexpect;
      }
    }
    do {
      wg= where(abs(pdiffs(..,1) - pdmed) > outlier);
      wl= where(abs(pdiffs(..,1) - pdmed) <= outlier);
      outlier*= 2;
    } while (numberof(wl) < 2);
    pdavg= pdiffs(..,1)(wl)(avg);
    if(!is_null(wg)) pdiffs(wg,1)= pdavg;
  }

  _mkMapData= img_new(pdiffs(,1,3), pdiffs(1,,2), pdiffs(..,1));
  
  if(type == "polysurf") {
    af= fitsurf(pdiffs(..,1), pdiffs(..,2), pdiffs(..,3), degree= 1);
    return dvv_mkSurf(img, af, type= "polysurf");
  } else if (type == "tl2rat") {
    pdtab= tl2cub( , pdiffs(,1,3), pdiffs(1,,2), pdiffs(..,1));
    return dvv_mkSurf(img, pdtab, type= "tl2rat");
  } else if (type == "tl2cub") {
    pdtab= tl2rat( , pdiffs(,1,3), pdiffs(1,,2), pdiffs(..,1));
    return dvv_mkSurf(img, pdtab, type= "tl2cub");
  }
}

func dvv_mkSurf(img, af, type=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  s= img_copy(img);
  xg= (*img.xscale)(,-:1:img.ny);
  yg= (*img.yscale)(-:1:img.nx,);
  if (type == "polysurf") s.data= &polysurf(af, yg, xg);
  else if (type == "tl2cub") s.data= &tl2cub(af, xg, yg);
  else if (type == "tl2rat") s.data= &tl2rat(af, xg, yg);
  else error, "type= must be 'polysurf', 'tl2cub' or 'tl2rat'";
  return s;
}

func dvv_FitPhase(x, a, &g, deriv=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _na, _nb;
  xg= x(..,1);
  yg= x(..,2);
  na2= _na*_na;
  nb2= _nb*_nb;
  bf= reform(a(1:nb2),   [2,_nb,_nb]);
  af= reform(a(nb2+1:0), [2,_na,_na]);
  if(deriv) {
    p= polysurf(af,yg,xg,agrad,deriv=1);
    z= polysurf(bf,yg,xg,zgrad,deriv=1);
    //    ga= -sin(p)(..,-:1:na2)*transpose(agrad,0);
    //    gb=  cos(p)(..,-:1:nb2)*transpose(zgrad,0);
    ga= cos(p)(..,-:1:na2)*transpose(agrad,0);
    gb= sin(p)(..,-:1:nb2)*transpose(zgrad,0);
    g= _(gb, ga);
  } else {
    p= polysurf(af,yg,xg);
    z= polysurf(bf,yg,xg);
  }
  //  return z*cos(p);
  return z*sin(p);
}

func dvv_fitProfile(img, &af, degree=, select_box=, select_circle=, exclude_box=)
/* DOCUMENT dvv_fitProfile, img, af, degree=. select_circle=, exclude_box=

   KEYWORDS:
     degree=
     select_box= select data withing this box: [x1,x2,y1,y2]
     select_circle= select data within this circle: [x,y,radius]
     exclude_box=
     
   SEE ALSO:
 */
{
  if(is_void(degree)) degree= 2;
  xg= (*img.xscale)(,-:1:img.ny);
  yg= (*img.yscale)(-:1:img.nx,);
  zg= (*img.data);
  if(!is_void(select_circle)) {
    circ= select_circle;
    wi= where((xg - circ(1))^2 + (yg - circ(2))^2 < circ(3)^2);
    if(numberof(wi) < 2*degree+1) error, "Insufficient data selected for fit!";
    xg1= xg(wi); yg1= yg(wi); zg1= zg(wi);
  } else if(!is_void(select_box)) {
    box= select_box;
    if (struct_type(box) == "REGION") {
      x1= box.x1; x2= box.x2; y1= box.y1; y2= box.y2;
    } else if (numberof(box) == 4) {
      x1= box(1); x2= box(2); y1= box(3); y2= box(4);
    }
    wi= where(xg >= x1 & xg <= x2 & yg >= y1 & yg <= y2); 
    if(numberof(wi) < 2*degree+1) error, "Insufficient data selected for fit!";
    xg1= xg(wi); yg1= yg(wi); zg1= zg(wi);
  } else {
    xg1= xg; yg1= yg; zg1= zg;
  }
  
  if(!is_void(exclude_rg)) {
    rgx= exclude_region;
    wx= where(xg1 < rgx.x1 | xg1 > rgx.x2 | yg1 < rgx.y1 | yg1 > rgs.y2);
    xg1= xg(wx); yg1= yg(wx); zg1= zg(wx);
  }
  af= fitsurf(zg1, yg1, xg1, degree= degree);
  s= img_copy(img);
  s.data= &polysurf(af, yg, xg);
  return s;
}

func dvv_Vel(img, scale=, shift=)
/* DOCUMENT dvv_Vel, img, shift=

   Returns a velocity image based on the unwrapped phase image, using
   the VPF of the interferometer.  Natural scale is units of km/s

   KEYWORDS:
     scale=  scale factor: 1 for km/s, 1000 for m/s
     shift=  additive fringe shift constant
   
   SEE ALSO:
 */
{
  extern _VPF;

  if(is_void(scale)) scale= 1000.0;
  if(is_void(shift)) n= 0; else n= shift*2*pi;
  f= img_copy(img);
  f.data= &(((*f.data)/2/pi + n)*_VPF*scale);
  if (scale == 1.0) {
    f.z_unit = "km/s";
  } else if (scale == 1000.0) {
    f.z_unit= "m/s";
  }
  f.z_label = "Velocity";
  return f;
}

func dvv_ShiftShock(V, nfringes, crcl, bkbox, scale=)
/* DOCUMENT dvv_ShiftShock, V, nfringes, crcl, bkbox, scale=

   Given an input velocity map return a velocity map that is shifted
   by nfringes fringes.  The shift is applied to a circular region
   defined by the parameter circl=[[x0,y0],[r,r]].  If a background
   region is defined the average phase in the background region is
   subtracted from the velocity field before the shift is applied.
   
   SEE ALSO:
     dvv_Vel
 */
{
  extern _VPF;
  
  Vs= img_copy(V);
  if(is_void(scale)) scale= 1.0;
  if(is_void(bkbox)) Vbkg= 0.0;
  else {Vx= img_extract(V, bkbox); Vbkg= (*Vx.data)(*)(avg);}
  if(is_void(crcl)) {x0= y0= 0.0; r= 350.0;}
  else {x0= crcl(1,1); y0= crcl(2,1); r= crcl(1,2);}
  xg= (*V.xscale)(,-:1:V.ny);
  yg= (*V.yscale)(-:1:V.nx,);
  ws= where((xg - x0)^2 + (yg - y0)^2 < r*r);
  (*Vs.data)(ws)= (*V.data)(ws) - Vbkg + nfringes*_VPF*scale;
  return Vs;
}

func dvv_ImgPhase(img)
{
  x= img_copy(img);
  xg= (*x.xscale)(,-:1:x.ny);
  yg= (*x.yscale)(-:1:x.nx,);
  //md= exp(-(0.+1.i)*2*pi*(xg*_FRINGE_MODE(1) + yg*_FRINGE_MODE(2)));
  //x.data= &((*x.data)*md);
  x.data= &atan((*x.data).im, (*x.data).re);
  return x;
}

func dvv_ImgMag(img)
{
  x= img_copy(img);
  x.data= &abs(*x.data);
  return x;
}

func dvv_FieldShift(img, dz)
{
  extern _LAMBDA, _DVV_FS_FWD, _DVV_FS_REV;

  //Temporary copies
  ft= x= img_copy(img);
  dx= (*x.xscale)(dif)(avg);
  dy= (*x.yscale)(dif)(avg);

  //Take dft
  if(is_void(_DVV_FS_FWD)) _DVV_FS_FWD= fftw_plan(dimsof(*x.data), +1);
  ft.data= &roll(fftw(*x.data, _DVV_FS_FWD));
  ft.xscale= &roll(fftw_indgen(ft.nx)/double(ft.nx)/dx);
  ft.yscale= &roll(fftw_indgen(ft.ny)/double(ft.ny)/dy);
  (*ft.xscale)(1)= - abs((*ft.xscale)(1));
  (*ft.yscale)(1)= - abs((*ft.yscale)(1));

  //Construct Feit-Fleck propagator (2 pi factors are mostly canceled out)
  xg= (*ft.xscale)(,-:1:ft.ny);
  yg= (*ft.yscale)(-:1:ft.nx,);
  xg2= xg*xg;
  yg2= yg*yg;
  k= 1./_LAMBDA;
  k2= k*k;
  FFprop= exp(-(0.+1.i)*pi*2*dz*(xg2 + yg2)/(sqrt(k2 - xg2 - yg2) + k));
  ft.data= &roll(((*ft.data)*FFprop));

  //Invert dft
  if(is_void(_DVV_FS_REV)) _DVV_FS_REV= fftw_plan(dimsof(*ft.data), -1);
  x.data= &(fftw(*ft.data, _DVV_FS_REV)/ft.nx/ft.ny);
  return x;
}

func dvv_optContrast(x, a)
/* DOCUMENT dvv_optContrast, x, a
     

   Used for optimizing the lissajous parameters by examining the PSD
   of the fringe contrast map at the fundamental and 2nd harmonic of
   the background fringe frequency.

   ARGUMENTS:
     x - an oxy_object containing the following members:
           frdat - a FrameSet containing fringe data
           rbkg - a background fringe phase image
           fm - the fringe mode of the fringe pattern
           reg - a subregion in the dataset
           m1 - fringe mode mask array - fundamental
           m2 - fringe mode mask array - 2nd harmonic
      a - the Lissajous parameters: [phi, a, xoff, yoff]
     
   SEE ALSO:
    dvv_optPhase, dvv_setPA_Maps, dvv_MapPhaseDiff
 */
{
  restore, x, frdat, fm, reg, m1, m2, np;
  c= dvv_Contrast(frdat, a, reg= reg, method= [], median_filter= 1);
  psd= dvv_PSD(c, wndw= "hanning");
  d1= dvv_PSD_rms(psd, mask= m1);
  d2= dvv_PSD_rms(psd, mask= m2);
  return array(d1*d1 + d2*d2, np);
}

func dvv_optPhase(x, a)
/* DOCUMENT dvv_optPhase, x, a

   Used for optimizing the lissajous parameters by examining the PSD of
   the unwrapped phase at the fundamental and 2nd harmonic of the
   background fringe frequency.

   ARGUMENTS:
     x - an oxy_object containing the following members:
           frdat - a FrameSet containing fringe data
           rbkg - a background fringe phase image
           fm - the fringe mode of the fringe pattern
           reg - a subregion in the dataset
           m1 - fringe mode mask array - fundamental
           m2 - fringe mode mask array - 2nd harmonic
      a - the Lissajous parameters: [phi, a, xoff, yoff]
     
   SEE ALSO:
    dvv_optContrast, dvv_setPA_Maps, dvv_MapPhaseDiff
 */
{
  restore, x, frdat, rbkg, fm, reg, m1, m2, np;
  rfp= dvv_Phase(frdat, a, reg= reg, method= [], median_filter= 1);
  //r= unwrap2d(img_sub(rfp, rbkg));
  r= unwrap2d(img_sub(rfp, rbkg));
  psd= dvv_PSD(r, wndw= "hanning");
  d1= dvv_PSD_rms(psd, mask= m1);
  d2= dvv_PSD_rms(psd, mask= m2);
  return array(d1*d1 + d2*d2, np);
}

func dvv_optLissajous(x, a, &d1P, &d2P, &d1C, &d2C)
/* DOCUMENT dvv_optLissajous, x, a

   Used for optimizing the lissajous parameters by examining the PSD of
   the unwrapped phase at the fundamental and 2nd harmonic of the
   background fringe frequency.

   ARGUMENTS:
     x - an oxy_object containing the following members:
           frdat - a FrameSet containing fringe data
           rbkg - a background fringe phase image
           fm - the fringe mode of the fringe pattern
           reg - a subregion in the dataset
           m1 - fringe mode mask array - fundamental
           m2 - fringe mode mask array - 2nd harmonic
      a - the Lissajous parameters: [phi, a, xoff, yoff]
     
   SEE ALSO:
    dvv_optContrast, dvv_setPA_Maps, dvv_MapPhaseDiff
 */
{
  restore, x, frdat, rbkg, fm, reg, m1, m2, np;
  rfp= dvv_Phase(frdat, a, reg= reg);
  psdP= dvv_PSD(unwrap2d(img_sub(rfp, rbkg)), wndw= "hanning");
  psdC= dvv_PSD(dvv_Contrast(frdat, a, reg= reg, method= [], median_filter= 1), wndw= "hanning");

  d1P= dvv_PSD_rms(psdP, mask= m1);
  d2P= dvv_PSD_rms(psdP, mask= m2);
  d1C= dvv_PSD_rms(psdC, mask= m1);
  d2C= dvv_PSD_rms(psdC, mask= m2);
  
  return array(d1P*d1P + d2P*d2P + d1C*d1C + d2C*d2C, np);
}

func dvv_optPA(x, a)
/* DOCUMENT dvv_optPA, x, a

   Used for optimizing the phase and amplitude by examining the PSD of
   the unwrapped phase at the 2nd harmonic of the background fringe
   frequency.

   ARGUMENTS:
     x - an oxy_object containing the following members:
           frdat - a FrameSet containing fringe data
           rbkg - a background fringe phase image
           fm - the fringe mode of the fringe pattern
           reg - a subregion in the dataset
           m1 - fringe mode mask array - fundamental
           m2 - fringe mode mask array - 2nd harmonic
      a - the Lissajous parameters: [phi, a, xoff, yoff]
     
   SEE ALSO:
    dvv_optContrast, dvv_setPA_Maps, dvv_MapPhaseDiff
 */
{
  restore, x, frdat, rbkg, fm, reg, m1, m2, np;
  rfp= dvv_Phase(frdat, a, reg= reg, method= [], median_filter= 1);
  //r= unwrap2d(img_sub(rfp, rbkg));
  r= unwrap2d(img_sub(rfp, rbkg));
  psd= dvv_PSD(r, wndw= "hanning");
  //d1= dvv_PSD_rms(psd, mask= m1);
  d2= dvv_PSD_rms(psd, mask= m2);
  return array(d2*d2, np);
}

func dvv_optXY(x, a)
/* DOCUMENT dvv_optXY, x, a

   Used for optimizing the offset parameters by examining the PSD of
   the unwrapped phase at the background fringe frequency.

   ARGUMENTS:
     x - an oxy_object containing the following members:
           frdat - a FrameSet containing fringe data
           rbkg - a background fringe phase image
           fm - the fringe mode of the fringe pattern
           reg - a subregion in the dataset
           m1 - fringe mode mask array - fundamental
           m2 - fringe mode mask array - 2nd harmonic
      a - the Lissajous parameters: [phi, a, xoff, yoff]
     
   SEE ALSO:
    dvv_optContrast, dvv_setPA_Maps, dvv_MapPhaseDiff
 */
{
  restore, x, frdat, rbkg, fm, reg, m1, m2, np;
  rfp= dvv_Phase(frdat, a, reg= reg, method= [], median_filter= 1);
  //r= unwrap2d(img_sub(rfp, rbkg));
  r= unwrap2d(img_sub(rfp, rbkg));
  psd= dvv_PSD(r, wndw= "hanning");
  d1= dvv_PSD_rms(psd, mask= m1);
  //d2= dvv_PSD_rms(psd, mask= m2);
  return array(d1*d1, np);
}

func dvv_plPhasePoints(i)
{
  extern _c1s, _c2s, _c1p, _c2p;
  if(i > numberof(_c1s(,2))) return  0;
  av=  [_c1s(i,2), _c2s(i,2), _c1p(i,2), _c2p(i,2)](avg);
  plmk, _c2p(i,2)/av, [_c2s(i,2),_c1s(i,2)](avg)/av, marker= 5, msize= 0.3, width= 10, color= "red";
  plmk, _c1p(i,2)/av, [_c2s(i,2),_c1s(i,2)](avg)/av, marker= 5, msize= 0.3, width= 10, color= "red";
  plmk, [_c2p(i,2),_c1p(i,2)](avg)/av, _c2s(i,2)/av,  marker= 5, msize= 0.3, width= 10, color= "blue";
  plmk, [_c2p(i,2),_c1p(i,2)](avg)/av, _c1s(i,2)/av,  marker= 5, msize= 0.3, width= 10, color= "blue";
  
  return 1;
}

func dvv_plQuadrature(i)
{
  extern c1s, c2s, c1p, c2p;
  if(i > numberof(img_data(c1s))) return  0;
  plmk, [img_data(c2p)(i), img_data(c1p)(i)], [img_data(c2s)(i), img_data(c1s)(i)], width= 10, msize= 0.3, color= "blue";
  plg, [img_data(c2p)(i), img_data(c1p)(i)], [img_data(c2s)(i), img_data(c1s)(i)], width= 3, marks= 0, color= "blue";
  plmk, [img_data(c2p)(i), img_data(c1p)(i)], [img_data(c1s)(i), img_data(c2s)(i)], width= 10, msize= 0.3, color= "blue";
  plg, [img_data(c2p)(i), img_data(c1p)(i)], [img_data(c1s)(i), img_data(c2s)(i)], width= 3, marks= 0, color= "blue";

}

func dvv_jumpList(pts, img, side=)
/* DOCUMENT dvv_jumpMask, pts, img, side=

   Generates a list of pixel coordinates that can be used to apply a
   fringe jump across a cutline defined by a series of points.  The
   cutline cannot be a closed loop, but must vary monotonically along
   either the x-coordinate or the y-coordinate from the beginning
   point to the end point - e.g. a mostly vertical line which divides
   the image into +x and -x segments, or a mostly horizontal line
   which divides the image into +y and -y segments.

   The points can be input in three possible forms as
     (a) an array of POINT(x=x1,y=x1) structs
     (b) 2 x n array: [[x1,y1],[x2,y2],...]
     (c) a 1 x m array: [x1,y1,x2,y2,x3,y3,...]

   KEYWORDS:
     side= "+x" - selects all points on the +x side of the line
           "-x","+y","-y" similarly
   
   SEE ALSO:
 */
{
  xg= img_grid(img, 1);
  yg= img_grid(img, 2);
  if (struct_type(pts) == "POINT") {
    yp= pts.y; xp= pts.x;
  } else if (struct_type(pts) == "double") {
    if (rankof(pts) == 2) {
      xp= pts(1,);
      yp= pts(2,);
    } else {
      pts= reform(pts, [2,2,numberof(pts)/2]);
      xp= pts(1,);
      yp= pts(2,);
    }
  }

  wr= [];
  if (side == "+x" || side == "-x") {
    if (!is_monotonic(yp))
      error, "y coordinates must be montonic for +/- x mask";
    if (yp(dif)(avg) < 0) {
      yp= yp(::-1);
      xp= xp(::-1);
    }
    for(j= 2; j<= numberof(yp); j++) {
      if (xp(j) == xp(j-1)) {
        test = (side == "+x" ? yg >= yp(j-1) & yg < yp(j) & xg > xp(j) :
                yg >= yp(j-1) & yg < yp(j) & xg <= xp(j));
      } else {
        m= (yp(j) - yp(j-1))/(xp(j) - xp(j-1));
        b= (xp(j)*yp(j-1) - xp(j-1)*yp(j))/(xp(j) - xp(j-1));
        test= (side == "+x" ? yg >= yp(j-1) & yg < yp(j) & xg > (yg - b)/m :
               yg >= yp(j-1) & yg < yp(j) & xg <= (yg - b)/m);
      }
      wr= _(wr, where(test));
    }
    if (is_array(wr)) return wr;
  } else if (side == "+y" || side == "-y") {
    if (!is_monotonic(xp))
      error, "x coordinates must be montonic for +/- y mask";
    if (xp(dif)(avg) < 0) {
      yp= yp(::-1);
      xp= xp(::-1);
    }
    for(j= 2; j<= numberof(yp); j++) {
      if (yp(j) == yp(j-1)) {
        test = (side == "+y" ? xg >= xp(j-1) & xg < xp(j) & yg > yp(j) :
                xg >= xp(j-1) & xg < xp(j) & yg <= yp(j));
      } else {
        m= (yp(j) - yp(j-1))/(xp(j) - xp(j-1));
        b= (xp(j)*yp(j-1) - xp(j-1)*yp(j))/(xp(j) - xp(j-1));
        test= (side == "+y" ? xg >= xp(j-1) & xg < xp(j) & yg > m * xg + b :
               xg >= xp(j-1) & xg < xp(j) & yg <= m * xg + b);
      }
      wr= _(wr, where(test));
    }
    if (is_array(wr)) return wr;
  } else {
    error, "No mask applied";
  } 
}

func dvv_adjustPhase(pUW, p, bkg1, bkg2)
{
  if (is_void(bkg2))
    pX= img_wrap(img_sub(p, pUW, bkg1));
  else
    pX= img_wrap(img_sub(p, pUW, bkg1, bkg2));
  return img_add(pUW, pX);
}


func dvv_delPhi(phi, fr, delta, ar, &result, case=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  if (is_void(case)) case= 1;
  sd= img_data(img_sub(fr.ch2S, fr.ch1S));
  pd= img_data(img_sub(fr.ch2P, fr.ch1P));
  phi= img_data(phi);
  
  if (struct_type(ar) == "IMG_DAT") {
    ar= img_data(ar);
  } else if (numberof(ar) > 1) {
      xg= (*sdiff.xscale)(,-:1:sdiff.ny);
      yg= (*sdiff.yscale)(-:1:sdiff.nx,);
      ar= polysurf(amp_ratio, yg, xg);
  }


  if (struct_type(delta) == "IMG_DAT") {
    delta= img_data(delta);
  } else if (numberof(delta) > 1) {
      xg= (*sdiff.xscale)(,-:1:sdiff.ny);
      yg= (*sdiff.yscale)(-:1:sdiff.nx,);
      delta= polysurf(delta, yg, xg);
  }

  den= sqrt(pd*pd + ar*ar*sd*sd - 2*ar*pd*sd*cos(delta));
  num= ar*sd*sin(delta - phi) + pd*sin(phi);
  arg= num/den;
  wp= where(arg > 1.0); if (is_array(wp)) arg(wp)= 1.0;
  wm= where(arg < -1.0); if (is_array(wm)) arg(wm)= -1.0;
  result= save(num= num, den= den, arg= arg);
  if (case == 1) {
    return img_copy(fr.ch1S, data= -acos(-arg), title= "del-phi");
  } else if (case == 2) {
    return img_copy(fr.ch1S, data= +acos(-arg), title= "del-phi");
  } else if (case == 3) {
    return img_copy(fr.ch1S, data= -acos(+arg), title= "del-phi");
  } else if (case == 4) {
    return img_copy(fr.ch1S, data= +acos(+arg), title= "del-phi");
  } else error, "Invalid case.";
}

func dvv_PhaseAmplScan(fr, &ooc, reg=, xoff=, yoff=)
/* DOCUMENT rpa= dvv_PhaseAmpleScan(fr, ooc, reg=, xoff=, yoff=)

   Performs a 2D scan of the phase and amplitude lissajous parameters
   and examines the unwrapped phase energy in the fundamental and 2nd
   harmonics.

   Example usage:

   > rpa= dvv_PhaseAmplScan(rf, ooc, xoff= -0.005, yoff= 0.005)
   > dvv_DisplayScan, rpa, w= 20
     
   SEE ALSO:
 */
{
  if (is_void(reg)) reg= [-150,150,-150,150];
  if (is_void(xoff)) xoff= 0.0;
  if (is_void(yoff)) yoff= 0.0;
  
  tot= dvv_frameSum(fr, stot, ptot, average= 1);
  ssig= img_zclip(img_div(img_sub(fr.ch2S, fr.ch1S), tot), [-1.5, 1.5]);
  psig= img_zclip(img_div(img_sub(fr.ch2P, fr.ch1P), tot), [-1.5, 1.5]);
  fm= dvv_getFringeMode(ssig);
  mode1= swrite(format="mode_select %7.5f %7.5f %7.5f", fm(1), fm(2), 0.007);
  mode2= swrite(format="mode_select %7.5f %7.5f %7.5f", 2*fm(1), 2*fm(2), 0.007);
  m1= dvv_GenFilterMask(img_extract(fr.ch1S, reg), mode1);
  m2= dvv_GenFilterMask(img_extract(fr.ch1S, reg), mode2);
  rbkg= img_extract(dvv_bkgPhase(fr), reg);
  ooc= save(frdat= fr, fm= fm, reg= reg, m1= m1, m2= m2, rbkg= rbkg, np= 10);
  
  pscl= span(1.45, 1.60, 10);
  ascl= span(0.95, 1.1, 10);
  pgr= pscl(,-:1:10);
  agr= ascl(-:1:10,);sb, 300
  p1gr= p2gr= c1gr= c2gr= array(0.0, 10, 10);
  for (i=1; i<= 10; i++) {
    for (j= 1; j<= 10; j++) {
      dvv_optLissajous, ooc, [pgr(i,j), agr(i,j), xoff, yoff], p1, p2, c1, c2;
      p1gr(i,j)= p1; p2gr(i,j)= p2;
      c1gr(i,j)= c1; c2gr(i,j)= c2;
    }
  }
  c1im= img_new(pscl, ascl, c1gr, title= "Contrast - fundamental", xlabel= "Phase", ylabel= "Amplitude", xunit= "rad", yunit="nil");
  c2im= img_new(pscl, ascl, c2gr, title= "Contrast - 2nd harmonic", xlabel= "Phase", ylabel= "Amplitude", xunit= "rad", yunit="nil");
  p1im= img_new(pscl, ascl, p1gr, title= "Phase - fundamental", xlabel= "Phase", ylabel= "Amplitude", xunit= "rad", yunit="nil");
  p2im= img_new(pscl, ascl, p2gr, title= "Phase - 2nd harmonic", xlabel= "Phase", ylabel= "Amplitude", xunit= "rad", yunit="nil");
  return save(c1im, c2im, p1im, p2im, type= "Phase/ampl scan");
}

func dvv_OffsetScan(fr, &ooc, reg=, phs=, ampl=)
/* DOCUMENT rpa= dvv_OffsetScan(fr, ooc, reg=, phs=, ampl=)

   Performs a 2D scan of the x- and y-offset lissajous parameters and
   examines the unwrapped phase energy in the fundamental and 2nd
   harmonics.

   Example usage:

   > offs= dvv_PhaseAmplScan(rf, ooc, phs= 1.55, ampl= 1.04)
   > dvv_DisplayScan, offs, w= 20
     
   SEE ALSO:
 */
{
  if (is_void(reg)) reg= [-150,150,-150,150];
  if (is_void(phs)) phs= 1.5;
  if (is_void(ampl)) ampl= 1.02;
  
  tot= dvv_frameSum(fr, stot, ptot, average= 1);
  ssig= img_zclip(img_div(img_sub(fr.ch2S, fr.ch1S), tot), [-1.5, 1.5]);
  psig= img_zclip(img_div(img_sub(fr.ch2P, fr.ch1P), tot), [-1.5, 1.5]);
  fm= dvv_getFringeMode(ssig);
  mode1= swrite(format="mode_select %7.5f %7.5f %7.5f", fm(1), fm(2), 0.007);
  mode2= swrite(format="mode_select %7.5f %7.5f %7.5f", 2*fm(1), 2*fm(2), 0.007);
  m1= dvv_GenFilterMask(img_extract(fr.ch1S, reg), mode1);
  m2= dvv_GenFilterMask(img_extract(fr.ch1S, reg), mode2);
  rbkg= img_extract(dvv_bkgPhase(fr), reg);
  ooc= save(frdat= fr, fm= fm, reg= reg, m1= m1, m2= m2, rbkg= rbkg, np= 10);
  
  xscl= span(-0.1, 0.1, 10);
  yscl= span(-0.1, 0.1, 10);
  xgr= xscl(,-:1:10);
  ygr= yscl(-:1:10,);
  p1gr= p2gr= c1gr= c2gr= array(0.0, 10, 10);
  for (i=1; i<= 10; i++) {
    for (j= 1; j<= 10; j++) {
      dvv_optLissajous, ooc, [phs, ampl, xgr(i,j), ygr(i,j)], p1, p2, c1, c2;
      p1gr(i,j)= p1; p2gr(i,j)= p2;
      c1gr(i,j)= c1; c2gr(i,j)= c2;
    }
  }
  c1im= img_new(xscl, yscl, c1gr, title= "Contrast - fundamental", xlabel= "x-offset", ylabel= "y-offet", xunit= "nil", yunit="nil");
  c2im= img_new(xscl, yscl, c2gr, title= "Contrast - 2nd harmonic", xlabel= "x-offset", ylabel= "y-offset", xunit= "nil", yunit="nil");
  p1im= img_new(xscl, yscl, p1gr, title= "Phase - fundamental", xlabel= "x-offset", ylabel= "y-offset", xunit= "rad", yunit="nil");
  p2im= img_new(xscl, yscl, p2gr, title= "Phase - 2nd harmonic", xlabel= "x-offset", ylabel= "y-offset", xunit= "rad", yunit="nil");
  return save(c1im, c2im, p1im, p2im, type="Offset scan");
}
