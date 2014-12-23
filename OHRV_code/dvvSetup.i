/*
  DVVSETUP.I

  Definitions to setup a dVV/OHRV data set for analysis

  P.M. Celliers
  LLNL
  December, 2007

*/

/*

  dvv_Setup - driver to set up the reference data (frame mapping & dewarping)
  dvv_GenerateFrameSet - generates a four-frame set from raw image data
  dvv_DefineBkg - defines the CCD background image _DVV_BKG
  dvv_SetVISAR - sets the VPF
  dvv_setMapping - sets the mapping/dewarping coefficients
  dvv_ApplyFilter - apply a filter to a frameset
  dvv_ApplyWarp - apply the warping transformation
  dvv_polyDewarp - applies a 2D polynomial dewarping transformation to a frame
  dvv_SetWarp - determines the dewarp reference data and defines the dewarping coefficients
  dvv_setPolyWarp - determines the dewarping coefficients for 2D polynomial dewarping
  dvv_setTPSWarp - determines the dewarping coefficients for thin plate spline dewarping
  dvv_FixWarpSet - detect & repair outliers in the dewarping reference data
  dvv_processWarpSet - apply a polynomial fitting procedure to repair outliers in the dewarping data
  dvv_bicubicDewarp - dewarp with a 2D cubic spine inteprpolation 
  dvv_InputMapping - determine a rough mapping from raw pixels coordinates to image coordinates
  dvv_GetGridCoords - get raw image coordinates from user input
  dvv_AngleFromGrid - estimate rotation angle correction for raw frame
  dvv_MagFromGrid - estimate magnification from raw frame
  dvv_Smooth - smooth an image by convolving with a gaussian
  dvv_xtrFrame - extract a frame from a raw image
  dvv_SetFrameMap - define the frame mapping
  dvv_MapFrame - map the frame from the raw data set
  dvv_ReverseMap - compute a reverse mapping from grid-coordinates to raw coordinates
  dvv_mkCross - make a synthetic grid cross element
  dvv_mkGaussian - make a gaussian psf-type intensity image
  dvv_mkBox - make a box
  dvv_Correlate - compute a correlation map of two images
  dvv_CoordMin - finds local grid intersection point
  dvv_SetChFilters - select parasitic mode on each channel for defining a filter
  dvv_SetupCommonSumFilter - define common mode summation filter
  dvv_RegionList - format a region list output string for defining batch processing input
  dvv_ParseRegionList - parse a formatted region list string
  dvv_FilterArg - formats a list of filter definitions for defining batch processing input
  dvv_ListFilters - generates a formatted list of filter definitions
  dvv_SelectMode - use user input to select a mode
  dvv_GetVelBkg - get user-driven mouse input to identify background regions
  dvv_getAnalysisBoxCenter - get user-driven mouse input to identify the analysis region
  dvv_getBreakoutCenter - get user-driven mouse input to identify the center of breakout
  dvv_setOption - generate analysis option list
  dvv_getFringeMode - extract the fringe mode from an data image
  dvv_bkgPhase - 
  dvv_imgOffset - apply an image cross-correlation technique to determine local image offsets
  dvv_shiftFrame - shift an OHRV frame set by a small offset (through resmapling the image)
  dvv_frameOffsets - determine frame-to-frame offsets based on cross-correlation measurements
  dvv_reportOffsets - report frame-to-frame offsets
  dvv_diffOffsets
  dvv_frameSum - comput sums of the frames
  dvv_region - return a region vector from a point and delta input
  dvv_tpsDewarp - apply thin plate spline dewarping
  dvv_tpsGen - generate a thin plate spline dewarping map
  dvv_tpsEval - evaluate a thin plate spline dewarping function
  dvv_framesetOffsets - generate a report of frame-to-frame offset statistics
  dvv_plotOffsets - plot frame-to-frame offsets
  dvv_analyzeOffsets - 
  dvv_fitOffsets
  dvv_mapQuality - evaluate frame-to-frame registration quality
  dvv_corrScan - perform fine-grained correlation scan of two images
  dvv_corrScanSet - generate a set of correlation scans for a frame set
  
 */

func dvv_Setup(shot, base=, getpoints=, other_map=, dwmethod=,
               prefix=, ff_files=, bk_files=, ssdewarp=, phiexpect=,
               swmethod=, swsmooth=, magfactor=, pamethod=)
/* DOCUMENT dvv_Setup, shot, base=, getpoints=, other_map=,
                prefix=, ff_files=, bk_files=, ssdewarp=, phiexpect=,
                swmethod=, swsmooth=, dwmethod=, magfactor=, pamethod=

   Generates reference data for importing accurately registered and
   dewarped 4-frame data sets.  For setting up the analysis of a day's
   campaign this routine should be called with getpoints= 1 to drive
   mouse-driven user input to locate the center and [+/-300, +/-300]
   grid points that are used to identify the approximate location of
   each frame on the image (to ~2 microns precision).  This is accurate
   enough to execute more accurate dewarping/import mapping functions
   that can register each frame to sub-pixel precision.

   KEYWORDS:
     base=   look for "baseline-reference.pdb" file to identify the grid
     other_map= import the dewarping map from another shot reference file
     getpoints=  get mouse-driven user input to identify the grid
     prefix=  prefix string for the file name
     ff_files=  list of data files to assemble the flat field
     bk_files=  list of data files to assemble the background
     ssdewarp=  use single step dewarp method of dewarping
     phiexpect=  expected phase
     swmethod= 1, use correlation with ideal cross image
               2, <default> use cross-correlations between channels to derive small offsets
     swsmooth=  set warp smoothing <default is 0 or void>
     dwmethod=  dewarp method <default = 2>
     magfactor= magnification factor
     
   SEE ALSO:
     dvv_setMapping, dvv_SetWarp
 */
{
  extern dt, rf;
  extern _FRINGE_MODE, _DVV_MODE;
  extern _DVV_PSF_RF, _DVV_PSF_FF, _DVV_PSF_DT;
  extern _DVV_FF, _DVV_REF, _DVV_MAPPING;

  if (is_void(swmethod)) swmethod= 2;
  if (is_void(dwmethod)) dwmethod= "polywarp";
  if (is_void(prefix)) prefix= "";
  
  if(is_void(_DVV_REF)) {
    if (base) {
      if (other_map && other_map != "") {
        write, format="dvv_Setup:: defining reference using other shot reference file - %s\n", other_map;
        rs= dvv_RestoreReferenceData(other_map);
      } else {
        if (open("baseline-reference.pdb", "r", 1))
          dvv_RestoreReferenceData, "baseline-reference.pdb";
        else if (open("../baseline/baseline-reference.pdb", "r", 1))
          dvv_RestoreReferenceData, "../baseline/baseline-reference.pdb";
        _DVV_FF=[];
        print, "dVVprocess: Redefining the grid reference ... ";
        //grd= hdf_SDget(prefix+shot+".hdf", setname="References", plane= 2, setdouble= 1);
        //bkg= hdf_SDget(prefix+shot+".hdf", setname="Streak_array", plane= 2, setdouble= 1);
        print, "dvv_Setup: _DVV_REF = ", _DVV_REF;
        dvv_setMapping, prefix+shot+".hdf", getpoints= getpoints, swmethod= swmethod,
                             swsmooth= swsmooth, magfactor= magfactor;
      }

      // if (!getpoints) {
      //   if (open("baseline-reference.pdb", "r", 1))
      //     dvv_RestoreReferenceData, "baseline-reference.pdb";
      //   else if (open("../baseline/baseline-reference.pdb", "r", 1))
      //     dvv_RestoreReferenceData, "../baseline/baseline-reference.pdb";
      //   rs= dvv_SetReference(grd, dwmethod= dwmethod, dwsmooth= dwsmooth);
      // } else {
      //   grd= hdf_SDget(prefix+shot+".hdf", setname="References", plane= 2);
      //   print, "dvv_Setup: _DVV_REF = ", _DVV_REF;
      //   rs= dvv_SetReference(grd, getpoints= getpoints, dwmethod= dwmethod, dwsmooth= dwsmooth);
      // }

      dvv_SaveReferenceData, shot+"-reference.pdb";
      if(!is_void(bk_files)) dvv_ImportBkg, bk_files; else dvv_ImportBkg, prefix+shot+".hdf";

      _DVV_PSF_FF= 0.0;
      //if (_DVV_PSF_FF) psf= dvv_mkGaussian(rs.ch1S, _DVV_PSF_FF);
      if (!is_void(ff_files)) dvv_ImportFlatField, ff_files, dwmethod= dwmethod, ssdewarp= ssdewarp; else
        dvv_ImportFlatField, prefix+shot+".hdf", dwmethod= dwmethod, ssdewarp= ssdewarp;

      rf= dvv_ImportDataSet(prefix+shot+".hdf", set= "r", ssdewarp= ssdewarp, dwmethod= dwmethod, psf= 5.0);
      _FRINGE_MODE= dvv_getFringeMode(img_sub(rf.ch1S, rf.ch2S));
      _DVV_MODE= dvv_GenFilterMask(_DVV_FF.ch1S, _(_FRINGE_MODE, 0.025), type= "mode_pass");
      //if (pamaps) dvv_SetPA_Maps, rf, phiexpect= phiexpect;
      //dvv_SaveReferenceData, shot+"-reference.pdb";
      //dvv_SetPA_Maps, rf, phiexpect= phiexpect, method= pamethod;
      //dvv_SetPA_Maps, rf, phiexpect= phiexpect, method= 3;
      _DVV_PSF_RF= 0.0;
      _DVV_PSF_DT= 0.0;

      dvv_SaveReferenceData, shot+"-reference.pdb";
    } else dvv_RestoreReferenceData, shot+"-reference.pdb";
  }

  //   if(is_void(label)) label= "";
  //   if(is_void(rf)) rf= dvv_ImportDataSet(shot+"_ref.tif", noff= noff);
  //   //window, 1, legends= 0;
  //   dvv_DisplayImg, dvv_Ampl(rf, shot= shot+"REF-"+label), w= 1, 50, 15000, lgscl= 1;
  //   palette, "stern.gp";
  //   wprt, 1, shot+"-REF-ampl-"+label+".ps";
  //   //window, 2, legends= 0;
  //   dvv_DisplayImg, dvv_Phase(rf, shot= shot+"REF-"+label), w= 2;
  //   palette, "gray.gp";
  //   wprt, 2, shot+"-REF-phase-"+label+".ps";
  
  //   dt= dvv_ImportDataSet(shot+"_dat.tif", noff= noff);
  //   //window, 3, legends= 0;
  //   dvv_DisplayImg, dvv_Ampl(dt, shot= shot+"DAT-"+label), w= 3, 50, 15000, lgscl= 1;
  //   palette, "stern.gp";
  //   wprt, 3, shot+"-DAT-ampl-"+label+".ps";
  //   //window, 4, legends= 0;
  //   dvv_DisplayImg, w=4, dvv_Phase(dt, shot= shot+"DAT-"+label);
  //   palette, "gray.gp";
  //   wprt, 4, shot+"-DAT-phase-"+label+".ps";
}

//func dvv_GenerateFrameSet(img, rs, smooth=, filter=, dwmethod=, noresmpl=, ssdewarp=, binning=)
func dvv_GenerateFrameSet(img, map, smooth=, filter=, dwmethod=, noresmpl=, ssdewarp=, binning=, dbg=)
/* DOCUMENT fs= GenerateFrameSet(img, map, smooth=, filter=, dwmethod=, noresmpl=, ssdewarp=, binning=, dbg=)

    Generates a frame set consisting of a set of four properly
    registered frames from a raw image using the registration data
    stored in the reference set rs.

   ARGUMENTS:
     img - IMG_DAT structure containing the raw data
     map  - FrameMapping structure containing the frame registration map
     
   KEYWORDS:
     smooth=  gives a psf function to use to smooth (low pass filter) the data;
              smoothing is applyed after the data are mapped & dewarped
     filter=  applies a set of filter mask specifications and any gaussian smoothing
              convolutions to each frame data after all other transformations are
              completed
     dwmethod=  --void-- apply the basic rotation, translation and scaling
                         operations, but no non-linear warping
                1 | "polywarp"
                2 | "bicubic"
                3 | "thin plate spline"
     ssdewarp= use the newer single step dewarping algorithm instead of the older
               multi-step resampling method
     binning= binning factor
     dbg=
     
     
   SEE ALSO:
     dvv_FrameSet, dvv_setMapping, dvv_FrameMapping, IMG_DAT
*/
{
  extern _pcounter;
  extern _DVV_FIELD, _DVV_MAGFACTOR;
  if (is_void(_DVV_MAGFACTOR)) _DVV_MAGFACTOR= 1.;
  if (is_void(_DVV_FIELD)) _DVV_FIELD= 400.0 * 1.1125;
  
  if (dbg) print, "GenerateFrameSet: Generate frame set, pc=",_pcounter," ... ";
  pc= _pcounter;
  tic, _pcounter++;

  fs= dvv_FrameSet();  //new fram set
  print, "dvv_GenerateFrameSet: copying map ...";
  fs.map= dvv_SetMappingCopy(map);  //copy over the mapping data
  //fs.map= rs.map;
  //fs.raw= img;
  //mapping= rs.map;
  mapping= map;

  if(!is_null(binning)) fs.binning= binning;
  if(is_void(ssdewarp)) ssdewarp= 1;  //Default is to use single step dewarp

  ch1S= dvv_MapFrame(img, mapping.ch1S);
  ch1P= dvv_MapFrame(img, mapping.ch1P);
  ch2S= dvv_MapFrame(img, mapping.ch2S);
  ch2P= dvv_MapFrame(img, mapping.ch2P);

  /*
    ch1S= dvv_mapFrame(img, *mapping.ch1S.gridList,
    mapping.ch1S.sclx, mapping.ch1S.scly);
    ch1P= dvv_mapFrame(img, *mapping.ch1P.gridList,
    mapping.ch1P.sclx, mapping.ch1P.scly);
    ch2S= dvv_mapFrame(img, *mapping.ch2S.gridList,
    mapping.ch2S.sclx, mapping.ch2S.scly);
    ch2P= dvv_mapFrame(img, *mapping.ch2P.gridList,
    mapping.ch2P.sclx, mapping.ch2P.scly);
  */

  //if (is_void(noresmpl)) {
  if (!noresmpl) {
    xscl= yscl= span(-_DVV_FIELD, _DVV_FIELD, 2100);
    fs.ch1S= img_resample(ch1S, xscl, yscl);
    fs.ch1P= img_resample(ch1P, xscl, yscl);
    fs.ch2S= img_resample(ch2S, xscl, yscl);
    fs.ch2P= img_resample(ch2P, xscl, yscl);
  } else {
    reg= REGION(x1= -_DVV_FIELD, x2= _DVV_FIELD,
                y1= -_DVV_FIELD, y2= _DVV_FIELD);
    fs.ch1S= img_extract(ch1S, reg);
    fs.ch1P= img_extract(ch1P, reg);
    fs.ch2S= img_extract(ch2S, reg);
    fs.ch2P= img_extract(ch2P, reg);
  }
    
  if (!ssdewarp) {
    
    //Default is to apply the dewarping algorithm
    if (dwmethod) {
      fs.dwmethod= dvv_dwmethod(dwmethod);
      dvv_ApplyWarp, fs, dbg= dbg, dwmethod= dwmethod;
    } else {
      fs.dwmethod= "--none--";
      print, "No dewarping applied.";
    }
    //if (is_void(nowarp)) dvv_ApplyWarp, fs, dbg= dbg;
    fs.ssdewarp= 0;

  } else {

    fs.ssdewarp= ssdewarp;
    // xscl= yscl= span(-_DVV_FIELD, _DVV_FIELD, 2100);
    
    // frm= IMG_DAT();
    // frm.xscale= &xscl;
    // frm.yscale= &yscl;
    // frm.ny= frm.nx= 2100;
    // frm.y_label= frm.x_label= "Position";
    // frm.y_unit= frm.x_unit= "!mm";
    // frm.z_label= "Intensity";
    // frm.z_unit= "ADU";
    // //frm.shotid= fs.raw.shotid;
    // frm.shotid= img.shotid;
    // frm.data= &array(0.0, frm.nx, frm.ny);
    // if (fs.binning) frm= img_rebin(frm, fs.binning);
    
    // fs.ch1S= img_copy(frm);
    // fs.ch2S= img_copy(frm);
    // fs.ch1P= img_copy(frm);
    // fs.ch2P= img_copy(frm);
    
    //    if(fs.binning) dvv_ApplyWarp, fs, raw= img_rebin(fs.raw, fs.binning);
    //else dvv_ApplyWarp, fs, raw= fs.raw;
    //Don't dewarp if dwmethod= 0 or void
    if (!is_void(dwmethod) && dwmethod != "") {
      // if (fs.binning) {
      //   dvv_ApplyWarp, fs, raw=img_rebin(img, fs.binning), dbg= dbg, dwmethod= dwmethod;
      // } else {
        dvv_ApplyWarp, fs, raw= img, dbg= dbg, dwmethod= dwmethod;
        fs.dwmethod= dvv_dwmethod(dwmethod);
       //}
    } else {
      fs.dwmethod= "--none--";
      print, "No dewarping applied.";
    }
  }
  
  fs.ch1S.shotid+= " - ch1S";
  fs.ch2S.shotid+= " - ch2S";
  fs.ch1P.shotid+= " - ch1P";
  fs.ch2P.shotid+= " - ch2P";

  if (_DVV_MAGFACTOR != 1.0) {
    (*fs.ch1S.xscale)*= _DVV_MAGFACTOR;
    (*fs.ch1S.yscale)*= _DVV_MAGFACTOR;
    (*fs.ch2S.xscale)*= _DVV_MAGFACTOR;
    (*fs.ch2S.yscale)*= _DVV_MAGFACTOR;
    (*fs.ch1P.xscale)*= _DVV_MAGFACTOR;
    (*fs.ch1P.yscale)*= _DVV_MAGFACTOR;
    (*fs.ch2P.xscale)*= _DVV_MAGFACTOR;
    (*fs.ch2P.yscale)*= _DVV_MAGFACTOR;
    //_DVV_FIELD= (*fs.ch1S.xscale)(0);
    fs.magfactor= _DVV_MAGFACTOR;
  } else {
    fs.magfactor= 1;
  }
          
  //Apply binning if specified - note that this is applied after the
  //frame mapping and dewarping
  if(fs.binning > 1) {
    fs.ch1S= img_rebin(fs.ch1S, fs.binning);
    fs.ch2S= img_rebin(fs.ch2S, fs.binning);
    fs.ch1P= img_rebin(fs.ch1P, fs.binning);
    fs.ch2P= img_rebin(fs.ch2P, fs.binning);
  }

  //Smoothing with a guassian psf
  if(!is_void(smooth)  && smooth > 0.0) {
    if (dbg) print, "Smoothing with the specified psf function";
    fs.psf= smooth;
    fs.ch1S= dvv_Smooth(fs.ch1S, smooth);
    fs.ch2S= dvv_Smooth(fs.ch2S, smooth);
    fs.ch1P= dvv_Smooth(fs.ch1P, smooth);
    fs.ch2P= dvv_Smooth(fs.ch2P, smooth);
  }

  //Additional filtering
  if(!is_void(filter)) {
    //    fs.import_filter= filter;
    fs.filter= filter;
    fset= parse_line(filter, ";");
    m= [];

    //Apply all binary masking filters
    for(i= 1; i<= numberof(fset); i++) {
      type= parse_line(fset(i),",")(1);
      if (type != "gaussian_psf") {
        if (dbg) print, "Generating mask for: ", fset(i);
        m= _(m, &dvv_GenFilterMask(fs.ch1S, fset(i)));
      }
    }

    if (numberof(m) > 0) {
      m= dvv_CombineMasks(m);
      fs.ch1S= dvv_Filter(fs.ch1S, m);
      fs.ch2S= dvv_Filter(fs.ch2S, m);
      fs.ch1P= dvv_Filter(fs.ch1P, m);
      fs.ch2P= dvv_Filter(fs.ch2P, m);
    }

    //Apply any gaussian convolution, if specified
    for(i= 1; i<= numberof(fset); i++) {
      type= parse_line(fset(i),",")(1);
      if (type == "gaussian_psf") {
        sm= 0.0;
        sread, format="%e", parse_line(fset(i),",")(2), sm;
        if (dbg) print, "Smoothing with: ", fset(i), "smoothing parameter ", sm;
        psf= dvv_mkGaussian(fs.ch1S, sm);
        fs.ch1S= dvv_Smooth(fs.ch1S, psf);
        fs.ch2S= dvv_Smooth(fs.ch2S, psf);
        fs.ch1P= dvv_Smooth(fs.ch1P, psf);
        fs.ch2P= dvv_Smooth(fs.ch2P, psf);
      }
    }
  } else fs.filter= "";

  if (dbg) print, "GenerateFrameSet: ... done",tac(pc);
  _pcounter--;
  return fs;
}

func dvv_DefineBkg(bkg)
/* DOCUMENT dvv_DefineBkg(fs)

   Defines a CCD background image to be subtracted from
   all imported data files.

   ARGUMENTS:
     bkg  - IMG_DAT structure containing the background image

   KEYWORDS:
     
     
   SEE ALSO:
     dvv_ImportFlatField, dvv_ImportData
*/
{
  extern _DVV_BKG;
  _DVV_BKG= img_copy(bkg);
  return;
}

func dvv_SetVISAR(tau, index=)
/* DOCUMENT dvv_SetVISAR, tau, index=

   Sets the VISAR constants coressponding to the given etalon
   delay.

   ARGUMENTS:
     tau - etalon delay in ps
     
   SEE ALSO:
 */
{
  extern _WAV, _TAU, _ETALON_DELTA, _VPF;

  _WAV= 0.395;
  if(!is_void(tau)) _TAU= tau/1000.0;
  _ETALON_DELTA= silica_delta(_WAV)(1);
  _VPF= _WAV/2./_TAU/(1.0 + _ETALON_DELTA);
  if(!is_void(index)) _VPF/= index;
}

func dvv_dwmethod(dwmethod)
{
  _DWMETHODS=["polywarp", "bicubic", "TPS"];
  if ((typeof(dwmethod) == "long" || typeof(dwmethod) == "int") &&
      dwmethod > 0 &&
      dwmethod < 4) return _DWMETHODS(dwmethod);
  else if (typeof(dwmethod) == "string") return dwmethod;
  else error, print("Invalid dwmethod: ", dwmethod);
}

// func dvv_SetReference(img, noresmpl=, smooth=, order=, lims=, getpoints=,
//                       swmethod=, swsmooth=, dbg=, magfactor=, case=, ng= )
// /* DOCUMENT dvv_SetReference, img, noresmpl=, smooth=, order=, lims=,
//                       getpoints=, swmethod=, swsmooth=, dbg=, magfactor=, case=, ng=
//
func dvv_setMapping(file, noresmpl=, smooth=, order=, lims=, getpoints=,
                    swmethod=, swsmooth=, dbg=, magfactor=, case=, ng=, focustest= )
/* DOCUMENT dvv_setMapping, img, noresmpl=, smooth=, order=, lims=,
                      getpoints=, swmethod=, swsmooth=, dbg=, magfactor=, case=, ng=

   Defines the mapping transformations needed for importing a data
   image into a frame set.  The mapping as well as the raw image file
   are returned in a "reference" frame set.  The mapping is also used
   to define the global structure variable _DVV_MAPPING, which is
   automatically used by dvv_ImportDataSet when a data set is loaded
   from the data file.

   KEYWORDS:
     noresmpl=
     smooth=
     order=
     lims=
     getpoints=
     swmethod= 1, use correlation with ideal cross image
               2, <default> use cross-correlations between channels to derive small offsets
     swsmooth=
     dbg=
     magfactor=   magnification factor
   
   SEE ALSO:
     dvv_SetWarp, dvv_mapFrame
 */
{
  extern _DVV_REF, _DVV_MAPPING, _DVV_MAGFACTOR, _DVV_FIELD;

  if (focustest)  {
    //grd= hdf_SDget(file, setname="Streak_array", plane= 1, setdouble= 1);
    grd= hdf_SDget(file, setname="Streak_array", setdouble= 1);
  } else {
    grd= hdf_SDget(file, setname="References", plane= 2, setdouble= 1);
  }
  //bkg= hdf_SDget(file, setname="Streak_array", plane= 2, setdouble= 1);
  //img= img_sub(grd, bkg);
  img= img_copy(grd);
  
  if (is_void(magfactor)) magfactor= 1.;
  _DVV_MAGFACTOR= magfactor;

  if (is_void(_DVV_MAPPING) && !is_void(_DVV_REF)) _DVV_MAPPING= _DVV_REF.map;
  if (is_void(_DVV_REF)) _DVV_REF= dvv_FrameSet();

  if(getpoints || is_void(_DVV_MAPPING)) {
    if (dbg) print, "dvv_setMapping: getpoints= ", getpoints, "_DVV_MAPPING= ", _DVV_MAPPING;
    mapping= dvv_InputMapping(img, lims= lims);
    print, mapping;
  } else mapping= _DVV_MAPPING;

  if (struct_compare(_DVV_REF.map, dvv_SetMappingCopy(mapping)) < 0)  _DVV_REF= dvv_FrameSet();
  _DVV_REF.map= dvv_SetMappingCopy(mapping);
  
  //if(is_void(smooth))_DVV_REF.raw= img;
  //else _DVV_REF.raw= dvv_Smooth(img, mkGaussian(img, smooth));
  
  ch1S= dvv_MapFrame(img, mapping.ch1S);
  ch1P= dvv_MapFrame(img, mapping.ch1P);
  ch2S= dvv_MapFrame(img, mapping.ch2S);
  ch2P= dvv_MapFrame(img, mapping.ch2P);
  
//   ch2S= dvv_mapFrame(img, *mapping.ch2S.gridList,
//                  mapping.ch2S.sclx, mapping.ch2S.scly);
//   ch2P= dvv_mapFrame(img, *mapping.ch2P.gridList,
//                  mapping.ch2P.sclx, mapping.ch2P.scly);

  _DVV_FIELD= 400.0 * 1.1125;
  
  //if (is_void(noresmpl)) {
  if (!noresmpl) {
    xscl= yscl= span(-_DVV_FIELD, _DVV_FIELD, 2100);
    _DVV_REF.ch1S= img_resample(ch1S, xscl, yscl);
    _DVV_REF.ch1P= img_resample(ch1P, xscl, yscl);
    _DVV_REF.ch2S= img_resample(ch2S, xscl, yscl);
    _DVV_REF.ch2P= img_resample(ch2P, xscl, yscl);
  } else {
    reg= REGION(x1= -_DVV_FIELD, x2= _DVV_FIELD,
                y1= -_DVV_FIELD, y2= _DVV_FIELD);
    _DVV_REF.ch1S= img_extract(ch1S, reg);
    _DVV_REF.ch1P= img_extract(ch1P, reg);
    _DVV_REF.ch2S= img_extract(ch2S, reg);
    _DVV_REF.ch2P= img_extract(ch2P, reg);
  }
  
  _DVV_REF.ch1S.shotid+= "-ch1S";
  _DVV_REF.ch2S.shotid+= "-ch2S";
  _DVV_REF.ch1P.shotid+= "-ch1P";
  _DVV_REF.ch2P.shotid+= "-ch2P";

  if(is_void(order)) order= 2;
  //Set the initial warp correction using 1D fix method, and ideal grid 
  dvv_SetWarp, _DVV_REF, order= 1, fix= 1, fixorder= 0, swmethod= 2, rmsf= 2.5, 
    swsmooth= swsmooth, fixall= 1, case= case, dbg= dbg, offset= 0.0, ng= ng;
  _DVV_MAPPING= _DVV_REF.map;
  h5save, shot+"_mapping_0.h5", map0= _DVV_MAPPING;

  //Iterative sequence to converge to an accurate map
  if (!focustest) {
    //First pass of correction
    f= dvv_ImportDataSet(file, set= "g", dwmethod= 1, ssdewarp= 1, binning= 1, noff= 1);
    corr_map= dvv_SetWarp(f, fix= 1, order= 1, fixorder= 1, swmethod= 2,
                          fixall= 1, offset= 0.0, rmsf= 2.5);
    dvv_correctMap, _DVV_MAPPING, corr_map, corr_dw= 1;
    h5save, shot+"_mapping_1.h5", map1= _DVV_MAPPING, corr1= corr_map;

    //Second pass of correction - 
    f2= dvv_ImportDataSet(file, set= "g", dwmethod= 1, ssdewarp= 1);
    corr_map2= dvv_SetWarp(f2, fix= 1, order= 2, fixorder= 1, swmethod= 2,
                           fixall= 1, case= "svsolve", offset= 0.0, rmsf= 2.75);
    dvv_correctMap, _DVV_MAPPING, corr_map2, corr_dw= 1;
    h5save, shot+"_mapping_2.h5", map2= _DVV_MAPPING, corr2= corr_map2;
    
    //Third pass of correction
    f3= dvv_ImportDataSet(file, set= "g", dwmethod= 3, ssdewarp= 1);
    corr_map3= dvv_SetWarp(f3, fix= 1, order= 2, fixorder= 1, swmethod= 2,
                           fixall= 1, case= "svsolve", offset= 0.0, rmsf= 3.0);
    dvv_correctMap, _DVV_MAPPING, corr_map3, corr_dw= 0;
    h5save, shot+"_mapping_3.h5", map3= _DVV_MAPPING, corr3= corr_map3;
    
    //Fourth pass of correction - final pass
    f4= dvv_ImportDataSet(file, set= "g", dwmethod= 3, ssdewarp= 1);
    corr_map4= dvv_SetWarp(f4, fix= 1, order= order, fixorder= order, swmethod= 2,
                           fixall= 1, case= "svsolve", offset= 0.0, rmsf= 3.5);
    dvv_correctMap, _DVV_MAPPING, corr_map4, corr_dw= 0;
    h5save, shot+"_mapping_4.h5", map4= _DVV_MAPPING, corr4= corr_map4;
    
    _DVV_REF.map= _DVV_MAPPING;
  }
  if (!am_subroutine()) return _DVV_MAPPING;
}

func dvv_correctMap(&map, map1, corr_dw=, order=, dbg=)
/* DOCUMENT dvv_correctMap, fs, map1, corr_dw=, order=, dbg=

   Applies a second order correction to a dewarping map based on an
   updated map determined from warp-corrected data.

   Example usage:
     ff= dvv_ImportDataSet("sXXXXX.hdf", set= "f", dwmethod= 1, ssdewarp= 1);
     corr_map= dvv_SetWarp(ff, fix= 1, fixorder= 3, swmethod= 2, case= "svsolve", offset= 0.0)
     dvv_correctMap, _DVV_MAPPING, corr_map, corr_dw= 1

   KEYWORDS:
     corr_dw=  specified the dewarp method used to produce the data set
               used to generate the correction
     
   SEE ALSO:
 */
{
  refpts= (*map1.warp_refpts);
  x= refpts.x; y= refpts.y;
  nmap= dvv_SetMappingCopy(map);

  //Determine offsets of correction map at each grid point
  d1sx= (*map1.ch1S.warp_pts).x - x;
  d1sy= (*map1.ch1S.warp_pts).y - y;
  d2sx= (*map1.ch2S.warp_pts).x - x;
  d2sy= (*map1.ch2S.warp_pts).y - y;
  d1px= (*map1.ch1P.warp_pts).x - x;
  d1py= (*map1.ch1P.warp_pts).y - y;
  d2px= (*map1.ch2P.warp_pts).x - x;
  d2py= (*map1.ch2P.warp_pts).y - y;

  if (is_void(corr_dw) || corr_dw == 1 || corr_dw == "polywarp") {
    //Apply the offsets relative to the polywarp map -- makes sense if
    //the correction map was generated against a data set that was
    //dewarped with the polywarp method (corr_dw == 1)
    x1s= x2s= x1p= x2p= x;
    y1s= y2s= y1p= y2p= y;
    _dewarp, (*map.ch1S.kx), (*map.ch1S.ky), x1s, y1s;
    _dewarp, (*map.ch2S.kx), (*map.ch2S.ky), x2s, y2s;
    _dewarp, (*map.ch1P.kx), (*map.ch1P.ky), x1p, y1p;
    _dewarp, (*map.ch2P.kx), (*map.ch2P.ky), x2p, y2p;
    (*nmap.ch1S.warp_pts).x= x1s + d1sx;
    (*nmap.ch1S.warp_pts).y= y1s + d1sy;
    (*nmap.ch2S.warp_pts).x= x2s + d2sx;
    (*nmap.ch2S.warp_pts).y= y2s + d2sy;
    (*nmap.ch1P.warp_pts).x= x1p + d1px;
    (*nmap.ch1P.warp_pts).y= y1p + d1py;
    (*nmap.ch2P.warp_pts).x= x2p + d2px;
    (*nmap.ch2P.warp_pts).y= y2p + d2py;
    
  } else {
    //Apply the offsets relative to the previous grid points; the
    //correction map was generated against a data set that was
    //dewarped with the other interpolated methods (cubic spline or
    //thin plate spline)
    (*nmap.ch1S.warp_pts).x+= d1sx;
    (*nmap.ch1S.warp_pts).y+= d1sy;
    (*nmap.ch2S.warp_pts).x+= d2sx;
    (*nmap.ch2S.warp_pts).y+= d2sy;
    (*nmap.ch1P.warp_pts).x+= d1px;
    (*nmap.ch1P.warp_pts).y+= d1py;
    (*nmap.ch2P.warp_pts).x+= d2px;
    (*nmap.ch2P.warp_pts).y+= d2py;
  }

  if (is_void(order)) order= dimsof(*nmap.ch1S.kx)(2)-1;

  //Redefine the polywarp map using the new mapping points
  dvv_setPolyWarp, nmap, order;

  //Redefine the TPS dewarping coefficients
  dvv_setTPSWarp, nmap, dbg= dbg;

  nmap.note+=" -- correction applied: "+map1.note;
  
  if (!am_subroutine()) return nmap; else map= nmap;
}

func dvv_ApplyFilter(fs, mask)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  print, "dvv_ApplyFilter - fscopy";
  newfs= dvv_FSCopy(fs);
  
  newfs.ch1S= dvv_Filter(fs.ch1S, mask);
  newfs.ch2S= dvv_Filter(fs.ch2S, mask);
  newfs.ch1P= dvv_Filter(fs.ch1P, mask);
  newfs.ch2P= dvv_Filter(fs.ch2P, mask);
  return newfs;
}

func dvv_ApplyWarp(&fs, dwmethod=, force=, dev=, raw=, dbg=)
/* DOCUMENT dvv_ApplyWarp, fs, dwmethod=, force=, dev=, raw=, dbg=

   Applies a warping transformation to a frame set, and sets
   an internal flag indicating the data has been dewarped.

   ARGUMENTS:
     fs  - FrameSet structure containing data

   KEYWORDS:
     dwmethod=  "polywarp" | 1 - uses IDL-like polywarp algorithm
                "bicubic" | 2 (default) - uses bicubic interpolation over grid segments
                "thin plate spline" | 3 - thin plate spline dewarping

     force=  Set = 1 to force dewarping transformation. Normally
             this function will aborted if a warping
             transformation has already been applied.
             
     dev=    use deviations rather than absolute coordinates for the polywarp
     raw=    Perform a single-step mapping and dewarping operation by interpolating
             from the dewarping function directly into the raw CCD image
     
   SEE ALSO:
     dvv_GenerateFrameSet
*/
{

  extern _pcounter;
  //if (dbg) print, "ApplyWarp: pc=",_pcounter," ... ";
  print, "ApplyWarp: pc=",_pcounter,"dwmethod= ", dwmethod, " ... ";
  pc= _pcounter;
  tic, _pcounter++;

  //if (dbg) print, "dvv_ApplyWarp ... "; tic;
  if(is_void(force)) force= 0;
  if(is_void(dwmethod)) dwmethod= "bicubic";
  
  if(fs.warped && !force) {
    error, "Warping transformation has already been applied.";
  } else {
    if(fs.warped) print, "WARNING (dvv_ApplyWarp) Warping transformation previously "+
                    "applied, dewarping is forced ...";
    mapping= fs.map;

    if (dwmethod == "polywarp"  || dwmethod == 1) {
      if (dbg) print, "(dvv_ApplyWarp) Applying polynomial dewarp ...";
      // dw1s= dewarp(fs.ch1S, *mapping.ch1S.kx, *mapping.ch1S.ky, dev= dev);
      // dw1p= dewarp(fs.ch1P, *mapping.ch1P.kx, *mapping.ch1P.ky, dev= dev);
      // dw2s= dewarp(fs.ch2S, *mapping.ch2S.kx, *mapping.ch2S.ky, dev= dev);
      // dw2p= dewarp(fs.ch2P, *mapping.ch2P.kx, *mapping.ch2P.ky, dev= dev);
      dw1s= dvv_polyDewarp(fs.ch1S, mapping.ch1S, raw= raw, dev= dev);
      dw1p= dvv_polyDewarp(fs.ch1P, mapping.ch1P, raw= raw, dev= dev);
      dw2s= dvv_polyDewarp(fs.ch2S, mapping.ch2S, raw= raw, dev= dev);
      dw2p= dvv_polyDewarp(fs.ch2P, mapping.ch2P, raw= raw, dev= dev);
    } else if (dwmethod == "tiled" || dwmethod == "bicubic" || dwmethod == 2) {
      if (dbg) print, "(dvv_ApplyWarp) Applying bicubic dewarp ...";
      ref= *fs.map.warp_refpts;
      dw1s= dvv_bicubicDewarp(fs.ch1S, fs.map.ch1S, raw= raw, dbg= dbg);
      dw1p= dvv_bicubicDewarp(fs.ch1P, fs.map.ch1P, raw= raw, dbg= dbg);
      dw2s= dvv_bicubicDewarp(fs.ch2S, fs.map.ch2S, raw= raw, dbg= dbg);
      dw2p= dvv_bicubicDewarp(fs.ch2P, fs.map.ch2P, raw= raw, dbg= dbg);
    } else if (dwmethod == "thin plate spline" || dwmethod == "TPS" || dwmethod == 3) {
      if (dbg) print, "(dvv_ApplyWarp) Applying thin plate spline dewarp ...";
      ref= *fs.map.warp_refpts;
      dw1s= dvv_tpsDewarp(fs.ch1S, fs.map.ch1S, fs.map, raw= raw);
      dw1p= dvv_tpsDewarp(fs.ch1P, fs.map.ch1P, fs.map, raw= raw);
      dw2s= dvv_tpsDewarp(fs.ch2S, fs.map.ch2S, fs.map, raw= raw);
      dw2p= dvv_tpsDewarp(fs.ch2P, fs.map.ch2P, fs.map, raw= raw);
    } else {
      error, "Unrecognized dewarp method.";
    }
    
    fs.ch1S= dw1s;
    fs.ch1P= dw1p;
    fs.ch2S= dw2s;
    fs.ch2P= dw2p;

    fs.warped= 1L;
    if (dbg) print, "ApplyWarp ... done: ", tac(pc);
    _pcounter--;
    return;
  }
}

func dvv_polyDewarp(img, map, raw=, dev=, dbg=)
/* DOCUMENT dvv_polyDewarp, img, map, raw=, dev=, dbg= 

   Applies the 2D polynomial dewarping model to dewarp a frame, using
   the map data provided in the second argument.

   KEYWORDS:
     dbg=
     raw=
     dev=
     
   SEE ALSO:
     dvv_bicubicDewarp, dvv_tpsDewarp, dvv_ApplyWarp
 */
{
  
  extern _pcounter;
  if (dbg) print, "polyDewarp: pc=",_pcounter," ... ";
  pc= _pcounter;
  tic, _pcounter++;

  if (raw) {
    //Apply the dewarp to x- and y- grids
    xg= img_copy(img, data= img_grid(img, 1));
    yg= img_copy(img, data= img_grid(img, 2));
    xx= img_data(dewarp(xg, *map.kx, *map.ky, dev= dev));
    yy= img_data(dewarp(yg, *map.kx, *map.ky, dev= dev));
    dvv_ReverseMap, map, yy, xx;
    dw= img_resample(img_extract(raw, dvv_xtrFrame(raw, *map.gridList)),
                     xx, yy, *img.xscale, *img.yscale, meth= 2);
    dw.shotid= img.shotid;
    dw.x_label= img.x_label;
    dw.y_label= img.y_label;
    dw.z_label= img.z_label;
    dw.x_unit= img.x_unit;
    dw.y_unit= img.y_unit;
    dw.z_unit= img.z_unit;
    dw.zrange= img.zrange;
  } else {
    dw= dewarp(img, *map.kx, *map.ky, dev= dev);
  }

  if (dbg) print, "polyDewarp ... done", tac(pc);
  _pcounter-= 1;
  return dw;  
}

func dvv_SetWarp(&fs, mask=, order=, dev=, fix=, fixorder=, fixall=, rmsf=,
                 swmethod=, dbg=, ofact=, swsmooth=, ng=, offset=, case=)
/* DOCUMENT dvv_SetWarp, fs, mask=, order=, dev=, fix=, fixorder=, fixall=, rmsf=,
                 swmethod=, dbg=, ofact=, swsmooth=, ng=, offset=, case=

   Performs a correlation scan of a series of grid reference images to
   locate the coodinates of each grid point in order to generate a
   data set that can be used for later dewarping of the image.

   The images contain patterns showing a regular grid with 100 microns
   spacings in the x- and y-directions, over total size of 800 x 800
   microns.  The image is cross-correlated with a synthetically
   generated image of a cross.

   KEYWORDS:
     mask=  set to mask off corner and edge points from the fit
     order= polynomial order used for the dewarp transformation (when dwmethod == 1)
     dev=
     dbg=
     fix=   fill in "missing" or outlier grid points using polynomial models
            1 uses a single 2D polynomial fit to detect outliers
            2 use set of 1D polynomials to detect outliers
     fixorder= poynomial order to use for "fixing" the data
     fixall= fix all grid points, don't try to detect "bad" fixes
     rmsf=  rms factor to use for rejecting outliers
     ofact= parameter that controls the offset correction magnitude
     swmethod= 1, use correlation with ideal cross image
               2, <default> use cross-correlations between channels to derived small offsets
     swsmooth=
     ng=
     offset=
     case=
     
   SEE ALSO:
     dvv_setMapping
*/
{
  extern w1s, w1p, w2s, w2p;
  extern m1s, m1p, m2s, m2p;
  extern d1s, d2s, d1p, d2p;
  extern _DVV_CORR_OFFSET;

  //Set defaults
  if (is_void(case)) case= "old_method";
  if (is_void(fix)) fix= 1;
  if (is_void(swmethod)) swmethod= 2;
  if (is_void(order)) order= 1;
  if (is_void(fixorder)) fixorder= order;

  //offset correction factor for x-correlation (method 2)
  if (is_void(ofact)) {if (case == "old_method") ofact= 4./5.; else ofact= 1.0; }

  if (is_void(swsmooth)) swsmooth= 0;     //smooth the dewarping offsets
  if (is_void(ng)) ng= 9;
  if (!is_void(offset)) {
    _DVV_CORR_OFFSET= offset*(*fs.ch1S.xscale)(dif)(avg)/2.;
  } else {
    _DVV_CORR_OFFSET= (*fs.ch1S.xscale)(dif)(avg)/2.;
  }
  
  mapping= fs.map;
  if (dbg) print, "dvv_SetWarp ....";
  if (dbg) print, "mapping [init] = ", mapping;
  
  pts= ptsf= array(POINT, ng*ng);
  pts.x= reform((span(-4, 4, ng)*100)(,-:1:ng), [1,ng*ng]);
  pts.y= reform((span(-4, 4, ng)*100)(-:1:ng,), [1,ng*ng]);
  if(!is_void(mask)) {
    print, "Applying mask ... ";
    //     wcorner= _(where(abs(pts.x) == 400 & abs(pts.y) == 400), 
    //                where(abs(pts.x) == 400 & abs(pts.y) == 200), 
    //                where(abs(pts.x) == 200 & abs(pts.y) == 400), 
    //                where(abs(pts.x) == 400 & abs(pts.y) == 300), 
    //                where(abs(pts.x) == 300 & abs(pts.y) == 400));
    wcorner= _(where(abs(pts.x) == 400 & abs(pts.y) == 400)); 
    mask= array(1L, ng*ng);
    mask(wcorner)= 0;
    pts= pts(where(mask != 0));
  }
  //window, 0;
  //plmk, pts.y, pts.x;
  npts= numberof(pts);
  mapping.warp_refpts= &pts;

  // Correlation aginst an ideal cross used to locate grid intersections
  if (swmethod == 1) {
    m1s= m2s= m1p= m2p= [];
    d1s= d2s= d1p= d2p= [];
    w1s= w2s= w1p= w2p= [];
    cx= dvv_mkCross(fs.ch1S, 3.25, 40.0);
    cx= img_sub(cx, img_avg(cx));

    // Find correlation offset from auto correlation
    // This is needed because correlation result is always
    // 0.5 to 1.5 pixels off from the nominal centering
    // so instead we find the coordinate of the correlation
    // maximum and subtract this value
    cxc= dvv_Correlate(cx, cx);
    ww= where(*cxc.data == max(*cxc.data));
    xoff= (*cxc.xscale)((ww-1)/cxc.nx+1);
    yoff= (*cxc.yscale)((ww-1)%cxc.nx+1);
    write, format= "Correlation offset is [%f, %f].\n", xoff, yoff;

    msk= img_genFilterMask(fs.ch1S, "high_pass 0.01");
    c1s= img_filter(dvv_Correlate(img_median(fs.ch1S), cx), msk);
    c1p= img_filter(dvv_Correlate(img_median(fs.ch1P), cx), msk);
    c2s= img_filter(dvv_Correlate(img_median(fs.ch2S), cx), msk);
    c2p= img_filter(dvv_Correlate(img_median(fs.ch2P), cx), msk);

    pts1S= pts1P= pts2S= pts2P= pts;
    for(i= 1; i<= npts; i++) pts1S(i)= dvv_CoordMin(c1s, pts(i), 20.0, 2.0, offset= [xoff, yoff]);
    for(i= 1; i<= npts; i++) pts1P(i)= dvv_CoordMin(c1p, pts(i), 20.0, 2.0, offset= [xoff, yoff]);
    for(i= 1; i<= npts; i++) pts2S(i)= dvv_CoordMin(c2s, pts(i), 20.0, 2.0, offset= [xoff, yoff]);
    for(i= 1; i<= npts; i++) pts2P(i)= dvv_CoordMin(c2p, pts(i), 20.0, 2.0, offset= [xoff, yoff]);

  // Cross-correlation comparing frame-to-frame at points localized at the grid nodes
  } else if (swmethod == 2) {
    
    if (dbg) print, "dvv_SetWarp: cross-correlation dewarping";
    pts1S= pts1P= pts2S= pts2P= pts;
    m1s= m2s= m1p= m2p= array(0,numberof(pts));
    w1s= w2s= w1p= w2p= array(0.0,numberof(pts));
    for(i= 1; i<= npts; i++) {
      offsets= dvv_frameOffsets(fs, pts(i), case= case);
      pts1S(i).x= pts(i).x + ofact*offsets(1,1); 
      pts1P(i).x= pts(i).x + ofact*offsets(1,2);
      pts2S(i).x= pts(i).x + ofact*offsets(1,3);
      pts2P(i).x= pts(i).x + ofact*offsets(1,4);
      pts1S(i).y= pts(i).y + ofact*offsets(2,1);
      pts1P(i).y= pts(i).y + ofact*offsets(2,2);
      pts2S(i).y= pts(i).y + ofact*offsets(2,3);
      pts2P(i).y= pts(i).y + ofact*offsets(2,4);
      m1s(i)= m1s(i) | offsets(3,1) != 0.0;
      m1p(i)= m1p(i) | offsets(3,2) != 0.0;
      m2s(i)= m2s(i) | offsets(3,3) != 0.0;
      m2p(i)= m2p(i) | offsets(3,4) != 0.0;  
      w1s(i)= offsets(4,1);
      w1p(i)= offsets(4,2);
      w2s(i)= offsets(4,3);
      w2p(i)= offsets(4,4);
      if (dbg) print, "i= ", i, pts1S(i),pts1P(i),pts2S(i),pts2P(i);
      if (dbg) print, "i= ", i, m1s(i),m1p(i),m2s(i),m2p(i);
      if (dbg) print, "i= ", i, w1s(i),w1p(i),w2s(i),w2p(i);
    } 
  } else error, swrite(format= "Invalid method: %g specified.", method);

  //Initial archive of the dewarping reference points;
  mapping.ch1S.warp_pts= &pts1S; wp1S= pts1S;
  mapping.ch1P.warp_pts= &pts1P; wp1P= pts1P;
  mapping.ch2S.warp_pts= &pts2S; wp2S= pts2S;
  mapping.ch2P.warp_pts= &pts2P; wp2P= pts2P;
  if (dbg) print, "mapping [1] = ", mapping;

  //Process the warp set to remove outliers and/or apply
  //smoothing fits
  if (fix == 1) {
    dvv_FixWarpSet, wp1S, pts, "ch1S", rp1S, order= fixorder, fixall= fixall, mask= m1s, smooth= swsmooth, rmsf= rmsf;
    if (numberof(rp1S)) mapping.ch1S.repairedList= &rp1S;
    if (allof(wp1S.x == pts1S.x) && allof(wp1S.y == pts1S.y)) print, "ch1S::Nothing fixed ...";
    
    dvv_FixWarpSet, wp1P, pts, "ch1P", rp1P, order= fixorder, fixall= fixall, mask= m1p, smooth= swsmooth, rmsf= rmsf;
    if (numberof(rp1P)) mapping.ch1P.repairedList= &rp1P;
    if (allof(wp1P.x == pts1P.x) && allof(wp1P.y == pts1P.y)) print, "ch1P::Nothing fixed ...";
    
    dvv_FixWarpSet, wp2S, pts, "ch2S", rp2S, order= fixorder, fixall= fixall, mask= m2s, smooth= swsmooth, rmsf= rmsf;
    if (numberof(rp2S)) mapping.ch2S.repairedList= &rp2S;
    if (allof(wp2S.x == pts2S.x) && allof(wp2S.y == pts2S.y)) print, "ch2S::Nothing fixed ...";
    
    dvv_FixWarpSet, wp2P, pts, "ch2P", rp2P, order= fixorder, fixall= fixall, mask= m2p, smooth= swsmooth, rmsf= rmsf;
    if (numberof(rp2P)) mapping.ch2P.repairedList= &rp2P;
    if (allof(wp2P.x == pts2P.x) && allof(wp2P.y == pts2P.y)) print, "ch2P::Nothing fixed ...";
    
  } else if (fix == 2) {
    dvv_processWarpSet, wp1S, pts, d1s, ch= "ch1S", order= fixorder, mask= m1s, weight= w1s;
    dvv_processWarpSet, wp1P, pts, d1p, ch= "ch1P", order= fixorder, mask= m1p, weight= w1p;
    dvv_processWarpSet, wp2S, pts, d2s, ch= "ch2S", order= fixorder, mask= m2s, weight= w2s;
    dvv_processWarpSet, wp2P, pts, d2p, ch= "ch2P", order= fixorder, mask= m2p, weight= w2p;
  }
  
  //Archive the dewarping reference points;
  mapping.ch1S.warp_pts= &wp1S;
  mapping.ch1P.warp_pts= &wp1P;
  mapping.ch2S.warp_pts= &wp2S;
  mapping.ch2P.warp_pts= &wp2P;
  if (dbg) print, "mapping [2] = ", mapping;

  dvv_setPolyWarp, mapping, order;
  dvv_setTPSWarp, mapping, dbg= dbg;
  
  if (dbg) print, "mapping [3] = ", mapping;
  if (swmethod == 1) {
    mapping.note=swrite(format="Initial- sw:%d,order:%d,fix:%d,fixorder:%d,swsmooth:%d,ng:%d",
                        swmethod,order,fix,fixorder,swsmooth,ng);
  } else if (swmethod == 2) {
    mapping.note=swrite(format="Initial- sw:%d,case:%s,offset:%4.2f,ofact:%4.2f,fix:%d,fixorder:%d,swsmooth:%d,ng:%d",
                        swmethod,case,offset,ofact,fix,fixorder,swsmooth,ng);
  }  else {
    mapping.note=swrite(format="Initial- sw:%d,order:%d,fix:%d,fixorder:%d,swsmooth:%d,ng:%d",
                        swmethod,order,fix,fixorder,swsmooth,ng);
  }

  if (am_subroutine()) {
    fs.map= mapping; 
    if (dbg) print, "fs.mapping [4] = ", mapping;
    return;
  } else {
    return mapping;
  }
}

func dvv_setPolyWarp(&mapping, order)
/* DOCUMENT dvv_setPolyWarp, mapping

   Sets the polywarp coefficients for a dewarp mapping.  The argument
   <order> sets the polynomial order.
     
   SEE ALSO:
     dvv_SetWarp
 */
{
  local pts, pts1S, pts1P, pts2S, pts2P;
  if (is_void(order)) order= 2;
  if (order < 1) error, "polynomial dewarping order must be >= 1";
  
  //Define the dewarping transformation from the grid reference (to be
  //used with the IDL-style polywarp function)
  pts= *mapping.warp_refpts;
  pts1S= *mapping.ch1S.warp_pts;
  pts2S= *mapping.ch2S.warp_pts;
  pts1P= *mapping.ch1P.warp_pts;
  pts2P= *mapping.ch2P.warp_pts;
  
  ww= where(pts.x < 400 & pts.x > -400 & pts.y < 400 & pts.y > -400);

  //Deterimine the polywarp coefficients, all 4 channels
  if (dev) {
    xi= pts1S.x - pts.x;
    yi= pts1S.y - pts.y;
  } else {
    xi= pts1S.x;
    yi= pts1S.y;
  }
  polywarp, xi(ww), yi(ww), pts(ww).x, pts(ww).y, order, kx1s, ky1s;

  if (dev) {
    xi= pts1P.x - pts.x;
    yi= pts1P.y - pts.y;
  } else {
    xi= pts1P.x;
    yi= pts1P.y;
  }
  polywarp, xi(ww), yi(ww), pts(ww).x, pts(ww).y, order, kx1p, ky1p;

  if (dev) {
    xi= pts2S.x - pts.x;
    yi= pts2S.y - pts.y;
  } else {
    xi= pts2S.x;
    yi= pts2S.y;
  }
  polywarp, xi(ww), yi(ww), pts(ww).x, pts(ww).y, order, kx2s, ky2s;

  if (dev) {
    xi= pts2P.x - pts.x;
    yi= pts2P.y - pts.y;
  } else {
    xi= pts2P.x;
    yi= pts2P.y;
  }
  polywarp, xi(ww), yi(ww), pts(ww).x, pts(ww).y, order, kx2p, ky2p;
  
  //Update the transformation coefficients
  mapping.ch1S.kx= &kx1s;
  mapping.ch1S.ky= &ky1s;
  mapping.ch1P.kx= &kx1p;
  mapping.ch1P.ky= &ky1p;
  mapping.ch2S.kx= &kx2s;
  mapping.ch2S.ky= &ky2s;
  mapping.ch2P.kx= &kx2p;
  mapping.ch2P.ky= &ky2p;
}

func dvv_setTPSWarp(&mapping, dbg=, dlim=)
/* DOCUMENT dvv_setTPSWarp, mapping, dbg=, dlim=

   Sets the thin plate spline coefficients for dewarping using the
   thin plate spline interpolation method.
     
   SEE ALSO:
     dvv_SetWarp, dvv_setPolyWarp
 */
{
  local pts, pts1S, pts1P, pts2S, pts2P;
  //if (dbg) print, "dvv_setTPSWarp ... ";
  print, "dvv_setTPSWarp ... ";

  if (is_void(dlim)) dlim= 425.0;
  
  //Define the dewarping transformation from the grid reference (to be
  //used with the IDL-style polywarp function)
  pts= *mapping.warp_refpts;
  pts1S= *mapping.ch1S.warp_pts;
  pts2S= *mapping.ch2S.warp_pts;
  pts1P= *mapping.ch1P.warp_pts;
  pts2P= *mapping.ch2P.warp_pts;
  
  ww= where(pts.x < dlim & pts.x > -dlim & pts.y < dlim & pts.y > -dlim);

  //Set up the reference points for the mapping of all frames
  xc= pts(ww).x;
  yc= pts(ww).y;
  mapping.xc= &xc;
  mapping.yc= &yc;

  //Deterimine the thin plate spline coefficients, each frame
  dvv_tpsGen, xc, yc, pts1S(ww).x, pts1S(ww).y, aX, aY;
  mapping.ch1S.aX= &aX;
  mapping.ch1S.aY= &aY;

  dvv_tpsGen, xc, yc, pts2S(ww).x, pts2S(ww).y, aX, aY;
  mapping.ch2S.aX= &aX;
  mapping.ch2S.aY= &aY;

  dvv_tpsGen, xc, yc, pts1P(ww).x, pts1P(ww).y, aX, aY;
  mapping.ch1P.aX= &aX;
  mapping.ch1P.aY= &aY;

  dvv_tpsGen, xc, yc, pts2P(ww).x, pts2P(ww).y, aX, aY;
  mapping.ch2P.aX= &aX;
  mapping.ch2P.aY= &aY;
}

func dvv_FixWarpSet(&pts, ref, ch, &wout, order=, rmsf=, mask=, smooth=, fixall=, dbg=)
/* DOCUMENT dvv_FixWarpSet, pts, ref, ch, order=, rmsf=, smooth=, fixall=, mask=, dbg=

     Detects and repairs outliers found in the warp set.  An outlier
     is a an offset deviates more than rmsf standard deviations away
     from the mean deviation.

     KEYWORDS:
       order=
       rmsf=   rms factor for identifying outliers, default is 2 std dev
       mask=   a bitwise mask vector for masking out points as determined
               by an external criterion
       smooth= smooth the data set by returning the points evaluated
               on the 2D polynomial fit
       fixall= repair all outliers (don't try to detect "bad" fixes)
       dbg=

     SEE ALSO:
       dvv_SetWarp, dvv_processWarpSet
     
*/
{
  if (is_void(order)) order= 1;
  if (is_void(rmsf)) rmsf= 2.0;
  opts= pts;
  dy= ref.y - pts.y;
  dx= ref.x - pts.x;
  
  //Find obvious outliers, all points should be within 10 units of nominal
  dev0= abs(dy) < 10.0 & abs(dx) < 10.0;  w0= where(dev0);
  if (dbg) print, "Raw dewarping point differences: channel", ch;
  if (dbg) print, "dy : ", dy;
  if (dbg) print, "dx : ", dx;
  if (dbg) print, "medians: dx, dy ", median(dx), median(dy);
  if (dbg) print, "rms: dx, dy ", dx(rms), dy(rms);
  //dev= (abs(dy - dy(avg)) > 2.0) | (abs(dx - dx(avg)) > 2.0);
  //dev= (abs(dy - dy(avg)) > 0.5*dy(rms)) | (abs(dx - dx(avg)) > 0.5*dx(rms));
  //dev= (abs(dy - dy(avg)) > rmsf*dy(rms)) | (abs(dx - dx(avg)) > rmsf*dx(rms));
  dev= (abs(dy - median(dy(w0))) > rmsf*dy(w0)(rms)) | (abs(dx - median(dx(w0))) > rmsf*dx(w0)(rms));
    
  if (dbg) print, "dev: ", dev;
  if (dbg) print, "=======================================================";
  if (dbg) print, "dx-rms = ", dx(w0)(rms), "dy-rms = ", dy(w0)(rms);
  
  //If a mask vector is supplied include the masked out points
  if (!is_void(mask)) {
    dev= dev | mask;
    w= where(mask);
    if (dbg) {
      if (numberof(w) > 0) {
        write, format= "Channel %s: There are %d masked off points: \n", ch, numberof(w);
        for(j= 1; j<= numberof(w); j++) write, format= " ==> [%g, %g]\n", ref(w(j)).x, ref(w(j)).y;
      }
    }
  }
  
  //Find all the outliers
  dev1= dev; NIT= 0;
  do {
    dev= dev1;
    wout= where(dev);
    wn= where(!dev);
    gpts= pts(wn);
    gpts.x= pts(wn).x - ref(wn).x;
    gpts.y= pts(wn).y - ref(wn).y;
    rpts= ref(wn);
    ax= fitsurf(gpts.x, rpts.y, rpts.x, degree= order);
    ay= fitsurf(gpts.y, rpts.y, rpts.x, degree= order);
    xfit= polysurf(ax, ref.y, ref.x);
    yfit= polysurf(ay, ref.y, ref.x);
    xdif= xfit - dx;
    ydif= yfit - dy;
    dev1= abs(xdif - median(xdif)) > rmsf*xdif(wn)(rms) | abs(ydif - median(ydif)) > rmsf*ydif(wn)(rms);
    if (!is_void(mask)) dev1= dev1 | mask;
    if (dbg) print, "numberof dx : ", numberof(dx), "numberof outliers", numberof(wout), dev1 != dev;
    NIT++;
  } while (anyof(dev1 != dev) && NIT < 10);

  if (dbg) print, "Outlier search completed, NIT=", NIT;
  
  //Good data points
  /*
  dxavg_good= (gpts.x - rpts.x)(avg);
  dxrms_good= (gpts.x - rpts.x)(rms);
  dyavg_good= (gpts.y - rpts.y)(avg);
  dyrms_good= (gpts.y - rpts.y)(rms);
  */
  dxavg_good= (gpts.x)(avg);
  dxrms_good= (gpts.x)(rms);
  dyavg_good= (gpts.y)(avg);
  dyrms_good= (gpts.y)(rms);
  
  if (numberof(wout) > 0) {

    //Generate a 2D polynomial fit to the good data
    //ax= fitsurf(gpts.x, rpts.y, rpts.x, degree= order);
    //ay= fitsurf(gpts.y, rpts.y, rpts.x, degree= order);
    
    //Interpolate and/or extrapolate the fit to fix the outliers
    //xfit= polysurf(ax, ref(wout).y, ref(wout).x);
    //yfit= polysurf(ay, ref(wout).y, ref(wout).x);
    //print, "xfit = ", xfit;
    //print, "yfit = ", yfit;

    //Any of the repairs still bad, should be within 1.5 microns (~3 pixels)
    //bad= abs(yfix - ref(wout).y - dyavg_good) > max(rmsf*dyrms_good, 1.5) |
    //  abs(xfix - ref(wout).x - dxavg_good) > max(rmsf*dxrms_good, 1.5);
    /*
    bad= abs(yfix - ref(wout).y - dyavg_good) > 1.5 |
      abs(xfix - ref(wout).x - dxavg_good) > 1.5;
    */
    bad= abs(yfit - dyavg_good) > 1.5 | abs(xfit - dxavg_good) > 1.5;

    //Print out a list of the fixes
    if (dbg) print, "Outliers, "+ch;
    if (dbg) print, "================";
    for(j= 1; j<= numberof(wout); j++) {
      if (!bad(j) || fixall) {
        if (dbg) write, format= "REPAIRED GRID POINT [%f,%f]:   [%f, %f] --> [%f, %f]\n",
          ref(wout(j)).x, ref(wout(j)).y, pts(wout(j)).x, pts(wout(j)).y,
          xfit(wout(j))+ref(wout(j)).x, yfit(wout(j))+ref(wout(j)).y;
        if (dbg) write, format= "*** yfit = %f, dy = %f, dyavg = %f dyrms = %f; xfit = %f, dx = %f, dxavg = %f, dxrms = %f\n",
          yfit(wout(j)) + ref(wout(j)).y, yfit(wout(j)), dyavg_good, dyrms_good,
          xfit(wout(j)) + ref(wout(j)).x, xfit(wout(j)), dxavg_good, dxrms_good;
      } else {
        if (dbg) write, format= "*** UNREPAIRABLE GRID POINT [%f,%f]: [%f, %f] **reverting to nominal\n",
          ref(wout(j)).x, ref(wout(j)).y, pts(wout(j)).x, pts(wout(j)).y;
        if (dbg) write, format= "*** yfit = %f, dy = %f, dyavg = %f dyrms = %f; xfit = %f, dx = %f, dxavg = %f, dxrms = %f\n",
          yfit(wout(j)) + ref(wout(j)).y, yfit(wout(j)), dyavg_good, dyrms_good,
          xfit(wout(j)) + ref(wout(j)).x, xfit(wout(j)), dxavg_good, dxrms_good;
        /*
        xfit(wout(j))= ref(wout(j)).x;
        yfit(wout(j))= ref(wout(j)).y;
        */
        xfit(wout(j))= 0.0;
        yfit(wout(j))= 0.0;
      }
    }
    
   //Apply the fixes
    print, "Fixing these points:", wout;
    pts(wout).x= xfit(wout) + ref(wout).x;
    pts(wout).y= yfit(wout) + ref(wout).y;
  }

  dy= ref.y - pts.y;
  dx= ref.x - pts.x;
  if (dbg) print, "Second pass:: Raw dewarping point differences: channel", ch;
  if (dbg) print, "dy : ", dy;
  if (dbg) print, "dx : ", dx;
  if (dbg) print, "medians: dx, dy ", median(dx), median(dy);
  if (dbg) print, "rms: dx, dy ", dx(rms), dy(rms);
  dev= (abs(dy - median(dy)) > rmsf*dy(rms)) | (abs(dx - median(dx)) > rmsf*dx(rms));
  if (dbg) print, "dev: ", dev;
  if (dbg) print, "=========================================================";
  
  //Return a smoothed dewarping set using a 2D polynomial fit to the
  //dewarping offsets.
  if (smooth) {    
    print, "dvv_FixWarpSet: applying polynomial smoothing";
    //Generate a 2D polynomial fit to the good data
    ax= fitsurf(gpts.x-rpts.x, rpts.y, rpts.x, degree= order);
    ay= fitsurf(gpts.y-rpts.y, rpts.y, rpts.x, degree= order);
    //Interpolate and/or extrapolate the fit to fix the outliers
    pts.x= ref.x + polysurf(ax, ref.y, ref.x);
    pts.y= ref.y + polysurf(ay, ref.y, ref.x);
  }
}

func dvv_processWarpSet(&pts, ref, &dpts, ch=, order=, mask=, weight=, dbg=)
/* DOCUMENT dvv_processWarpSet, pts, ref, dpts, ch=, order=, mask=, weight=, dbg=

     Processes a warp set matrix by weaving a series of 1d fits along
     each grid line in both the x- and y-directions together to get an
     average result for the dewarping offsets at each grid point.
     Replaces the dewarping offset array with the processed points.

     KEYWORDS:
       ch=
       order=
       mask=   a bitwise mask vector for masking out points as determined
               by an external criterion
       weight= weight to be used when applying the fit
       dbg=

     SEE ALSO:
       dvv_SetWarp, dvv_bicubicDewarp
     
*/
{
  if (is_void(order)) order= 2;

  dy= ref.y - pts.y;
  dx= ref.x - pts.x;
  dxx= dyx= dxy= dyy= array(0.0, dimsof(dy));

  ng= sqrt(numberof(pts));
  v= (ref.x)(1:int(ng));
  if (is_void(weight)) weight= array(1.0, dimsof(dy));
  if (is_void(mask)) mask= array(0, dimsof(dy));

  //Apply a series of fits along the x-direction (fixed y values)
  for(i= 1; i<= ng; i++) {
    selym= where(ref.y == v(i) & mask == 0);
    sely= where(ref.y == v(i));
    err= 1./weight(selym);
    xfit= poly1_fit(dx(selym), ref(selym).x, order, 1./(err*err));
    yfit= poly1_fit(dy(selym), ref(selym).x, order, 1./(err*err));
    dxx(sely)= poly1(ref(sely).x, xfit);
    dyx(sely)= poly1(ref(sely).x, yfit);
  }
  
  //Apply a series of fits along the y-direction (fixed x values)
  for(i= 1; i<= ng; i++) {
    selxm= where(ref.x == v(i) & mask == 0.0);
    selx= where(ref.x == v(i));
    err= 1./weight(selxm);
    xfit= poly1_fit(dx(selxm), ref(selxm).y, order, 1./(err*err));
    yfit= poly1_fit(dy(selxm), ref(selxm).y, order, 1./(err*err));
    dxy(selx)= poly1(ref(selx).y, xfit);
    dyy(selx)= poly1(ref(selx).y, yfit);
  }
  
  dxavg= (dxx + dxy)/2.;
  dxdif= (dxx - dxy);
  dyavg= (dyx + dyy)/2.;
  dydif= (dyx - dyy);

  dpts= pts;
  dpts.x= dxdif;
  dpts.y= dydif;

  if (am_subroutine()) {
    pts.x= ref.x - dxavg;
    pts.y= ref.y - dyavg;
  } else {
    npts= pts;
    npts.x= ref.x - dxavg;
    npts.y= ref.y - dyavg;
    return npts;
  }
}

// func dvv_FixWarpSet(&pts, ref, ch, order=)
// /* DOCUMENT dvv_FixWarpSet, pts, ref, ch, order=

//      Detects and repairs outliers found in the warp set

//      KEYWORDS:

//      SEE ALSO:
     
// */
// {
//   if(is_void(order)) order= 1;
//   dy= ref.y - pts.y;
//   dx= ref.x - pts.x;
//   dev= (abs(dy - dy(avg)) > 2.0) | (abs(dx - dx(avg)) > 2.0);
//   wout= where(dev);
//   wn= where(!dev);
//   gpts= pts(wn);
//   rpts= ref(wn);
//   if(numberof(wout) > 0) {
//     print, "Outliers, "+ch;
//     print, "================";
//     for(j= 1; j<= numberof(wout); j++) {
//       print, "ORIGINAL: ", ref(wout(j)).x, "-->", pts(wout(j)).x,
//              ref(wout(j)).y, "-->", pts(wout(j)).y;
      
//       //Generate repaired x-coordinate
//       wy= where(ref(wout(j)).x == rpts.x); 
//       Y= gpts(wy).x;
//       X= gpts(wy).y;
//       //print, "X- repair:"; 
//       //print, "Y -", Y; 
//       //print, "X -", X; 
//       newx= poly2(ref(wout(j)).y, fitpoly(order, X, Y));
      
//       //Generate repaired y-coordinate
//       wx= where(ref(wout(j)).y == rpts.y);
//       Y= gpts(wx).y;
//       X= gpts(wx).x;
//       //print, "Y- repair:"; 
//       //print, "Y -", Y; 
//       //print, "X -", X; 
//       newy= poly2(ref(wout(j)).x, fitpoly(order, X, Y));
      
//       //Reset fixed points
//       pts(wout(j)).x= newx;
//       pts(wout(j)).y= newy;
//       print, "REPAIRED: ", ref(wout(j)).x, "-->", pts(wout(j)).x,
//              ref(wout(j)).y, "-->", pts(wout(j)).y;
//     }
//   }
//   return;
// }


func dvv_bicubicDewarp(img, mapping, dbg=, raw=)
/* DOCUMENT dvv_bicubicDewarp, img, mapping, dbg=, raw=

   KEYWORDS:
     dbg= set to turn on debugging information
     raw= Use a reverse mapping transformation to interpolate the dewarped
          data directly from the raw CCD image
          
   SEE ALSO:
     dvv_SetWarp, dvv_FixWarpSet
 */
{
  extern _pcounter, _xx, _yy, ref;
  
  if (dbg) print, "bicubicDewarp: dbg= ", dbg, ", generate frame set, pc=",_pcounter," ... ";
  pc= _pcounter;
  tic, _pcounter++;

  //Dewarping points
  pts= *mapping.warp_pts;
  ng= int(sqrt(numberof(pts)));
  scl=(ref.x)(1:ng);
  ds= scl(dif)(avg);

  //Generate bicubic spline fit to the dewarping mesh
  xscl= yscl= span(scl(1)-ds, scl(0)+ds, ng+2);
  //yscl= span(-5, 5, 11)*100.0;
  //xwrp= ywrp= array(0.0, 11, 11);
  xwrp= ywrp= array(0.0, ng+2, ng+2);

  //Fill the center of the warp mesh with the data
  //xwrp(2:10, 2:10)=  reform(pts.x,9,9);
  //ywrp(2:10, 2:10)=  reform(pts.y,9,9);
  xwrp(2:ng+1, 2:ng+1)=  reform(pts.x,ng,ng);
  ywrp(2:ng+1, 2:ng+1)=  reform(pts.y,ng,ng);
  
  //Pad out the edges
  // xwrp(1, 2:10)= xwrp(2, 2:10) - 100.0;
  // xwrp(11, 2:10)= xwrp(10, 2:10) + 100.0;
  // xwrp(2:10, 1)= xwrp(2:10, 2);
  // xwrp(2:10, 11)= xwrp(2:10, 10);
  // xwrp(1,1)= xwrp(2,2) - 100.0;
  // xwrp(11,11)= xwrp(10,10) + 100.0;
  // xwrp(1,11)= xwrp(2,10) - 100.0;
  // xwrp(11,1)= xwrp(10,2) + 100.0;
  
  // ywrp(1, 2:10)= ywrp(2, 2:10);
  // ywrp(11, 2:10)= ywrp(10, 2:10);
  // ywrp(2:10, 1)= ywrp(2:10, 2) - 100.0;
  // ywrp(2:10, 11)= ywrp(2:10, 10) + 100.0;
  // ywrp(1,1)= ywrp(2,2) - 100.0;
  // ywrp(11,11)= ywrp(10,10) + 100.0;
  // ywrp(1,11)= ywrp(2,10) + 100.0;
  // ywrp(11,1)= ywrp(10,2) - 100.0;

  xwrp(1, 2:ng+1)= xwrp(2, 2:ng+1) - ds;
  xwrp(ng+2, 2:ng+1)= xwrp(ng+1, 2:ng+1) + ds;
  xwrp(2:ng+1, 1)= xwrp(2:ng+1, 2);
  xwrp(2:ng+1, ng+2)= xwrp(2:ng+1, ng+1);
  xwrp(1,1)= xwrp(2,2) - ds;
  xwrp(ng+2,ng+2)= xwrp(ng+1,ng+1) + ds;
  xwrp(1,ng+2)= xwrp(2,ng+1) - ds;
  xwrp(ng+2,1)= xwrp(ng+1,2) + ds;
  if (dbg) print, "xwrp = ", xwrp;
  
  ywrp(1, 2:ng+1)= ywrp(2, 2:ng+1);
  ywrp(ng+2, 2:ng+1)= ywrp(ng+1, 2:ng+1);
  ywrp(2:ng+1, 1)= ywrp(2:ng+1, 2) - ds;
  ywrp(2:ng+1, ng+2)= ywrp(2:ng+1, ng+1) + ds;
  ywrp(1,1)= ywrp(2,2) - ds;
  ywrp(ng+2,ng+2)= ywrp(ng+1,ng+1) + ds;
  ywrp(1,ng+2)= ywrp(2,ng+1) + ds;
  ywrp(ng+2,1)= ywrp(ng+1,2) - ds;
  if (dbg) print, "ywrp = ", ywrp;
  
  if(dbg) {
    print, "xscl = ", xscl;
    print, "yscl = ", yscl;
    window, 6; pli, xwrp - xscl(,-:1:ng+2);
    window, 7; pli, ywrp - yscl(-:1:ng+2,);
  }
  xtab= tl2cub(,xscl,yscl,xwrp);
  ytab= tl2cub(,xscl,yscl,ywrp);

  //Generate grid of interpolated dewarping data over
  //entire image
  xgrid= (*img.xscale)(,-:1:img.ny);
  ygrid= (*img.yscale)(-:1:img.nx,);
  xx= tl2cub(xtab, xgrid, ygrid);
  yy= tl2cub(ytab, xgrid, ygrid);

  if (dbg) {
    _xx= xx - xgrid;
    _yy= yy - ygrid;
    window, 0; pli, xx - xgrid;
    window, 1; pli, yy - ygrid;
  }
  
  //Resample the original image using the dewarping mesh

  // -- using bilinear interpolation
  //dw= img_copy(img);
  //dw.data= &interp2d(yy, xx, *img.data, *img.yscale, *img.xscale, fill= 0.0);

  // -- using bicubic interpolation
  if(!raw) {
    if (dbg) print, "Using original method ...";
    dw= img_resample(img, xx, yy, meth= 2);
  } else {
    if (dbg) print, "Dewarping with reverse map.";
    dvv_ReverseMap, mapping, yy, xx;
    dw= img_resample(img_extract(raw, dvv_xtrFrame(raw, *mapping.gridList)),
                  xx, yy, *img.xscale, *img.yscale, meth= 2);
    dw.shotid= img.shotid;
    dw.x_label= img.x_label;
    dw.y_label= img.y_label;
    dw.z_label= img.z_label;
    dw.x_unit= img.x_unit;
    dw.y_unit= img.y_unit;
    dw.z_unit= img.z_unit;
    dw.zrange= img.zrange;
  }
  
  if (dbg) print, "bicubicDewarp ... done", tac(pc);
  _pcounter-= 1;
  return dw;  
}

func dvv_InputMapping(img, lims=)
/* DOCUMENT dvv_InputMapping, img

   Driver function to prompt the user to select a series of grid
   vertices that are used to define the grid mapping function for
   the data set.  For each frame the pixel coordinates of the following
   set of vertices are returned (labeled in object coordinates below):
   
      Center vertex at x= 0, y= 0
      Lower left vertex at x= -300, y= -300
      Upper left vertex at x= -300, y= +300
      Upper right vertex at x= +300, y= +300
      Lower right vertex at x= + 300, y= -300

   KEYWORDS:
     lims=
     
   SEE ALSO:
     dvv_setMapping, dvv_GetGridCoords
 */
{
  print, "dvv_InputMapping ...";

  set= dvv_SetMapping();

  if(is_void(lims)) {
    ch1slims= [img.nx/2, img.nx, 1, img.ny/2, 0.0];
    ch1plims= [1, img.nx/2,      1, img.ny/2, 0.0];
    ch2slims= [img.nx/2, img.nx, img.ny/2, img.ny, 0.0];
    ch2plims= [1, img.nx/2,      img.ny/2, img.ny, 0.0];
  }
  
  set.ch1S= dvv_SetFrameMap(img, dvv_GetGridCoords(img, "ch1S", lims= ch1slims));
  print, set.ch1S;

  set.ch1P= dvv_SetFrameMap(img, dvv_GetGridCoords(img, "ch1P", lims= ch1plims));
  print, set.ch1P;

  set.ch2S= dvv_SetFrameMap(img, dvv_GetGridCoords(img, "ch2S", lims= ch2slims));
  print, set.ch2S;

  set.ch2P= dvv_SetFrameMap(img, dvv_GetGridCoords(img, "ch2P", lims= ch2plims));
  print, set.ch2P;


  //Flip the frame scales appropriately for mirror reflections
  set.ch1S.sclx*=  1.0;
  set.ch1S.scly*=  1.0;

  set.ch2S.sclx*=  1.0;
  set.ch2S.scly*= -1.0;

  set.ch1P.sclx*= -1.0;
  set.ch1P.scly*=  1.0;

  set.ch2P.sclx*= -1.0; 
  set.ch2P.scly*= -1.0;
  
  print, set;
  return set;
}

func dvv_GetGridCoords(g, prompt, lims=)
/* DOCUMENT getGridCoords, g, prompt, lims=

   Prompts the user to select a sequence of vertices for a gridded
   subframe of a data set.  The user must select the center vertex
   followed by the lower left, upper left, upper right and lower right
   at x= +/- 300 and y= +/- 300 microns.

   KEYWORDS:
     lims=  Specifies the extent of the frame, used to zero in on each
            vertex quickly

   SEE ALSO:
     dvv_InputMapping

*/
{
  extern _BINNING;
  if(is_void(_BINNING)) _BINNING= 1;
  
  wrst, 0;
  plist= [];
  palette, "stern.gp";
  //gd= img_data(g);
  //gd= abs(gd(dif,dif)(pcen,pcen));
  //sh, 0, g, lgscl= 1;
  sh, 0, g;
  if(is_void(lims)) {
    loc= 1;
    lims= limits();
  } else {
    loc= 0;
    limits, lims;
    xavg= lims(1:2)(avg);
    yavg= lims(3:4)(avg);
    xdif= lims(1:2)(dif);
    ydif= lims(3:4)(dif);
    print, "xavg= ", xavg, "yavg= ", yavg;
    print, "xdif= ", xavg, "ydif= ", ydif;
    dx= dy= 300/_BINNING;
    //ox= oy= 720/_BINNING;
    ox= oy= 660/_BINNING;
  }
  //if(loc)  {
    r= greg(0, prompt= prompt+": locate the center",count=1);
    limits, r.x1, r.x2, r.y1, r.y2;
  //} else limits, [xavg-dx, xavg+dx, yavg-dy, yavg+dy, 0.0];
  pts= gpts(0, prompt="Identify the center point",count=1);
  limits, lims;
  grow, plist, pts;
  xavg= pts.x; yavg= pts.y;

  f1= 0.34;
  f2= 0.34;
  if(loc) {
    r= greg(0, prompt= "locate the Lower Left [-300, -300]",count=1);
    limits, r.x1, r.x2, r.y1, r.y2;
  } else limits, [xavg-ox-dx, xavg-ox+dx, yavg-oy-dy, yavg-oy+dy, 0.0];
  pts= gpts(0, prompt="Identify the center point",count=1);
  limits, lims;
  grow, plist, pts;
  if(loc) {
    r= greg(0, prompt= "locate the Upper Left [-300, 300]",count=1);
    limits, r.x1, r.x2, r.y1, r.y2;
  } else limits, [xavg-ox-dx, xavg-ox+dx, yavg+oy-dy, yavg+oy+dy, 0.0];
  pts= gpts(0, prompt="Identify the center point",count=1);
  limits, lims;
  grow, plist, pts;
  if(loc) {
    r= greg(0, prompt= "locate the Upper Right [300, 300]",count=1);
    limits, r.x1, r.x2, r.y1, r.y2;
  } else limits, [xavg+ox-dx, xavg+ox+dx, yavg+oy-dy, yavg+oy+dy, 0.0];
  pts= gpts(0, prompt="Identify the center point",count=1);
  limits, lims;
  grow, plist, pts;
  if(loc) {
    r= greg(0, prompt= "locate the Lower Right [300, -300]",count=1);
    limits, r.x1, r.x2, r.y1, r.y2;
  } else limits, [xavg+ox-dx, xavg+ox+dx, yavg-oy-dy, yavg-oy+dy, 0.0];
  pts= gpts(0, prompt="Identify the center point",count=1);
  limits, lims;
  grow, plist, pts;
  return plist;
}

func dvv_AngleFromGrid(clist)
/* DOCUMENT dvv_AngleFromGrid(clist)

   Determines an average rotation angle from a list of 5 coordinates
   indicating points on the mapping grid.  The points are assumed to
   be listed in the following order:
   
    (0,0), (-300, -300), (-300, 300), (300, 300), (300, -300)

   The set of vectors from the center point to the following four
   corner points should be [-135, 135, 45, -45] degrees.  Returns the
   small difference between the actual angle and the nominal values.
   This value is used to rotate the data into a horizontal/vertical
   orientation during mapping operations.
     
   SEE ALSO:
     dvv_MagFromGrid, dvv_SetFrameMap, dvv_MapFrame
 */
{
  xl= clist.x;
  yl= clist.y;

  dy= yl(2:) - yl(1);
  dx= xl(2:) - xl(1);

  aa= atan(dy, dx) - [-135.0, 135.0, 45.0, -45.0]*pi/180.0;
  return aa(avg);
}

func dvv_MagFromGrid(clist)
/* DOCUMENT dvv_MagFromGrid, clist

   Determines an average scaling factor from a list of 5 coordinates
   indicating points on the mapping grid.  The points are assumed to
   be listed in the following order:
   
    (0,0), (-300, -300), (-300, 300), (300, 300), (300, -300)

   The set of vectors from the center point to the following four
   corner points should be sqrt(2)*300.0 microns.  Using the grid data
   this function returns the scaling factor required to map pixels to
   microns (magnification factor).
     
   SEE ALSO:
     dvv_AngleFromGrid, dvv_SetFrameMap, dvv_MapFrame
 */
{
  
  xl= clist.x;
  yl= clist.y;

  dy= yl(2:) - yl(1);
  dx= xl(2:) - xl(1);

  rr= sqrt(dy*dy + dx*dx);

  return sqrt(2)*300.0/rr(avg);
}

func dvv_Smooth(img, psf)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  if ((struct_type(psf) == "double" ||
       struct_type(psf) == "float") &&
       numberof(psf) == 1) {
    psdat= img_data(dvv_mkGaussian(img, psf));
  } else if (struct_type(psf) == "IMG_DAT") {
    psdat= img_data(psf);
  }

  sm= img_copy(img);
  imdat= img_data(img);
    //psdat= *psf.data;
  sm.data= &fftw_convolve(double(imdat), double(psdat), 0);
  //sm.data= &convoln(imdat, psdat);
  return sm;
}

func dvv_xtrFrame(img, clist, binning=)
/* DOCUMENT dvv_xtrFrame, img, clist, binning=
     
   SEE ALSO:
 */
{
  if(is_void(binning)) binning= 1;
  
  NX= img.nx;
  NY= img.ny;
  iC= iPts(img, clist(1));
  x1= min(max(1.0, iC.ix - 1040/binning), NX);
  x2= min(max(1.0, iC.ix + 1040/binning), NX);
  y1= min(max(1.0, iC.iy - 1040/binning), NY);
  y2= min(max(1.0, iC.iy + 1040/binning), NY);
  return iREGION(ix1= x1, ix2= x2, iy1= y1, iy2= y2); 
}

func dvv_SetFrameMap(img, gc)
/* DOCUMENT dvv_SetFrameMap, img, gc

   Extracts a frame from a raw data set, then rotates it and applies a
   linear rescaling to produce a rough mapping onto the experiment
   coordinate system.

   returns a dvv_FrameMapping structure

   ARGUMENTS:
     img - raw data image containing the frame
     gc  - identifying points of the frame: gc(1) is the center
    
   SEE ALSO:
 */
{
  extern _BINNING;
  if(is_void(_BINNING)) _BINNING= 1;

  print, "dvv_SetFrameMap ...";
  print, "dvv_SetFrameMap, gc=", gc;

  mapping= dvv_FrameMapping();
  mapping.gridList= &gc;
  mapping.angle= dvv_AngleFromGrid(gc);
  mapping.sclx= mapping.scly= dvv_MagFromGrid(gc);

  return mapping;
}

func dvv_MapFrame(img, mapping)
/* DOCUMENT dvv_mapFrame, img, mapping

   Extracts a frame from a raw data set, then rotates it and applies a
   linear rescaling to produce a rough mapping onto the experiment
   coordinate system.

   ARGUMENTS:
     img - raw data image containing the frame
     mapping - dvv_FrameMapping structure
    
   SEE ALSO:
 */
{
  extern _BINNING;
  if(is_void(_BINNING)) _BINNING= 1;
  
  gc= *mapping.gridList;

  smx= SCLSET(axis= "x", unit= "!mm", label= "position");
  smy= SCLSET(axis= "y", unit= "!mm", label= "position");
  smx.oldOrigin= gc(1).x;
  smy.oldOrigin= gc(1).y;

  //Extract the quadrant containing the frame
  //ximg= img_extract(img, dvv_xtrFrame(img, gc, binning= _BINNING));
  ximg= img_extract(img, dvv_xtrFrame(img, gc));
    
  //Rotation operation will cause degradation of high frequency data
  //Need to use a higher order interpolation scheme here
  ximg= img_rotate(ximg, mapping.angle, pt=gc(1), meth= 2);
  smx.scale= mapping.sclx;
  smy.scale= mapping.scly;

  //Apply linear scaling (magnification)
  return remapi(ximg, smx, smy);
}

func dvv_ReverseMap(mapping, &yg, &xg)
/* DOCUMENT dvv_ReverseMap, mapping, yg, xg

   From a given mapping specification and a set of object plane points
   (i.e. coordinates within the field of view, -450 < xg,yg < 450
   microns), determine the corresponding pixel locations on the raw
   image.  For each frame the mapping is different, so the same set of
   coordinates can be separately reverse-mapped onto the raw image.

   ARGUMENTS:
     mapping - dvv_FrameMapping structure containing the mapping parameters
     xg - x-coordinates of the input/output points
     yg - y-coordinates of the input/output points

   EXAMPLE:
     xg= (*map.warp_pts).x;         //Coordinates of the dewarping grid ...
     yg= (*map.warp_pts).y;
     dvv_ReverseMap, map, yg, xg;  // ... mapped to the raw CCD image
     plmk, yg, xg;

   WARNING:
     This routine replaces the input values in xg, yg with the new
     values
     
   SEE ALSO:
     dvv_InputMapping, dvv_SetFrameMap, dvv_MapFrame
 */
{
  gc= *mapping.gridList;
  X0= gc(1).x;
  Y0= gc(1).y;
  xg= xg/mapping.sclx;
  yg= yg/mapping.scly;
  rg= sqrt(xg*xg + yg*yg);
  aa= atan(yg, xg);
  xg= rg*cos(aa + mapping.angle);
  yg= rg*sin(aa + mapping.angle);
  xg+=X0;
  yg+=Y0;
  return;
}

func dvv_mkCross(img, width, length, xoff, yoff, box=)
/* DOCUMENT dvv_mkCross, img, width, length, xoff, yoff, box=

   Generates a synthetic image mapped onto the same coordinates as the
   input "img".  The returned image contains a uniform background over
   which a pattern of a cross is placed.  The cross consists of arms
   of specified width and length, centered at location (xoff, yoff);
   within these features the intensity is set to zero.

   Arguments:
     img
     fwhm
     length
     xoff
     yoff

   KEYWORDS:
     box=   
     
   SEE ALSO:
    dvv_mkGaussian, dvv_mkBox
 */
{
  if(is_void(xoff)) xoff= 0.0;
  if(is_void(yoff)) yoff= 0.0;
  cx= img_copy(img);
  dat= array(1.0, cx.nx, cx.ny);
  hbary= where(*cx.yscale < yoff+width/2.0 & *cx.yscale > yoff-width/2.0);
  hbarx= where(*cx.xscale < xoff+length/2.0 & *cx.xscale > xoff-length/2.0);
  dat(hbarx, hbary)= 0.0;
  vbarx= where(*cx.xscale < xoff+width/2.0 & *cx.xscale > xoff-width/2.0);
  vbary= where(*cx.yscale < yoff+length/2.0 & *cx.yscale > yoff-length/2.0);
  dat(vbarx, vbary)= 0.0;
  if(!is_void(box)) {
    boxy= where(*cx.yscale <= yoff-length/2.0 | *cx.yscale >= yoff+length/2.0);
    boxx= where(*cx.xscale <= xoff-length/2.0 | *cx.xscale >= xoff+length/2.0);
    dat(boxx, )= 0.0;
    dat(,boxy)= 0.0;
  }
  cx.data= &dat;
  cx.shotid= "Cross";
  return cx;
}

func dvv_mkGaussian(img, fwhm)
/* DOCUMENT dvv_mkGaussian, img, fwhm

   Generates a synthetic image mapped onto the same coordinates as the
   input "img".  The returned image contains a gaussian intensity
   distribution with the specified fwhm. The distribution is scaled so
   that its integral (summed over all pixel elements) is unity.  This
   image type is designed to be used with dvv_Smooth to allow
   smoothing of images with FFT convolution.

   Arguments:
     img
     fwhm
     
   SEE ALSO:
     dvv_mkBox, dvv_mkCross, dvv_Smooth
 */
{
  cx= img_copy(img);
  xgrid= img_grid(img, 1);
  ygrid= img_grid(img, 2);
  r2= xgrid*xgrid + ygrid*ygrid;
  tau2= (fwhm/2.0)^2/log(2.0);
  dd= array(0.0, dimsof(r2));
  wn= where(abs(r2/tau2) <= 100.0);
  if(numberof(wn) > 0) dd(wn)= exp(-r2(wn)/tau2)/pi/tau2;
  dsum= dd(sum,sum);
  cx.data= &(dd/dsum);
  return cx;
}

func dvv_mkBox(img, fwhm)
/* DOCUMENT dvv_mkBox, img, fwhm

   Generates a synthetic image mapped onto the same coordinates
   as the input "img".  The returned image contains a square
   box intensity distribution with side set to fwhm

   The distribution is scaled so that its integral (summed over
   all pixel elements) is unity.

   Arguments:
     img
     fwhm
     
   SEE ALSO:
     dvv_mkGaussian, dvv_mkCross, dvv_Smooth
 */
{
  cx= img_copy(img);
  hwidth= fwhm/2.0;
  xb= where(*cx.xscale < hwidth & *cx.xscale > -hwidth);
  yb= where(*cx.yscale < hwidth & *cx.yscale > -hwidth);
  dd= array(0.0, img.nx, img.ny);
  dd(xb, yb)= 1.0;
  dsum= dd(sum,sum);
  cx.data= &(dd/dsum);
  return cx;
}

func dvv_Correlate(img1, img2, reg=)
/* DOCUMENT cimg= dvv_Correlate(img1, img2, reg=)

   Returns the 2D correlation img1 & img2
     
   SEE ALSO:
 */
{
  corr= img_correlate(img1, img2, reg= reg);
  corr.shotid= "Correlation: "+img1.shotid+" + "+img2.shotid;
  return corr;
}

func dvv_CoordMin(img, pt, radius, r2, dbg=, method=, offset=)
/* DOCUMENT pt= dvv_CoordMin(img, pt, radius, r2, dbg=, offset=)

   Returns a POINT data structure containing the coordinates of
   the point corresponding to local grid intersection point
   in an image frame.
     img - IMG_DAT structure containing input image
     pt - location to start looking
     radius - initial search radius
     r2 - refined (second stage) search radius
     
   

   KEYWORDS:
     dbg=    set to turn on debugging output
     method= minimization method for refined output
             "brute" - brute force evaluation of a large interpolated table
     offset= correlation offset (subtracted from final result)
     
   SEE ALSO:
 */
{

  if(is_void(offset)) offset= [0.0, 0.0]; 
  if(is_void(method)) method = "brute";
  ix= digitize(pt.x, *img.xscale);
  iy= digitize(pt.y, *img.yscale);
  xgrid= (*img.xscale)(,-:1:img.ny);
  ygrid= (*img.yscale)(-:1:img.nx,);
  dat= *img.data;

  // Find pixel of local maximum
  xrge= where(*img.xscale > pt.x - radius & *img.xscale < pt.x + radius);
  yrge= where(*img.yscale > pt.y - radius & *img.yscale < pt.y + radius);
  dd= dat(xrge, yrge);
  xg= xgrid(xrge, yrge);
  yg= ygrid(xrge, yrge);
  r= (xg - pt.x)^2 + (yg - pt.y)^2;
  w= where(r > radius^2);
  dd= dd - dd(avg,avg);
  dd(w)= dd(min,min);
  wx= where(dd == dd(max,max));
  //print, "wx= ", wx;
  nxx= numberof(xrge); nyy= numberof(yrge);
  iy= (wx-1)/nxx + 1;
  ix= (wx-1)%nxx + 1;
  x= ((*img.xscale)(xrge))(ix);
  y= ((*img.yscale)(yrge))(iy);
  //return POINT(x= x, y= y);
  
  // Refinement: integrate in x- and y- directions and locate
  // the maximum by interpolation of the difference
  xrge= where(*img.xscale > x - r2 & *img.xscale < x + r2);
  yrge= where(*img.yscale > y - r2 & *img.yscale < y + r2);
  xscl= (*img.xscale)(xrge);
  yscl= (*img.yscale)(yrge);

  //Generate interpolation
  dtab= tl2cub(,xscl, yscl, dat(xrge, yrge));
  
  if(method == "brute") {
    nx= numberof(xrge)*50;
    ny= numberof(yrge)*50;
    xg= span(xscl(1), xscl(0), nx)(,-:1:ny);
    yg= span(yscl(1), yscl(0), ny)(-:1:nx,);
    dd= tl2cub(dtab, xg, yg);
    wmax= where(dd == dd(max,max));
    ix= (wmax - 1)%nx + 1;
    iy= (wmax - 1)/nx + 1;
    if(numberof(ix) > 1) error, "maximum over more than one point.";
  } else if (method == "minimization") {
    ///extern dtab;
    error, "Currently not implemented ... ";
    //Implementation note: probably compute the auto correlation
    //of the cross function centered at (0,0); then use this
    //to compute an interpolation that can be re-centered
    //near the actual cross-point. Use lmfit to fit relative to
    //F= a0 * F_auto(x - x0, y - y0)
  }

  if(!is_void(dbg)) {
    extern _dd, _dat;
    //print, "xscl= ", xscl;
    //print, "xg(,1)= ", xg(,1);
    //print, "yg(1,)= ", yg(1,);
    print, "dd(max,max)= ", dd(max,max);
    print, "where dd(max,max)= ", where(dd == dd(max,max));
    print, "ix= ", ix, "iy= ", iy;
    _dd= dd;
    _dat= dat(xrge, yrge);
    info, dat(xrge, yrge);
    info, dd;
    imm= IMG_DAT(nx= nx, ny= ny, xscale= &xg(,1), yscale= &yg(1,));
    imm.data= &dd;
    sh, 0, imm;
    window, 1;
    pli, dat(xrge, yrge);
    print, "dat average", dat(xrge, yrge)(avg,avg); 
  }
  
//   xdiff= integ(dd, xg(,1), xg(0,1), 1)(dif);
//   ydiff= integ(dd, yg(1,), yg(1,0), 2)(dif);
//   ny= numberof(ydiff);
//   nx= numberof(xdiff);
//   xdiff= _(monotonic(xdiff(nx/2:1:-1), incr= 1)(0:1:-1),
//            monotonic(xdiff(nx/2+1:nx), incr= -1));
//   ydiff= _(monotonic(ydiff(ny/2:1:-1), incr= 1)(0:1:-1),
//            monotonic(ydiff(ny/2+1:ny), incr= -1));
  //d2y= where(ydiff(pcen)(dif) < 0);
  //d2x= where(xdiff(pcen)(dif) < 0);
  //if(!is_void(dbg)) {
    //plg, "xg= ", xg(,1);
    //plg, "yg= ", yg(1,);
    //window, 2;
    //plg, xdiff, xg(zcen,1), color= "red";
    //plg, ydiff, yg(1,zcen), color= "blue";
    //print, "xdiff = ", xdiff, xg(dif, 1);
    //print, "ydiff= ", ydiff, yg(1, dif);
  //}
  //   xx= integ(dd*xg, xg(,1), xg(0,1), 1)/integ(dd, xg(,1), xg(0,1), 1);
  //   yy= integ(dd*yg, yg(1,), yg(1,0), 2)/integ(dd, yg(1,), yg(1,0), 2);
  //   return POINT(x= xx(avg), y= yy(avg));
  //xx= interp(xg(zcen, 1), xdiff, 0.0);
  //yy= interp(yg(1, zcen), ydiff, 0.0);
  if(!is_void(dbg)) {
    window, 0;
    plmk, yg(1,iy) - offset(2), xg(ix,1) - offset(1);
  }
  return POINT(x= xg(ix,1) - offset(1), y= yg(1,iy) - offset(2));
}

// func dvv_CMin(dtab, a)
// /* DOCUMENT 
     
//    SEE ALSO:
//  */
// {
//   extern dtab;
//   return tl2cub(xtab, xvec, yvec);
// }

func dvv_SetChFilters(fs, rad, w=, common=, highfreq=)
/* DOCUMENT dvv_SetChFilters, fs, rad, w=, common=, highfreq=

   Drives a user-input sequence using mouse clicks applied to 2D power
   spectra displayed for each of the 4 channels, in order to isolate
   and filter out the high frequency parasitic ghost fringe content
   that is present in each OHRV channel.

   The user can input up to _DVV_CHFILTER_MAXCOUNT clicks per channel
     
   SEE ALSO:
 */
{
  extern _DVV_CH1S_FILTER, _DVV_CH2S_FILTER;
  extern _DVV_CH1P_FILTER, _DVV_CH2P_FILTER;
  extern _DVV_CMMN_FILTER, _DVV_CHFILTER_MAXCOUNT;
  extern _DVV_CHS_FILTER, _DVV_CHP_FILTER;

  if (is_void(w)) w= 0;
  if (is_void(rad)) rad= 0.01;
  window, w;
  fma;

  if (!common && !highfreq) {
    _DVV_CH1S_FILTER= [];
    _DVV_CH2S_FILTER= [];
    _DVV_CH1P_FILTER= [];
    _DVV_CH2P_FILTER= [];
    ch1S= img_copy(fs.ch1S);
    ch2S= img_copy(fs.ch2S);
    sdiff= img_div(img_sub(ch2S, ch1S), img_add(ch1S, ch2S, average= 1));
    dvv_DisplayImg, dvv_PSD(ch1S, wndw="hanning"), 1, 1e4, lgscl= 1, w= w+1, pltitl="Ch1S";
    plsys, 1;
    limits, -0.2, 0.2, -0.2, 0.2; limsReset;
    dvv_DisplayImg, dvv_PSD(sdiff, wndw="hanning"), 1e-8, 1e-5, lgscl= 1, w= w, pltitl="S-channels";
    plsys, 1;
    limits, -0.2, 0.2, -0.2, 0.2; limsReset;
    mode= gpts(w, prompt= "Select ch1s parasitic mode", count= _DVV_CHFILTER_MAXCOUNT, plcircl= rad);
    for(i= 1; i <= numberof(mode); i++)
      _DVV_CH1S_FILTER= _(_DVV_CH1S_FILTER,
                          swrite(format= "mode_reject %8.5f %8.5f %8.5f", mode(i).x, mode(i).y, rad));
    pause, 500;

    window, w+1; fma;
    dvv_DisplayImg, dvv_PSD(ch2S, wndw="hanning"), 1, 1e4, lgscl= 1,w= w+1, pltitl="Ch2S";
    plsys, 1;
    limits, -0.2, 0.2, -0.2, 0.2; limsReset;
    mode= gpts(w, prompt= "Select ch2s parasitic mode", count= _DVV_CHFILTER_MAXCOUNT, plcircl= rad);
    for(i= 1; i <= numberof(mode); i++)
      _DVV_CH2S_FILTER= _(_DVV_CH2S_FILTER,
                          swrite(format= "mode_reject %8.5f %8.5f %8.5f", mode(i).x, mode(i).y, rad));
    pause, 500;

    ch1P= img_copy(fs.ch1P);
    ch2P= img_copy(fs.ch2P);
    pdiff= img_div(img_sub(ch2P, ch1P), img_add(ch1P, ch2P, average= 1));
    window, w+1; fma;
    dvv_DisplayImg, dvv_PSD(ch1P, wndw="hanning"), 1, 1e4, lgscl= 1, w= w+1, pltitl="Ch1P";
    plsys, 1;
    limits, -0.2, 0.2, -0.2, 0.2; limsReset;
    dvv_DisplayImg, dvv_PSD(pdiff, wndw="hanning"), 1e-8, 1e-5, lgscl= 1, w= w, pltitl="P-channels";
    plsys, 1;
    limits, -0.2, 0.2, -0.2, 0.2; limsReset;
    mode= gpts(w, prompt= "Select ch1p parasitic mode", count= _DVV_CHFILTER_MAXCOUNT, plcircl= rad);
    for(i= 1; i <= numberof(mode); i++)
      _DVV_CH1P_FILTER= _(_DVV_CH1P_FILTER,
                          swrite(format= "mode_reject %8.5f %8.5f %8.5f", mode(i).x, mode(i).y, rad));
    pause, 500;

    window, w+1; fma;
    dvv_DisplayImg, dvv_PSD(ch2P, wndw="hanning"), 1, 1e4, lgscl= 1, w= w+1, pltitl="Ch2P";
    plsys, 1;
    limits, -0.2, 0.2, -0.2, 0.2; limsReset;
    mode= gpts(w, prompt= "Select ch2p parasitic mode", count= _DVV_CHFILTER_MAXCOUNT, plcircl= rad);
    for(i= 1; i <= numberof(mode); i++)
      _DVV_CH2P_FILTER= _(_DVV_CH2P_FILTER,
                          swrite(format= "mode_reject %8.5f %8.5f %8.5f", mode(i).x, mode(i).y, rad));
    pause, 500;

    print, "CH1S Filter = ", _DVV_CH1S_FILTER;
    print, "CH2S Filter = ", _DVV_CH2S_FILTER;
    print, "CH1P Filter = ", _DVV_CH1P_FILTER;
    print, "CH2P Filter = ", _DVV_CH2P_FILTER;
  } 

  if (highfreq) {
    _DVV_CHS_FILTER= [];
    _DVV_CHP_FILTER= [];
    ch1S= img_copy(fs.ch1S);
    ch2S= img_copy(fs.ch2S);
    sdiff= img_div(img_sub(ch2S, ch1S), img_add(ch1S, ch2S, average= 1));
    dvv_DisplayImg, dvv_PSD(sdiff, wndw="hanning"), 1e-9, 5e-6, lgscl= 1, w= w, pltitl="S-channels";
    plsys, 1;
    limits, -0.5, 0.5, -0.5, 0.5; limsReset;
    mode= gpts(w, prompt= "Select s-channel high freq mode", count= _DVV_CHFILTER_MAXCOUNT, plcircl= rad);
    for(i= 1; i <= numberof(mode); i++)
      _DVV_CHS_FILTER= _(_DVV_CHS_FILTER,
                         swrite(format= "mode_reject %8.5f %8.5f %8.5f", mode(i).x, mode(i).y, rad));
    pause, 500;

    ch1P= img_copy(fs.ch1P);
    ch2P= img_copy(fs.ch2P);
    pdiff= img_div(img_sub(ch2P, ch1P), img_add(ch1P, ch2P, average= 1));
    dvv_DisplayImg, dvv_PSD(pdiff, wndw="hanning"), 1e-9, 5e-6, lgscl= 1, w= w, pltitl="P-channels";
    plsys, 1;
    limits, -0.5, 0.5, -0.5, 0.5; limsReset;
    mode= gpts(w, prompt= "Select p-channel high freq mode", count= _DVV_CHFILTER_MAXCOUNT, plcircl= rad);
    for(i= 1; i <= numberof(mode); i++)
      _DVV_CHP_FILTER= _(_DVV_CHP_FILTER,
                         swrite(format= "mode_reject %8.5f %8.5f %8.5f", mode(i).x, mode(i).y, rad));

    print, "CHS Filter = ", _DVV_CHS_FILTER;
    print, "CHP Filter = ", _DVV_CHP_FILTER;
  }
  
  if (common) {
    _DVV_CMMN_FILTER= [];
    ch1P= img_copy(fs.ch1P);
    ch2P= img_copy(fs.ch2P);
    pdiff= img_div(img_sub(ch2P, ch1P), img_add(ch1P, ch2P, average= 1));
    window, w+1; fma;
    dvv_DisplayImg, dvv_PSD(ch1P, wndw="hanning"), 1, 1e4, lgscl= 1, w= w+1, pltitl="Ch1P";
    plsys, 1;
    limits, -0.1, 0.1, -0.1, 0.1; limsReset;
    dvv_DisplayImg, dvv_PSD(pdiff, wndw="hanning"), 1e-8, 1e-5, lgscl= 1, w= w, pltitl="P-channels";
    plsys, 1;
    limits, -0.1, 0.1, -0.1, 0.1; limsReset;
    mode= gpts(w, prompt= "Select common mode filter", count= _DVV_CHFILTER_MAXCOUNT, plcircl= rad);
    for(i= 1; i <= numberof(mode); i++)
      _DVV_CMMN_FILTER= _(_DVV_CMMN_FILTER,
                          swrite(format= "mode_reject %8.5f %8.5f %8.5f", mode(i).x, mode(i).y, rad));
    pause, 500;
    print, "Common Mode Filter = ", _DVV_CMMN_FILTER;
  }
}

func dvv_SetupCommonSumFilter(..)
{
  extern _DVV_CMMN_FILTER, _DVV_CMMN_SUM_FILTER;
  
  if(numberof(_DVV_CMMN_FILTER) != 2)
    error, "Common filter set requires two and only two modes.";
  
  m1= strread(parse_line(_DVV_CMMN_FILTER(1), " ", 4)(2:0), "%e");
  m2= strread(parse_line(_DVV_CMMN_FILTER(2), " ", 4)(2:0), "%e");
  //Need to sum the filter coordinates with the complex conjugate
  //of the opposing lobe
  mtot= [m1, -m2](,avg);
  _DVV_CMMN_SUM_FILTER= swrite(format= "mode_reject %8.5f %8.5f %8.5f", mtot(1), mtot(2), m1(3));
}

func dvv_RegionList(region_array)
/* DOCUMENT dvv_RegionList(region_array)

   Processes an array of regions to produce a string containing the
   region specifications as a single comma-delimited string that can
   be used in the command-line argument list of the mkOHRV command.

   e.g. if _DVV_VEL_BKG = [REGION(x1=12.3,x2=34.5,y1=78.9,y2=92.1),
                           REGION(x1=32.3,x2=38.5,y1=18.9,y2=29.1)]

   then dvv_RegionList(_DVV_VEL_BKG) returns:
   
   "region,12.3,34.5,78.9,92.1,region,32.3,38.5,18.9,29.1"

   To specify this region array and store it in _DVV_VEL_BKG in
   the command line use:

   mkOHRV -s <shot> -vbg region,12.3,34.5,78.9,92.1,region,32.3,38.5,18.9,29.1
     
   SEE ALSO:
     dvv_ParseRegionList
 */
{
  if (struct_type(region_array) != "REGION") error, "input argument must be of struct type REGION()";
  for(i= 1; i<= numberof(region_array); i++) {
    r= region_array(i);
    if (i== 1) rstr= swrite(format="region,%f,%f,%f,%f",r.x1,r.x2,r.y1,r.y2);
    else rstr+= swrite(format=",region,%f,%f,%f,%f",r.x1,r.x2,r.y1,r.y2);
  }
  return rstr;
}

func dvv_ParseRegionList(rstr)
/* DOCUMENT dvv_ParseRegionList(region_list)

   Processes a region list string to produce an array of REGION
   structures corresponding to the region list string.  The
   region_list is usually generated from dvv_RegionList, and entered
   as a string in the command line of mkOHRV.
   
   e.g. if region_list =
        "region,12.3,34.5,78.9,92.1,region,32.3,38.5,18.9,29.1"

   returns:
         [REGION(x1=12.3,x2=34.5,y1=78.9,y2=92.1),
          REGION(x1=32.3,x2=38.5,y1=18.9,y2=29.1)]
          
   SEE ALSO:
     dvv_RegionList
 */
{
  if(numberof(rstr) == 1) l= parse_line(rstr, ",");
  else l= rstr;
  rg= [];
  do {
    p= l(1:5);
    if(numberof(l) > 5) l= l(6:0);
    else l= [];
    if (p(1) != "region") {
      if(is_void(rg)) error, "apparently not a region list";
      else break;  // Have reached the end of the list
    }
    xy= strread(p(2:0), "%e");
    if(allof(xy) == 0.0) break;  //Can't accept all zeros, end of list?
    rg= _(rg, REGION(x1= xy(1), x2= xy(2), y1= xy(3), y2= xy(4)));
  } while (numberof(l) >= 5);
  return rg;
}

func dvv_FilterArg(filt)
/* DOCUMENT dvv_FilterArg(filt)

   Processes a string array containing filter array specification
   parameters and returns a single comma-delimited string that can be
   used in the command-line argument list of the mkOHRV command.  The
   filter array is stored in a format that dvv_GenFilterMaskSet can
   recognize; but this format has to be translated into a different
   text form for external input from the command line.

   e.g. if _DVV_CH1S_FILTER =
        ["mode_reject -0.08579  0.10258  0.01000",
         "mode_reject -0.07702  0.09632  0.01000",
         "mode_reject -0.08454  0.09256  0.01000"]
   which specifies a three-lobed filter, then dvv_FilterArg(_DVV_CH1S_FILTER)
   returns:
   
   "mode_reject,-0.08579,0.10258,0.01000,mode_reject,-0.07702,0.09632,
        0.01000,mode_reject,-0.08454,0.09256,0.01000"

   To specify this filter array and store it in _DVV_CH1S_FILTER in
   the command line use:

   mkOHRV -s <shot> -fc1s mode_reject,-0.08579,0.10258,0.01000,\
        mode_reject,-0.07702,0.09632,0.01000,mode_reject,-0.08454,\
        0.09256,0.01000
     
   SEE ALSO:
     dvv_SetChFilters, dvv_ChFilters, dvv_GenFilterMaskSet
 */
{
  pl= parse_line((filt+" ")(sum), " ", 20);
  pl= pl(where(pl!= ""));
  if(!is_void(pl)) return strpart((pl+",")(sum), :-1);
  else return "";
}

func dvv_ListFilters(f)
/* DOCUMENT dvv_ListFilters, f

   Prints out 4 lines giving the text output that parameterizes a set
   of channel filters.  The text output can be pasted into a mkOHRV
   configuration file such that this set of filter specifications will
   be applied during OHRV command-line processing (via a makefile).
     
   SEE ALSO:
     dvv_SetChFilters, dvv_FilterArg
 */
{
  fmt0="%s\n";
  if(!is_void(_DVV_CH1S_FILTER))
    write, f, format= fmt0, "FCH1S="+dvv_FilterArg(_DVV_CH1S_FILTER);
  if(!is_void(_DVV_CH2S_FILTER))
    write, f, format= fmt0, "FCH2S="+dvv_FilterArg(_DVV_CH2S_FILTER);
  if(!is_void(_DVV_CH1P_FILTER))
    write, f, format= fmt0, "FCH1P="+dvv_FilterArg(_DVV_CH1P_FILTER);
  if(!is_void(_DVV_CH2P_FILTER))
    write, f, format= fmt0, "FCH2P="+dvv_FilterArg(_DVV_CH2P_FILTER);
  if(!is_void(_DVV_CHS_FILTER))
    write, f, format= fmt0, "FCHS="+dvv_FilterArg(_DVV_CHS_FILTER);
  if(!is_void(_DVV_CHP_FILTER))
    write, f, format= fmt0, "FCHP="+dvv_FilterArg(_DVV_CHP_FILTER);
  if(!is_void(_DVV_CMMN_FILTER))
    write, f, format= fmt0, "FCMMN="+dvv_FilterArg(_DVV_CMMN_FILTER);
}

func dvv_SelectMode(dat, rad, zmin=, zmax=, w=, lims=)
/* DOCUMENT dvv_SelectMode, fs, rad, w=

   Drives a user-input sequence using mouse clicks applied to 2D power
   spectra displayed for a given data set to allow the user to identify
   a single mode to be selected for analysis.

   The user can input up to _DVV_CHFILTER_MAXCOUNT clicks per channel
     
   SEE ALSO:
 */
{
  if (is_void(w)) w= 0;
  if (is_void(rad)) rad= 0.01;
  window, w;
  fma;

  psd= dvv_PSD(dat, wndw= "hanning");
  dvv_DisplayImg, psd, zmin, zmax, w= w;
  plsys, 1;
  if(!is_void(lims)) limits, lims; else limits, -0.2, 0.2, -0.2, 0.2;
  limsReset;
  mode= gpts(w, prompt= "Select the data mode", count= 1, plcircl= rad);
  mode= swrite(format= "mode_select %8.5f %8.5f %8.5f", mode(i).x, mode(i).y, rad);
  write, format="%s\n", "The mode is: "+mode;
  mf= dvv_GenFilterMask(dat, mode);
  mode_rms= dvv_PSD_rms(psd, mask= mf);
  write, format="%s\n", "The mode rms, amplitude and p-t-p is: ",
    mode_rms, mode_rms*sqrt(2), mode_rms*2*sqrt(2);  
}

func dvv_getVelBkg(img, w=)
/* DOCUMENT dvv_getVelBkg, img, w=

   Obtains mouse-directed input from the user to identify a list of regions
   that will be used for background velocity correction.
     
   SEE ALSO:
     dvv_RegionList
 */
{
  if(is_void(w)) w= 0;
  winkill, w;
  yn= "";
  rgs= [];
  sh, w, img;
  do {
    limits;
    r= greg(w, prompt= "Isolate the area for selections ... >",count=1);
    limits, r.x1, r.x2, r.y1, r.y2;
    do {
      rgs= _(rgs, greg(w, ret, prompt="Identify a region to be added ... > ", plot_box= 1, count= 1));
    } while (!ret);
    read, yn, prompt= "Continue over another area?";
  } while (yn == "y");
  return rgs;
}

func dvv_getAnalysisBoxCenter(img, w=, size=, rge=)
/* DOCUMENT dvv_getAnalysisBox, img, w=, size=

   Obtains mouse-directed input from the user to identify a list of regions
   that will be used for background velocity correction.
     
   SEE ALSO:
     dvv_RegionList
 */
{
  if(is_void(w)) w= 0;
  if(is_void(rge)) {ll= []; uu= [];}
  else {ll= rge(1); uu= rge(2);}
  winkill, w;
  yn= "";
  rgs= [];
  dvv_DisplayImg, img, ll, uu, w= w;
  do {
    //limits;
    //r= greg(w, prompt= "Isolate the area for selections ... >",count=1);
    //limits, r.x1, r.x2, r.y1, r.y2;
    if(is_void(size)) {
      r= ""; read, r, prompt= "Enter box size (or 0 to cancel)? >";
      b= 0.0; nr= sread(r, b);
    } else b= size;
    ctr= gpts(w, prompt= "Identify the box center ... > ", count= 1);
    box= REGION(x1= ctr.x - b/2, x2= ctr.x+b/2, y1= ctr.y-b/2, y2= ctr.y+b/2);
    lh= LO_POS(or= "h", x1= ctr.y+1.0, x2= ctr.y-1.0);
    lv= LO_POS(or= "v", x1= ctr.x+1.0, x2= ctr.x-1.0);
    fma;
    dvv_DisplayImg, img, ll, uu, w= w, lh= lh, lv= lv, shlo= 1, regs= box;
    read, yn, prompt= "Satisfactory?";
  } while (yn == "n");
  return [ctr.x, ctr.y](*);
}

func dvv_getBreakoutCenter(img, w=, radius=, rge=)
/* DOCUMENT dvv_getBreakouCenter, img, w=, radius=, rge=

   Obtains mouse-directed input from the user to identify a list of regions
   that will be used for background velocity correction.
     
   SEE ALSO:
     dvv_RegionList
 */
{
  if(is_void(w)) w= 0;
  if(is_void(rge)) {ll= []; uu= [];}
  else {ll= rge(1); uu= rge(2);}
  winkill, w;
  yn= "";
  rgs= [];
  dvv_DisplayImg, img, ll, uu, w= w;
  do {
    //limits;
    //r= greg(w, prompt= "Isolate the area for selections ... >",count=1);
    //limits, r.x1, r.x2, r.y1, r.y2;
    if(is_void(radius)) {
      r= ""; read, r, prompt= "Enter circle radius (or 0 to cancel)? >";
      rad= 0.0; nr= sread(r, rad);
    } else rad= radius;
    ctr= gpts(w, prompt= "Identify the center ... > ", count= 1);
    lh= LO_POS(or= "h", x1= ctr.y+1.0, x2= ctr.y-1.0);
    lv= LO_POS(or= "v", x1= ctr.x+1.0, x2= ctr.x-1.0);
    fma;
    dvv_DisplayImg, img, ll, uu, w= w, lh= lh, lv= lv, shlo= 1;
    plsys, 1;
    dvv_plCircle, [ctr.x, ctr.y](*), rad, color= "magenta", width= 3, type= 2, marks= 0;
    read, yn, prompt= "Satisfactory?";
  } while (yn == "n");
  return [ctr.x, ctr.y](*);
}

func dvv_setOption(opt, str)
{
  extern shot;
  f= open(shot+".cfg", "rw");
  
  if (opt == "bxctr") swrite(format="BXCTR=%7.4f,%7.4f",ctr.x,ctr.y);
}

func dvv_getFringeMode(img)
/* DOCUMENT dvv_getFringeMode, img

   Analyzes a reference interferogram to extract the bakground fringe
   mode in frequency coordinates, i.e. (fx, fy) of the background
   fringe pattern.  The mode is found by searching in the lower
   half-plane of the PSD in the interval 0.0025 < fx < 1 for the
   highest peak.
     
   SEE ALSO:
     _FRINGE_MODE
*/
{
  psd= dvv_PSD(img, wndw= "hanning");
  xg= img_grid(psd, 1);
  yg= img_grid(psd, 2);
  half= where(xg > 0.0025);
  ix= img_data(psd)(half)(mxx);
  return [xg(half)(ix), yg(half)(ix)];
}

func dvv_bkgPhase(fs, wrap=, invert=)
/* DOCUMENT phs= dvv_bkgPhase(fs, wrap=, invert=)

   Returns a phase ramp that matches the dominant mode of the
   background fringe pattern found in an image.

   KEYWORDS:
     wrap= generate a wrapped phase verion of the ramp
     
   SEE ALSO:
     dvv_getFringeMode
 */
{
  if (!is_void(invert) && invert) {
    //print, "inverting the phase";
    pf= -1.0; 
  } else {
    pf= 1.0;
  }
 
  fm= dvv_getFringeMode(img_sub(fs.ch2S, fs.ch1S));
  xg= img_grid(fs.ch1S, 1);
  yg= img_grid(fs.ch1S, 2);
  p= pf*2*pi*(fm(1)*xg + fm(2)*yg);
  if (wrap) return img_wrap(img_copy(fs.ch1S, data= p)); else
    return img_copy(fs.ch1S, data= p);
}

func dvv_imgOffset(img1, img2, pt=, reg=, dbg=, method=)
/* DOCUMENT dvv_imgOffset, img1, img2, pt=, reg=, dbg=, method=

   Computes an offset between two nominally equivalent images (and/or
   subregions thereof), e.g. between the ch1S and ch2S frames of a
   grid image.  An image cross-correlation is applied to compare the
   relative positions of the two images, and the resulting
   cross-correlation peak is examined within a 20 micron square region
   near zero. Returned value is [delta-x, delta-y, err].  For a
   successful offset determination err= 0.

   The function uses non-linear fitting to determine the parameters of
   a 2D gaussian fitted onto the cross-correlation peak.  If the peak
   is ill-defined in the cross-correlation result then the fitting
   attempt fails and the 3rd element of the return vector is non-zero.

   KEYWORDS:
     pt=  [x,y] or POINT(x=X,y=Y), analysis centers on this point
     reg= a subregion of the image, should contain pt in the center
     dbg= 
     
   SEE ALSO:
     dvv_region, dvv_frameOffsets, dvv_reportOffsets
 */
{
  extern _DVV_CORR_OFFSET;
  extern _DVV_AC, _DVV_CC, _wgt;
  
  mask= 0.0;
  
  //reg= [-250, 250, -250, 250];
  if(is_void(pt)) {
    p= pt= [0.0, 0.0];
    off= [0.0, 0.0, 0.0, 0.0];
  } else  {
    if (struct_type(pt) == "POINT") {
      off= [pt.x, pt.x, pt.y, pt.y];
      p= [pt.x, pt.y];
    } else {
      p= [0.0, 0.0];
      off= [pt(1), pt(1), pt(2), pt(2)];
    }
  }

  dr= [-10, 10, -10, 10]*2;
  dr1= dr/10.0;
  
  if (is_void(method) || method == 1) {
    //ac= img_extract(img_correlate(img1, reg= reg, wndw= "hanning"), off+dr);
    //cc= img_extract(img_correlate(img1, img2, reg= reg, wndw= "hanning"), off+dr);
    //ac= img_extract(img_correlate(img1, reg= reg, wndw= "bartlett"), off+dr);
    //cc= img_extract(img_correlate(img1, img2, reg= reg, wndw= "bartlett"), off+dr);
    ac= img_extract(img_correlate(img1, reg= reg), off+dr);
    cc= img_extract(img_correlate(img1, img2, reg= reg), off+dr);
    _DVV_AC= ac; _DVV_CC= cc;
    //ac= img_extract(img_correlate(img1), [-10, 10, -10, 10]);
    //cc= img_extract(img_correlate(img1, img2), [-10, 10, -10, 10]);
    zr= img_data(cc)(*)(max);
    if (dbg) sh, 60, ac, zr/1.2, 1.2*zr;
    if (dbg) sh, 61, cc, zr/1.2, 1.2*zr;
    img1avg=img_data(img_extract(img1, off+dr1))(*)(avg);
    img2avg=img_data(img_extract(img2, off+dr1))(*)(avg);
    imgavg= (img1avg+img2avg)/2;

    //ix= ac.nx/2+1;
    //iy= ac.ny/2+1;
    //ixv= [-1,0,1,-1,1,-1,0,1];
    //iyv= [1,1,1,0,0,-1,-1,-1];
    ad= img_data(ac);  //ad(ix, iy)= median(ad(ix+ixv, iy+iyv)(*));
    //ac= img_copy(ac, data= ad);
    cd= img_data(cc);  //cd(ix, iy)= median(cd(ix+ixv, iy+iyv)(*));
    //cc= img_copy(cc, data= cd);

    amx= ad(*)(max);
    ija= ij(ad(*)(mxx), ad);
    xa= (*ac.xscale)(ija(1));
    ya= (*ac.yscale)(ija(2));
    if (dbg) print, "ija= ", ija, "[x,y] = ", [xa, ya];
    acx= img_extract(ac, dvv_region(POINT(x=xa,y=ya), 2.5));
    //aa= [amx, 0.0, 1.0, 0.0, 1.0, 0.5, amx/5.0];
    aa= [amx, xa, 2.5, ya, 2.5, amx/2.0];
    fit0= [1,2,4,6];
    if (dbg) print, "aa initial = ", aa;
    //r= lmfit(gaussian2D, [img_grid(ac,1)(*), img_grid(ac,2)(*)], aa, img_data(img_median(ac))(*));
    //r= lmfit(gaussian2D, [img_grid(ac,1)(*), img_grid(ac,2)(*)], aa, img_data(ac)(*));
    _wgt= clip(img_data(acx)(*)/amx, 0., 1.0);
    r= lmfit(gaussian2D, [img_grid(acx,1)(*), img_grid(acx,2)(*)], aa, img_data(acx)(*), _wgt, fit= fit0);
    if (dbg || r.niter >= 100) print, "A1: pt= ", p, "aa final = ", aa;
    //aa= _(aa, 0.1);
    //fit1= [1,2,4,6,7];
    //aa(7)= aa(7)%pi;
    //r= lmfit(gaussian2D, [img_grid(ac,1)(*), img_grid(ac,2)(*)], aa, img_data(img_median(ac))(*));
    //r= lmfit(gaussian2D, [img_grid(ac,1)(*), img_grid(ac,2)(*)], aa, img_data(ac)(*));
    //r= lmfit(gaussian2D, [img_grid(acx,1)(*), img_grid(acx,2)(*)], aa, img_data(acx)(*), fit= fit1);
    //if (dbg || r.niter >= 100) {
    //  print, "A2: pt= ", p, "aa final = ", aa;
    //  mask= 1.0;  //fit is no good set a mask value so we can ignore this one
    // }
    //aa(6)= aa(6)%pi;
    //lmfit, gaussian2D, [img_grid(ac,1)(*), img_grid(ac,2)(*)], aa, img_data(img_median(ac))(*);

    bmx= cd(*)(max);
    ijb= ij(cd(*)(mxx), cd);
    xb= (*cc.xscale)(ijb(1));
    yb= (*cc.yscale)(ijb(2));
    if (dbg) print, "ijb= ", ijb, "[x,y] = ", [xb, yb];
    ccx= img_extract(cc, dvv_region(POINT(x=xb,y=yb), 2.5));
    //bb= [bmx, 0.0, 1.0, 0.0, 1.0, 0.5, bmx/5.0];
    bb= [bmx, xb, 2.5, yb, 2.5, bmx/2.0];
    if (dbg) print, "bb initial = ", bb;
    //r= lmfit(gaussian2D, [img_grid(cc,1)(*), img_grid(cc,2)(*)], bb, img_data(img_median(cc))(*));
    //r= lmfit(gaussian2D, [img_grid(cc,1)(*), img_grid(cc,2)(*)], bb, img_data(cc)(*));
    _wgt= clip(img_data(ccx)(*)/bmx, 0., 1.0);
    r= lmfit(gaussian2D, [img_grid(ccx,1)(*), img_grid(ccx,2)(*)], bb, img_data(ccx)(*), _wgt, fit= fit0);
    if (dbg || r.niter >= 100) print, "B1: pt= ", p, "bb final = ", bb;
    //bb= _(bb, 0.1);
    //bb(7)= bb(7)%pi;
    //r= lmfit(gaussian2D, [img_grid(cc,1)(*), img_grid(cc,2)(*)], bb, img_data(img_median(cc))(*));
    //r= lmfit(gaussian2D, [img_grid(cc,1)(*), img_grid(cc,2)(*)], bb, img_data(cc)(*));
    //r= lmfit(gaussian2D, [img_grid(ccx,1)(*), img_grid(ccx,2)(*)], bb, img_data(ccx)(*), fit= fit1);
    //if (dbg || r.niter >= 100) {
    //  print, "B2: pt= ", p, "bb final = ", bb;
    //  mask= 1.0;  //fit is no good set a mask value so we can ignore this one
    // }

    //bb(6)= bb(6)%pi;
    //lmfit, gaussian2D, [img_grid(cc,1)(*), img_grid(cc,2)(*)], bb, img_data(img_median(cc))(*);
    if (dbg) {
      print, "aa= ", aa;
      print, "bb= ", bb;
    }
    //offset= [_DVV_CORR_OFFSET, _DVV_CORR_OFFSET];
    if (dbg) {
      window, 60; plmk, ya, xa, msize= 0.75, color= "blue", marker= 5, width= 3;
      window, 61; plmk, yb, xb, msize= 0.75, color= "blue", marker= 5, width= 3;
      window, 60; plmk, aa(4), aa(2), msize= 0.5, color= "red", marker= 5, width= 10;
      window, 61; plmk, bb(4), bb(2), msize= 0.5, color= "red", marker= 5, width= 10;
    }
    //offset= [_DVV_CORR_OFFSET, _DVV_CORR_OFFSET];
    //return _(bb([2,4]) - aa([2,4]) - offset, mask, imgavg);
    return _(bb([2,4]) - aa([2,4]), mask, imgavg);
    //return _(bb([2,4]), mask, imgavg);
    
  } else if (method == 2) {
    img1x= img_extract(img1,reg);
    imr1= img_resample(img1x, img1x.nx*3, img1x.ny*3);
    imr2= img_resample(img_extract(img2, reg), img1x.nx*3, img1x.ny*3);
    ac= img_correlate(imr1);
    cc= img_correlate(imr1, imr2);
    zr= img_data(cc)(*)(max);
    img_max, ac, ac_coords;
    img_max, cc, cc_coords;
    if (dbg) {
      print, "cc_coords= ",cc_coords, "ac_coords= ", ac_coords;
      sh, 60, ac, zr/1.2, 1.2*zr, c= 1;
      plmk, ac_coords(2), ac_coords(1), marker= 5, color= "red", msize= 0.5, width= 10;
      sh, 61, cc, zr/1.2, 1.2*zr, c= 1;
      plmk, cc_coords(2), cc_coords(1), marker= 5, color= "red", msize= 0.5, width= 10;
    }
    return _(cc_coords - ac_coords, 0, (img_avg(img1, reg= reg) + img_avg(img2, reg= reg))/2.);
    
  } else if (method == 3) {
    ac= img_correlate(img1, reg= reg);
    cc= img_correlate(img1, img2, reg= reg);
    img_max, ac, ac_coords;
    img_max, cc, cc_coords;
    offac= [ac_coords(1), ac_coords(1), ac_coords(2), ac_coords(2)];
    offcc= [cc_coords(1), cc_coords(1), cc_coords(2), cc_coords(2)];
    dr= [-2, 2, -2, 2];
    _DVV_AC= ac; _DVV_CC= cc;
    zr= img_data(cc)(*)(max);
    acr= img_resample(img_extract(ac, offac+dr), 5);
    ccr= img_resample(img_extract(cc, offcc+dr), 5);
    img_max, acr, ac_coords;
    img_max, ccr, cc_coords;
    if (dbg) {
      print, "cc_coords= ",cc_coords, "ac_coords= ", ac_coords;
      sh, 60, ac, zr/1.2, 1.2*zr, c= 1;
      plmk, ac_coords(2), ac_coords(1), marker= 5, color= "red", msize= 0.5, width= 10;
      sh, 61, cc, zr/1.2, 1.2*zr, c= 1;
      plmk, cc_coords(2), cc_coords(1), marker= 5, color= "red", msize= 0.5, width= 10;
    }
    return _(cc_coords - ac_coords, 0, (img_avg(img1, reg= reg) + img_avg(img2, reg= reg))/2.);
  }
}

func dvv_shiftFrame(&fr, offset)
/* DOCUMENT dvv_shiftFrame, fr, offset

   Shifts an OHRV frame set by a small offset = [delta-x, delta-y].
   The offset is usually obtained using dvv_imgOffset.  This is useful
   for shifting a data image to agree with the reference image.  Owing
   to vibration, probably in the TIM telescope, the images are all at
   slightly shift relative to each other, with the shift being
   predominantly in the y-direction.

   Typical usage would be:
     > rtot= dvv_frameSum(rf);
     > dtot= dvv_frameSum(dt);
     > off= dvv_imgOffset(dtot, rtot);
     > dt_shifted= dvv_shiftFrame(dt, off);
     
   SEE ALSO:
     dvv_imgOffset, dvv_frameSum
 */
{
  xg= img_grid(fr.ch1S,1)+offset(1);
  yg= img_grid(fr.ch1S,2)+offset(2);
  if(!am_subroutine()) {
    print, "dvv_ShiftFrame - fscopy";
    nfr= dvv_FSCopy(fr);
    nfr.ch1S= img_resample(fr.ch1S, xg, yg);
    nfr.ch2S= img_resample(fr.ch2S, xg, yg);
    nfr.ch1P= img_resample(fr.ch1P, xg, yg);
    nfr.ch2P= img_resample(fr.ch2P, xg, yg);
    return nfr;
  } else {
    fr.ch1S= img_resample(fr.ch1S, xg, yg);
    fr.ch2S= img_resample(fr.ch2S, xg, yg);
    fr.ch1P= img_resample(fr.ch1P, xg, yg);
    fr.ch2P= img_resample(fr.ch2P, xg, yg);
  }
}

func dvv_frameOffsets(fr, pt, dbg=, check=, case=)
/* DOCUMENT dvv_frameOffset, fr, pt, dbg=, check=, case=

   Operates on a 4-frame set: For a given coordinate point pt: ( [x,y]
   or POINT(x=,y=) ) compute a set of frame-to-frame registration
   offsets for each frame in a data set.  These offsets can then be
   used to construct a dewarping map that can be used to "dewarp" each
   frame such that the frames are overlapped globally to within a
   fraction of a pixel.  The set of offsets represent small local
   shifts that need to be applied to the coordinates of each frame to
   achieve the registration.

   Returns a 3x4 array of offsets, [[delta-x, delta-y, err, weight], ...] for
   ch1s, ch2s, ch1p, ch2p as measured near the point pt.  The 3rd
   element is non-zero if <dvv_imgOffset> signals an error in determining
   the offset value (i.e. to enable discarding this point in later
   processing).

   This function is intended to be applied to the grid image.

   NOTE:
     This function does not produce sensible results for fringing
     data sets (i.e. the reference and data images frame sets).  It
     should only be used for non-fringing data such as the grid 
     or the flat field frame sets.

   KEYWORDS:
     case=  "svsolve"
            "old_method"
            "sp_separate"
            "ch1s_master"
            
   SEE ALSO:
     dvv_imgOffset, dvv_reportOffsets, dvv_ApplyDewarp
     
 */
{
  reg_size= 30.0;
  if (case == "svsolve") {
    
    // There are 6 comparisons among the 4 images which means we have a
    // system of six equations and 4 unknowns.  There is no unique solution.
    o1= dvv_imgOffset(fr.ch1S, fr.ch2S, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o2= dvv_imgOffset(fr.ch1S, fr.ch2P, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o3= dvv_imgOffset(fr.ch1S, fr.ch1P, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o4= dvv_imgOffset(fr.ch1P, fr.ch2S, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o5= dvv_imgOffset(fr.ch1P, fr.ch2P, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o6= dvv_imgOffset(fr.ch2P, fr.ch2S, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg);

    //Use SVsolve (singular value decomposition) to find the best
    //solution set. Row vector is --> [1s, 1p, 2s, 2p]
    aa=transpose([[1,0,-1,0],
                  [1,0,0,-1],
                  [1,-1,0,0],
                  [0,1,-1,0],
                  [0,1,0,-1],
                  [0,0,-1,1]]);
    b= transpose([o1,o2,o3,o4,o5,o6]);
    oo= SVsolve(aa, b);

    oo(,3)= b(avg,3);
    oo(,4)= b(avg,4);

    if (check) {
      write, format="o1 = [%f, %f], diff= [%f,%f]\n", o1(1), o1(2), o1(1) - (oo(1,1)-oo(3,1)), o1(2) - (oo(1,2)-oo(3,2));
      write, format="o2 = [%f, %f], diff= [%f,%f]\n", o2(1), o2(2), o2(1) - (oo(1,1)-oo(4,1)), o2(2) - (oo(1,2)-oo(4,2));
      write, format="o3 = [%f, %f], diff= [%f,%f]\n", o3(1), o3(2), o3(1) - (oo(1,1)-oo(2,1)), o3(2) - (oo(1,2)-oo(2,2));
      write, format="o4 = [%f, %f], diff= [%f,%f]\n", o4(1), o4(2), o4(1) - (oo(2,1)-oo(3,1)), o4(2) - (oo(2,2)-oo(3,2));
      write, format="o5 = [%f, %f], diff= [%f,%f]\n", o5(1), o5(2), o5(1) - (oo(2,1)-oo(4,1)), o5(2) - (oo(2,2)-oo(4,2));
      write, format="o6 = [%f, %f], diff= [%f,%f]\n", o6(1), o6(2), o6(1) - (oo(4,1)-oo(3,1)), o6(2) - (oo(4,2)-oo(3,2));
    }

    //Returns the data
    return transpose(oo);
    
  } else if (case == "old_method") {
    
    //Method developed in 2012 (incorrect)
    //Algorithm preserved here to compare with newer versions
    // There are 6 comparisons among the 4 images
    o1= dvv_imgOffset(fr.ch1S, fr.ch2S, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o2= dvv_imgOffset(fr.ch1S, fr.ch2P, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o3= dvv_imgOffset(fr.ch1S, fr.ch1P, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o4= dvv_imgOffset(fr.ch1P, fr.ch2S, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o5= dvv_imgOffset(fr.ch1P, fr.ch2P, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o6= dvv_imgOffset(fr.ch2P, fr.ch2S, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg);

    o1s= o1p= o2s= o2p= [0.0, 0.0, 0.0, 0.0];

    //This algebra balances the corrections indicated by the 6
    //measurements
    o1s(1:2)= (o1+o2+o3)(1:2)/3.;
    o1p(1:2)= (-o3+o4+o5)(1:2)/3.;
    o2s(1:2)= (-o1-o4-o6)(1:2)/3.;
    o2p(1:2)= (-o2-o5+o6)(1:2)/3.;

    //Accumulate any possible errors in the 3rd element
    //and signal levels (weights) in the 4th element
    o1s(3:4)= (o1+o2+o3)(3:4)/3.;
    o1p(3:4)= (o1+o2+o3)(3:4)/3.;
    o2s(3:4)= (o1+o2+o3)(3:4)/3.;
    o2p(3:4)= (o1+o2+o3)(3:4)/3.;
  
    return [o1s, o1p, o2s, o2p];

  } else if (case == "sp_separate") {

    //Determine offsets for S- and P- frames independently
    o1= dvv_imgOffset(fr.ch1S, fr.ch2S, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o5= dvv_imgOffset(fr.ch1P, fr.ch2P, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o3= dvv_imgOffset(fr.ch1S, fr.ch1P, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o6= dvv_imgOffset(fr.ch2P, fr.ch2S, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg);

    o1s= o1p= o2s= o2p= [0.0, 0.0, 0.0, 0.0];

    o1s(1:2)= o1(1:2)/2.;  o1s(3:4)= o1(3:4);
    o2s(1:2)= -o1(1:2)/2.; o2s(3:4)= o1(3:4);
    o1p(1:2)= o5(1:2)/2.;  o1p(3:4)= o5(3:4);
    o2p(1:2)= -o5(1:2)/2.; o2p(3:4)= o5(3:4);
   
    return [o1s, o1p, o2s, o2p];

  } else if (case == "ch1s_master") {

    // Updated algebra
    // Use ch1S as a master reference
    o1= dvv_imgOffset(fr.ch1S, fr.ch2S, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o2= dvv_imgOffset(fr.ch1S, fr.ch2P, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 
    o3= dvv_imgOffset(fr.ch1S, fr.ch1P, pt= pt, reg= dvv_region(pt, reg_size), dbg= dbg); 

    o1s= o1p= o2s= o2p= [0.0, 0.0, 0.0, 0.0];

    // Ch1S is left intact, the others are warped to match
    o1s(1:2)= [0.0, 0.0];
    o1p(1:2)= -o3(1:2);
    o2s(1:2)= -o1(1:2);
    o2p(1:2)= -o2(1:2);

    //Accumulate any possible errors in the 3rd element
    //and signal levels (weights) in the 4th element
    o1s(3:4)= (o1+o2+o3)(3:4)/3.;
    o1p(3:4)= (o1+o2+o3)(3:4)/3.;
    o2s(3:4)= (o1+o2+o3)(3:4)/3.;
    o2p(3:4)= (o1+o2+o3)(3:4)/3.;
  
    return [o1s, o1p, o2s, o2p];
  } else {
    error, "No case specified.";
  }
}

func dvv_reportOffsets(fr, fr1, rmsf=, pts=)
/* DOCUMENT dvv_reportOffsets, fr

   Using frame-to-frame correlations report statistics on the
   frame-to-frame registration accuracy as evaluated globally over all
   of the node points of the alignment grid.  Reported offsets are in
   units of microns.  For a properly dewarped set of images the x- and
   y- average offsets (magnitudes) are usually < 0.030 um, and the rms
   variation is < 0.075 um.  Since each pixel is 0.45 um, this
   alignment quality is 1/10 to 2/10 of a pixel.

   Example usage:
   
   > dvv_reportOffsets, gr, rmsf= 3
   *** Warning: LMFIT reached maximum number of iterations (100).
   "B2: pt= "  [200,0]  "bb final = "  [5.31906e+07,200.55,3.04283,0.37276,
   2.50503,1.50516e+08,-0.806243]
             Ch1S          Ch1P           Ch2S           Ch2P
   number:    77            80             74             78
   x - avg: -0.018         0.006          0.002          0.004
   x - rms:  0.059         0.059          0.037          0.038
   y - avg: -0.009         0.011          0.007         -0.003
   y - rms:  0.062         0.076          0.036          0.038


   KEYWORDS:
     rmsf=  sigma factor for rejecting outliers (e.g. the default is rmsf= 2 or 2-sigma)
     
   SEE ALSO:
     dvv_frameOffsets
 */
{
  extern _DVV_CORR_OFFSET;
  _DVV_CORR_OFFSET= (*fr.ch1S.xscale)(dif)(avg)/2.;
  
  if(is_void(rmsf)) rmsf= 2.0;
  if (is_void(pts)) {
    refpts= (*fr.map.warp_refpts);
  } else {
    refpts= pts;
  }
  
  n= numberof(refpts);
  
  if (is_void(fr1)) {
    offsets= array(double, [3,4,4,n]);
    for(i= 1; i<= n; i++) {
      offsets(..,i)= dvv_frameOffsets(fr, refpts(i), case= "svsolve");
    }
  } else {
    offsets= dvv_diffOffsets(fr, fr1, refpts);
  }

  // Remove outliers
  if (rmsf) {
    wch1so= where(abs(offsets(1,1,..) - median(offsets(1,1,..))) < rmsf*offsets(1,1,..)(*)(rms) &
                  abs(offsets(2,1,..) - median(offsets(2,1,..))) < rmsf*offsets(2,1,..)(*)(rms) &
                  offsets(3,1,..) == 0.0);
    wch1po= where(abs(offsets(1,2,..) - median(offsets(1,2,..))) < rmsf*offsets(1,3,..)(*)(rms) &
                  abs(offsets(2,2,..) - median(offsets(2,2,..))) < rmsf*offsets(2,3,..)(*)(rms) &
                  offsets(3,2,..) == 0.0);
    wch2so= where(abs(offsets(1,3,..) - median(offsets(1,3,..))) < rmsf*offsets(1,2,..)(*)(rms) &
                  abs(offsets(2,3,..) - median(offsets(2,3,..))) < rmsf*offsets(2,2,..)(*)(rms) &
                  offsets(3,3,..) == 0.0);
    wch2po= where(abs(offsets(1,4,..) - median(offsets(1,4,..))) < rmsf*offsets(1,4,..)(*)(rms) &
                  abs(offsets(2,4,..) - median(offsets(2,4,..))) < rmsf*offsets(2,4,..)(*)(rms) &
                  offsets(3,4,..) == 0.0);
  } else {
    wch1so= where(offsets(3,1,..) == 0.0);
    wch1po= where(offsets(3,2,..) == 0.0);
    wch2so= where(offsets(3,3,..) == 0.0);
    wch2po= where(offsets(3,4,..) == 0.0);
  }

  write, format= "%s\n", "          Ch1S          Ch1P           Ch2S           Ch2P";
  write, format= "number:    %2d            %2d             %2d             %2d\n",
    numberof(wch1so), numberof(wch2so), numberof(wch1po), numberof(wch2po);
  write, format= "x - avg: %6.3f        %6.3f         %6.3f         %6.3f\n",
    offsets(1,1,wch1so)(*)(avg), offsets(1,2,wch2so)(*)(avg),
    offsets(1,3,wch1po)(*)(avg), offsets(1,4,wch2po)(*)(avg);
  write, format= "x - rms: %6.3f        %6.3f         %6.3f         %6.3f\n",
    offsets(1,1,wch1so)(*)(rms), offsets(1,2,wch2so)(*)(rms),
    offsets(1,3,wch1po)(*)(rms), offsets(1,4,wch2po)(*)(rms);
  write, format= "y - avg: %6.3f        %6.3f         %6.3f         %6.3f\n",
    offsets(2,1,wch1so)(*)(avg), offsets(2,2,wch2so)(*)(avg),
    offsets(2,3,wch1po)(*)(avg), offsets(2,4,wch2po)(*)(avg);
  write, format= "y - rms: %6.3f        %6.3f         %6.3f         %6.3f\n",
    offsets(2,1,wch1so)(*)(rms), offsets(2,2,wch2so)(*)(rms),
    offsets(2,3,wch1po)(*)(rms), offsets(2,4,wch2po)(*)(rms);

  if (n == 81) {
    window, next_window(), legends= 0, style= "square2x2img.gs";
    plsys, 1; plfc, reform(offsets(1,1,),9,9), reform(refpts.y,9,9), reform(refpts.x,9,9),levs= span(-0.4,0.4,11);
    plsys, 2; plfc, reform(offsets(1,3,),9,9), reform(refpts.y,9,9), reform(refpts.x,9,9),levs= span(-0.4,0.4,11);
    plsys, 3; plfc, reform(offsets(1,2,),9,9), reform(refpts.y,9,9), reform(refpts.x,9,9),levs= span(-0.4,0.4,11);
    plsys, 4; plfc, reform(offsets(1,4,),9,9), reform(refpts.y,9,9), reform(refpts.x,9,9),levs= span(-0.4,0.4,11);
    plt, "X-offsets", 0.1, 0.96, tosys= 0, height= 18, font= "helveticaB";
    plsys, 1; plt, "Ch1s", -350, 350, tosys= 1, height= 18, font= "helveticaB", color= "red";
    plsys, 2; plt, "Ch2s", -350, 350, tosys= 1, height= 18, font= "helveticaB", color= "red";
    plsys, 3; plt, "Ch1p", -350, 350, tosys= 1, height= 18, font= "helveticaB", color= "red";
    plsys, 4; plt, "Ch2p", -350, 350, tosys= 1, height= 18, font= "helveticaB", color= "red";
    color_bar1, span(-0.4, 0.4, 81), labs= 20, edges= 0;
    plt, shot+" - frame-to-frame registration", 0.1, 1.01, justify="LC", height= 18, font= "helveticaB", tosys= 0;
    y= 0.195; dy= 0.03;
    plt, swrite(format= "%s", "          Ch1S          Ch1P           Ch2S           Ch2P"), 0.1, y, font= "courier"; y-=dy;
    plt, swrite(format= "number:    %2d            %2d             %2d             %2d",
                numberof(wch1so), numberof(wch2so), numberof(wch1po), numberof(wch2po)), 0.1, y, font= "courier"; y-=dy;
    plt, swrite(format= "x - avg: %6.3f        %6.3f         %6.3f         %6.3f",
                offsets(1,1,wch1so)(*)(avg), offsets(1,2,wch2so)(*)(avg),
                offsets(1,3,wch1po)(*)(avg), offsets(1,4,wch2po)(*)(avg)), 0.1, y, font= "courier"; y-=dy;
    plt, swrite(format= "x - rms: %6.3f        %6.3f         %6.3f         %6.3f",
                offsets(1,1,wch1so)(*)(rms), offsets(1,2,wch2so)(*)(rms),
                offsets(1,3,wch1po)(*)(rms), offsets(1,4,wch2po)(*)(rms)), 0.1, y, font= "courier"; y-=dy;

    window, next_window(), legends= 0, style= "square2x2img.gs";
    plsys, 1; plfc, reform(offsets(2,1,),9,9), reform(refpts.y,9,9), reform(refpts.x,9,9),levs= span(-0.4,0.4,11);
    plsys, 2; plfc, reform(offsets(2,3,),9,9), reform(refpts.y,9,9), reform(refpts.x,9,9),levs= span(-0.4,0.4,11);
    plsys, 3; plfc, reform(offsets(2,2,),9,9), reform(refpts.y,9,9), reform(refpts.x,9,9),levs= span(-0.4,0.4,11);
    plsys, 4; plfc, reform(offsets(2,4,),9,9), reform(refpts.y,9,9), reform(refpts.x,9,9),levs= span(-0.4,0.4,11);
    plt, "Y-offsets", 0.1, 0.96, tosys= 0, height= 18, font= "helveticaB";
    plsys, 1; plt, "Ch1s", -350, 350, tosys= 1, height= 18, font= "helveticaB", color= "red";
    plsys, 2; plt, "Ch2s", -350, 350, tosys= 1, height= 18, font= "helveticaB", color= "red";
    plsys, 3; plt, "Ch1p", -350, 350, tosys= 1, height= 18, font= "helveticaB", color= "red";
    plsys, 4; plt, "Ch2p", -350, 350, tosys= 1, height= 18, font= "helveticaB", color= "red";
    color_bar1, span(-0.4, 0.4, 81), labs= 20, edges= 0;
    plt, shot+" - frame-to-frame registration", 0.1, 1.01, justify="LC", height= 18, font= "helveticaB", tosys= 0;
    y= 0.195; dy= 0.03;
    plt, swrite(format= "%s", "          Ch1S          Ch1P           Ch2S           Ch2P"), 0.1, y, font= "courier"; y-=dy;
    plt, swrite(format= "number:    %2d            %2d             %2d             %2d",
                numberof(wch1so), numberof(wch2so), numberof(wch1po), numberof(wch2po)), 0.1, y, font= "courier"; y-=dy;
    plt, swrite(format= "y - avg: %6.3f        %6.3f         %6.3f         %6.3f",
                offsets(2,1,wch1so)(*)(avg), offsets(2,2,wch2so)(*)(avg), 
                offsets(2,3,wch1po)(*)(avg), offsets(2,4,wch2po)(*)(avg)), 0.1, y, font= "courier"; y-=dy;
    plt, swrite(format= "y - rms: %6.3f        %6.3f         %6.3f         %6.3f",
                offsets(2,1,wch1so)(*)(rms), offsets(2,2,wch2so)(*)(rms),
                offsets(2,3,wch1po)(*)(rms), offsets(2,4,wch2po)(*)(rms)), 0.1, y, font= "courier"; y-=dy;
  }
}

func dvv_diffOffsets(fs1, fs2, pts)
{
  if(struct_type(pts) == "POINT") pts= transpose([pts.x, pts.y],0);
  rk= 1+dimsof(pts);
  pdim= dimsof(pts)(3:);
  np= numberof(pts)/2;
  pts= reform(pts, [2,2,np]);
  of1= of2= array(0.0, [3,4,4,np]);
  for(j= 1; j<= np; j++) of1(..,j)= dvv_frameOffsets(fs1, pts(,j));
  for(j= 1; j<= np; j++) of2(..,j)= dvv_frameOffsets(fs2, pts(,j));
  fdim= _([3,4], pdim);
  rk= numberof(fdim);
  return reform(of1 - of2, _(rk, fdim));
}

func dvv_frameSum(fr, &stot, &ptot, average=)
/* DOCUMENT dvv_frameSum, fr, stot, ptot, average=

   Returns the summation of all 4 frames in a data set.  Also returns
   the S- and P-channel sums in the 2nd and 3rd arguments.
     
   SEE ALSO:
 */
{
  stot= img_add(fr.ch1S, fr.ch2S, average= average);
  ptot= img_add(fr.ch1P, fr.ch2P, average= average);
  if (!am_subroutine()) return img_add(stot, ptot, average= average);
}

func dvv_region(pt, delta)
/* DOCUMENT dvv_region, pt, delta

   Returns a 4-element vector defining a local region centered on
   <pt>: ( [x,y] pair or POINT(x=,y=) ) extending distance +/- delta
   in x- and y-extent from the point, i.e. [x-delta, x+delta, y-delta,
   y+delta].
     
   SEE ALSO:
     dvv_imgOffset
 */
{
  if(struct_type(pt) == "POINT") {
    return [pt.x-delta, pt.x+delta, pt.y-delta, pt.y+delta]; 
  } else if (struct_type(pt) == "double" || struct_type(pt) == "long") {
    return [pt(1)-delta, pt(1)+delta, pt(2)-delta, pt(2)+delta];
  }
}

func dvv_tpsDewarp(img, mapping, set_mapping, dy=, dx=, dbg=, raw= )
/* DOCUMENT dvv_tpsDewarp, img, mapping, dy=, dx=, dbg=, raw=

   Applies a thin plate spline dewarping transformation to a raw data
   image.  This transformation requires the thin plate spline
   coefficients be evaluated prior to the dewarp using dvv_tpsGen.

   KEYWORDS:
     dy=
     dx=
     dbg=
     raw=
     
   SEE ALSO:
     dvv_tpsGen, dvv_tpsEval, dvv_tpsDeltax, dvv_tpsDeltay
 */
{
  if (dbg) print, "tpsDewarp: dbg= ", dbg, ", generate frame set, pc=",_pcounter," ... ";
  pc= _pcounter;
  tic, _pcounter++;
  //restore, tps, aX, aY, xc, yc;

  xc= *set_mapping.xc;
  yc= *set_mapping.yc;
  aX= *mapping.aX;
  aY= *mapping.aY;
  
  xg= is_void(dx) ? img_grid(img,1) : img_grid(img,1) + dx; 
  yg= is_void(dy) ? img_grid(img,2) : img_grid(img,2) + dy;

  xx= xg - dvv_tpsEval(yg, xg, aX, yc, xc);
  yy= yg - dvv_tpsEval(yg, xg, aY, yc, xc);

  // -- using bicubic interpolation
  if(!raw) {
    if (dbg) print, "Resampling from mapping frame ...";
    dw= img_resample(img, xx, yy, meth= 2);
  } else {
    if (dbg) print, "Dewarping with reverse map ...";
    dvv_ReverseMap, mapping, yy, xx;
    dw= img_resample(img_extract(raw, dvv_xtrFrame(raw, *mapping.gridList)),
                  xx, yy, *img.xscale, *img.yscale, meth= 2);
    dw.shotid= img.shotid;
    dw.x_label= img.x_label;
    dw.y_label= img.y_label;
    dw.z_label= img.z_label;
    dw.x_unit= img.x_unit;
    dw.y_unit= img.y_unit;
    dw.z_unit= img.z_unit;
    dw.zrange= img.zrange;
  }
  
  if (dbg) print, "tpsDewarp ... done", tac(pc);
  _pcounter-= 1;
  return dw;
}

func dvv_tpsGen(xc, yc, xd, yd, &aX, &aY)
/* DOCUMENT dvv_tpsGen, xc, yc, xd, yd, &aX, &aY

   Returns the coefficients of a thin plate spline dewarping
   transformation that can be applied to an input image.  The thin
   plate spline dewarping method is based on the following references:

   [1] I. Barrodale et al., Pattern Recognition, vol. 26 pp 375-376 (1993)
   [2] C.A. Zala et al., Advances in Computational Mathematics,
       vol. 11 pp 211-227 (1999)

   The input parameters are a list of control points (xc, yc) and a
   list of associated "data" coordinates (xd, yd).  The data
   coordinates are determined from a calibration image, and the
   control points (xc, yc) are the corresponding coordinates of the
   data points in a distortion-corrected image.  The transformation
   maps points (xd, yd) -> (xc, yc) over the entire image using a thin
   plate spline interpolation method.  The spline is fit to the deltas
   dx= (xc - xd) and dy = (yc - yd).

   The returned data are:
      aX - coefficients of the spline fit for the x-coordinate deviations
      aY - coefficients of the spline fit for the Y-coordinate deviations

   The fitting algorithm dictates that for N control points, there are
   N+3 spline fitting coefficients, so aX and aY always have 3 more
   elements than xc and yc.
     
   SEE ALSO:
     dvv_tpsEval, r2logr2
   
 */
{
  N= numberof(xd);

  //TPS matrix
  xx= xc(-,) - xc(,-);
  yy= yc(-,) - yc(,-);
  rij2= xx*xx + yy*yy;
  rMat= array(0.0, N, N);
  //offdiag= where(indgen(N)(-,) != indgen(N)(,-));
  offdiag= where(rij2 != 0.0);
  rMat(offdiag)= rij2(offdiag)*log(rij2(offdiag));
  tMat= array(0.0, N+3, N+3);
  tMat(4:,4:)= rMat;
  tMat(1,4:)= tMat(4:,1)= array(1.0, N);
  tMat(2,4:)= tMat(4:,2)= xc;
  tMat(3,4:)= tMat(4:,3)= yc;

  //X coefficients
  bX= _([0.,0.,0.], xc - xd);
  aX= SVsolve(tMat, bX);

  //Y coefficients
  bY= _([0.,0.,0.], yc - yd);
  aY= SVsolve(tMat, bY);

  //return save("xc", xc, "yc", yc, aX, aY);
}

func dvv_tpsEval(y, x, coeff, yc, xc)
/* DOCUMENT dvv_tpsEval, y, x, coeff, yc, xc

   Evaluate a thin plate spline at the points (x,y) where the spline
   is defined by the coefficient list <coeff> and the list of control
   points xc, yc.

   The thin plate spline dewarping method is based on the following
   references:

   [1] I. Barrodale et al., Pattern Recognition, vol. 26 pp 375-376 (1993)
   [2] C.A. Zala et al., Advances in Computational Mathematics,
       vol. 11 pp 211-227 (1999)

   SEE ALSO:
     dvv_tpsGen, r2logr2
 */
{
  extern _tps_meth;
  if (is_void(_tps_meth)) _tps_meth= 2;

  z= array(0.0, dimsof(y));
  z= coeff(1) + coeff(2)*x + coeff(3)*y;
  if (_tps_meth == 1) {
    for (i= 1; i<= numberof(yc); i++) z+= coeff(i+3)*r2logr2(y, x, yc(i), xc(i));
  } else if (_tps_meth == 2) {
    z+= r2logr2(y, x, yc, xc, coeff(4:));
  } else error, "Invalid setting for _tps_meth";
  return z;
}

func dvv_framesetOffsets(file, binning=, ssdewarp=, dbg=)
/* DOCUMENT dvv_framesetOffsets, file, binning=, ssdewarp=, dbg=

   Reports a set of image-to-image offset statistics to try to help
   determine movement or offsets of the image scene from one image to
   the next.  The offsets are computed for relative movements between
   the reference, flatfield, data and grid images.  The relative
   movement can be up to 5 or 6 microns and is probably caused by
   vibrations in the TIM.

   Example usage:
   
    > dvv_framesetOffsets, "s65159.hdf", binning= 2, ssdewarp= 1
    Offset data --> ref:   [-0.0496664, 2.44614]
    Offset grid --> ref:   [-0.0882079, 5.18211]
    Offset ff   --> ref:   [0.085693, -0.430211]
    Offset data --> ff:    [0.170907, 2.88854]
    Offset data --> grid:  [0.148394, -2.58164]
    Offset ff   --> grid:  [0.249909, -5.65182]

     
   SEE ALSO:
     dvv_imgOffset
 */
{
  extern ff, rf, dt, gr;
  
  ff= dvv_ImportDataSet(file, set= "f", ssdewarp= ssdewarp, binning= binning, dbg= dbg);
  rf= dvv_ImportDataSet(file, set= "r", ssdewarp= ssdewarp, binning= binning, dbg= dbg);
  dt= dvv_ImportDataSet(file, set= "d", ssdewarp= ssdewarp, binning= binning, dbg= dbg);
  gr= dvv_ImportDataSet(file, set= "g", ssdewarp= ssdewarp, binning= binning, dbg= dbg);

  rtot= dvv_frameSum(rf);
  ftot= dvv_frameSum(ff);
  gtot= dvv_frameSum(gr);
  dtot= dvv_frameSum(dt);

  offdr= dvv_imgOffset(dtot, rtot);
  offgr= dvv_imgOffset(gtot, rtot);
  offfr= dvv_imgOffset(ftot, rtot);
  offdf= dvv_imgOffset(dtot, ftot);
  offdg= dvv_imgOffset(dtot, gtot);
  offfg= dvv_imgOffset(ftot, gtot);
  
  write, format= "Offset data --> ref:   [%g, %g]\n", offdr(1), offdr(2);
  write, format= "Offset grid --> ref:   [%g, %g]\n", offgr(1), offgr(2);
  write, format= "Offset ff   --> ref:   [%g, %g]\n", offfr(1), offfr(2);
  write, format= "Offset data --> ff:    [%g, %g]\n", offdf(1), offdf(2);
  write, format= "Offset data --> grid:  [%g, %g]\n", offdg(1), offdg(2);
  write, format= "Offset ff   --> grid:  [%g, %g]\n", offfg(1), offfg(2);
}

func dvv_plotOffsets(offs, vec, order=, type=, marker=, rge=, w=, errscl=, resid=)
/* DOCUMENT dvv_plotOffsets, offs, vec, order=, type=, marker=, rge=, w=, errscl=, resid=
     
   SEE ALSO:
     dvv_analyzeOffsets
 */
{
  if(is_void(errscl)) scl= 1.0; else scl= errscl;
  if(is_void(order)) order= 2;
  if (is_void(stype)) stype= 1;
  if(is_void(rge)) rge= [-2,2];
  if(is_void(w)) {
    window, next_window(), style= "square2x2img.gs";
    get_style, ld, sy, lg, cl;
    sy.ticks.frame= 1;
    sy.ticks.horiz.flags=0x008+0x003;
    sy.ticks.vert.flags=0x008+0x003;
    sy.ticks.vert.nMajor=4.5;
    sy(3).ticks.horiz.flags+=0x020;
    sy(3).ticks.vert.flags+=0x020;
    set_style, ld, sy, lg, cl;
  } else {
    window, w;
  }
  
  plsys, 1;
  gd= where(offs(3,1,) == 0.0);
  //err= scl/sqrt(abs(offs(4,1,gd)));
  err= scl/abs(offs(4,1,gd));
  xfit= poly1_fit(offs(1,1,gd), vec(gd), order, 1./(err*err));
  yfit= poly1_fit(offs(2,1,gd), vec(gd), order, 1./(err*err));
  //xfit1= dvv_fitOffsets(offs(1,1,), vec, err, order);
  //yfit1= dvv_fitOffsets(offs(2,1,), vec, err, order);
  if (resid) {
    pleb, poly1(vec(gd), xfit) - offs(1,1,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "red";
    pleb, poly1(vec(gd), yfit) - offs(2,1,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "blue";
  } else {
    pleb, offs(1,1,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "red";
    pleb, offs(2,1,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "blue";
    plg, poly1(vec, xfit), vec, width= 3, type= type, color= "red";
    plg, poly1(vec, yfit), vec, width= 3, type= type, color= "blue";
    //plg, poly1(vec, xfit1), vec, width= 3, type= type, color= "magenta";
    //plg, poly1(vec, yfit1), vec, width= 3, type= type, color= "cyan";
  }
  range, rge(1), rge(2);

  plsys, 2;
  gd= where(offs(3,3,) == 0.0);
  //err= scl/sqrt(abs(offs(4,3,gd)));
  err= scl/abs(offs(4,3,gd));
  xfit= poly1_fit(offs(1,3,gd), vec(gd), order, 1./(err*err));
  yfit= poly1_fit(offs(2,3,gd), vec(gd), order, 1./(err*err));
  if (resid) {
    pleb, poly1(vec(gd), xfit) - offs(1,3,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "red";
    pleb, poly1(vec(gd), yfit) - offs(2,3,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "blue";
  } else {
    pleb, offs(1,3,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "red";
    pleb, offs(2,3,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "blue";
    plg, poly1(vec, xfit), vec, width= 3, type= type, color= "red";
    plg, poly1(vec, yfit), vec, width= 3, type= type, color= "blue";
  }
  range, rge(1), rge(2);

  plsys, 3;
  gd= where(offs(3,2,) == 0.0);
  //err= scl/sqrt(abs(offs(4,2,gd)));
  err= scl/abs(offs(4,2,gd));
  xfit= poly1_fit(offs(1,2,gd), vec(gd), order, 1./(err*err));
  yfit= poly1_fit(offs(2,2,gd), vec(gd), order, 1./(err*err));
  if (resid) {
    pleb,  poly1(vec(gd), xfit) - offs(1,2,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "red";
    pleb,  poly1(vec(gd), yfit) - offs(2,2,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "blue";
  } else {
    pleb, offs(1,2,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "red";
    pleb, offs(2,2,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "blue";
    plg, poly1(vec, xfit), vec, width= 3, type= type, color= "red";
    plg, poly1(vec, yfit), vec, width= 3, type= type, color= "blue";
  }
  range, rge(1), rge(2);

  plsys, 4;
  gd= where(offs(3,4,) == 0.0);
  //err= scl/sqrt(abs(offs(4,4,gd)));
  err= scl/abs(offs(4,4,gd));
  xfit= poly1_fit(offs(1,4,gd), vec(gd), order, 1./(err*err));
  yfit= poly1_fit(offs(2,4,gd), vec(gd), order, 1./(err*err));
  if (resid) {
    pleb,  poly1(vec(gd), xfit) - offs(1,4,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "red";
    pleb,  poly1(vec(gd), yfit) - offs(2,4,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "blue";
  } else {
    pleb, offs(1,4,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "red";
    pleb, offs(2,4,gd), vec(gd), dy= err, msize= 0.3, width= 3, marker= marker, color= "blue";
    plg, poly1(vec, xfit), vec, width= 3, type= type, color= "red";
    plg, poly1(vec, yfit), vec, width= 3, type= type, color= "blue";
  }
  range, rge(1), rge(2);
}

func dvv_analyzeOffsets(file, n=, resid=, rge=, errscl=)
/* DOCUMENT dvv_analyzeOffsets, file, n=, resid=, rge=, errscl=

   Imports the flat field and grid reference images and determines
   the frame offsets along the central x- and y-axes of the image.
     
   SEE ALSO:
 */
{
  if (is_void(n)) n= 50;
  gr= dvv_ImportDataSet(file, set= "g", noff= 1, dwmethod= 1, ssdewarp= 1);
  ff= dvv_ImportDataSet(file, set= "f", noff= 1, dwmethod= 1, ssdewarp= 1);
  vecx= transpose([span(-400,400,n), array(0.0,n)]);
  vecy= vecx([2,1],);
  ffox= ffoy= grox= groy= array(0.0, [3,4,4,n]);
  for(j= 1; j<= n; j++) ffox(..,j)= dvv_frameOffsets(ff, vecx(,j), case= "svsolve");
  for(j= 1; j<= n; j++) ffoy(..,j)= dvv_frameOffsets(ff, vecy(,j), case= "svsolve");
  for(j= 1; j<= n; j++) grox(..,j)= dvv_frameOffsets(gr, vecx(,j), case= "svsolve");
  for(j= 1; j<= n; j++) groy(..,j)= dvv_frameOffsets(gr, vecy(,j), case= "svsolve");
  
  dvv_plotOffsets, ffox, vecx(1,), resid= resid, rge= rge, errscl= errscl;
  dvv_plotOffsets, grox, vecx(1,), w= current_window(), type= 2, resid= resid, rge= rge, errscl= errscl; 
  plt, strtr(file,'_','-')+": offsets @ y=0", 0.41, 0.97, tosys= 0, justify="CC";    
  if (resid) {
    wpdf, current_window(), strpart(file,1:-4)+":offsets@y=0,resid";
  } else {
    wpdf, current_window(), strpart(file,1:-4)+":offsets@y=0";
  }

  dvv_plotOffsets, ffoy, vecy(2,), resid= resid, rge= rge, errscl= errscl;
  dvv_plotOffsets, groy, vecy(2,), w= current_window(), type= 2, resid= resid, rge= rge, errscl= errscl;
  plt, strtr(file,'_','-')+": offsets @ x=0", 0.41, 0.97, tosys= 0, justify="CC";
  if (resid) {
    wpdf, current_window(), strpart(file,1:-4)+":offsets@x=0,resid";
  } else {
    wpdf, current_window(), strpart(file,1:-4)+":offsets@x=0";
  }

}

func dvv_fitOffsets(y,x,err,order,thresh=)
/* DOCUMENT dvv_fitOffsets, x, y, err, order, thresh=

   Performs an iterative fit to an offset data set and attempts to
   eliminate outliers from the fit in order to improve the quality of
   the fit.
     
   SEE ALSO:
 */
{
  if (is_void(thresh)) thresh= 2.0;
  if (is_void(order)) order= 2;
   v= indgen(numberof(y));
   ff= 1; n= 0;
   fit= array(0.0, order+1);
   do {
     ofit= fit;
     fit= poly1_fit(y(v), x(v), order, 1./(err(v)*err(v)));
     residv= y(v) - poly1(x(v), fit);
     residmed= median(residv);
     out= abs(residv - residmed) > exp(ff)*thresh*err(v);
     wo= where(out);
     residt= y - poly1(x, fit);
     residtmed= median(residt);
     v= where(abs(residt - residtmed) <= exp(ff)*thresh*err);
     ff/=1.5; n++;
     print, "n= ", n, "nv= ", numberof(v), "out = ", wo;
   } while(is_array(wo) && !allof(fit == ofit));
   print, "Threshold = ", exp(ff)*thresh;
   out= abs(residt) > thresh*err;
   return fit;
}

func dvv_mapQuality(f1, ..)
/* DOCUMENT dvv_mapQuality, f1

   Evaluates and reports on the quality of the frame-to-frame
   registration (the frame map) by examining the rms deviation of the
   ratio of pairs of the frame set normalized by the total frame
   intensity.  Only makes sense for an imported flat field or imported
   grid frame set.

   Exxample usage:
     ff= dvv_ImportDataSet("sxxxxx.hdf", set= "f", dwmethod= 1, ssdewarp= 1);
     dvv_mapQuality, ff;

   SEE ALSO:
     dvv_corrScan, dvv_corrScanSet
 */
{
  f1tot= dvv_frameSum(f1, average= 1);
  c1s2s= img_rms(img_div(img_sub(f1.ch1S, f1.ch2S), f1tot), reg= [-250,250,250,250]);
  c1s1p= img_rms(img_div(img_sub(f1.ch1S, f1.ch1P), f1tot), reg= [-250,250,250,250]);
  c1s2p= img_rms(img_div(img_sub(f1.ch1S, f1.ch2P), f1tot), reg= [-250,250,250,250]);
  c2s1p= img_rms(img_div(img_sub(f1.ch2S, f1.ch1P), f1tot), reg= [-250,250,250,250]);
  c2s2p= img_rms(img_div(img_sub(f1.ch2S, f1.ch2P), f1tot), reg= [-250,250,250,250]);
  c1p2p= img_rms(img_div(img_sub(f1.ch1P, f1.ch2P), f1tot), reg= [-250,250,250,250]);

  write, format="%s\n%s\n", "Case    rms", "=============";
  write, format="1S->2S   %5.3f\n", c1s2s;
  write, format="1S->1P   %5.3f\n", c1s1p;
  write, format="1S->2P   %5.3f\n", c1s2p;
  write, format="2S->1P   %5.3f\n", c2s1p;
  write, format="2S->2P   %5.3f\n", c2s2p;
  write, format="1P->2P   %5.3f\n", c1p2p;
}

func dvv_corrScan(img1, img2, or=, w=, comment=, method=)
/* DOCUMENT dvv_corrScan, img1, img2, or=, w=, comment=, method=

   Performs a fine-grained scan of two images that are nominally the
   same; e.g. the flat field images of ch1S and ch2S of a frameset.
   Plots the results of 5 scans within the center +/- 100 microns
   of the field either
     
   SEE ALSO:
 */
{
  if (is_void(comment)) comment= "";
  if(is_void(w)) w= next_window();
  window, w;
  plsplit, 1, 2, dpi= 100;
  pts= array(POINT(), 400)(-:1:5,);
  if(or =="h" || or == 1) {
    pts.x= span(-400,400,400)(-:1:5,);
    pts.y= span(-100,100,5)(,-:1:400);
  } else {
    pts.y= span(-400,400,400)(-:1:5,);
    pts.x= span(-100,100,5)(,-:1:400);
  }
  op= array(double,4,5,400);
  for(j= 1; j<= 5; j++) {
    for(k= 1; k<= 400; k++) {
      op(,j,k)= dvv_imgOffset(img1, img2, pt= pts(j,k), reg= dvv_region(pts(j,k), 40), method= method);
    }
  }
  cst= 45;
  ol= olegend(new,);
  if (or == "h" || or == 1) {
    for(j= 1; j<= 5; j++) {
      wg= where(op(3,j,) == 0);
      plsys, 1;
      plg, op(1,j,wg), pts(j,wg).x, width= 3,  color=  _rgb(cst+15*j), marks= 0;
      plsys, 2;
      plg, op(2,j,wg), pts(j,wg).x, type= 1, width= 3,  color=  _rgb(cst+15*j), marks= 0;
      ol, add, swrite(format="y=%3.0f", pts(j,1).y), clr= _rgb(cst+15*j), type= 1;
    }
    w300= where(op(3,..) == 0 & pts.x >= -300 & pts.x <= 300);
    avgrms1= swrite(format="average,rms [-300,300]: %5.3f,%5.3f",op(1,w300)(*)(avg),op(1,w300)(*)(rms));
    avgrms2= swrite(format="average,rms [-300,300]: %5.3f,%5.3f",op(2,w300)(*)(avg),op(2,w300)(*)(rms));
    plsys, 1;
    xytitles, "Position (!mm)", "X-Offset (!mm)";
    plsys, 2;
    xytitles, "Position (!mm)", "Y-Offset (!mm)";
    plt, "Horizontal correlation scan\n"+comment, 0.2, 0.90, tosys= 0, height=12, font="helveticaB";
    plt, avgrms1, 0.42, 0.38, tosys= 0, height= 10, font="helvetica";             
    plt, avgrms2, 0.42, 0.36, tosys= 0, height= 10, font="helvetica";             

  } else {
    for(j= 1; j<= 5; j++) {
      wg= where(op(3,j,) == 0);
      plsys, 1;
      plg, op(1,j,wg), pts(j,wg).y, width= 3, color=  _rgb(cst+15*j), marks= 0;
      plsys, 2;
      plg, op(2,j,wg), pts(j,wg).y, type= 1, width= 3, color=  _rgb(cst+15*j), marks= 0;
      ol, add, swrite(format="x=%3.0f", pts(j,1).x), clr= _rgb(cst+15*j), type= 1;
    }
    w300= where(op(3,..) == 0 & pts.y >= -300 & pts.y <= 300);
    avgrms1= swrite(format="average,rms [-300,300]: %5.3f,%5.3f",op(1,w300)(*)(avg),op(1,w300)(*)(rms));
    avgrms2= swrite(format="average,rms [-300,300]: %5.3f,%5.3f",op(2,w300)(*)(avg),op(2,w300)(*)(rms));
    plsys, 1;
    xytitles, "Position (!mm)", "X-Offset (!mm)";
    plsys, 2;
    xytitles, "Position (!mm)", "Y-Offset (!mm)";
    plt, "Vertical correlation scan\n"+comment, 0.2, 0.90, tosys= 0, height=12, font="helveticaB";
    plt, avgrms1, 0.42, 0.38, tosys= 0, height= 10, font="helvetica";             
    plt, avgrms2, 0.42, 0.36, tosys= 0, height= 10, font="helvetica";             
  }
  ol, set, 0.20, 0.38;
  ol, plot;
  plsys, 1; range, -0.25, 0.25; plsys,2; range, -0.25, 0.25;
  return save(pts= pts, op= op, or= or);
}

func dvv_corrScanSet(ff, comment=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  csg= dvv_corrScan(ff.ch1S, ff.ch2S, or= 1, comment= "ch1S,ch2S,"+comment, method= 1, w= next_window());
  csg= dvv_corrScan(ff.ch1S, ff.ch2S, or= 2, comment= "ch1S,ch2S,"+comment, method= 1, w= next_window());
  csg= dvv_corrScan(ff.ch1P, ff.ch2P, or= 1, comment= "ch1P,ch2P,"+comment, method= 1, w= next_window());
  csg= dvv_corrScan(ff.ch1P, ff.ch2P, or= 2, comment= "ch1P,ch2P,"+comment, method= 1, w= next_window());
}
