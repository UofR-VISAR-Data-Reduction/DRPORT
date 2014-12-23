/*
  DVVFILE.I
  
  Definitions for I/O related routines for processing dVV/OHRV
  data sets

  P.M. Celliers
  LLNL
  December, 2007

*/

/*
  
  dvv_ImportDataSet
  dvv_ImportBkg
  dvv_SaveBkg
  dvv_RestoreBkg
  
  dvv_SaveReferenceData
  dvv_RestoreReferenceData
  dvv_SavePhaseData
  dvv_RestorePhaseData
  dvv_SavePAmaps
  dvv_RestorePAmaps
  dvv_SaveFrameSet
  dvv_RestoreFrameSet
  
  dvv_Marinak
  dvv_Report
  dvv_writeMap
  dvv_AssembleSet
  dvv_PAmap_average
  dvv_exportData
  
 */

func dvv_ImportDataSet(file, &img, set=, noff=, ffmethod=, dwmethod=, nobkgsub=, psf=, ssdewarp=, binning=,
                       filter=, median_filter=, dbg= )
/* DOCUMENT fs= dvv_ImportDataSet(file, &img, set= noff=, ffmethod=, dwmethod=, nobkgsub=, psf=, ssdewarp=,
                       binning=, filter=, median_filter=, dbg= )

   Returns a frameset consisting of the 4 frames contained in the
   raw data record contained in <file>.

   KEYWORDS:
     set=  if undefined (default) import the on-shot data image
          -or- specify "d", "D" or "dat" for the on-shot data image
                       "g", "G" or "grid" for the grid reference
                       "r", "R" or "ref" for the background fringe reference
                       "f", "F" or "ff" for the flatfield reference
                       "b", "B" or "bkg" for the background reference
     noff= set this to omit flat field correction
     ffmethod= leave undefined or 0 for the original 4-frame flat field method
               1 apply flat field before frame separation
     nobkgsub=
     dwmethod= leave undefined or zero to omit detailed dewarping
               1 - use original polywarp method
               2 - use bicubic dewarp
               3 - use fine scale TPS dewarp
     psf=      apply filter
     ssdewarp= use single step dewarp (if not specified assumed as default)
     binning=  set binning
     filter=
     median_filter= apply a median filter to raw data upon import
     
   SEE ALSO:
     dvv_ImportBkg, dvv_ImportFlatField
 */
{
  extern _DVV_BKG, _DVV_REF, _DVV_FF, _DVV_MAGFACTOR;
  extern _DVV_FF_IMG;
  extern _DVV_MAPPING;

  if(is_void(set)) {plane= 1; setname= "Streak_array";}
  else if (set == "g" || set == "G" || set == "grid") {plane= 2; setname= "References";}
  else if (set == "r" || set == "R" || set == "ref") {plane= 3; setname= "References";}
  else if (set == "f" || set == "F" || set == "ff") {plane= 1; setname= "References";}
  else if (set == "b" || set == "B" || set == "bkg") {plane= 2; setname= "Streak_array";}
  else if (set == "d" || set == "D" || set == "dat") {plane= 1; setname= "Streak_array";}
  else error, "Invalid specification for set keyword.";

  if (is_void(nobkgsub) && is_void(_DVV_BKG) && (set != "b" && set != "B")) {
    print, "(dvv_ImportDataSet) WARNING: using internal background subtration";
    img= img_sub(hdf_SDget(file, setname= setname, plane= plane, setdouble= 1),
                 hdf_SDget(file, setname= "Streak_array", plane= 2, setdouble= 1));

  } else if (is_void(nobkgsub) && !is_void(_DVV_BKG)) {
    img= img_sub(hdf_SDget(file, setname= setname, plane= plane, setdouble= 1), _DVV_BKG);

  } else {
    print, "(dvv_ImportDataSet) WARNING: not doing background subtration";
    img= hdf_SDget(file, setname= setname, plane= plane, setdouble= 1);
  }

  if (!is_void(median_filter)) {
    print, "(dvv_ImportDataSet) Applying median filter on import";
    //Deal with the 4 frames separately because the median filter function chokes on large
    //data sets (4k x 4k is too big)
    nx= img.nx; ny= img.ny
    imx1= img_median(img_extract(img, xrange=1:nx/2, yrange=1:ny/2), median_filter);
    imx2= img_median(img_extract(img, xrange=nx/2+1:nx, yrange=1:ny/2), median_filter);
    imx3= img_median(img_extract(img, xrange=nx/2+1:nx, yrange=ny/2+1:ny), median_filter);
    imx4= img_median(img_extract(img, xrange=1:nx/2, yrange=ny/2+1:ny), median_filter);
    
    img= img_insert(img,, img_data(imx1), xrange=1:nx/2, yrange=1:ny/2);
    img= img_insert(img,, img_data(imx2), xrange=nx/2+1:nx, yrange=1:ny/2);
    img= img_insert(img,, img_data(imx3), xrange=nx/2+1:nx, yrange=ny/2+1:ny);
    img= img_insert(img,, img_data(imx4), xrange=1:nx/2, yrange=ny/2+1:ny);
  }
  if (dbg) print, "Mapping note: ", _DVV_MAPPING.note;

  if (!am_subroutine()) {
    if (is_void(noff)) {
      //if(is_void(_DVV_REF)) error, "Can't import data, no reference set defined";
      if(is_void(_DVV_MAPPING)) error, "Can't import data, no mapping defined";
      if(is_void(_DVV_FF)) {
        print, "(dvv_ImportDataSet) WARNING: no flat field reference defined, not doing flat fielding";
        //fs= dvv_GenerateFrameSet(img, _DVV_REF, smooth= psf, dwmethod= dwmethod,
        fs= dvv_GenerateFrameSet(img, _DVV_MAPPING, smooth= psf, dwmethod= dwmethod,
                                 ssdewarp= ssdewarp, binning= binning, filter= filter, dbg= dbg);
        //} else fs= dvv_ApplyFlatField(dvv_GenerateFrameSet(img, _DVV_REF, smooth= psf,
      } else {
        if (is_void(ffmethod) || ffmethod == 0) {
          print, "(dvv_ImportDataSet) INFO: applying 4-frame flat field method";
          fs= dvv_ApplyFlatField(
                                 dvv_GenerateFrameSet(img, _DVV_MAPPING, smooth= psf,
                                                      dwmethod= dwmethod, ssdewarp= ssdewarp, binning= binning,
                                                      filter= filter, dbg= dbg));
        } else if (ffmethod == 1) {
          print, "(dvv_ImportDataSet) INFO: applying full image flat field";
          fs= dvv_GenerateFrameSet(img_div(img, _DVV_FF_IMG), _DVV_MAPPING, smooth= psf,
                                   dwmethod= dwmethod, ssdewarp= ssdewarp, binning= binning,
                                   filter= filter, dbg= dbg);
        }
      }
    } else {
      //fs= dvv_GenerateFrameSet(img, _DVV_REF, smooth= psf, dwmethod= dwmethod,
      print, "(dvv_ImportDataSet) WARNING: no flat field reference defined, not doing flat fielding";
      fs= dvv_GenerateFrameSet(img, _DVV_MAPPING, smooth= psf, dwmethod= dwmethod,
                               ssdewarp= ssdewarp, binning= binning, filter= filter, dbg= dbg);
    }
    if (structof(psf) == double) fs.psf= psf;
    return fs;
  }
}

func dvv_ImportBkg(file, ..)
/* DOCUMENT     dvv_ImportBkg, <file1>, <file2>, ...
           -or- dvv_ImportBkg, [<file1>, <file2>, ...]

   Reads background reference files into memory. Multiple background
   files may be listed.  The background is defined as the average of
   the data contained in the file set. The file can be sepcified
   sequentially in the argument list, or supplied in a string array.
     
   SEE ALSO:
     dvv_SaveBkg, dvv_RestoreBkg, dvv_ImportFlatField, dvv_ImportDataSet
 */
{
  extern _DVV_BKG, _DVV_REF, _DVV_FF, _DVV_BKG_FILES;

  if(numberof(file) == 1) {
    //img= getTIFF(file);
    //img= hdf_SDget(file, setname="Streak_array", plane= 2, setdouble= 1);
    _DVV_BKG_FILES= file;
    fstring=  file;
    nn= 1;
    while(more_args()) {
      file= next_arg();
      //img= img_add(img, hdf_SDget(file, setname="Streak_array", plane= 2, setdouble= 1));
      _DVV_BKG_FILES= _(_DVV_BKG_FILES, file);
      fstring+="+"+file;
      nn++;
    }
  } else {
    nn= 1;
    //img= getTIFF(file(1));
    //img= hdf_SDget(file(1), setname="Streak_array", plane= 2, setdouble= 1);
    _DVV_BKG_FILES= file(1);
    fstring=  file(1);
    for (i= 2; i<= numberof(file); i++) {
      if(strlen(file(i)) > 0) {
        //img= img_add(img, hdf_SDget(file(i), setname="Streak_array", plane= 2, setdouble= 1));
        fstring+="+"+file(i);
        nn++;
        _DVV_BKG_FILES= _(_DVV_BKG_FILES, file(i));
      }
    }
  }

  if (dvv_RestoreBkg()) {
    if (!am_subroutine()) return _DVV_BKG;
  } else {
    //nn= numberof(_DVV_BKG_FILES);
    for(j= 1; j<= nn; j++) {
      if (j == 1)
        img= hdf_SDget(_DVV_BKG_FILES(1), setname="Streak_array", plane= 2, setdouble= 1);
      else
        img= img_add(img, hdf_SDget(_DVV_BKG_FILES(j), setname="Streak_array", plane= 2, setdouble= 1));
    }
    img.data= &((*img.data)/nn);
    img.shotid= fstring;
    _DVV_BKG= img_copy(img);
    dvv_SaveBkg;
    if (!am_subroutine()) return img;
  }
}

func dvv_SaveBkg(..)
/* DOCUMENT dvv_SaveBkg

   Checks the current CCD background data set (file list) against the
   contents of the file "../baseline/background.pdb" if that file
   exists.  If the file list differs then the
   "../baseline/background.pdb" file is updated with the current
   background data set.  if "../baseline/background/pdb" does not
   exist, then a new copy is created.
     
   SEE ALSO:
     dvv_RestoreBkg, dvv_ImportBkg
 */
{
  extern _DVV_BKG, _DVV_BKG_FILES;
  if (open("../baseline/background.pdb","rb",1)) {
    f= openb("../baseline/background.pdb");
    vars= *get_vars(f)(1);
    if (is_member(vars, "_DVV_BKG_FILES")) restore, f, "_DVV_BKG_FILES", flist;
    close, f;
    if (numberof(flist) == numberof(_DVV_BKG_FILES) &&
        allof(flist == _DVV_BKG_FILES)) {
      return;
    } 
  }
  f= createb("../baseline/background.pdb");
  _p= POINT();
  save, f, _p;
  if(!is_void(_DVV_BKG_FILES)) save, f, _DVV_BKG_FILES;
  if(!is_void(_DVV_BKG)) save, f, _DVV_BKG= img_floatData(_DVV_BKG);
  close, f;
}

func dvv_RestoreBkg(..)
/* DOCUMENT dvv_RestoreBkg

   Checks the current CCD background data set (file list) against the
   contents of the file "../baseline/background.pdb" if that file
   exists.  If the file list matches the list stored in
   "../baseline/background.pdb" then the background data set is restored
   from that file and this function returns 1. If there is a mismatch
   nothing is imported and this function returns 0.

   Returns 1 - successful restore
           0 - no restore took place
     
   SEE ALSO:
     dvv_SaveBkg, dvv_ImportBkg
     
*/
{
  extern _DVV_BKG, _DVV_BKG_FILES;
  if (open("../baseline/background.pdb","rb",1)) {
    f= openb("../baseline/background.pdb");
    vars= *get_vars(f)(1);
    if (is_member(vars, "_DVV_BKG_FILES")) restore, f, "_DVV_BKG_FILES", flist;
    if (numberof(flist) == numberof(_DVV_BKG_FILES) &&
        allof(flist == _DVV_BKG_FILES)) {
      restore, f, _DVV_BKG;
      close, f;
      return 1;
    } else {
      close, f;
      return 0;
    }
  } else {
    return 0;
  }
}

func dvv_SaveReferenceData(file, comment=)
/* DOCUMENT dvv_SaveReferenceData, file, comment=

   ARGUMENTS:
     file  - file name

   KEYWORDS:
     comment= annotations, entered as a string, or array of strings;

     
   SEE ALSO:
     dvv_RestoreReferenceData
*/
{
  extern _DVV_REF, _DVV_FF, _DVV_BKG, _DVV_COMMENT;
  extern _DVV_FF_IMG;
  extern _DVV_PHASE_MAP, _DVV_AMPL_MAP, _DVV_PMAP, _DVV_AMAP;
  extern _DVV_XOFF_MAP, _DVV_YOFF_MAP;
  extern _DVV_PHASE_DATA, _DVV_AMPL_DATA;
  extern _DVV_XOFF_DATA, _DVV_YOFF_DATA, _DVV_MASK_DATA;
  extern _DVV_MODE, _FRINGE_MODE;
  extern _WAV, _TAU, _ETALON_DELTA, _VPF;
  extern _DVV_FF_MASK, _BINNING;
  extern _DVV_PSF_FF, _DVV_PSF_RF, _DVV_PSF_DT;
  extern _DVV_FFM1, _DVV_FFM2;
  extern _DVV_CH1S_FILTER, _DVV_CH2S_FILTER;
  extern _DVV_CH1P_FILTER, _DVV_CH2P_FILTER;
  extern _DVV_MAPPING;
  extern _DVV_MAGFACTOR;
  extern _DVV_FIELD;
//   extern _DVV_LISS_PARAMS, _DVV_GHOST_PARAMS;
//   extern _DVV_LISS_PARAMS_C, _DVV_LISS_PARAMS_P;

  if(is_void(comment)) cmnt= timestamp()+"\n";
  else cmnt= timestamp()+": "+comment+"\n";

  if(is_void(_DVV_COMMENT)) _DVV_COMMENT= cmnt;
  else _DVV_COMMENT+= cmnt;
    
  _p= POINT();
  _fm= dvv_FrameMapping();
  _sm= dvv_SetMapping();
  
  f= createb(file);
  save, f, _p;
  save, f, _fm;
  save, f, _sm;
  save, f, _DVV_COMMENT;
  if(!is_void(_DVV_BKG_FILES)) save, f, _DVV_BKG_FILES;
  if(!is_void(_DVV_FF_FILES)) save, f, _DVV_FF_FILES;

  if(!is_void(_DVV_PSF_FF)) save, f, _DVV_PSF_FF;
  if(!is_void(_DVV_PSF_RF)) save, f, _DVV_PSF_RF;
  if(!is_void(_DVV_PSF_DT)) save, f, _DVV_PSF_DT;
  
  if(!is_void(_DVV_CH1S_FILTER)) save, f, _DVV_CH1S_FILTER;
  if(!is_void(_DVV_CH2S_FILTER)) save, f, _DVV_CH2S_FILTER;
  if(!is_void(_DVV_CH1P_FILTER)) save, f, _DVV_CH1P_FILTER;
  if(!is_void(_DVV_CH2P_FILTER)) save, f, _DVV_CH2P_FILTER;
  if(!is_void(_DVV_CMMN_FILTER)) save, f, _DVV_CMMN_FILTER;
  
  if(!is_void(_FRINGE_MODE)) save, f, _FRINGE_MODE;
  if(!is_void(_DVV_MODE)) save, f, _DVV_MODE;
  if(!is_void(_WAV)) save, f, _WAV;
  if(!is_void(_TAU)) save, f, _TAU;
  if(!is_void(_ETALON_DELTA)) save, f, _ETALON_DELTA;
  if(!is_void(_VPF)) save, f, _VPF;
  if(!is_void(_BINNING)) save, f, _BINNING;
  //if(!is_void(_DVV_BKG)) save, f, _DVV_BKG= img_floatData(_DVV_BKG);
  //if(!is_void(_DVV_REF)) save, f, _DVV_REF;
  if(!is_void(_DVV_MAPPING)) save, f, _DVV_MAPPING;
  if(!is_void(_DVV_MAGFACTOR)) save, f, _DVV_MAGFACTOR;
  if(!is_void(_DVV_DEWARP_METHOD)) save, f, _DVV_DEWARP_METHOD;
  if(!is_void(_DVV_FIELD)) save, f, _DVV_FIELD;
  // if(!is_void(_DVV_FF)) save, f, _DVV_FF= dvv_FSCopy(_DVV_FF, floatData= 1);
  // if(!is_void(_DVV_FF_IMG)) save, f, _DVV_FF_IMG= img_floatData(_DVV_FF_IMG);
  // if(!is_void(_DVV_PMAP)) save, f, _DVV_PMAP= img_floatData(_DVV_PMAP);
  // if(!is_void(_DVV_AMAP)) save, f, _DVV_AMAP= img_floatData(_DVV_AMAP);
  // if(!is_void(_DVV_PHASE_MAP)) save, f, _DVV_PHASE_MAP= img_floatData(_DVV_PHASE_MAP);
  // if(!is_void(_DVV_AMPL_MAP)) save, f, _DVV_AMPL_MAP= img_floatData(_DVV_AMPL_MAP);
  // if(!is_void(_DVV_XOFF_MAP)) save, f, _DVV_XOFF_MAP= img_floatData(_DVV_XOFF_MAP);
  // if(!is_void(_DVV_YOFF_MAP)) save, f, _DVV_YOFF_MAP= img_floatData(_DVV_YOFF_MAP);
  // if(!is_void(_DVV_PHASE_DATA)) save, f, _DVV_PHASE_DATA= img_floatData(_DVV_PHASE_DATA);
  // if(!is_void(_DVV_AMPL_DATA)) save, f, _DVV_AMPL_DATA= img_floatData(_DVV_AMPL_DATA);
  // if(!is_void(_DVV_XOFF_DATA)) save, f, _DVV_XOFF_DATA= img_floatData(_DVV_XOFF_DATA);
  // if(!is_void(_DVV_YOFF_DATA)) save, f, _DVV_YOFF_DATA= img_floatData(_DVV_YOFF_DATA);
  // if(!is_void(_DVV_MASK_DATA)) save, f, _DVV_MASK_DATA= img_floatData(_DVV_MASK_DATA);
  if(!is_void(_DVV_FF_MASK)) save, f, _DVV_FF_MASK;
  if(!is_void(_DVV_FFM1)) save, f, _DVV_FFM1;
  if(!is_void(_DVV_FFM2)) save, f, _DVV_FFM2;
  // if(!is_void(_DVV_REF_LISS_PARAMS)) save, f, _DVV_REF_LISS_PARAMS;
  // if(!is_void(_DVV_REF_LISS_PARAMS_C)) save, f, _DVV_REF_LISS_PARAMS_C;
  // if(!is_void(_DVV_REF_LISS_PARAMS_P)) save, f, _DVV_REF_LISS_PARAMS_P;
  // if(!is_void(_DVV_DAT_LISS_PARAMS)) save, f, _DVV_DAT_LISS_PARAMS;
  // if(!is_void(_DVV_DAT_LISS_PARAMS_C)) save, f, _DVV_DAT_LISS_PARAMS_C;
  // if(!is_void(_DVV_DAT_LISS_PARAMS_P)) save, f, _DVV_DAT_LISS_PARAMS_P;
//   if(!is_void(_DVV_GHOST_PARAMS)) save, f, _DVV_GHOST_PARAMS;
  close, f;
}

func dvv_RestoreReferenceData(file)
/* DOCUMENT dvv_RestoreReferenceData, file

   ARGUMENTS:
     file  - file name

   KEYWORDS:
     
   SEE ALSO:
     dvv_SaveReferenceData
*/
{
  extern _DVV_MAPPING;
  extern _DVV_MAGFACTOR;
  extern _DVV_FIELD;
  //extern _DVV_REF, _DVV_FF, _DVV_BKG;
  //extern _DVV_FF_IMG;
  //extern _DVV_PMAP, _DVV_AMAP;
  //extern _DVV_PHASE_MAP, _DVV_AMPL_MAP;
  //extern _DVV_XOFF_MAP, _DVV_YOFF_MAP;
  //extern _DVV_PHASE_DATA, _DVV_AMPL_DATA;
  //extern _DVV_XOFF_DATA, _DVV_YOFF_DATA, _DVV_MASK_DATA;
  extern _DVV_COMMENT;
  extern _DVV_MODE, _FRINGE_MODE;
  extern _WAV, _TAU, _ETALON_DELTA, _VPF;
  extern _DVV_FF_MASK, _BINNING;
  extern _DVV_PSF_FF, _DVV_PSF_RF, _DVV_PSF_DT;
  extern _DVV_FFM1, _DVV_FFM2;
  extern _DVV_CH1S_FILTER, _DVV_CH2S_FILTER;
  extern _DVV_CH1P_FILTER, _DVV_CH2P_FILTER;
  // extern _DVV_REF_LISS_PARAMS_C, _DVV_REF_LISS_PARAMS_P;
  // extern _DVV_REF_LISS_PARAMS;
  // extern _DVV_DAT_LISS_PARAMS_C, _DVV_DAT_LISS_PARAMS_P;
  // extern _DVV_DAT_LISS_PARAMS;
  extern _DVV_GHOST_PARAMS;
  extern _DVV_DEWARP_METHOD;

  sdw= _DVV_DEWARP_METHOD;

  print, "dvv_RestoreReferenceData: restoring ", file;
  f= openb(file);
  restore, f;
  if(is_void(_DVV_MAPPING) && !is_void(_DVV_REF)) _DVV_MAPPING= _DVV_REF.map;
  close, f;
  _DVV_DEWARP_METHOD= sdw;
  if(!am_subroutine()) return _DVV_REF;
}

func dvv_SavePhaseData(file, comment=)
/* DOCUMENT dvv_SavePhaseData, file, comment=

   ARGUMENTS:
     file  - file name

   KEYWORDS:
     comment= annotations, entered as a string, or array of strings;

     
   SEE ALSO:
     dvv_RestorePhaseData
*/
{
  extern _DVV_REF, _DVV_FF, _DVV_BKG, _DVV_COMMENT;
  extern _DVV_PHASE_MAP, _DVV_AMPL_MAP, _DVV_PMAP, _DVV_AMAP;
  extern _DVV_XOFF_MAP, _DVV_YOFF_MAP;
  extern _DVV_PHASE_DATA, _DVV_AMPL_DATA;
  extern _DVV_XOFF_DATA, _DVV_YOFF_DATA, _DVV_MASK_DATA;
  extern _DVV_MODE, _FRINGE_MODE;
  extern _DVV_MAGFACTOR, _DVV_FIELD;
  extern _WAV, _TAU, _ETALON_DELTA, _VPF;
  extern _DVV_FF_MASK, _BINNING;
  extern _DVV_REF_LISS_PARAMS_C, _DVV_REF_LISS_PARAMS_P;
  extern _DVV_REF_LISS_PARAMS;
  extern _DVV_DAT_LISS_PARAMS_C, _DVV_DAT_LISS_PARAMS_P;
  extern _DVV_DAT_LISS_PARAMS;
  extern _DVV_GHOST_CORRECTION;
  extern _DVV_DEWARP_METHOD;

  extern dtp, dtpUW, dtpUWs, dtpUWf, ampl, vel;
  extern rfp, rfpUW, rfpUWs;
  extern Sa, Ra, Savg;
  extern dtV, rfV, dtVaz, rfVaz;
  extern MLO, HLO, VLO;

  if(is_void(comment)) cmnt= timestamp()+"\n";
  else cmnt= timestamp()+": "+comment+"\n";

  if(is_void(_DVV_COMMENT)) _DVV_COMMENT= cmnt;
  else _DVV_COMMENT+= cmnt;
    
  _p= POINT();
  _img= IMG_DAT();
  f= createb(file);
  save, f, _p;
  save, f, _img;
  save, f, _DVV_COMMENT;

  if(!is_void(_PSF_FF)) save, f, _PSF_FF;
  if(!is_void(_PSF_RF)) save, f, _PSF_RF;
  if(!is_void(_PSF_DT)) save, f, _PSF_DT;
  
  if(!is_void(_DVV_CH1S_FILTER)) save, f, _DVV_CH1S_FILTER;
  if(!is_void(_DVV_CH2S_FILTER)) save, f, _DVV_CH2S_FILTER;
  if(!is_void(_DVV_CH1P_FILTER)) save, f, _DVV_CH1P_FILTER;
  if(!is_void(_DVV_CH2P_FILTER)) save, f, _DVV_CH2P_FILTER;
  if(!is_void(_DVV_CMMN_FILTER)) save, f, _DVV_CMMN_FILTER;
 
  if(!is_void(_DVV_MODE)) save, f, _DVV_MODE;
  if(!is_void(_FRINGE_MODE)) save, f, _FRINGE_MODE;
  if(!is_void(_DVV_MAGFACTOR)) save, f, _DVV_MAGFACTOR;
  if(!is_void(_DVV_DEWARP_METHOD)) save, f, _DVV_DEWARP_METHOD;
  if(!is_void(_DVV_FIELD)) save, f, _DVV_FIELD;
  if(!is_void(_WAV)) save, f, _WAV;
  if(!is_void(_TAU)) save, f, _TAU;
  if(!is_void(_ETALON_DELTA)) save, f, _ETALON_DELTA;
  if(!is_void(_VPF)) save, f, _VPF;
  if(!is_void(_DVV_PMAP)) save, f, _DVV_PMAP= img_floatData(_DVV_PMAP);
  if(!is_void(_DVV_AMAP)) save, f, _DVV_AMAP= img_floatData(_DVV_AMAP);
  if(!is_void(_DVV_PHASE_MAP)) save, f, _DVV_PHASE_MAP= img_floatData(_DVV_PHASE_MAP);
  if(!is_void(_DVV_AMPL_MAP)) save, f, _DVV_AMPL_MAP= img_floatData(_DVV_AMPL_MAP);
  if(!is_void(_DVV_XOFF_MAP)) save, f, _DVV_XOFF_MAP= img_floatData(_DVV_XOFF_MAP);
  if(!is_void(_DVV_YOFF_MAP)) save, f, _DVV_YOFF_MAP= img_floatData(_DVV_YOFF_MAP);
  if(!is_void(_DVV_PHASE_DATA)) save, f, _DVV_PHASE_DATA= img_floatData(_DVV_PHASE_DATA);
  if(!is_void(_DVV_AMPL_DATA)) save, f, _DVV_AMPL_DATA= img_floatData(_DVV_AMPL_DATA);
  if(!is_void(_DVV_XOFF_DATA)) save, f, _DVV_XOFF_DATA= img_floatData(_DVV_XOFF_DATA);
  if(!is_void(_DVV_YOFF_DATA)) save, f, _DVV_YOFF_DATA= img_floatData(_DVV_YOFF_DATA);
  if(!is_void(_DVV_MASK_DATA)) save, f, _DVV_MASK_DATA= img_floatData(_DVV_MASK_DATA);
  if(!is_void(_DVV_LISS_PARAMS)) save, f, _DVV_LISS_PARAMS;
  if(!is_void(_DVV_LISS_PARAMS_C)) save, f, _DVV_LISS_PARAMS_C;
  if(!is_void(_DVV_LISS_PARAMS_P)) save, f, _DVV_LISS_PARAMS_P;
  if(!is_void(_DVV_GHOST_CORRECTION)) save, f, _DVV_GHOST_CORRECTION;
  if(!is_void(MLO)) save, f, MLO;
  if(!is_void(HLO)) save, f, HLO;
  if(!is_void(VLO)) save, f, VLO;
  if(!is_void(dtp)) save, f, dtp= img_floatData(dtp);
  if(!is_void(dtpUW)) save, f, dtpUW= img_floatData(dtpUW);
  if(!is_void(rfp)) save, f, rfp= img_floatData(rfp);
  if(!is_void(rfpUW)) save, f, rfpUW= img_floatData(rfpUW);
  if(!is_void(rfBkg)) save, f, rfBkg= img_floatData(rfBkg);
  if(!is_void(ampl)) save, f, ampl= img_floatData(ampl);
  if(!is_void(vel)) save, f, vel= img_floatData(vel);
  if(!is_void(dtV)) save, f, dtV= img_floatData(dtV);
  if(!is_void(rfV)) save, f, rfV= img_floatData(rfV);
  if(!is_void(dtVaz)) save, f, dtVaz;
  if(!is_void(rfVaz)) save, f, rfVaz;

  if(!is_void(rfpUWs)) save, f, rfpUWs= img_floatData(rfpUWs);
  if(!is_void(dtpUWs)) save, f, dtpUWs= img_floatData(dtpUWs);
  if(!is_void(dtpUWf)) save, f, dtpUWf= img_floatData(dtpUWf);
  if(!is_void(bkgPhase)) save, f, bkgPhase= img_floatData(bkgPhase);
  if(!is_void(Sa)) save, f, Sa;
  if(!is_void(Ra)) save, f, Ra;
  if(!is_void(Savg)) save, f, Savg;
  
  close, f;
}

func dvv_RestorePhaseData(file)
/* DOCUMENT dvv_RestorePhaseData, file

   ARGUMENTS:
     file  - file name

   KEYWORDS:
     
   SEE ALSO:
     dvv_SavePhaseData
*/
{
  extern _DVV_REF, _DVV_FF, _DVV_BKG;
  extern _DVV_PHASE_MAP, _DVV_AMPL_MAP, _DVV_PMAP, _DVV_AMAP;
  extern _DVV_XOFF_MAP, _DVV_YOFF_MAP;
  extern _DVV_PHASE_DATA, _DVV_AMPL_DATA;
  extern _DVV_XOFF_DATA, _DVV_YOFF_DATA, _DVV_MASK_DATA;
  extern _DVV_COMMENT;
  extern _DVV_MODE, _FRINGE_MODE;
  extern _DVV_FIELD, _DVV_MAGFACTOR;
  extern _WAV, _TAU, _ETALON_DELTA, _VPF;
  extern _DVV_FF_MASK, _BINNING;
  extern _PSF_FF, PSF_DT, _PSF_RF;
  extern _DVV_LISS_PARAMS_C, _DVV_LISS_PARAMS_P;
  extern _DVV_LISS_PARAMS, _DVV_GHOST_PARAMS;
  extern _DVV_DEWARP_METHOD;

  extern dtp, dtpUW, ampl, vel;
  extern rfp, rfpUW;
  extern dtV, rfV, dtVaz, rfVaz;
  extern MLO, HLO, VLO;

  f= openb(file);
  restore, f;
  close, f;
}

func dvv_SavePAmaps(file)
/* DOCUMENT dvv_SavePAmaps, file

   ARGUMENTS:
     file  - file name

   KEYWORDS:
     
   SEE ALSO:
     dvv_SavePhaseData
*/
{
  extern _DVV_PHASE_MAP, _DVV_AMPL_MAP, _DVV_PMAP, _DVV_AMAP;
  extern _DVV_XOFF_MAP, _DVV_YOFF_MAP;
  extern _DVV_PHASE_DATA, _DVV_AMPL_DATA;
  extern _DVV_XOFF_DATA, _DVV_YOFF_DATA, _DVV_MASK_DATA;
  extern _DVV_LISS_PARAMS_C, _DVV_LISS_PARAMS_P;
  extern _DVV_LISS_PARAMS;

  f= createb(file);
  if (!is_void(_DVV_PHASE_MAP)) save, f, _DVV_PHASE_MAP;
  if (!is_void(_DVV_AMPL_MAP)) save, f, _DVV_AMPL_MAP;
  if (!is_void(_DVV_XOFF_MAP)) save, f, _DVV_XOFF_MAP;
  if (!is_void(_DVV_YOFF_MAP)) save, f, _DVV_YOFF_MAP;
  if (!is_void(_DVV_PHASE_DATA)) save, f, _DVV_PHASE_DATA;
  if (!is_void(_DVV_AMPL_DATA)) save, f, _DVV_AMPL_DATA;
  if (!is_void(_DVV_XOFF_DATA)) save, f, _DVV_XOFF_DATA;
  if (!is_void(_DVV_YOFF_DATA)) save, f, _DVV_YOFF_DATA;
  if (!is_void(_DVV_MASK_DATA)) save, f, _DVV_MASK_DATA;
  if (!is_void(_DVV_LISS_PARAMS)) save, f, _DVV_LISS_PARAMS;
  if (!is_void(_DVV_LISS_PARAMS_C)) save, f, _DVV_LISS_PARAMS_C;
  if (!is_void(_DVV_LISS_PARAMS_P)) save, f, _DVV_LISS_PARAMS_P;
  close, f;
}

func dvv_RestorePAmaps(file)
/* DOCUMENT dvv_RestorePAmaps, file

   Restores phase extraction calibration maps from a phase file or a
   PAmaps file.

   ARGUMENTS:
     file  - file name

   KEYWORDS:
     
   SEE ALSO:
     dvv_SavePhaseData
*/
{
  extern _DVV_PHASE_MAP, _DVV_AMPL_MAP, _DVV_PMAP, _DVV_AMAP;
  extern _DVV_XOFF_MAP, _DVV_YOFF_MAP;
  extern _DVV_PHASE_DATA, _DVV_AMPL_DATA;
  extern _DVV_XOFF_DATA, _DVV_YOFF_DATA, _DVV_MASK_DATA;
  extern _DVV_LISS_PARAMS_C, _DVV_LISS_PARAMS_P;
  extern _DVV_LISS_PARAMS;

  f= openb(file);
  vars= *get_vars(f)(1);
  if (is_member(vars, "_DVV_PHASE_MAP")) restore, f, _DVV_PHASE_MAP;
  if (is_member(vars, "_DVV_AMPL_MAP")) restore, f, _DVV_AMPL_MAP;
  if (is_member(vars, "_DVV_XOFF_MAP")) restore, f, _DVV_XOFF_MAP;
  if (is_member(vars, "_DVV_YOFF_MAP")) restore, f, _DVV_YOFF_MAP;
  if (is_member(vars, "_DVV_PHASE_DATA")) restore, f, _DVV_PHASE_DATA;
  if (is_member(vars, "_DVV_AMPL_DATA")) restore, f, _DVV_AMPL_DATA;
  if (is_member(vars, "_DVV_XOFF_DATA")) restore, f, _DVV_XOFF_DATA;
  if (is_member(vars, "_DVV_YOFF_DATA")) restore, f, _DVV_YOFF_DATA;
  if (is_member(vars, "_DVV_MASK_DATA")) restore, f, _DVV_MASK_DATA;
  if (is_member(vars, "_DVV_LISS_PARAMS")) restore, f, _DVV_LISS_PARAMS;
  if (is_member(vars, "_DVV_LISS_PARAMS_C")) restore, f, _DVV_LISS_PARAMS_C;
  if (is_member(vars, "_DVV_LISS_PARAMS_P")) restore, f, _DVV_LISS_PARAMS_P;
  close, f;
}

func dvv_SaveFrameSet(fset, file, comment=)
/* DOCUMENT dvv_SaveFrameSet, fset, file, comment=

   ARGUMENTS:
     fset - frameset name
     file  - file name

   KEYWORDS:
     comment= annotations, entered as a string, or array of strings;

     
   SEE ALSO:
     dvv_RestoreFrameSet
*/
{
  extern _DVV_COMMENT, _DVV_FSET;
  extern _DVV_FIELD, _DVV_MAGFACTOR;
  extern _DVV_DEWARP_METHOD;

  if(is_void(comment)) cmnt= timestamp()+"\n";
  else cmnt= timestamp()+": "+comment+"\n";

  if(is_void(_DVV_COMMENT)) _DVV_COMMENT= cmnt;
  else _DVV_COMMENT+= cmnt;
    
  _p= POINT();
  _fm= dvv_FrameMapping();
  _sm= dvv_SetMapping();
  f= createb(file);
  print, "dvv_SaveFrameSet - fscopy";
  _DVV_FSET= dvv_FSCopy(fset, floatData= 1);
  save, f, _p;
  save, f, _fm;
  save, f, _sm;
  save, f, _DVV_COMMENT;
  save, f, _DVV_FSET;
  if(!is_void(_DVV_MAGFACTOR)) save, f, _DVV_MAGFACTOR;
  if(!is_void(_DVV_DEWARP_METHOD)) save, f, _DVV_DEWARP_METHOD;
  if(!is_void(_DVV_FIELD)) save, f, _DVV_FIELD;
  close, f;
}

func dvv_RestoreFrameSet(file)
/* DOCUMENT dvv_SaveFrameSet, file

   ARGUMENTS:
     fset - string containing the frameset name
     file  - file name

   SEE ALSO:
     dvv_SaveFrameSet
*/
{
  extern _DVV_COMMENT;
  extern _DVV_FIELD, _DVV_MAGFACTOR;
  extern _DVV_DEWARP_METHOD;

  f= openb(file);
  restore, f;
  close, f;

  return _DVV_FSET;
}

func dvv_Marinak(.., tiled=)
/* DOCUMENT dvv_Marinak

   Returns Marty Marinak's Be velocity data set
     
   SEE ALSO:
*/
{
  pth= "~/Experiments/LLE/CAPSEEDS/Calculations/marinak/";
  f= openb(pth+"zdot8395.pdb");
  restore, f;
  close, f;

  mmA= IMG_DAT();
  mmA.nx= mmA.ny= 301;
  mmA.xscale= &x(,1);
  mmA.yscale= &y(1,);
  mmA.x_label= "X position";
  mmA.y_label= "Y position";
  mmA.z_label= "Velocity";
  *mmA.xscale *= 1e4;
  *mmA.yscale *= 1e4;
  mmA.x_unit= "!mm";
  mmA.y_unit= "!mm";
  mmA.z_unit= "m/s";
  mmA.data= &array(0.0, dimsof(zd));
  *mmA.data= 1e4*zd;
  if(tiled) return img_tile(mmA);
  else return mmA;
} 

func dvv_Report(.., intf=, ff=, bf=, chf=, maps=, all=, txtfile=)
/* DOCUMENT dvv_Report intf=, ff=, bf=, chf=, maps=, all=

   Generates a text report of the settings applied in the analysis of
   the current data set.

   KEYWORDS:
   SEE ALSO:
 */
{
  extern _DVV_REF;
  extern _DVV_BKG_FILES, _DVV_FF_FILES;

  if (txtfile && is_void(f)) f= open(shot+"-report.txt", "w");

  if(!is_void(all)) {intf= 1; ff= 1; bf= 1; chf= 1; maps= 1;}
 
  fmt0="%s\n";
  fmt2="  %s  --> (%7.4f, %7.4f)\n"; 
  fmt3="  %s  --> %5.4f\n";
  fmt4="  %s  --> (%7.4f, %7.4f, %8.5f, %8.5f)\n";
  

  //Interferometer setup;
  if (intf) {
    write, f, format= fmt0, "Interferometer parameters";
    write, f, format= fmt0, "=========================";
    write, f, format= fmt3, "WAVELENGTH (um):   ", _WAV;
    write, f, format= fmt3, "TAU (ps):          ", _TAU*1000.0;
    write, f, format= fmt3, "ETALON_DELTA:      ", _ETALON_DELTA;
    write, f, format= fmt3, "VPF (km/s/fringe): ", _VPF;
    write, f, format= fmt2, "FRINGE MODE:       ", _FRINGE_MODE(1), _FRINGE_MODE(2);
    if (!is_void(_DVV_LISS_PARAMS)) write, f, format= fmt4, "LISSAJOUS PARAMS:  ",
     _DVV_LISS_PARAMS(1),_DVV_LISS_PARAMS(2),_DVV_LISS_PARAMS(3),_DVV_LISS_PARAMS(4), ;
    if (!is_void(_DVV_LISS_PARAMS_C)) write, f, format= fmt4, "LISSAJOUS PARAMS C:",
     _DVV_LISS_PARAMS_C(1),_DVV_LISS_PARAMS_C(2),_DVV_LISS_PARAMS_C(3),_DVV_LISS_PARAMS_C(4), ;
    if (!is_void(_DVV_LISS_PARAMS_P)) write, f, format= fmt4, "LISSAJOUS PARAMS P:",
     _DVV_LISS_PARAMS_P(1),_DVV_LISS_PARAMS_P(2),_DVV_LISS_PARAMS_P(3),_DVV_LISS_PARAMS_P(4), ;
    write, f, format= fmt0, "";
  }

  //Flat fielding and background files;
  if (ff) {
    write, f, format= fmt0, "Flat fielding files";
    write, f, format= fmt0, "===================";
    for(i= 1; i<= numberof(_DVV_FF_FILES); i++) 
      write, f, format= fmt0, _DVV_FF_FILES(i);
    write, format= fmt0, "";
  }
  
  //Background reference files
  if (bf) {
    write, f, format= fmt0, "CCD background files";
    write, f, format= fmt0, "====================";
    for(i= 1; i<= numberof(_DVV_BKG_FILES); i++) 
      write, f, format= fmt0, _DVV_BKG_FILES(i);
    write, format= fmt0, "";
  }
  
  //Filter lists
  if (chf) {
    if(!is_void(_DVV_CH1S_FILTER) ||
       !is_void(_DVV_CH2S_FILTER) ||
       !is_void(_DVV_CH1P_FILTER) ||
       !is_void(_DVV_CH2P_FILTER) ||
       !is_void(_DVV_CMMN_FILTER)) {
      write, f, format= fmt0, "Filter lists for each channel";
      write, f, format= fmt0, "=============================";
      if (!is_void(_DVV_CH1S_FILTER)) {
        write, f, format= fmt0, "   Ch1S filter: ";
        if (numberof(_DVV_CH1S_FILTER) > 1) for (i= 1; i<= numberof(_DVV_CH1S_FILTER); i++)
          write, f, format= fmt0, "      "+_DVV_CH1S_FILTER(i);
        else
          write, f, format= fmt0, "      "+_DVV_CH1S_FILTER;
      } else write, f, format= fmt0, "   No Ch1S filter ";
  
      if (!is_void(_DVV_CH2S_FILTER)) {
        write, f, format= fmt0, "   Ch2S filter: ";
        if (numberof(_DVV_CH2S_FILTER) > 1) for (i= 1; i<= numberof(_DVV_CH2S_FILTER); i++)
          write, f, format= fmt0, "      "+_DVV_CH2S_FILTER(i);
        else
          write, f, format= fmt0, "      "+_DVV_CH2S_FILTER;
      } else write, f, format= fmt0, "   No Ch2S filter ";

      if (!is_void(_DVV_CH1P_FILTER)) {
        write, f, format= fmt0, "   Ch1P filter: ";
        if (numberof(_DVV_CH1P_FILTER) > 1) for (i= 1; i<= numberof(_DVV_CH1P_FILTER); i++)
          write, f, format= fmt0, "      "+_DVV_CH1P_FILTER(i);
        else
          write, f, format= fmt0, "      "+_DVV_CH1P_FILTER;
      } else write, f, format= fmt0, "   No Ch1P filter ";

      if (!is_void(_DVV_CH2P_FILTER)) {
        write, f, format= fmt0, "   Ch2P filter: ";
        if (numberof(_DVV_CH2P_FILTER) > 1) for (i= 1; i<= numberof(_DVV_CH2P_FILTER); i++)
          write, f, format= fmt0, "      "+_DVV_CH2P_FILTER(i);
        else
          write, f, format= fmt0, "      "+_DVV_CH2P_FILTER;
      } else write, f, format= fmt0, "   No Ch2P filter ";
      if(!is_void(_DVV_CMMN_FILTER)) {
        write, f, format= fmt0, "   Common mode filter: ";
        if (numberof(_DVV_CMMN_FILTER) > 1) for (i= 1; i<= numberof(_DVV_CMMN_FILTER); i++)
          write, f, format= fmt0, "      "+_DVV_CMMN_FILTER(i);
        else
          write, f, format= fmt0, "      "+_DVV_CMMN_FILTER;
      } else write, f, format= fmt0, "   No common mode filter ";
    } else {
      write, f, format= fmt0, "No filters specified any channel";
      write, f, format= fmt0, "================================";
    }
    write, format= fmt0, "";
  }
  
  //List the coordinate mapping data
  if (maps) {
    //map= _DVV_REF.map;
    map= _DVV_MAPPING;
    write, f, format= fmt0, "Ch1S mapping parameters:";
    write, f, format= fmt0, "========================";
    //dvv_writeMap,_DVV_REF.map.ch1S, _DVV_REF.map.warp_refpts;
    dvv_writeMap,_DVV_MAPPING.ch1S, _DVV_MAPPING.warp_refpts;
    write, f, format= fmt0, "Ch2S mapping parameters:";
    write, f, format= fmt0, "========================";
    //dvv_writeMap, _DVV_REF.map.ch2S, _DVV_REF.map.warp_refpts;
    dvv_writeMap, _DVV_MAPPING.ch2S, _DVV_MAPPING.warp_refpts;
    write, f, format= fmt0, "Ch1P mapping parameters:";
    write, f, format= fmt0, "========================";
    //dvv_writeMap, _DVV_REF.map.ch1P, _DVV_REF.map.warp_refpts;
    dvv_writeMap, _DVV_MAPPING.ch1P, _DVV_MAPPING.warp_refpts;
    write, f, format= fmt0, "Ch2P mapping parameters:";
    write, f, format= fmt0, "========================";
    //dvv_writeMap, _DVV_REF.map.ch2P, _DVV_REF.map.warp_refpts;
    dvv_writeMap, _DVV_MAPPING.ch2P, _DVV_MAPPING.warp_refpts;
  }
  close, f;
}

func dvv_writeMap(map, ref)
/* DOCUMENT dvv_writeMap, map, ref

   Writes out the dewarp mapping points (grid points) that are used to
   map the raw pixel coordinates onto the target plane coordinate
   system.

   SEE ALSO:
     dvv_SetWarp, dvv_FixWarp, dvv_bicubicDewarp, dvv_Report
 */
{
  extern f;
  
  fmt0="%s\n";
  fmt1="  (%6.1f, %6.1f) --> (%6.1f, %6.1f), dx,dy: (%6.3f, %6.3f)\n";
  fmt2="  %s:  --> (%7.2f, %7.2f)\n";
  fmt3="  %s:  --> %5.4f\n";
  //write, format= fmt0, "[";
  g= *map.gridList;
  write, f, format= fmt2, "CENTER     :  ", g(1).x, g(1).y;
  write, f, format= fmt2, "BOTTOM LEFT:  ", g(2).x, g(2).y;
  write, f, format= fmt2, "TOP LEFT:     ", g(3).x, g(3).y;
  write, f, format= fmt2, "TOP RIGHT:    ", g(4).x, g(4).y;
  write, f, format= fmt2, "BOTTOM RIGHT: ", g(5).x, g(5).y;
  write, f, format= fmt0, "";
  write, f, format= fmt3, "ANGLE (DEG):  ", map.angle * 180.0/pi;
  write, f, format= fmt3, "X-SCALE:      ", map.sclx;
  write, f, format= fmt3, "Y-SCALE:      ", map.scly;
  write, f, format= fmt0, "";
  pts= *map.warp_pts;
  rpt= *ref
  for(i= 1; i<= numberof(pts); i++) {
    write, f, format= fmt1, rpt(i).x, rpt(i).y,  pts(i).x, pts(i).y,
    pts(i).x - rpt(i).x, pts(i).y - rpt(i).y;
  }
  write, f, format= fmt0, "";
}

func dvv_AssembleSet(file, prefix=, suffix=)
/* DOCUMENT dvv_AssembleSet, file, prefix=, suffix=

   Assembles an OHRV data set from a collection of TIF files that are
   named around the base filename <file> and writes them out to a
   single HDF file with the same structure as the hdf files generated
   by the LLE OHRV diagnostic.

   KEYWORDS:
     prefix=  prefix string applied to the base file name (default is "interferometer_")
     suffix=  suffix string applied to the base file name (default = ".TIF")
     
   SEE ALSO:
 */
{
  
  if(is_void(prefix)) prefix= "interferometer_";
  if(is_void(suffix)) suffix= ".TIF";
  
  bkg= getTIFF(prefix+file+"_bkg"+suffix);  //CCD dark current background
  dat= getTIFF(prefix+file+"_dat"+suffix);  //shot data
  flt= getTIFF(prefix+file+"_flt"+suffix);  //flat field data
  if (open(prefix+file+"_grid"+suffix, "r", 1))
    grid= getTIFF(prefix+file+"_grid"+suffix);//grid data
  else if (open(prefix+file+"_grd"+suffix, "r", 1))
    grid= getTIFF(prefix+file+"_grd"+suffix);//grid data
  ref= getTIFF(prefix+file+"_ref"+suffix);  //reference data

  strk= _(img_data(dat)(..,-), img_data(bkg));
  refr= _(img_data(flt)(..,-), img_data(grid), img_data(ref));
  
  fid= SDstart(prefix+file+".hdf", "w");

  if(typeof(strk) == "short") dt= _DFNT_UINT16; else dt= _hdf_YtoHtype(typeof(strk));
  sds= SDcreate(fid, "Streak_array", dt, dimsof(strk));
  SDwritedata, sds, strk, ignore_type= 1;

  if(typeof(refr) == "short") dt= _DFNT_UINT16; else dt= _hdf_YtoHtype(typeof(strk));
  sds= SDcreate(fid, "References", dt, dimsof(refr));
  SDwritedata, sds, refr, ignore_type= 1;
  SDend, fid;
}

func dvv_PAmap_average(shots)
{
  n= numberof(shots);
  mampl= mphase= mxoff= myoff= array(IMG_DAT(), n);
  dampl= dphase= dxoff= dyoff= dmsk= array(IMG_DAT(), n);
  
  for (i= 1; i <= n; i++) {
    ff= shots(i)+"/psf/"+shots(i)+"-phase-results.pdb";
    if (open(ff, "r", 1))
      f= openb(ff);
    else
      f= openb("../"+ff);
    
    restore, f, _DVV_AMPL_MAP, _DVV_PHASE_MAP, _DVV_XOFF_MAP, _DVV_YOFF_MAP,
      _DVV_AMPL_DATA, _DVV_PHASE_DATA, _DVV_XOFF_DATA, _DVV_YOFF_DATA,
      _DVV_MASK_DATA;
    
    mampl(i)= _DVV_AMPL_MAP;
    mphase(i)= _DVV_PHASE_MAP;
    mxoff(i)= _DVV_XOFF_MAP;
    myoff(i)= _DVV_YOFF_MAP;
    dampl(i)= _DVV_AMPL_DATA;
    dphase(i)= _DVV_PHASE_DATA;
    dxoff(i)= _DVV_XOFF_DATA;
    dyoff(i)= _DVV_YOFF_DATA;
    dmsk(i)= _DVV_MASK_DATA;
    close, f;
    
    if (i == 1) {
      amplS= img_copy(mampl(i), data= 0.0);
      phaseS= img_copy(mampl(i), data= 0.0);
      xoffS= img_copy(mampl(i), data= 0.0);
      yoffS= img_copy(mampl(i), data= 0.0);
      nx= amplS.nx; ny= amplS.ny;
      
      damplS= img_copy(dampl(i), data= 0.0);
      dphaseS= img_copy(dampl(i), data= 0.0);
      dxoffS= img_copy(dampl(i), data= 0.0);
      dyoffS= img_copy(dampl(i), data= 0.0);
      dmskS= img_copy(dampl(i), data= 0.0);
    }
    
    amplS= img_add(amplS, img_resample(mampl(i), nx, ny));
    phaseS= img_add(phaseS, img_resample(mphase(i), nx, ny));
    xoffS= img_add(xoffS, img_resample(mxoff(i), nx, ny));
    yoffS= img_add(yoffS, img_resample(myoff(i), nx, ny));
    
    damplS= img_add(damplS, dampl(i));
    dphaseS= img_add(dphaseS, dphase(i));
    dxoffS= img_add(dxoffS, dxoff(i));
    dyoffS= img_add(dyoffS, dyoff(i));
    dmskS= img_add(dmskS, dmsk(i));
  }
  
  amplS= img_div(amplS, n);
  phaseS= img_div(phaseS, n);
  xoffS= img_div(xoffS, n);
  yoffS= img_div(yoffS, n);
  
  damplS= img_div(damplS, n);
  dphaseS= img_div(dphaseS, n);
  dxoffS= img_div(dxoffS, n);
  dyoffS= img_div(dyoffS, n);
  dmskS= img_div(dmskS, n);
  
  obj= save(_DVV_AMPL_MAP= amplS, _DVV_PHASE_MAP= phaseS,
            _DVV_XOFF_MAP= xoffS, _DVV_YOFF_MAP= yoffS,
            _DVV_AMPL_DATA= damplS, _DVV_PHASE_DATA= dphaseS,
            _DVV_XOFF_DATA= dxoffS, _DVV_YOFF_DATA= dyoffS,
            _DVV_MAKS_DATA= dmskS, _DVV_PAmap_AVG_SHOTS= shots);
    
  return obj;
}

func dvv_exportData(..)
/* DOCUMENT dvv_exportData

   ARGUMENTS:

   KEYWORDS:

     
   SEE ALSO:

*/
{
  extern rf, dt;
  extern dtpUWf, rfpUW;
  extern _BINNING;
  
  file= shot+"-export-data.pdb";
 
  f= createb(file);
  vdat= dvv_Vel(dtpUWf);
  vbkg= dvv_Vel(rfpUW);
  rftot= dvv_frameSum(rf, average= 1);
  dttot= dvv_frameSum(dt, average= 1);
  ch1s= img_copy(dt.ch1S);

  if (img_dims(vdat)(2) == 2100) img_rebin, vdat, 2, average= 1;
  if (img_dims(vbkg)(2) == 2100) img_rebin, vbkg, 2, average= 1;
  if (img_dims(rftot)(2) == 2100) img_rebin, rftot, 2, average= 1;
  if (img_dims(dttot)(2) == 2100) img_rebin, dttot, 2, average= 1;
  if (img_dims(ch1s)(2) == 2100) img_rebin, ch1s, 2, average= 1;

  save, f, shot= shot;
  save, f, xscale= *vdat.xscale, yscale= *vdat.yscale;
  save, f, vdat= float(*vdat.data);
  save, f, vbkg= float(*vbkg.data);
  save, f, sdat= float(*dttot.data);
  save, f, sbkg= float(*rftot.data);
  save, f, ch1s= float(*ch1s.data);
  save, f, vunit= vdat.z_unit, sunit= dttot.z_unit,
    xunit= vdat.x_unit, yunit=vdat.y_unit;
  close, f;
}
