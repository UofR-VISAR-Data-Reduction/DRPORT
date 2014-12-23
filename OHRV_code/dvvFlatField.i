/*
  DVVFLATFIELD.I
  
  Definitions for I/O related routines for dVV/OHRV flat field data
  sets

  P.M. Celliers
  LLNL
  December, 2007

*/

/*
  
  dvv_ImportFlatField
  dvv_FilterFlatField
  dvv_ApplyFlatField
  
  dvv_SaveFFimg
  dvv_RestoreFFimg

  dvv_SaveFF_fset
  dvv_RestoreFF_fset

 */

func dvv_ImportFlatField(file, .., mask=, psf=, allplanes=, dwmethod=, ssdewarp=, noavg=)
/* DOCUMENT      dvv_ImportFlatField, <file1>, <file2>, ..., mask=,
                    psf=, allplanes=, ssdewarp=, dwmethod=, noavg=
            -or- dvv_ImportFlatField, [<file1>, <file2>, ...], mask=,
                    psf=, allplanes=, ssdewarp=, dwmethod=, noavg=
     
   Reads flat fielding reference data from shot files into
   memory. Multiple files may be listed either by passing a string
   array in the first argument, or by using multiple arguments, each
   specifying a single file name.  The flat field is defined as the
   average of the data contained in all the file sets.

   KEYWORDS:
     mask=
     psf=
     allplanes=
     dwmethod=
     ssdewarp=
     
   SEE ALSO:
 */
{
  extern _DVV_BKG, _DVV_REF, _DVV_FF;
  extern _DVV_FF_FILES, _DVV_MAPPING;
  extern _DVV_FF_IMG, _DVV_FF_ALLPLANES;
  extern _DVV_BKG_FILES;
  
  if(is_void(_DVV_BKG)) dvv_ImportBkg, _DVV_BKG_FILES;
  if(is_void(_DVV_MAPPING)) error, "Can't import flat field, no frame mapping defined";
  if (allplanes) _DVV_FF_ALLPLANES= 1;

  if(numberof(file) == 1) {
    //img= img_sub(getTIFF(file), _DVV_BKG);
    _DVV_FF_FILES= file;
    // if(allplanes) { 
    //   for(j= 1; j<=3; j++)
    //     img= img_add(img, img_sub(hdf_SDget(file, setname="References", plane= j, setdouble= 1), _DVV_BKG));
    //   img= img_add(img, img_sub(hdf_SDget(file, setname="Streak_array", plane= 1, setdouble= 1), _DVV_BKG));
    //   nn= 4;
    // } else {
    //   img= img_sub(hdf_SDget(file, setname="References", plane= 1, setdouble= 1), _DVV_BKG);
    //   nn= 1;
    // }
    fstring=  file;
    while(more_args()) {
      file= next_arg();
      // if(allplanes) {
      //   for(j= 1; j<=3; j++)
      //     img= img_add(img, img_sub(hdf_SDget(file, setname="References", plane= j, setdouble= 1), _DVV_BKG));
      //   img= img_add(img, img_sub(hdf_SDget(file, setname="Streak_array", plane= 1, setdouble= 1), _DVV_BKG));
      //   nn+=4;
      // } else {
      //   img= img_add(img, img_sub(hdf_SDget(file, setname="References", plane= 1, setdouble= 1), _DVV_BKG));
      //   nn++;
      // }
      _DVV_FF_FILES= _(_DVV_FF_FILES, file);
      fstring+="+"+file;
    }
  } else {
    fstring=  file(1);
    _DVV_FF_FILES= file(1);
    // if (allplanes) { 
    //   for(j= 1; j<=3; j++)
    //     img= img_add(img, img_sub(hdf_SDget(file(1), setname="References", plane= j, setdouble= 1), _DVV_BKG));
    //   img= img_add(img, img_sub(hdf_SDget(file(1), setname="Streak_array", plane= 1, setdouble= 1), _DVV_BKG));
    //   nn= 4;
    // } else {
    //   img= img_sub(hdf_SDget(file(1), setname="References", plane= 1, setdouble= 1), _DVV_BKG);
    //   nn= 1;
    // }
    for (i= 2; i<= numberof(file); i++) {
      if(strlen(file(i)) > 0) {
        // if(allplanes) {
        //   for(j= 1; j<=3; j++)
        //     img= img_add(img, img_sub(hdf_SDget(file(i), setname="References", plane= j, setdouble= 1), _DVV_BKG));
        //   img= img_add(img, img_sub(hdf_SDget(file(i), setname="Streak_array", plane= 1, setdouble= 1), _DVV_BKG));
        //   nn+=4;
        // } else {
        //   img= img_add(img, img_sub(hdf_SDget(file(i), setname="References", plane= 1, setdouble= 1), _DVV_BKG));
        //   nn++;
        // }
        fstring+="+"+file(i);
        _DVV_FF_FILES= _(_DVV_FF_FILES, file(i));
      }
    }
    //_DVV_FF_FILES= fstring;
  }
  
  if (!dvv_RestoreFFimg(allplanes= allplanes)) {
    //nn= numberof(_DVV_BKG_FILES);
    print, "dvv_ImportFlatField - importing files";
    nn= 0;
    img= img_copy(_DVV_BKG, data= 0.0);
    for(k= 1; k<= numberof(_DVV_FF_FILES); k++) {
      file= _DVV_FF_FILES(k);
      // if (k == 1) {
      //   if (allplanes) {
      //     for(j= 1; j<=3; j++)
      //       img= img_add(img, img_sub(hdf_SDget(file, setname="References", plane= j, setdouble= 1), _DVV_BKG));
      //     img= img_add(img, img_sub(hdf_SDget(file, setname="Streak_array", plane= 1, setdouble= 1), _DVV_BKG));
      //     nn+= 4;
      //   } else {
      //     img= img_add(img, img_sub(hdf_SDget(file, setname="References", plane= 1, setdouble= 1), _DVV_BKG));
      //     nn+=1;
      //   }
      // } else {
        if (allplanes) {
          for(j= 1; j<=3; j++)
            img= img_add(img, img_sub(hdf_SDget(file, setname="References", plane= j, setdouble= 1), _DVV_BKG));
          img= img_add(img, img_sub(hdf_SDget(file, setname="Streak_array", plane= 1, setdouble= 1), _DVV_BKG));
          nn+= 4;
        } else {
          img= img_add(img, img_sub(hdf_SDget(file, setname="References", plane= 1, setdouble= 1), _DVV_BKG));
          nn+= 1;
        }
      // }
    }
    //img.data= &((*img.data)/nn);
    *img.data/= nn;
    img.shotid= fstring;
    imgavg= img_data(img)(*)(avg);
    _DVV_FF_IMG= img_copy(img, data= float(img_data(img)/imgavg));
    dvv_SaveFFimg, allplanes= allplanes;
  }
  
  //fs= dvv_GenerateFrameSet(img, _DVV_REF, ssdewarp= ssdewarp);
  print, "dvv_ImportFlatField - generating frame set";
  fs= dvv_GenerateFrameSet(_DVV_FF_IMG, _DVV_MAPPING, ssdewarp= ssdewarp, dwmethod= dwmethod);
  chavg= img_add(fs.ch1P, fs.ch2P, fs.ch1S, fs.ch2S, average= 1);
  print, "dvv_ImportFlatField - fscopy";
  _DVV_FF= dvv_FSCopy(fs);
  _DVV_FF.ch1S= img_div(fs.ch1S, chavg);
  _DVV_FF.ch2S= img_div(fs.ch2S, chavg);
  _DVV_FF.ch1P= img_div(fs.ch1P, chavg);
  _DVV_FF.ch2P= img_div(fs.ch2P, chavg);
}

func dvv_FilterFlatField(.., mask=, psf=)
/* DOCUMENT      ff= dvv_FilterFlatField, mask=, psf=
            -or- ff_filtered= dvv_FilterFlatField(mask=, psf=)

   Defines a flat field reference using the data contained in a frame
   set containing flat field data.  The flat field reference is stored
   in the global variable _DVV_FF, and also returned if this function
   is called as a function.

   KEYWORDS:
     mask=
     psf=

   SEE ALSO:
     dvv_ImportFlatField, dvv_GenerateFrameSet
*/
{
  extern _DVV_FF, _DVV_FF_MASK;

  print, "dvv_FilterFlatField - fscopy";
  FF_filtered= dvv_FSCopy(_DVV_FF);
    
  if (!is_void(mask)) {
    if(typeof(mask) == "pointer") {
      for(i= 1; i<= numberof(mask); i++) {
        write, format="Applying FF filter mask #%d to FF images ... \n", i;
        FF_filtered.ch1S= dvv_Filter(FF_filtered.ch1S, *mask(i));
        FF_filtered.ch2S= dvv_Filter(FF_filtered.ch2S, *mask(i));
        FF_filtered.ch1P= dvv_Filter(FF_filtered.ch1P, *mask(i));
        FF_filtered.ch2P= dvv_Filter(FF_filtered.ch2P, *mask(i));
      }
    } else {
      print, "Applying FF filter mask to FF images ... ";
      FF_filtered.ch1S= dvv_Filter(FF_filtered.ch1S, mask);
      FF_filtered.ch2S= dvv_Filter(FF_filtered.ch2S, mask);
      FF_filtered.ch1P= dvv_Filter(FF_filtered.ch1P, mask);
      FF_filtered.ch2P= dvv_Filter(FF_filtered.ch2P, mask);
    }
    _DVV_FF_MASK= mask;
  }

  if(!is_void(psf)) {
    info, psf;
    if (structof(psf) == double) {
      FF_filtered.psf= psf;
      psf= dvv_mkGaussian(FF_filtered.ch1S, psf);
    }
    print, "Applying FF smoothing PSF to FF images ... ";
    FF_filtered.ch1S= dvv_Smooth(FF_filtered.ch1S, psf);
    print, "ch1S done ... ";
    FF_filtered.ch2S= dvv_Smooth(FF_filtered.ch2S, psf);
    print, "Ch2S done ... ";
    FF_filtered.ch1P= dvv_Smooth(FF_filtered.ch1P, psf);
    print, "Ch1P done ... ";
    FF_filtered.ch2P= dvv_Smooth(FF_filtered.ch2P, psf);
    print, "Ch2P done ... ";
  }
  print, "Done FF filtering ... ";
  print, "dvv_FilterFlatField b - fscopy";
  if (!am_subroutine()) return FF_filtered; else _DVV_FF= dvv_FSCopy(FF_filtered);
}

func dvv_ApplyFlatField(fs, ff=)
/* DOCUMENT ff= dvv_ApplyFlatField(fs)

   Applies flat field reference data to an input frameset

   ARGUMENTS:
     fs  - FrameSet structure containing the flat field data

   KEYWORDS:
     ff=   Define this to use a flat field reference
           other than the default (which is _DVV_FF)
     
   SEE ALSO:
     dvv_GenerateFrameSet
*/
{
  extern _DVV_FF;
  
  if (is_void(ff)) ff= _DVV_FF;
  if (is_void(ff)) error, "undefined flat field data";
    
  print, "dvv_ApplyFlatField - fscopy";
  newfs= dvv_FSCopy(fs);

  //If the data set has been rebinned, need to rebin the flat
  //field equivalently (raw flat field data is not rebinned
  //when it is created)
  if (fs.binning > 1) {
    newfs.ch1S= img_div(fs.ch1S, img_rebin(ff.ch1S, fs.binning));
    newfs.ch2S= img_div(fs.ch2S, img_rebin(ff.ch2S, fs.binning));
    newfs.ch1P= img_div(fs.ch1P, img_rebin(ff.ch1P, fs.binning));
    newfs.ch2P= img_div(fs.ch2P, img_rebin(ff.ch2P, fs.binning));
  } else {
    newfs.ch1S= img_div(fs.ch1S, ff.ch1S);
    newfs.ch2S= img_div(fs.ch2S, ff.ch2S);
    newfs.ch1P= img_div(fs.ch1P, ff.ch1P);
    newfs.ch2P= img_div(fs.ch2P, ff.ch2P);
  }
  return newfs;
}

func dvv_SaveFFimg(.., allplanes=)
/* DOCUMENT dvv_SaveFFimg

   Checks the current CCD flat_field data set (file list) against the
   contents of the file "../baseline/flat_field.pdb" if that file
   exists.  If the file list and/or _DVV_FF_ALLPLANES differ then the
   "../baseline/flat_field.pdb" file is updated with the current
   flat_field data set.  if "../baseline/flat_field/pdb" does not
   exist, then a new copy is created.
     
   SEE ALSO:
     dvv_RestoreFFimg, dvv_ImportFlatField
 */
{
  extern _DVV_FF_IMG, _DVV_FF_FILES, _DVV_FF_ALLPLANES;
  if(is_void(allplanes)) allplanes= 0;
  if (open("../baseline/flat_field.pdb","rb",1)) {
    f= openb("../baseline/flat_field.pdb");
    vars= *get_vars(f)(1);
    if (is_member(vars, "_DVV_FF_FILES")) restore, f, "_DVV_FF_FILES", flist;
    if (is_member(vars, "_DVV_FF_ALLPLANES")) restore, f, "_DVV_FF_ALLPLANES", apl;
    close, f;
    if (numberof(flist) == numberof(_DVV_FF_FILES) &&
        allof(flist == _DVV_FF_FILES) &&
        allplanes == apl) {
      return;
    } 
  }
  f= createb("../baseline/flat_field.pdb");
  _p= POINT();
  save, f, _p;
  if(!is_void(_DVV_FF_FILES)) save, f, _DVV_FF_FILES;
  if(!is_void(_DVV_FF_IMG)) save, f, _DVV_FF_IMG= img_floatData(_DVV_FF_IMG);
  save, f, _DVV_FF_ALLPLANES= allplanes;
  close, f;
}

func dvv_RestoreFFimg(.., allplanes=)
/* DOCUMENT dvv_RestoreFFimg

   Checks the current CCD flat_field data set (file list) against the
   contents of the file "../baseline/flat_field.pdb" if that file
   exists.  If the file list matches the list stored in
   "../baseline/flat_field.pdb" then the flat_field data set is restored
   from that file and this function returns 1. If there is a mismatch
   nothing is imported and this function returns 0.

   Returns 1 - successful restore
           0 - no restore took place
     
   SEE ALSO:
     dvv_SaveFFimg, dvv_ImportFlatField
     
*/
{
  extern _DVV_FF_IMG, _DVV_FF_FILES, _DVV_FF_ALLPLANES;
  if(is_void(allplanes)) allplanes= 0;
  if (open("../baseline/flat_field.pdb","rb",1)) {
    f= openb("../baseline/flat_field.pdb");
    vars= *get_vars(f)(1);
    if (is_member(vars, "_DVV_FF_FILES")) restore, f, "_DVV_FF_FILES", flist;
    if (is_member(vars, "_DVV_FF_ALLPLANES")) restore, f, "_DVV_FF_ALLPLANES", apl;
    print, "apl= ", apl, "flist= ", flist;
    if (numberof(flist) == numberof(_DVV_FF_FILES) &&
        allof(flist == _DVV_FF_FILES) &&
        allplanes == apl) {
      restore, f, _DVV_FF_IMG;
      restore, f, _DVV_FF_ALLPLANES;
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

func dvv_SaveFF_fset(path)
/* DOCUMENT dvv_SaveFF_fset

   Checks the current CCD flat_field data set (file list) against the
   contents of the file "../baseline/flat_field.pdb" if that file
   exists.  If the file list and/or _DVV_FF_ALLPLANES differ then the
   "../baseline/flat_field.pdb" file is updated with the current
   flat_field data set.  if "../baseline/flat_field/pdb" does not
   exist, then a new copy is created.
     
   SEE ALSO:
     dvv_RestoreFFimg, dvv_ImportFlatField
 */
{
  extern _DVV_COMMENT, _DVV_FSET;
  extern _DVV_FIELD, _DVV_MAGFACTOR;
  extern _DVV_DEWARP_METHOD;
  extern _DVV_FF, _DVV_FF_FILES;
  extern _DVV_PSF_FF, _DVV_FF_MASK;
  extern _DVV_FFM1, _DVV_FFM2;

  print, "dvv_SaveFF_fset: saving to: "+path+"-fset-ff.pdb";
  dvv_SaveFrameSet(_DVV_FF, path+"-fset-ff.pdb");
  
}

func dvv_RestoreFF_fset(path)
/* DOCUMENT dvv_RestoreReferenceData, path

   ARGUMENTS:
     file  - file name

   KEYWORDS:
     
   SEE ALSO:
     dvv_SaveReferenceData
*/
{
  extern _DVV_FF;

  if (!open(path+"-fset-ff.pdb", "rb", 1)) {
    return 0;
  } else {
    print, "dvv_RestoreFF_fset: restoring from: "+path+"-fset-ff.pdb";
    _DVV_FF= dvv_RestoreFrameSet(path+"-fset-ff.pdb");
    return 1;
  }
}
