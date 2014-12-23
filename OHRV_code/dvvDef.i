/*
  DVVDEF.I
  Definitions for routines to process quadrature interferometer CCD data

  $Id$

  P.M. Celliers
  Lawrence Livermore National Laboratory
  December, 2007

*/

/*
  dvv_FrameMapping
  dvv_SetMapping
  dvv_FrameSet
  dvv_PhaseSet
  dvv_FocusSet
  dvv_FSCopy

 */

//#include "mpeg.i"
//#include "unwrap0.i"
#include "unwrap2d.i"
#include "yutil.i"
//#include "tiff_1.i"
#include "hdf.i"
//#include "yeti_fftw.i"
#include "fft_utils.i"
#include "yfftw.i"
#include "image.i"
#include "lineout.i"
#include "points.i"
#include "phase_fft.i"
#include "color.i"
#include "eos.i"
#include "polywarp.i"
#include "fitpoly.i"
#include "poly.i"
#include "fitsurf.i"
#include "movie.i"
#include "convol.i"
#include "util.i"
#include "lmfit.i"
#include "limutils.i"
#include "etalon-dispersion.i"
#include "random.i"
#include "gaussfit.i"

_red= to_rgb([0.0, 0.86, 0.69]);
_green= to_rgb([112.0, 0.86, 0.69]);
_blue= to_rgb([234, 0.86, 0.69]);
_magenta= to_rgb([289.0, 0.86, 0.69]);
_cyan= to_rgb([186.0, 0.86, 0.69]);
_orange= to_rgb([40.0, 0.98, 0.79]);
_black= to_rgb([0.0, 1.0, 0.0]);
_grey= to_rgb([50.0, 70.0, 30.0]);
_colorset=[_black, _red, _blue, _magenta, _green, _cyan, _orange];
_nclrs= 7;

/*
_DVV_COMMENT=[];
_DVV_BKG=[];
_DVV_REF=[];
_DVV_FF=[];
_DVV_FF_MASK=[];
*/

_Xcoord= _Ycoord= 0.0;
_Xrange= _Yrange= 145.0;

_REG= REGION(x1= -400, x2= 400, y1= -400, y2= 400);
_lh= LO_POS(x1= -0.001, x2= 0.001, or= "h");
_lv= LO_POS(x1= -0.001, x2= 0.001, or= "v");

_FrameCount= 0;
_FrameCountMax= 100;
_Zmin= 500;
_Zmax= 20000;

_DVV_PSF_FF= [];
_DVV_PSF_RF= [];
_DVV_PSF_DT= [];
_DVV_BINNING= [];
_LAMBDA= 0.395;

_DVV_CHFILTER_MAXCOUNT= 5;
_DVV_CH1S_FILTER=[];
_DVV_CH2S_FILTER=[];
_DVV_CH1P_FILTER=[];
_DVV_CH2P_FILTER=[];
_DVV_CMMN_FILTER=[];

_pcounter= 0;
_NO_FFTW=1;

struct dvv_FrameMapping
/* DOCUMENT dvv_FrameMapping()

     Struct definition contains basic data required for
     performing image rotation, translation, scaling and dewarping
     using an image of a grid as the basic data input.

     Elements are:
       pointer gridList;  - list of user-selected grid points
       double  angle;     - rotation correction
       double  sclx;      - pixel-to-micron x-scale factor
       double  scly;      - pixel-to-micron y-scale factor
       pointer kx;        - polynomial warp x-coefficients
       pointer ky;        - polynomial warp y-coefficients
       pointer warp_pts;  - warping grid grid coordinates
       pointer aX;        - thin plate spline warp x-coefficients
       pointer aY;        - thin plate spline warp y-coefficients
     
   SEE ALSO:
     dvv_SetMapping, dvv_FrameSet
*/
{
  pointer gridList;
  double  angle;
  double  sclx;
  double  scly;
  double  mag;
  pointer kx;
  pointer ky;
  pointer warp_pts;
  pointer aX;
  pointer aY;
  pointer repairedList;
}

func dvv_copyFrameMapping(fmap)
/* DOCUMENT dvv_copyFrameMapping, fmap

   Copies a frame mapping struct, and allows for the original to use
   an older definition.
     
   SEE ALSO:
     dvv_FrameMapping, dvv_SetMappingCopy
 */
{
  
  elems= struct_element(fmap);

  nfmap= dvv_FrameMapping();
  nfmap.gridList= fmap.gridList;
  nfmap.angle= fmap.angle;
  nfmap.sclx= fmap.sclx;
  nfmap.scly= fmap.scly;
  nfmap.kx= fmap.kx;
  nfmap.ky= fmap.ky;
  nfmap.warp_pts= fmap.warp_pts;
  //print, "nfmap= ", nfmap;
  if (anyof(elems == "aX")) nfmap.aX= &(*fmap.aX);
  if (anyof(elems == "aY")) nfmap.aY= &(*fmap.aY);
  if (anyof(elems == "repairedList")) nfmap.repairedList= &(*fmap.repairedList);
  return nfmap;
}
  
struct dvv_SetMapping
/* DOCUMENT dvv_SetMapping()

     Struct definition for the grid mapping data to map each of
     the four frames onto a common coordinate system.  This struct
     contains four instances of dvv_FrameMapping structs.
     
     Elements are:
       dvv_FrameMapping ch1S;  - channel 1, s-polarization frame
       dvv_FrameMapping ch1P;  - channel 1, p-polarization frame
       dvv_FrameMapping ch2S;  - channel 2, s-polarization frame
       dvv_FrameMapping ch2P;  - channel 2, p-polarization frame
     
   SEE ALSO:
     dvv_FrameMapping, dvv_FrameSet
*/
{
  string note;
  pointer warp_refpts;
  pointer xc;
  pointer yc;
  dvv_FrameMapping ch1S;
  dvv_FrameMapping ch1P;
  dvv_FrameMapping ch2S;
  dvv_FrameMapping ch2P;
}

func dvv_SetMappingCopy(map)
/* DOCUMENT dvv_SetMappingCopy()

   Copies a set mapping
     
   SEE ALSO:
 */
{
  nmap= dvv_SetMapping();
  elems= struct_element(map);
  if (anyof(elems == "note")) nmap.note= map.note;
  if (anyof(elems == "xc")) {xc= *map.xc; nmap.xc= &xc;}
  if (anyof(elems == "yc")) {yc= *map.yc; nmap.yc= &yc;}
  nmap.warp_refpts= map.warp_refpts;
  nmap.ch1S= dvv_copyFrameMapping(map.ch1S);
  nmap.ch2S= dvv_copyFrameMapping(map.ch2S);
  nmap.ch1P= dvv_copyFrameMapping(map.ch1P);
  nmap.ch2P= dvv_copyFrameMapping(map.ch2P);
  return nmap;
}

struct dvv_FrameSet
/* DOCUMENT dvv_FrameSet()

     Struct definition to contain all the data in a single set
     of frames (a single CCD exposure).  The individual channels
     have been mapped onto a common coordinate system using the
     mapping data contained 
     
     Elements are:
       long   warped; - indicates whether the warping transformation has been
                        applied to the data set
       long   binning; - binning factor used when reading data files
       double psf; - gaussian FWHM of smoothing function applied to each frame
       string filter; - text string containing filter definitions applied to each frame
       long   median_filter; - if non-zero gives median filter extent (pixels)
       string     dwmethod; - method used for dewarp
       long       ssdewarp; - single step or two-step dewarp
       double     magfactor; - magnification factor
       dvv_SetMapping   map;   - mapping data for registration
       IMG_DAT      raw;   - raw 4-frame image data
       IMG_DAT      ch1S;  - channel 1, s-polarization frame after mapping
       IMG_DAT      ch1P;  - channel 1, p-polarization frame after mapping
       IMG_DAT      ch2S;  - channel 2, s-polarization frame after mapping
       IMG_DAT      ch2P;  - channel 2, p-polarization frame after mapping
     
   SEE ALSO:
     dvv_SetMapping, dvv_FrameMapping, IMG_DAT
*/
{
  long       warped;
  long       binning;
  double     psf;
  string     filter;
  long       median_filter;
  string     dwmethod;
  long       ssdewarp;
  double     magfactor;
  dvv_SetMapping map;
  //IMG_DAT    raw;
  IMG_DAT    ch1S;
  IMG_DAT    ch1P;
  IMG_DAT    ch2S;
  IMG_DAT    ch2P;
}

func dvv_FSCopy(fs, reg=, floatData=)
/* DOCUMENT dvv_FSCopy(fs, reg=, floatData=)

     Generats a new copy of a dvv_FrameSet data object, including
     copies of each internal element.

   KEYWORDS:
     reg=  copies only the data within reg
     
   SEE ALSO:
     dvv_FrameSet, dvv_SetMapping, dvv_FrameMapping, IMG_DAT
*/
{
  newfs= dvv_FrameSet();
  newfs.warped= fs.warped;
  newfs.binning= fs.binning;
  if (is_member(fs, "psf")) newfs.psf= fs.psf;
  if (is_member(fs, "filter")) newfs.filter= fs.filter;
  if (is_member(fs, "median_filter")) newfs.median_filter= fs.median_filter;
  if (is_member(fs, "ssdewarp")) newfs.ssdewarp= fs.ssdewarp;
  if (is_member(fs, "dwmethod")) newfs.dwmethod= fs.dwmethod;
  if (is_member(fs, "magfactor")) newfs.magfactor= fs.magfactor;
  newfs.map= dvv_SetMappingCopy(fs.map);
  //newfs.raw= img_copy(fs.raw);
  if (is_void(reg)) {
    newfs.ch1S= img_copy(fs.ch1S);
    newfs.ch1P= img_copy(fs.ch1P);
    newfs.ch2S= img_copy(fs.ch2S);
    newfs.ch2P= img_copy(fs.ch2P);
  } else {
    newfs.ch1S= img_extract(fs.ch1S, reg);
    newfs.ch1P= img_extract(fs.ch1P, reg);
    newfs.ch2S= img_extract(fs.ch2S, reg);
    newfs.ch2P= img_extract(fs.ch2P, reg);
  }
  
  if (floatData) {
    newfs.ch1S= img_floatData(newfs.ch1S);
    newfs.ch2S= img_floatData(newfs.ch2S);
    newfs.ch1P= img_floatData(newfs.ch1P);
    newfs.ch2P= img_floatData(newfs.ch2P);
  }
  
  return newfs;
}

struct dvv_PhaseSet
/* DOCUMENT dvv_PhaseSet()

     Struct definition to contain all the reduced data from a single
     data set
     
     Elements are:
       pointer      MLO;   - Lineout through pre-imposed modulation
       pointer      HLO;   - Lineout through pre-imposed modulation
       pointer      VLO;   - Lineout through pre-imposed modulation
       pointer      dtVaz; - azimuthally averaged velocitty spectrum of the data 
       pointer      rfVaz; - azimuthally averaged velocitty spectrum of the reference 
       IMG_DAT      dtp;   - extracted phase (wrapped) from data frame
       IMG_DAT      dtpUW; - extracted phase (unwrapped) from data frame
       IMG_DAT      rfp;   - extracted phase (wrapped) from reference frame
       IMG_DAT      rfpUW; - extracted phase (unwrapped) from reference frame
       IMG_DAT      ampl;  - amplitude (non-fringing signal)
       IMG_DAT      vel;   - velocity (after applying _VPF)
       IMG_DAT      dtV;   - velocity spectrum
     
   SEE ALSO:
     dvv_SavePhaseData, dvv_RestorePhaseData
*/
{
  pointer MLO;
  pointer HLO;
  pointer VLO;
  pointer dtVaz;
  pointer rfVaz;
  IMG_DAT dtp;
  IMG_DAT dtpUW;
  IMG_DAT rfp;
  IMG_DAT rfpUW;
  IMG_DAT ampl;
  IMG_DAT vel;
  IMG_DAT dtV;
}

struct dvv_FocusSet
/* DOCUMENT dvv_FocusSet
     
   SEE ALSO:
 */
{
  string file;
  double pos;
  pointer ch1Sv;
  pointer ch2Sv;
  pointer ch1Pv;
  pointer ch2Pv;
  pointer ch1Sh;
  pointer ch2Sh;
  pointer ch1Ph;
  pointer ch2Ph;
}

include, "dvvFile.i";
include, "dvvFlatField.i";
include, "dvvMisc.i";
include, "dvvSetup.i";
include, "dvvShow.i";
include, "dvvSpectrum.i";
include, "dvvShow.i";
include, "dvvProcess.i";
require, "rgb4.i";

_rgb= _rgb_obj();
