/*
  DVVSHOW.I
  
  Definitions for displaying dVV/OHRV results from analysis of dVV/OHRV
  data sets

  P.M. Celliers
  LLNL
  December, 2007

*/

/*
  dvv_GridMovie
  dvv_GridAvg
  dvv_DisplayGridavg
  dvv_GridFrame
  dvv_DrawFrame
  dvv_ShowLineouts
  dvv_FocusLineout
  dvv_plFocus
  dvv_mkFocusSpect
  dvv_ShowFocusSpect
  dvv_ShowWarpSet
  dvv_ShowWarpDeviations
  dvv_DisplayImg
  dvv_DisplaySet
  dvv_DisplayPhaseMapping
  dvv_DisplayMapping
  dvv_Liss
  dvv_plCircle
  dvv_plsp
  dvv_ShowUWPhase
  dvv_ShowUWVelocity
  dvv_ShowUWVel3D
  dvv_ShowUWVelContour
  dvv_ShowSpectV
  dvv_ShowIContour
  dvv_CheckNorm
  dvv_ShowChannels
  dvv_ImprintPlots
  dvv_window
  
 */

func dvv_GridMovie(file, fs, lgscl=, warp_pts=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _currentSet, _FrameCount, _FrameCountMax;
  extern _Zmin, _Zmax;
  extern _Xrange, _Yrange;
  extern _lgscl, _warp_pts;

  _lgscl= lgscl;
  _currentSet=fs;
  _FrameCount= 0;
  _FrameCountMax= 81*4;
  _warp_pts= warp_pts;
  
  mpeg_movie, file, dvv_GridFrame, 1000.0, mpg_params=[1000000, 24, 1, 1];
}

func dvv_GridAvg(fs)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _Xrange, _Yrange;

  reg= REGION(x1= -_Xrange/2, x2= _Xrange/2,
              y1= -_Yrange/2, y2= _Yrange/2);
  
  //wref= (*_DVV_REF.map.warp_refpts);
  wref= (*_DVV_MAPPING.warp_refpts);
  
  ch1Sa= img_extract(fs.ch1S, reg);
  ch1Sa.data= &array(0.0, dimsof(*ch1Sa.data));
  ch1Pa= img_copy(ch1Sa); ch1Pa.shotid= "ch1P";
  ch2Sa= img_copy(ch1Sa); ch2Sa.shotid= "ch2S";
  ch2Pa= img_copy(ch1Sa); ch2Pa.shotid= "ch2P";

  //Generating interpolations
  tl1s= tl2cub(,*fs.ch1S.xscale, *fs.ch1S.yscale, *fs.ch1S.data);
  tl2s= tl2cub(,*fs.ch2S.xscale, *fs.ch2S.yscale, *fs.ch2S.data);
  tl1p= tl2cub(,*fs.ch1P.xscale, *fs.ch1P.yscale, *fs.ch1P.data);
  tl2p= tl2cub(,*fs.ch2P.xscale, *fs.ch2P.yscale, *fs.ch2P.data);
  
  xg= (*ch1Sa.xscale)(,-:1:ch1Sa.ny);
  yg= (*ch1Sa.yscale)(-:1:ch1Sa.nx,);

  nw= numberof(wref);
  nmid= nw/2+1;
  if(!fs.warped) {
    print, "unwarped data ...";
//     Xc1so= (*_DVV_REF.map.ch1S.warp_pts)(nmid).x;
//     Yc1so= (*_DVV_REF.map.ch1S.warp_pts)(nmid).y;
//     Xc2so= (*_DVV_REF.map.ch2S.warp_pts)(nmid).x;
//     Yc2so= (*_DVV_REF.map.ch2S.warp_pts)(nmid).y;
//     Xc1po= (*_DVV_REF.map.ch1P.warp_pts)(nmid).x;
//     Yc1po= (*_DVV_REF.map.ch1P.warp_pts)(nmid).y;
//     Xc2po= (*_DVV_REF.map.ch2P.warp_pts)(nmid).x;
//     Yc2po= (*_DVV_REF.map.ch2P.warp_pts)(nmid).y;
    Xc1so= (*_DVV_MAPPING.ch1S.warp_pts)(nmid).x;
    Yc1so= (*_DVV_MAPPING.ch1S.warp_pts)(nmid).y;
    Xc2so= (*_DVV_MAPPING.ch2S.warp_pts)(nmid).x;
    Yc2so= (*_DVV_MAPPING.ch2S.warp_pts)(nmid).y;
    Xc1po= (*_DVV_MAPPING.ch1P.warp_pts)(nmid).x;
    Yc1po= (*_DVV_MAPPING.ch1P.warp_pts)(nmid).y;
    Xc2po= (*_DVV_MAPPING.ch2P.warp_pts)(nmid).x;
    Yc2po= (*_DVV_MAPPING.ch2P.warp_pts)(nmid).y;
  } else {
    Xc1so= Yc1so= Xc2so= Yc2so= 0.0;
    Xc1po= Yc1po= Xc2po= Yc2po= 0.0;
  }

  //Sum the grid image  over all the nodes
  for(kk= 1; kk<= nw; kk++) {
    if (!fs.warped) {
//       Xc1s= (*_DVV_REF.map.ch1S.warp_pts)(kk).x - Xc1so;
//       Yc1s= (*_DVV_REF.map.ch1S.warp_pts)(kk).y - Yc1so;
//       Xc2s= (*_DVV_REF.map.ch2S.warp_pts)(kk).x - Xc2so;
//       Yc2s= (*_DVV_REF.map.ch2S.warp_pts)(kk).y - Yc2so;
//       Xc1p= (*_DVV_REF.map.ch1P.warp_pts)(kk).x - Xc1po;
//       Yc1p= (*_DVV_REF.map.ch1P.warp_pts)(kk).y - Yc1po;
//       Xc2p= (*_DVV_REF.map.ch2P.warp_pts)(kk).x - Xc2po;
//       Yc2p= (*_DVV_REF.map.ch2P.warp_pts)(kk).y - Yc2po;
      Xc1s= (*_DVV_MAPPING.ch1S.warp_pts)(kk).x - Xc1so;
      Yc1s= (*_DVV_MAPPING.ch1S.warp_pts)(kk).y - Yc1so;
      Xc2s= (*_DVV_MAPPING.ch2S.warp_pts)(kk).x - Xc2so;
      Yc2s= (*_DVV_MAPPING.ch2S.warp_pts)(kk).y - Yc2so;
      Xc1p= (*_DVV_MAPPING.ch1P.warp_pts)(kk).x - Xc1po;
      Yc1p= (*_DVV_MAPPING.ch1P.warp_pts)(kk).y - Yc1po;
      Xc2p= (*_DVV_MAPPING.ch2P.warp_pts)(kk).x - Xc2po;
      Yc2p= (*_DVV_MAPPING.ch2P.warp_pts)(kk).y - Yc2po;
    } else {
      Xc1s= Xc2s= Xc1p= Xc2p= wref(kk).x;
      Yc1s= Yc2s= Yc1p= Yc2p= wref(kk).y;
    }
    (*ch1Sa.data) += tl2cub(tl1s, xg-Xc1s, yg-Yc1s); 
    (*ch2Sa.data) += tl2cub(tl2s, xg-Xc2s, yg-Yc2s); 
    (*ch1Pa.data) += tl2cub(tl1p, xg-Xc1p, yg-Yc1p); 
    (*ch2Pa.data) += tl2cub(tl2p, xg-Xc2p, yg-Yc2p); 
  }
  (*ch1Sa.data)/= nw; 
  (*ch2Sa.data)/= nw;
  (*ch1Pa.data)/= nw;
  (*ch2Pa.data)/= nw;
  return [ch1Sa, ch2Sa, ch1Pa, ch2Pa];
}

func dvv_DisplayGridavg(ga, zmn, zmx)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  _lh= [LO_POS(x1= -20, x2= -10, or= "h"), LO_POS(x1= 10, x2= 20, or= "h")];
  _lv= [LO_POS(x1= -20, x2= -10, or= "v"), LO_POS(x1= 10, x2= 20, or= "v")];

  for(i= 1; i<= 4; i++) dvv_DisplayImg, ga(i), w= i, zmn, zmx;
}

func dvv_GridFrame(k)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _currentSet, _Xcoord, _Ycoord, _type;
  extern _FrameCount, _FrameCountMax;
  extern _Zmin, _Zmax;
  extern _Xrange, _Yrange;
  extern _lgscl, _warp_pts;

  k= (k-1)%(81*4)+ 1;
  kk= (k-1)/4 + 1;
  l= (k-1)%4 + 1;
  i= (kk-1)/9 + 1;
  j= (kk-1)%9 + 1;

  if(_warp_pts) {
    if(l == 1) {
      print, "warp_pts: 1";
//       _Xcoord= (*_DVV_REF.map.ch1S.warp_pts)(kk).x;
//       _Ycoord= (*_DVV_REF.map.ch1S.warp_pts)(kk).y;
      _Xcoord= (*_DVV_MAPPING.ch1S.warp_pts)(kk).x;
      _Ycoord= (*_DVV_MAPPING.ch1S.warp_pts)(kk).y;
    } else if (l == 2) {
      print, "warp_pts: 2";
//       _Xcoord= (*_DVV_REF.map.ch1P.warp_pts)(kk).x;
//       _Ycoord= (*_DVV_REF.map.ch1P.warp_pts)(kk).y;
      _Xcoord= (*_DVV_MAPPING.ch1P.warp_pts)(kk).x;
      _Ycoord= (*_DVV_MAPPING.ch1P.warp_pts)(kk).y;
    } else if (l == 3) {
      print, "warp_pts: 3";
//       _Xcoord= (*_DVV_REF.map.ch2S.warp_pts)(kk).x;
//       _Ycoord= (*_DVV_REF.map.ch2S.warp_pts)(kk).y;
      _Xcoord= (*_DVV_MAPPING.ch2S.warp_pts)(kk).x;
      _Ycoord= (*_DVV_MAPPING.ch2S.warp_pts)(kk).y;
    } else if (l == 4) {
      print, "warp_pts: 4";
//       _Xcoord= (*_DVV_REF.map.ch2P.warp_pts)(kk).x;
//       _Ycoord= (*_DVV_REF.map.ch2P.warp_pts)(kk).y;
      _Xcoord= (*_DVV_MAPPING.ch2P.warp_pts)(kk).x;
      _Ycoord= (*_DVV_MAPPING.ch2P.warp_pts)(kk).y;
    }
  } else {
    _Xcoord= -400 + (i-1)*100;
    _Ycoord= -400 + (j-1)*100;
  }
  
  ret= dvv_DrawFrame(l, type= _type, lgscl= _lgscl);
  print, "k= ", k, "l= ", l, "i= ",i, "j= ", "ret= ", ret;
  return ret;
}

func dvv_DrawFrame(i, type=, lgscl=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _currentSet;
  extern _Xcoord, _Ycoord;
  extern _Xrange, _Yrange;
  extern _Zmin, _Zmax;
  extern _FrameCount, _FrameCountMax;
  extern _shot, _warp_pts;

  if(is_void(_shot)) _shot= "";
  _FrameCount++;
  ii= (i-1)%4 + 1;
  if(is_void(type)) {
    if(ii==1) {
      sh, ,_currentSet.ch1S, _Zmin, _Zmax, notitle= 1, lgscl= _lgscl;
      pltitle, swrite(format= "%s:ch1S (%4.0f, %4.0f)", _shot, double(_Xcoord), double(_Ycoord));
    }
    if(ii==3) {
      sh, ,_currentSet.ch2S, _Zmin, _Zmax, notitle= 1,  lgscl= _lgscl;
      pltitle, swrite(format= "%s:ch2S (%4.0f, %4.0f)", _shot, double(_Xcoord), double(_Ycoord));
    }
    if(ii==2) {
      sh, ,_currentSet.ch1P, _Zmin, _Zmax, notitle= 1, lgscl= _lgscl;
      pltitle, swrite(format= "%s:ch1P (%4.0f, %4.0f)", _shot, double(_Xcoord), double(_Ycoord));
    }
    if(ii==4) {
      sh, ,_currentSet.ch2P, _Zmin, _Zmax, notitle= 1, lgscl= _lgscl;
      pltitle, swrite(format= "%s:ch2P (%4.0f, %4.0f)", _shot, double(_Xcoord), double(_Ycoord));
    } 
    if(_warp_pts) plmk, _Ycoord, _Xcoord, marker= 4, msize= 0.5, width= 10, color= "red";
  } else if (type == "dif") {
    if(ii==1) {
      sh, , img_sub(_currentSet.ch1S, _currentSet.ch2S), _Zmin, _Zmax, notitle= 1,  lgscl= _lgscl;
      pltitle, "difference: ch1S, ch2S";
    }
    if(ii==3) {
      sh, , img_sub(_currentSet.ch1S, _currentSet.ch2P), _Zmin, _Zmax, notitle= 1,  lgscl= _lgscl;
      pltitle, "difference: ch1S, ch2P";
    }
    if(ii==2) {
      sh, , img_sub(_currentSet.ch1P, _currentSet.ch2P), _Zmin, _Zmax, notitle= 1,  lgscl= _lgscl;
      pltitle, "difference: ch1P, ch2P";
    }
    if(ii==4) {
      sh, , img_sub(_currentSet.ch1P, _currentSet.ch2S), _Zmin, _Zmax, notitle= 1, lgscl= _lgscl;
      pltitle, "difference: ch1P, ch2S";
    }
  }
  limits, _Xcoord - _Xrange/2., _Xcoord + _Xrange/2.0,
    _Ycoord - _Yrange/2., _Ycoord + _Yrange/2.0;
  if(_FrameCount > _FrameCountMax) return 0;
  else return 1;
}

func dvv_ShowLineouts(fs, lo, w=, dofft=, lgscl=, rge=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  crvs=[];
  if(!is_void(w)) window, w;
  nl= numberof(lo);
  if(is_void(dofft)) {
    for(i= 1; i<= nl; i++) {
      t= (i-1)%5 + 1;
      crv1s= lineout1(fs.ch1S, lo(i));
      crv2s= lineout1(fs.ch2S, lo(i));
      crv1p= lineout1(fs.ch1P, lo(i));
      crv2p= lineout1(fs.ch2P, lo(i));
      plg, crv1s(,2), crv1s(,1), color= "red", width= 2, type= t;
      plg, crv2s(,2), crv2s(,1), color= "blue", width= 2, type= t;
      plg, crv1p(,2), crv1p(,1), color= "green", width= 2, type= t;
      plg, crv2p(,2), crv2p(,1), color= "magenta", width= 2, type= t;
      crvs= _(crvs, [crv1s, crv2s, crv1p, crv2p]);
    }
    dims= dimsof(crvs);
    dims(1)+=1; dims(0)/=2; dims= _(dims, 2);
    crvs= reform(crvs, dims);
  } else {
    ft1s= dvv_PSD(fs.ch1S, wndw= "hanning");
    ft2s= dvv_PSD(fs.ch2S, wndw= "hanning");
    ft1p= dvv_PSD(fs.ch1P, wndw= "hanning");
    ft2p= dvv_PSD(fs.ch2P, wndw= "hanning");
    for(i= 1; i<= nl; i++) {
      t= (i-1)%5 + 1;
      print, "i= ", i, "lineout = ", lo(i)
      crv1s= lineout1(ft1s, lo(i));
      crv2s= lineout1(ft2s, lo(i));
      crv1p= lineout1(ft1p, lo(i));
      crv2p= lineout1(ft2p, lo(i));
      plg, crv1s(,2), crv1s(,1), color= "red", width= 2, type= t;
      plg, crv2s(,2), crv2s(,1), color= "blue", width= 2, type= t;
      plg, crv1p(,2), crv1p(,1), color= "green", width= 2, type= t;
      plg, crv2p(,2), crv2p(,1), color= "magenta", width= 2, type= t;
      crvs= _(crvs, [crv1s, crv2s, crv1p, crv2p]);
    }
    dims= dimsof(crvs);
    dims(1)+=1; dims(0)/=2; dims= _(dims, 2);
    crvs= reform(crvs, dims);
  }
  if (!is_void(rge)) range, rge(1), rge(2);
  if (lgscl) logxy, 0, 1;
  return crvs;
}

func dvv_FocusLineout(file, pos, noff=, dwmethod=, dofft=, lgscl=, w=, ssdewarp=, doref=)
/* DOCUMENT dvv_FocusLineouts, file, pos, noff=, nowarp=, dofft=, lgscl=, w=. ssdewarp=, doref=
     
   SEE ALSO:
 */
{
  extern _lh, _lv;
  //if(!is_void(doref)) dvv_SetReference, hdf_SDget(file);
  if(!is_void(doref)) dvv_SetMapping, hdf_SDget(file);
  gr= dvv_ImportDataSet(file, noff= noff, dwmethod= dwmethod, ssdewarp= ssdewarp, nobkgsub= 1);
  crvs= dvv_ShowLineouts(gr, [_lh, _lv], dofft= dofft, rge= rge, lgscl= lgscl, w= w);
  pltitle, file;
  ff= dvv_FocusSet(file= file);
  ff.pos= pos;
  ff.ch1Sv= &crvs(,,1,2);
  ff.ch2Sv= &crvs(,,2,2);
  ff.ch1Pv= &crvs(,,3,2);
  ff.ch2Pv= &crvs(,,4,2);
  
  ff.ch1Sh= &crvs(,,1,1);
  ff.ch2Sh= &crvs(,,2,1);
  ff.ch1Ph= &crvs(,,3,1);
  ff.ch2Ph= &crvs(,,4,1);
  return ff;
}

func dvv_plFocus(focus, freq)
/* DOCUMENT dvv_plFocus, focus, freq


   SEE ALSO:
 */
{
  print, "dvv_plFocus ... ";
  nf= numberof(focus);
  nq= numberof(freq);
  print, "nf, nq=, freq=", nf, nq, freq;
  ch1sv= ch2sv= ch1pv= ch2pv= array(0.0, [2, nf, nq]);
  ch1sh= ch2sh= ch1ph= ch2ph= array(0.0, [2, nf, nq]);
  s= sort(focus.pos);
  for(i= 1; i<= nf; i++) {
    print, "interp=", interp((*focus(s(i)).ch1Sv)(,2), (*focus(s(i)).ch1Sv)(,1), freq);
    ch1sv(i,)= interp((*focus(s(i)).ch1Sv)(,2), (*focus(s(i)).ch1Sv)(,1), freq)/
      interp((*focus(s(i)).ch1Sv)(,2), (*focus(s(i)).ch1Sv)(,1), 0.0);
    ch2sv(i,)= interp((*focus(s(i)).ch2Sv)(,2), (*focus(s(i)).ch2Sv)(,1), freq)/
      interp((*focus(s(i)).ch2Sv)(,2), (*focus(s(i)).ch2Sv)(,1), 0.0);
    ch1pv(i,)= interp((*focus(s(i)).ch1Pv)(,2), (*focus(s(i)).ch1Pv)(,1), freq)/
      interp((*focus(s(i)).ch1Pv)(,2), (*focus(s(i)).ch1Pv)(,1), 0.0);
    ch2pv(i,)= interp((*focus(s(i)).ch2Pv)(,2), (*focus(s(i)).ch2Pv)(,1), freq)/
      interp((*focus(s(i)).ch2Pv)(,2), (*focus(s(i)).ch2Pv)(,1), 0.0);
    ch1sh(i,)= interp((*focus(s(i)).ch1Sh)(,2), (*focus(s(i)).ch1Sh)(,1), freq)/
      interp((*focus(s(i)).ch1Sh)(,2), (*focus(s(i)).ch1Sh)(,1), 0.0);
    ch2sh(i,)= interp((*focus(s(i)).ch2Sh)(,2), (*focus(s(i)).ch2Sh)(,1), freq)/
      interp((*focus(s(i)).ch2Sh)(,2), (*focus(s(i)).ch2Sh)(,1), 0.0);
    ch1ph(i,)= interp((*focus(s(i)).ch1Ph)(,2), (*focus(s(i)).ch1Ph)(,1), freq)/
      interp((*focus(s(i)).ch1Ph)(,2), (*focus(s(i)).ch1Ph)(,1), 0.0);
    ch2ph(i,)= interp((*focus(s(i)).ch2Ph)(,2), (*focus(s(i)).ch2Ph)(,1), freq)/
      interp((*focus(s(i)).ch2Ph)(,2), (*focus(s(i)).ch2Ph)(,1), 0.0);
  }
  fp= span(-2.25, 2.25, 200);
  for(i= 1; i<= nq; i++) {
    window, i-1, style= "rect1x2vg.gs";
    plsys, 1;
    pltitle, swrite(format="Frequency: %4.2f", freq(i));
    plmk, ch1sv(,i), focus(s).pos, color= "red", marker= i, width= 10, msize= 0.5;
    plmk, ch2sv(,i), focus(s).pos, color= "blue", marker= i, width= 10, msize= 0.5;
    plmk, ch1pv(,i), focus(s).pos, color= "green", marker= i, width= 10, msize= 0.5;
    plmk, ch2pv(,i), focus(s).pos, color= "magenta", marker= i, width= 10, msize= 0.5;
    plsys, 2;
    plmk, ch1sh(,i), focus(s).pos, color= "red", marker= i+1, width= 10, msize= 0.5;
    plmk, ch2sh(,i), focus(s).pos, color= "blue", marker= i+1, width= 10, msize= 0.5;
    plmk, ch1ph(,i), focus(s).pos, color= "green", marker= i+1, width= 10, msize= 0.5;
    plmk, ch2ph(,i), focus(s).pos, color= "magenta", marker= i+1, width= 10, msize= 0.5;
    //logxy, 0, 1;
    plsys, 1;
    plg, spline(ch1sv(,i), focus(s).pos, fp), fp, color= "red";
    plg, spline(ch2sv(,i), focus(s).pos, fp), fp, color= "blue";
    plg, spline(ch1pv(,i), focus(s).pos, fp), fp, color= "green";
    plg, spline(ch2pv(,i), focus(s).pos, fp), fp, color= "magenta";
    plsys, 2;
    plg, spline(ch1sh(,i), focus(s).pos, fp), fp, color= "red", type= 2;
    plg, spline(ch2sh(,i), focus(s).pos, fp), fp, color= "blue", type= 2;
    plg, spline(ch1ph(,i), focus(s).pos, fp), fp, color= "green", type= 2;
    plg, spline(ch2ph(,i), focus(s).pos, fp), fp, color= "magenta", type= 2;
  }
}

func dvv_mkFocusSpect(focus, ch, or)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  np= numberof(focus);
  nf= numberof((*focus(1).ch1Sv)(,1));
  yg=(*focus(1).ch1Sv)(,1);
  xsc= span(0.0, 1.17, 118);
  zg= xsc(,-:1:np);
  fsi= IMG_DAT();
  s= sort(focus.pos);
  ysc= focus(s).pos;
  if (or == 1 || or == "h") {
    if (ch == "ch1S" || ch == 1) {
      fsi.shotid= "ch1S, horizontal";
      for (i= 1; i<= np; i++)
        zg(,i)= interp((*focus(s(i)).ch1Sh)(,2)/(*focus(s(i)).ch1Sh)(*)(avg),yg,xsc);
    }
    if (ch == "ch2S" || ch == 2)  {
      fsi.shotid= "ch2S, horizontal";
      for (i= 1; i<= np; i++)
        zg(,i)= interp((*focus(s(i)).ch2Sh)(,2)/(*focus(s(i)).ch2Sh)(*)(avg),yg,xsc);
    }
    if (ch == "ch1P" || ch == 3)  {
      fsi.shotid= "ch1P, horizontal";
      for (i= 1; i<= np; i++)
        zg(,i)= interp((*focus(s(i)).ch1Ph)(,2)/(*focus(s(i)).ch1Ph)(*)(avg),yg,xsc);
    }
    if (ch == "ch2P" || ch == 4)  {
      fsi.shotid= "ch2P, horizontal";
      for (i= 1; i<= np; i++)
        zg(,i)= interp((*focus(s(i)).ch2Ph)(,2)/(*focus(s(i)).ch2Ph)(*)(avg),yg,xsc);
    }
  } else if (or == 2 || or == "v") {
    if (ch == "ch1S" || ch == 1)  {
      fsi.shotid= "ch1S, vertical";
      for (i= 1; i<= np; i++)
        zg(,i)= interp((*focus(s(i)).ch1Sv)(,2)/(*focus(s(i)).ch1Sv)(*)(avg),yg,xsc);
    }
    if (ch == "ch2S" || ch == 2)  {
      fsi.shotid= "ch2S, vertical";
      for (i= 1; i<= np; i++)
        zg(,i)= interp((*focus(s(i)).ch2Sv)(,2)/(*focus(s(i)).ch2Sv)(*)(avg),yg,xsc);
    }
    if (ch == "ch1P" || ch == 3)  {
      fsi.shotid= "ch1P, vertical";
      for (i= 1; i<= np; i++)
        zg(,i)= interp((*focus(s(i)).ch1Pv)(,2)/(*focus(s(i)).ch1Pv)(*)(avg),yg,xsc);
    }
    if (ch == "ch2P" || ch == 4)  {
      fsi.shotid= "ch2P, vertical";
      for (i= 1; i<= np; i++)
        zg(,i)= interp((*focus(s(i)).ch2Pv)(,2)/(*focus(s(i)).ch2Pv)(*)(avg),yg,xsc);
    }
  }
  fsi.nx= 118;
  fsi.ny= np;
  fsi.xscale= &xsc;
  fsi.yscale= &ysc;
  fsi.x_unit= "!mm^-1^";
  fsi.y_unit= "mm";
  fsi.x_label= "Focus position";
  fsi.y_label= "Spatial Frequency";
  fsi.z_label= "";
  fsi.z_unit= "";
  fsi.data= &zg;
  print, fsi.shotid;
  return fsi;
}

func dvv_ShowFocusSpect(focus, or, w=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  if(is_void(or)) or= 1;
  if(is_void(w)) w= 0;
 
  window, w, style= "square2x2vg.gs", width= 600, height= 700, legends= 0;
  for(i= 1; i<= 4; i++) {
    sh, w, dvv_mkFocusSpect(focus, i, or), 1e-5, 1.0, nolabels= 1,
      notitle= 1, lgscl= 1, plsy= i, c= 1, levs=[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0];
    limits, 0.0, 1.1;
    palette, "rainbow.gp";
  }
  if(or == 1) plt, "Horizontal gridlines", 0.4, 1.0, height= 18, justify="CC";
  else plt, "Vertical gridlines", 0.4, 1.0, height= 18, justify= "CC";
  plt, "ch1S", 0.22, 0.965, justify="CC", tosys= 0, height= 14;
  plt, "ch2S", 0.59, 0.965, justify="CC", tosys= 0, height= 14;
  plt, "ch1P", 0.22, 0.585, justify="CC", tosys= 0, height= 14;
  plt, "ch2P", 0.59, 0.585, justify="CC", tosys= 0, height= 14;
}

func dvv_ShowWarpSet(fs, fr, cmin, cmax, w=, nocx=)
/* DOCUMENT dvv_ShowWarpSet, fs, frame, w=, nocx=

   Displays the set of warp/dewarping data points used to
   generate the dewarping map for the data set

   fs - frame set containing the data
   fr - ={1,2,3 or 4}, identifies which of the 4 frames to show

   KEYOWRDS:
     w=  display on this window
     nocx= define this to show the raw data instead of the convolved data
     
   SEE ALSO:
 */
{
  if(is_void(w)) w= 0;
  cx= dvv_mkCross(fs.ch1S, 2.1, 40.0);
  window, w;
  if(fr == 1) {
    if(is_void(nocx)) sh, w, dvv_Correlate(fs.ch1S, cx);
    else sh, w, fs.ch1S, cmin, cmax;
    pts= *fs.map.ch1S.warp_pts;
  } else if(fr == 2) {
    if(is_void(nocx)) sh, w, dvv_Correlate(fs.ch1P, cx);
    else sh, w, fs.ch1P, cmin, cmax;
    pts= *fs.map.ch1P.warp_pts;
  } else if(fr == 3) {
    if(is_void(nocx)) sh, w, dvv_Correlate(fs.ch2S, cx);
    else sh, w, fs.ch2S, cmin, cmax;
    pts= *fs.map.ch2S.warp_pts;
  } else if(fr == 4) {
    if(is_void(nocx)) sh, w, dvv_Correlate(fs.ch2P, cx);
    else sh, w, fs.ch2P, cmin, cmax;
    pts= *fs.map.ch2P.warp_pts;
  } else error, "Invalid frame id.";
  dx= (*fs.ch1S.xscale)(dif)(avg);
  dy= (*fs.ch1S.yscale)(dif)(avg);
  plmk, pts.y, pts.x,
    msize= 0.5, color= "blue", width= 4, marker= 4;
}

func dvv_ShowWarpDeviations(fs, fr, w=, dev=, scale=, wfs=, pwarp=, nodc=)
/* DOCUMENT dvv_ShowWarpDeviations, fs, fr, w=, dev=, scale=

   Shows plots of the warp deviations of the grid data in each
   frame relative to the nominal grid values.

   ARGUMENTS:
     fs - frame set containing the warp mapping
     fr - select frame <1, 2, 3, or 4> for <ch1S, ch1P, ch2S, ch2P> frame

   KEYWORDS:
     w= specify window to plot to (plots to w and w+1)
     dev=
     scale=
     wfs=
     pwarp=
     nodc=
     
   SEE ALSO:
 */
{
  extern _colorset, _nclrs;
  extern _xpdev, _ypdev;
  extern xdev, ydev;

  olX= olegend(new,);
  olY= olegend(new,);

  if (is_void(wfs)) wfs= next_window();
  wfsk= window_exists(wfs);
  if (!wfsk) window, wfs;
  if (is_void(w)) {
    w= next_window();
    if (!wfsk) winkill, wfs;
  }
  if(is_void(scale)) scale= 100;
  ref= *fs.map.warp_refpts;
  window, w;
  if(fr == 1) {
    frm= "ch1S";
    pts= *fs.map.ch1S.warp_pts;
    if (pwarp) dewarp, fs.ch1S, *fs.map.ch1S.kx, *fs.map.ch1S.ky, shdev= 1, dev= dev;
  } else if(fr == 2) {
    frm= "ch2S";
    pts= *fs.map.ch2S.warp_pts;
    if (pwarp) dewarp, fs.ch2S, *fs.map.ch2S.kx, *fs.map.ch2S.ky, shdev= 1, dev= dev;
  } else if(fr == 3) {
    frm= "ch1P";
    pts= *fs.map.ch1P.warp_pts;
    if (pwarp) dewarp, fs.ch1P, *fs.map.ch1P.kx, *fs.map.ch1P.ky, shdev= 1, dev= dev;
  } else if(fr == 4) {
    frm= "ch2P";
    pts= *fs.map.ch2P.warp_pts;
    if (pwarp) dewarp, fs.ch2P, *fs.map.ch2P.kx, *fs.map.ch2P.ky, shdev= 1, dev= dev;
  } else error, "Invalid frame id.";

  ng= int(sqrt(numberof(ref)));
  cc= span(-4, 4, ng)*100.0;
  wout= where(abs(ref.y - pts.y) > 2.0 | abs(ref.x - pts.x) > 2.0);
  if(numberof(wout) > 0) {
    print, "Outliers: ";
    print, "==========";
    for(j= 1; j<= numberof(wout); j++)
      print, ref(wout(j)).x, "-->", pts(wout(j)).x,
        ref(wout(j)).y, "-->", pts(wout(j)).y;
  }

  xscl= (*fs.ch1S.xscale);
  yscl= (*fs.ch1S.yscale);
  
  if (pwarp) {
    //xdev= tl2cub(,xscl,yscl,_xpdev);
    //ydev= tl2cub(,xscl,yscl,_ypdev);
    xdev= img_copy(fs.ch1S, data= _xpdev);
    ydev= img_copy(fs.ch1S, data= _ypdev);
    if (dev) {
      xdev= img_sub(xdev, img_copy(xdev, data= img_grid(xdev, 1)));
      ydev= img_sub(ydev, img_copy(ydev, data= img_grid(ydev, 2)));
    }
  }
  
  for(i= 1; i<= ng; i++) {
    clr= _colorset(,i%_nclrs);
    mkr= (i-1)%7+1;

    ww= where(ref.y == cc(i));
    window, w;
    if(numberof(ww) > 0)
      plmk, (pts(ww).x - ref(ww).x), ref(ww).x, msize= 0.5, width= 11, color= clr, marker= mkr;
    if (pwarp) {
      //xdev_warp= tl2cub(xdev,xscl,cc(i));
      ll= lineout1(xdev, LO_POS(or= "h", x1= cc(i), x2= cc(i)));
      plg, ll(,2), ll(,1), color= clr, type= (i-1)%5+1, width= 3;
    } else {
      if(numberof(ww) > 1)
        plg, pts(ww).x - ref(ww).x, ref(ww).x, color= clr, type= (i-1)%5+1;
    }
    
    window, w+2;
    if(numberof(ww) > 0)
      plmk, (pts(ww).y - ref(ww).y), ref(ww).x, msize= 0.5, width= 11, color= clr, marker= mkr;
    if (pwarp) {
      //xdev_warp= tl2cub(xdev,xscl,cc(i));
      //plg, xdev_warp, xscl, color= clr, type= (i-1)%5+1, width= 3;
      ll= lineout1(ydev, LO_POS(or= "h", x1= cc(i), x2= cc(i)));
      plg, ll(,2), ll(,1), color= clr, type= (i-1)%5+1, width= 3;
    } else {
      if(numberof(ww) > 1)
        plg, pts(ww).y - ref(ww).y, ref(ww).x, color= clr, type= (i-1)%5+1;
    }
    olY, add, swrite(format="Y = %4.0f", cc(i)), clr= clr, type= (i-1)%5+1, width= 3, msize= 0.5, mwidth= 11, marker= mkr;

    ww= where(ref.x == cc(i));
    window, w+1;
    if(numberof(ww) > 0)
      plmk, (pts(ww).y - ref(ww).y), ref(ww).y, msize= 0.5, width= 11, color= clr, marker= mkr;
    if (pwarp) {
      //ydev_warp= tl2cub(ydev,xscl,cc(i));
      //plg, ydev_warp, xscl, color= clr, type= (i-1)%5+1, width= 3;
      ll= lineout1(ydev, LO_POS(or= "v", x1= cc(i), x2= cc(i)));
      plg, ll(,2), ll(,1), color= clr, type= (i-1)%5+1, width= 3;
    } else {
      if(numberof(ww) > 1)
        plg, pts(ww).y - ref(ww).y, ref(ww).y, color= clr, type= (i-1)%5+1;
    }

    window, w+3;
    if(numberof(ww) > 0)
      plmk, (pts(ww).x - ref(ww).x), ref(ww).y, msize= 0.5, width= 11, color= clr, marker= mkr;
    if (pwarp) {
      //ydev_warp= tl2cub(ydev,xscl,cc(i));
      //plg, ydev_warp, xscl, color= clr, type= (i-1)%5+1, width= 3;
      ll= lineout1(xdev, LO_POS(or= "v", x1= cc(i), x2= cc(i)));
      plg, ll(,2), ll(,1), color= clr, type= (i-1)%5+1, width= 3;
    } else {
      if(numberof(ww) > 1)
        plg, pts(ww).x - ref(ww).x, ref(ww).y, color= clr, type= (i-1)%5+1;
    }
    olX, add, swrite(format="X = %4.0f", cc(i)), clr= clr, type= (i-1)%5+1, width= 3, msize= 0.5, mwidth= 11, marker= mkr;
  }
  
  window, w;
  olY, set, 0.23, 0.33;
  olY, plot, w;
  pltitle, frm+": X-deviations versus X";
  xytitles, "X-position (!mm)", "Deviation (!mm)";
  aa= (pts.x - ref.x)(avg);
  range, aa-1, aa+1;
  //range, -1, 1;
  wpdf, w, shot+"-"+frm+"-Xdev_v_X";
  
  window, w+2;
  olY, set, 0.23, 0.33;
  olY, plot, w+2;
  pltitle, frm+": Y-deviations versus X";
  xytitles, "X-position (!mm)", "Deviation (!mm)";
  aa= (pts.y - ref.y)(avg);
  range, aa-1, aa+1;
  // range, -1, 1;
  wpdf, w+2, shot+"-"+frm+"-Ydev_v_X";
  
  window, w+1;
  olX, set, 0.23, 0.33;
  olX, plot, w+1;
  pltitle, frm+": Y-deviations versus Y";
  xytitles, "Y-position (!mm)", "Deviation (!mm)";
  aa= (pts.y - ref.y)(avg);
  range, aa-1, aa+1;
  //range, -1, 1;
  wpdf, w+1, shot+"-"+frm+"-Ydev_v_Y";

  window, w+3;
  olX, set, 0.23, 0.33;
  olX, plot, w+3;
  pltitle, frm+": X-deviations versus Y";
  xytitles, "Y-position (!mm)", "Deviation (!mm)";
  aa= (pts.x - ref.x)(avg);
  range, aa-1, aa+1;
  //range, -1, 1;
  wpdf, w+3, shot+"-"+frm+"-Xdev_v_Y";

  if (!window_exists(wfs)) window, wfs, legends= 0, style= "square2x2img.gs";
  window, wfs;
  plsys, fr;
  rr= reform(ref, [2,ng,ng]);
  pt= reform(pts, [2,ng,ng]);
  u= pt.x - rr.x;
  if (nodc) u-= u(*)(avg);
  v= pt.y - rr.y;
  if (nodc) v-= v(*)(avg);
  if (pwarp) plv, img_zvalue(ydev, rr) - _ypdev(*)(avg), img_zvalue(xdev, rr) - _xpdev(*)(avg), rr.y, rr.x, scale= 150, color= "red";
  plv, v, u, rr.y, rr.x, scale= scale;
  limits, -500, 500, -500, 500;
  if (fr == 4)  wpdf, wfs, shot+"-all-frames-dev";
}

func dvv_DisplayImg(img, zmin, zmax, lh=, lv=, w=, lgscl=, pltrms=, norm=, plsy=, style=,
                    pltitl=, box=, c=, levs=, rmslims=, rmsunit=, nodc=, shlo=, prlo=, regs=,
                    noresample=)
/* DOCUMENT dvv_DisplayImg, img, zmin, zmax, lh=, lv=, w=, lgscl=, pltrms=, plsy=, style=,
                    pltitl=, box=, c=, levs=, rmslims=, nodc=, shlo=, prlo=, regs=, noresample=
     
  SEE ALSO:
 */
{
  extern _DVV_RESAMPLE_DISPLAY;
  extern _lh, _lv;
  extern _PSF, _REFSUB;
  extern _VPF, _TAU;
  extern loh, lov;

  if(is_void(lh) && !is_void(_lh)) lh= _lh;
  if(is_void(lv) && !is_void(_lv)) lv= _lv;
  if(is_void(plsy)) plsy= 1;
  if(is_void(style)) style= "square+lo.gs";
  if(style == "square+lo.gs") {
    width= 600;
    height= 650;
  } else if (style == "beam_display.gs") {
    width= 800;
    height= 650;
  }

  if (_DVV_RESAMPLE_DISPLAY && !noresample) img= img_resample(img, _DVV_RESAMPLE_DISPLAY);

  x1= (*img.xscale)(1);
  x0= (*img.xscale)(0);
  y1= (*img.yscale)(1);
  y0= (*img.yscale)(0);

  if(is_void(w)) {w= 0; window, w, style= style, width= width, height= height, legends= 0;}
  else if (plsy > 1) window, w;
  else window, w, style= style, width= width, height= height, legends= 0;
  
  if(!is_void(box)) {
    print, "extracting data within box:", box;
    img2= img_extract(img, box);
  } else {
    img2= img_copy(img);
  }
  
  img1= img_copy(img);
  if (nodc) (*img1.data)-=(*img2.data)(*)(avg);
  
  if (style == "beam_display.gs") {
    sh, w, img1, zmin, zmax, lgscl= lgscl, notitle= 1, c= c,
      levs= levs, nolabels= 1, plsy= plsy;
    xytitles, img.x_label+"("+img.x_unit+")", "";
  } else {
    sh, w, img1, zmin, zmax, lgscl= lgscl, notitle= 1, c= c,
      levs= levs, plsy= plsy;
  }

  
  if(!is_void(lh)) {
    plsys, plsy+1;
    nh= numberof(lh);
    for(i= 1; i<= nh; i++) {
      loh= lineout1(img1, lh(i));
      plg, loh(,2), loh(,1), marks= 0, width= 2, color= _colorset(,(i-1)%_nclrs + 1), type= i;
      plt, swrite(format="Horizontal lineout at y=%3.1f",[lh(i).x1,lh(i).x2](avg)),
        0.099, 0.944, font= "courier", height= 10, tosys= 0;
      if(!is_void(shlo)) {
        plsys, 1;
        plg, [lh(i).x1, lh(i).x1], [x1,x0], color= "red", marks= 0;
        plg, [lh(i).x2, lh(i).x2], [x1,x0], color= "red", marks= 0;
        plsys, 2;
      }
      if(!is_void(prlo)) {
        fname= swrite(format="%s-hlo.%d.txt", prlo, i);
        fh= open(fname, "w");
        write, fh, format="# %s %d\n", prlo+" - Horizontal lineout number", i;
        write, fh, format="# VPF = %6.4f km/s/fringe, TAU = %7.3f \n", _VPF, _TAU;
        write, fh, format="# %s(%s) %s(%s) phase(mrad)\n",
          img.x_label, img.x_unit, img.z_label, img.z_unit;
        for(j= 1; j<= dimsof(loh)(2); j++) write, fh, loh(j,1), loh(j,2), loh(j,2)/_VPF*2*pi;
        close, fh;
      }
    }
    if(!(is_void(zmin) || is_void(zmax))) range, zmin, zmax;
    if(!is_void(lgscl)) logxy, 0, 1;
    if (plsy < 4) xytitles, "", img.z_label+"("+img.z_unit+")";
  }
  
  if(!is_void(lv)) {
    plsys, plsy+2;
    nv= numberof(lv);
    for(i= 1; i<= nv; i++) {
      lov= lineout1(img1, lv(i));
      plg, lov(,1), lov(,2), marks= 0, width= 2, color= _colorset(,(i-1)%_nclrs + 1), type= i;
      plt, swrite(format="Vertical lineout at x=%3.1f",[lv(i).x1,lv(i).x2](avg)),
        0.52, 0.73, font= "courier", height= 10, tosys= 0;
      if(!is_void(shlo)) {
        plsys, 1;
        plg, [y1,y0], [lv(i).x1, lv(i).x1], color= "red", marks= 0;
        plg, [y1,y0], [lv(i).x2, lv(i).x2], color= "red", marks= 0;
        plsys, 3;
      }
      if(!is_void(prlo)) {
        fname= swrite(format="%s-vlo.%d.txt", prlo, i);
        fh= open(fname, "w");
        write, fh, format="# %s %d\n", prlo+" - Vertical lineout number", i;
        write, fh, format="# VPF = %6.4f km/s/fringe, TAU = %7.3f \n", _VPF, _TAU;
        write, fh, format="# %s(%s) %s(%s) phase(mrad)\n",
          img.x_label, img.x_unit, img.z_label, img.z_unit;
        for(j= 1; j<= dimsof(lov)(2); j++) write, fh, lov(j,1), lov(j,2), lov(j,2)/_VPF*2*pi;
        close, fh;
      }
    }
    if(!(is_void(zmin) || is_void(zmax))) limits, zmin, zmax;
    if(!is_void(lgscl)) logxy, 1, 0;
    xytitles, img.z_label+"("+img.z_unit+")", "";
  }

  //Display region boxes
  if(!is_void(regs)) {
    plsys, 1;
    for (jj= 1; jj<= numberof(regs); jj++) plbox, regs(jj), color= "magenta", width=3;
  }

  //Display the integration box
  if(!is_void(box)) {
    plsys, 1;
    plg, [box.y2, box.y1, box.y1, box.y2, box.y2],
      [box.x2, box.x2, box.x1, box.x1, box.x2], width= 3, color= "red", marks= 0;
  }

  if (!is_void(pltrms)) {
    if(is_void(rmsunit)) rmsunit= img.z_unit;
    if(img.z_unit == "counts") fmt= "%8.0f";
    else if (img.z_unit == "m/s") fmt= "%6.2f";
    else fmt= "%12.5e";
    if(pltrms == 1) {
      if(is_void(rmslims)) {
        rmsval= dvv_rms(img2);
        if (!nodc) plt, swrite(format="Average: "+fmt+" %s", (*img2.data)(*)(avg), img.z_unit),
                     0.52, 0.92, tosys= 0, font= "courier", height= 12;
        if (box) plt, swrite(format="RMS: "+fmt+" %s within box.", rmsval, rmsunit),
                   0.52, 0.90, tosys= 0, font= "courier", height= 12;
        else plt, swrite(format="RMS: "+fmt+" %s", rmsval, rmsunit),
               0.52, 0.90, tosys= 0, font= "courier", height= 12;
        plt, swrite(format="%s", "Freq: all"),
          0.52, 0.88, tosys= 0, font= "courier", height= 12;
      } else {
        rmsval= dvv_rms(img2, rmslims);
        plt, swrite(format="RMS: %9.3e %s", rmsval, rmsunit),
          0.52, 0.92, tosys= 0, font= "courier", height= 12;
        plt, swrite(format="Freq: %4.2f-%4.2f %s", rmslims(1), rmslims(2), img.x_unit),
          0.52, 0.90, tosys= 0, font= "courier", height= 12;
      }
    
    } else if(pltrms == "psd-amplitude" || pltrms == "amplitude") {
      if(is_void(rmslims)) {
        rmsval= dvv_AMP_rms(img2, norm= norm);
        plt, swrite(format="RMS: %9.3e %s", rmsval, rmsunit),
          0.52, 0.92, tosys= 0, font= "courier", height= 12;
        plt, swrite(format="%s", "Freq: all"),
          0.52, 0.90, tosys= 0, font= "courier", height= 12;
      } else {
        rmsval= dvv_AMP_rms(img2, rmslims, norm= norm);
        plt, swrite(format="RMS: %9.3e %s", rmsval, rmsunit),
          0.52, 0.92, tosys= 0, font= "courier", height= 12;
        plt, swrite(format="Freq: %4.2f-%4.2f %s", rmslims(1), rmslims(2), img.x_unit),
          0.52, 0.90, tosys= 0, font= "courier", height= 12;
      }
  
    } else if(pltrms == "psd-power" || pltrms == "power") {
      if(is_void(rmslims)) {
        rmsval= dvv_PSD_rms(img2, norm= norm);
        plt, swrite(format="RMS: %9.3e %s", rmsval, rmsunit),
          0.52, 0.92, tosys= 0, font= "courier", height= 12;
        plt, swrite(format="%s", "Freq: all"),
          0.52, 0.90, tosys= 0, font= "courier", height= 12;
      } else {
        rmsval= dvv_PSD_rms(img2, rmslims, norm= norm);
        plt, swrite(format="RMS: %9.3e %s", rmsval, rmsunit),
          0.52, 0.92, tosys= 0, font= "courier", height= 12;
        plt, swrite(format="Freq: %4.2f-%4.2f %s", rmslims(1), rmslims(2), img.x_unit),
          0.52, 0.90, tosys= 0, font= "courier", height= 12;
      }
    }
  
    //     } else if(pltrms == "dft-ampl") {
    //       if(is_void(rmslims)) {
    //         rmsval= sqrt(((*img.data)^2)(*)(sum) - ((*img.data)(img.nx/2+1,img.ny/2+1))^2)/2.0;
    //         plt, swrite(format="RMS: %9.3e %s", rmsval, rmsunit),
    //           0.52, 0.92, tosys= 0, font= "courier", height= 12;
    //         plt, swrite(format="%s", "Freq: all"),
    //           0.52, 0.90, tosys= 0, font= "courier", height= 12;
    //       } else {
    //         xg= (*img.xscale)(,-:1:img.ny);
    //         yg= (*img.yscale)(-:1:img.nx,);
    //         r= sqrt(xg*xg + yg*yg);
    //         w= where(r >= rmslims(1) & r < rmslims(2));
    //         rmsval= sqrt(((*img.data)^2)(w)(sum))/2.0;
    //         plt, swrite(format="RMS: %9.3e %s", rmsval, rmsunit),
    //           0.52, 0.92, tosys= 0, font= "courier", height= 12;
    //         plt, swrite(format="Freq: %4.2f-%4.2f %s", rmslims(1), rmslims(2), img.x_unit),
    //           0.52, 0.90, tosys= 0, font= "courier", height= 12;
    
  }
  if(!is_void(pltitl)) plt, pltitl, 0.099, 0.98, tosys= 0, height=18;
}

func dvv_DisplaySet(fs, zmin, zmax, w=)
/* DOCUMENT dvv_DisplaySet, fs, zmin, zmax, w=

   Displays all four frames of a data set in a 2x2 grid
     
   SEE ALSO:
 */
{
  if (is_void(w)) w= next_window();
  if (is_void(zmin)) zmin= 0.9;
  if (is_void(zmax)) zmax= 1.1;
  window, w, legends= 0, style= "square2x2img.gs", width= 525, height=525 ;
  sh, w, fs.ch1S, zmin, zmax, plsy= 1, notitle= 1, nolabels= 1;
  sh, w, fs.ch2S, zmin, zmax, plsy= 2, notitle= 1, nolabels= 1;
  sh, w, fs.ch1P, zmin, zmax,  plsy= 3, notitle= 1, nolabels= 1;
  sh, w, fs.ch2P, zmin, zmax,  plsy= 4, notitle= 1, nolabels= 1;
}

func dvv_DisplayPhaseMapping(fs, zmin, zmax, w=)
{
  if (is_void(w)) w= 0;
  if (is_void(shot)) shot= "";
  window, w, legends= 0, style= "square2x2img.gs", height= 650, width= 600;
  get_style, ls, sy, lg, clg;
  //ls= 0;
  sy(1).ticks.horiz.flags= 0x003 + 0x008;
  sy(1).ticks.vert.flags= 0x003 + 0x008;
  sy(2).ticks.horiz.flags= 0x003 + 0x008;
  sy(2).ticks.vert.flags= 0x003 + 0x008;
  sy(3).ticks.horiz.flags= 0x003 + 0x028;
  sy(3).ticks.vert.flags= 0x003 + 0x028;
  sy(4).ticks.horiz.flags= 0x003 + 0x008;
  sy(4).ticks.vert.flags= 0x003 + 0x008;
  //sy(1).ticks.vert.labelOff= 0.005;
  //sy(1).viewport= [0.15,0.75,0.7,0.9];
  set_style, ls, sy, lg, clg;

  amp_avg= img_avg(_DVV_AMPL_DATA);
  phs_avg= img_avg(_DVV_PHASE_DATA);
  xoff_avg= img_avg(_DVV_XOFF_DATA);
  yoff_avg= img_avg(_DVV_YOFF_DATA);

  
  sh, w, _DVV_AMPL_DATA, amp_avg-0.1, amp_avg+0.1, plsy= 1, notitle= 1, nolabels= 1;
  sh, w, _DVV_PHASE_DATA, phs_avg-0.1, phs_avg+0.1, plsy= 2, notitle= 1, nolabels= 1;
  sh, w, _DVV_XOFF_DATA,  xoff_avg-0.1, xoff_avg+0.1,  plsy= 3, notitle= 1, nolabels= 1;
  sh, w, _DVV_YOFF_DATA,  yoff_avg-0.1, yoff_avg+0.1,  plsy= 4, notitle= 1, nolabels= 1;

  if (!is_void(_DVV_MASK_DATA)) {
    wm= where(img_data(_DVV_MASK_DATA));
    xgr= img_grid(_DVV_MASK_DATA, 1);
    ygr= img_grid(_DVV_MASK_DATA, 2);
  } else wm= [];
  if (is_array(wm)) {
    plsys, 1; plmk, ygr(wm), xgr(wm), marker= 6, msize= 0.5, color= "yellow", width= 10;
    plsys, 2; plmk, ygr(wm), xgr(wm), marker= 6, msize= 0.5, color= "yellow", width= 10;
    plsys, 3; plmk, ygr(wm), xgr(wm), marker= 6, msize= 0.5, color= "yellow", width= 10;
    plsys, 4; plmk, ygr(wm), xgr(wm), marker= 6, msize= 0.5, color= "yellow", width= 10;
  }

  if (is_void(_DVV_MASK_DATA)) {
    ww= indgen(numberof(img_data(_DVV_AMPL_DATA)));
  } else {
    ww= where(!img_data(_DVV_MASK_DATA) & xgr < 300.0 & xgr > -300.0 & ygr < 300.0 & ygr > -300.0);
  }
  plt, swrite(format= "AMPL - avg %5.3f, rms %5.3f",
              img_data(_DVV_AMPL_DATA)(ww)(avg),
              img_data(_DVV_AMPL_DATA)(ww)(rms)), 0.135, 0.88, tosys= 0, color= "red";
  plt, swrite(format= "PHASE - avg %5.3f, rms %5.3f",
              img_data(_DVV_PHASE_DATA)(ww)(avg),
              img_data(_DVV_PHASE_DATA)(ww)(rms)), 0.46, 0.88, tosys= 0, color= "red";
  plt, swrite(format= "XOFF - avg %5.3f, rms %5.3f",
              img_data(_DVV_XOFF_DATA)(ww)(avg),
              img_data(_DVV_XOFF_DATA)(ww)(rms)), 0.135, 0.57, tosys= 0, color= "red";
  plt, swrite(format= "YOFF - avg %5.3f, rms %5.3f",
              img_data(_DVV_YOFF_DATA)(ww)(avg),
              img_data(_DVV_YOFF_DATA)(ww)(rms)), 0.46, 0.57, tosys= 0, color= "red";
  plt, shot+" - Phase/amplitude mapping", 0.41, 0.97, tosys= 0, justify="CC",
    height= 16, font= "helveticaB";
  
}

func dvv_DisplayScan(sc, zmin, zmax, w=)
{
  if(is_void(w)) w= 0;
  window, w, legends= 0, style= "square2x2img.gs", height= 650, width= 600;
  get_style, ls, sy, lg, clg;
  //ls= 0;
  sy(1).ticks.horiz.flags= 0x003 + 0x008;
  sy(1).ticks.vert.flags= 0x003 + 0x008;
  sy(2).ticks.horiz.flags= 0x003 + 0x008;
  sy(2).ticks.vert.flags= 0x003 + 0x008;
  sy(3).ticks.horiz.flags= 0x003 + 0x028;
  sy(3).ticks.vert.flags= 0x003 + 0x028;
  sy(4).ticks.horiz.flags= 0x003 + 0x008;
  sy(4).ticks.vert.flags= 0x003 + 0x008;
  //sy(1).ticks.vert.labelOff= 0.005;
  //sy(1).viewport= [0.15,0.75,0.7,0.9];
  set_style, ls, sy, lg, clg;

  restore, sc;

  sh, w, p1im, plsy= 1, lgscl= 1, notitle= 1, nolabels= 1;
  sh, w, p2im, plsy= 2, lgscl= 1, notitle= 1, nolabels= 1;
  sh, w, c1im, plsy= 3, lgscl= 1, notitle= 1, nolabels= 1;
  sh, w, c2im, plsy= 4, lgscl= 1, notitle= 1, nolabels= 1;

  // plt, swrite(format= "AMPL - avg %5.3f, rms %5.3f",
  //             img_avg(_DVV_AMPL_DATA),
  //             img_rms(_DVV_AMPL_DATA)), 0.135, 0.88, tosys= 0, color= "red";
  // plt, swrite(format= "PHASE - avg %5.3f, rms %5.3f",
  //             img_avg(_DVV_PHASE_DATA),
  //             img_rms(_DVV_PHASE_DATA)), 0.46, 0.88, tosys= 0, color= "red";
  // plt, swrite(format= "XOFF - avg %5.3f, rms %5.3f",
  //             img_avg(_DVV_XOFF_DATA),
  //             img_rms(_DVV_XOFF_DATA)), 0.135, 0.57, tosys= 0, color= "red";
  // plt, swrite(format= "YOFF - avg %5.3f, rms %5.3f",
  //             img_avg(_DVV_YOFF_DATA),
  //             img_rms(_DVV_YOFF_DATA)), 0.46, 0.57, tosys= 0, color= "red";
}

func dvv_DisplayMapping(map, zmin, zmax, w=, color=, type=, comment=)
{

  if (is_void(comment)) comment="";
  if(is_void(w)) {
    w= next_window();
    window, w, legends= 0, style= "square2x2img.gs";
    w1= next_window();
    window, w1, legends= 0, style= "square2x2img.gs";
    w2= next_window();
    window, w2, legends= 0, style= "square2x2img.gs";
  } else {
    window, w, legends= 0, style= "square2x2img.gs";
    w1= next_window();
    window, w1, legends= 0, style= "square2x2img.gs";
    w2= next_window();
    window, w2, legends= 0, style= "square2x2img.gs";
  }

  if (is_void(zmin)) zmin= -10;
  if (is_void(zmax)) zmax= 10;
  if (is_member(map.ch1S, "repairedList")) {
    r1s= *map.ch1S.repairedList;
    r2s= *map.ch2S.repairedList;
    r1p= *map.ch1P.repairedList;
    r2p= *map.ch2P.repairedList;
  }

  window, w; plsys, 1;
  plg, (*map.ch1S.warp_pts).x - (*map.warp_refpts).x, color= color, type= 1, marks= 0, width= 3;
  plg, (*map.ch1S.warp_pts).y - (*map.warp_refpts).y, color= color, type= 2, marks= 0, width= 3;
  if (numberof(r1s)) {
    plmk, (*map.ch1S.warp_pts)(r1s).x - (*map.warp_refpts)(r1s).x, r1s, msize= 0.5, marker= 4, width= 13;
    plmk, (*map.ch1S.warp_pts)(r1s).y - (*map.warp_refpts)(r1s).y, r1s, msize= 0.5, marker= 4, width= 13; }        
  limits; range, zmin, zmax;

  window, w; plsys, 2;
  plg, (*map.ch2S.warp_pts).x - (*map.warp_refpts).x, color= color, type= 1, marks= 0, width= 3;
  plg, (*map.ch2S.warp_pts).y - (*map.warp_refpts).y, color= color, type= 2, marks= 0, width= 3;
  if (numberof(r2s)) {
    plmk, (*map.ch2S.warp_pts)(r2s).x - (*map.warp_refpts)(r2s).x, r2s, msize= 0.5, marker= 4, width= 13;
    plmk, (*map.ch2S.warp_pts)(r2s).y - (*map.warp_refpts)(r2s).y, r2s, msize= 0.5, marker= 4, width= 13; }        
  limits; range, zmin, zmax;

  window, w; plsys, 3;
  plg, (*map.ch1P.warp_pts).x - (*map.warp_refpts).x, color= color, type= 1, marks= 0, width= 3;
  plg, (*map.ch1P.warp_pts).y - (*map.warp_refpts).y, color= color, type= 2, marks= 0, width= 3;
  if (numberof(r1p)) {
    plmk, (*map.ch1P.warp_pts)(r1p).x - (*map.warp_refpts)(r1p).x, r1p, msize= 0.5, marker= 4, width= 13;
    plmk, (*map.ch1P.warp_pts)(r1p).y - (*map.warp_refpts)(r1p).y, r1p, msize= 0.5, marker= 4, width= 13; }        
  limits; range, zmin, zmax;

  window, w; plsys, 4;
  plg, (*map.ch2P.warp_pts).x - (*map.warp_refpts).x, color= color, type= 1, marks= 0, width= 3;
  plg, (*map.ch2P.warp_pts).y - (*map.warp_refpts).y, color= color, type= 2, marks= 0, width= 3;
  if (numberof(r2p)) {
    plmk, (*map.ch2P.warp_pts)(r2p).x - (*map.warp_refpts)(r2p).x, r2p, msize= 0.5, marker= 4, width= 13;
    plmk, (*map.ch2P.warp_pts)(r2p).y - (*map.warp_refpts)(r2p).y, r2p, msize= 0.5, marker= 4, width= 13; }        
  limits; range, zmin, zmax;

  xref= (*map.warp_refpts).x;
  yref= (*map.warp_refpts).y;

  window, w1;
  c1sdxi=img_new(span(-400,400,9),span(-400,400,9),reform((*map.ch1S.warp_pts).x - xref, 9,9));
  c1sdyi=img_new(span(-400,400,9),span(-400,400,9),reform((*map.ch1S.warp_pts).y - yref, 9,9));
  c1sdxi= img_sub(c1sdxi, img_data(c1sdxi)(*)(avg));
  c1sdyi= img_sub(c1sdyi, img_data(c1sdyi)(*)(avg));
  sh, w1, c1sdxi, zmin, zmax, plsy= 1, notitle= 1, nolabels= 1;
  if (numberof(r1s)) plmk, (*map.warp_refpts)(r1s).y, (*map.warp_refpts)(r1s).x, color= "red", msize= 0.5, marker= 4, width= 3;
  sh, w2, c1sdyi, zmin, zmax, plsy= 1, notitle= 1, nolabels= 1;
  if (numberof(r1s)) plmk, (*map.warp_refpts)(r1s).y, (*map.warp_refpts)(r1s).x, color= "red", msize= 0.5, marker= 4, width= 3;

  c2sdxi=img_new(span(-400,400,9),span(-400,400,9),reform((*map.ch2S.warp_pts).x - xref, 9,9));
  c2sdyi=img_new(span(-400,400,9),span(-400,400,9),reform((*map.ch2S.warp_pts).y - yref, 9,9));
  c2sdxi= img_sub(c2sdxi, img_data(c2sdxi)(*)(avg));
  c2sdyi= img_sub(c2sdyi, img_data(c2sdyi)(*)(avg));
  sh, w1, c2sdxi, zmin, zmax, plsy= 2, notitle= 1, nolabels= 1;
  if (numberof(r2s)) plmk, (*map.warp_refpts)(r2s).y, (*map.warp_refpts)(r2s).x, color= "red", msize= 0.5, marker= 4, width= 3;
  sh, w2, c2sdyi, zmin, zmax, plsy= 2, notitle= 1, nolabels= 1;
  if (numberof(r2s)) plmk, (*map.warp_refpts)(r2s).y, (*map.warp_refpts)(r2s).x, color= "red", msize= 0.5, marker= 4, width= 3;

  c1pdxi=img_new(span(-400,400,9),span(-400,400,9),reform((*map.ch1P.warp_pts).x - xref, 9,9));
  c1pdyi=img_new(span(-400,400,9),span(-400,400,9),reform((*map.ch1P.warp_pts).y - yref, 9,9));
  c1pdxi= img_sub(c1pdxi, img_data(c1pdxi)(*)(avg));
  c1pdyi= img_sub(c1pdyi, img_data(c1pdyi)(*)(avg));
  sh, w1, c1pdxi, zmin, zmax, plsy= 3, notitle= 1, nolabels= 1;
  if (numberof(r1p)) plmk, (*map.warp_refpts)(r1p).y, (*map.warp_refpts)(r1p).x, color= "red", msize= 0.5, marker= 4, width= 3;
  sh, w2, c1pdyi, zmin, zmax, plsy= 3, notitle= 1, nolabels= 1;
  if (numberof(r1p)) plmk, (*map.warp_refpts)(r1p).y, (*map.warp_refpts)(r1p).x, color= "red", msize= 0.5, marker= 4, width= 3;

  c2pdxi=img_new(span(-400,400,9),span(-400,400,9),reform((*map.ch2P.warp_pts).x - xref, 9,9));
  c2pdyi=img_new(span(-400,400,9),span(-400,400,9),reform((*map.ch2P.warp_pts).y - yref, 9,9));
  c2pdxi= img_sub(c2pdxi, img_data(c2pdxi)(*)(avg));
  c2pdyi= img_sub(c2pdyi, img_data(c2pdyi)(*)(avg));
  sh, w1, c2pdxi, zmin, zmax, plsy= 4, notitle= 1, nolabels= 1; 
  if (numberof(r2p)) plmk, (*map.warp_refpts)(r2p).y, (*map.warp_refpts)(r2p).x, color= "red", msize= 0.5, marker= 4, width= 3;
  sh, w2, c2pdyi, zmin, zmax, plsy= 4, notitle= 1, nolabels= 1;
  if (numberof(r2p)) plmk, (*map.warp_refpts)(r2p).y, (*map.warp_refpts)(r2p).x, color= "red", msize= 0.5, marker= 4, width= 3;

  window, w; plt, shot+": mapping - x, y deviations\n"+comment, 0.099, 0.955, tosys= 0, height=18;
  window, w1; plt, shot+": mapping - x deviations\n"+comment, 0.099, 0.955, tosys= 0, height=18;
  window, w2; plt, shot+": mapping - y deviations\n"+comment, 0.099, 0.955, tosys= 0, height=18;
}

func dvv_Liss(img1, img2, lh=, lv=, w=, reg=)
/* DOCUMENT dvv_liss

   img1, and img2 are the sin and cos images ...
     
   SEE ALSO:
 */
{
  extern _lh, _lv;

  theta= span(-pi, pi, 200);

  if(!is_void(reg)) {
    img1= img_extract(img1, reg);
    img2= img_extract(img2, reg);
  }

  if(is_void(lh) && !is_void(_lh)) lh= _lh;
  if(is_void(lv) && !is_void(_lv)) lv= _lv;

  if(is_void(w)) w= 0;
  window, w;
  if(!is_void(lh)) {
    loh1= lineout1(img1, lh);
    loh2= lineout1(img2, lh);
    plg, loh2(,2), loh1(,2);
  }
  if(!is_void(lv)) {
    lov1= lineout1(img1, lv);
    lov2= lineout1(img2, lv);
    plg, lov2(,2), lov1(,2), color= "red";
  }
  limits, -1, 1, -1, 1;
  gridxy, 1, 1;
  plg, sin(theta), cos(theta), color= "blue", width= 5, marks= 0;
}

func dvv_plCircle(ctr, radius, color=, width=, type=, marks=, nolimits=)
/* DOCUMENT dvv_plCircle, ctr, radius, color=, width=, type=, marks=
     
   SEE ALSO:
 */
{
  theta= span(-pi, pi, 200);
  plg, ctr(2)+radius*sin(theta), ctr(1)+radius*cos(theta),
    color= color, width= width, marks= marks, type= type;
  if (!nolimits) limits, ctr(1)-radius, ctr(1)+radius, ctr(2)-radius, ctr(2)+radius;
}

func dvv_plsp(ps, w=, color=, marks=, width=, legend=, sptype=)
/* DOCUMENT dvv_plsp, ps, w=, color=, marks=, width=, legend=, sptype=

   Plots a power spectrum of ...
     
   SEE ALSO:
 */
{
  if(is_void(w)) w= 0;
  window, w;
  if(is_void(sptype)) sptype= "azavg";

  if(sptype == "azavg") {
    dtVaz= *ps.dtVaz;
    mtf= (dvv_DiffMTF(*dtV.xscale, 1.0/(2.0*0.395*3.0)))(-(numberof(dtVaz(,2))-1):0);
    plg, dtVaz(,2)/mtf*1000.0, dtVaz(,1),
      color= color, width= width, marks= marks, legend= legend;
  } else if (sptype == "preimposed") {
    MLO= *ps.MLO;
    mtf= dvv_DiffMTF(*dtV.xscale, 1.0/(2.0*0.395*3.0));
    plg, MLO(,2)/mtf*1000.0, MLO(,1),
      color= color, width= width, marks= marks, legend= legend;
  } else error, "unknown sptype";

  limits, 0.02, 5, 0, 6;
  logxy, 1, 0;
}

func dvv_ShowUWPhase(pUW, w=, suffix=, lims=, pal=, rge=, nodc=, box=, lh=, lv=, crcl=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _DVV_RESAMPLE_DISPLAY;
  extern dir, shot, reg;
  extern ctr, radius;
  
  if(is_void(w)) w= 0;
  window, w, legends= 0;
  if(is_void(radius)) rr= 300.0;
  else rr= radius;

//   reg1= reg;
//   reg1.x1+=ctr(1);
//   reg1.x2+=ctr(1);
//   reg1.y1+=ctr(2);
//   reg1.y2+=ctr(2);

//   lh= LO_POS(or= "h", x1= ctr(2)+1.0, x2= ctr(2)-1.0);
//   lv= LO_POS(or= "v", x1= ctr(1)+1.0, x2= ctr(1)-1.0);

  if (numberof(rge) == 2) dvv_DisplayImg, pUW, rge(1), rge(2), w= w, nodc= nodc,
    pltrms= 1, box= box, lh= lh, lv= lv, shlo= 1, pltitl=shot+": "+suffix+" phase";
  else dvv_DisplayImg, pUW, w= w, nodc= nodc,
    pltrms= 1, box= box, lh= lh, lv= lv, shlo= 1, pltitl=shot+": "+suffix+" phase";
  
  plsys, 1;
  //if (is_void(lims)) limits, -rr+ctr(1), rr+ctr(1), -rr+ctr(2), rr+ctr(2);
  //else limits, lims(1)+ctr(1), lims(2)+ctr(1), lims(3)+ctr(2), lims(4)+ctr(2);
  limits;
  if (!is_void(crcl)) dvv_plCircle, crcl(,1), crcl(1,2), width= 3, type= 2, marks= 0, color= "red", nolimits= 1;
  if (is_void(pal)) pal= "earth.gp";
  palette, pal;
  limsReset;
  //wprt, w, dir+"/"+shot+"-"+suffix+"-phase.ps";
  wpdf, w, dir+"/"+shot+"-"+suffix+"-phase";
}

func dvv_ShowUWVelocity(pUW, &vel, w=, suffix=, lims=, rmslims=, pal=, rge=,
                           nodc=, box=, lh=, lv=, crcl=, scale=, nfringes=, bkbox=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _DVV_RESAMPLE_DISPLAY;
  extern dir, shot, reg;
  extern ctr, radius;
  
  if(is_void(w)) w= 0;
  if(is_void(suffix)) suffix= "";
  window, w, legends= 0;
  if(!is_void(radius)) rr= radius;

  vel= dvv_Vel(pUW, scale= scale);
  if (nfringes) vel= dvv_ShiftShock(vel, nfringes, [ctr, radius], bkbox);
  else (*vel.data) -= (*vel.data)(*)(avg);
  //if (!is_void(box)) vel= img_extract(vel, box);

  if (numberof(rge) == 2) dvv_DisplayImg, vel, rge(1), rge(2), w= w, nodc= nodc, rmslims= rmslims,
    pltrms= 1, box= box, lh= lh, lv= lv, shlo= 1, pltitl=shot+": "+suffix+" velocity";
  else dvv_DisplayImg, vel, w= w, nodc= nodc, rmslims= rmslims,
    pltrms= 1, box= box, lh= lh, lv= lv, shlo= 1, pltitl=shot+": "+suffix+" velocity";
  
  plsys, 1;
  //if(is_void(lims)) limits, -rr+ctr(1), rr+ctr(1), -rr+ctr(2), rr+ctr(2);
  //else limits, lims(1)+ctr(1), lims(2)+ctr(1), lims(3)+ctr(2), lims(4)+ctr(2);
  limits;
  if (!is_void(crcl)) dvv_plCircle, crcl(,1), crcl(1,2), width= 3, type= 2, marks= 0, color= "red", nolimits= 1;
  if (is_void(pal)) pal= "earth.gp";
  palette, pal;
  limsReset;
  //wprt, w, dir+"/"+shot+"-"+suffix+"-velocity.ps";
  wpdf, w, dir+"/"+shot+"-"+suffix+"-velocity";
}

func dvv_ShowUWVel3D(pUW, &vel, w=, suffix=, pal=, rge=, box=, alt=, az=, edges=, scale=, nfringes=, bkbox=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _DVV_RESAMPLE_DISPLAY;
  extern dir, shot, reg;
  extern ctr, radius;
  
  if (is_void(w)) w= 0;
  if (is_void(suffix)) suffix= "";
  window, w, legends= 0, wait= 1, style= "nobox.gs";
  if(is_void(pal)) pal= "earth.gp";
  palette, pal;
  if(!is_void(radius)) rr= radius;

  if (_DVV_RESAMPLE_DISPLAY) pUW= img_resample(pUW, _DVV_RESAMPLE_DISPLAY);

  vel= dvv_Vel(pUW, scale= scale);
  if (nfringes) vel= dvv_ShiftShock(vel, nfringes, [ctr, radius], bkbox);
  if (!is_void(box)) vel= img_extract(vel, box);
  if(!nfringes) (*vel.data) -= (*vel.data)(*)(avg);
  

  if (numberof(rge) == 2) pl3s, *vel.data, *vel.yscale, *vel.xscale, axis= 1, fill= 1,
                edges= edges, height= 12, alt= alt, az= az;
  //                edges= edges, height= 12, zrange= rge, alt= alt, az= az;
  else pl3s, *vel.data, *vel.yscale, *vel.xscale, axis= 1, fill= 1,
         edges= edges, height= 12, alt= alt, az= az;
  
  //wprt, w, dir+"/"+shot+"-"+suffix+"-vel3D.ps";
  wpdf, w, dir+"/"+shot+"-"+suffix+"-vel3D";
}

func dvv_ShowUWVelContour(pUW, &vel, w=, suffix=, pal=, rge=, box=, marks=,
                          nlevs=, ncolors=, color=, type=, scale=, nfringes=, bkbox=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _DVV_RESAMPLE_DISPLAY;
  extern dir, shot, reg;
  extern ctr, radius;
  
  if(is_void(w)) w= 0;
  if(is_void(suffix)) suffix= "";
  window, w, legends= 0;
  if(is_void(pal)) pal= "earth.gp";
  palette, pal;
  if(!is_void(radius)) rr= radius;
  if(!is_void(box)) vel= img_extract(pUW, box);
  else vel= img_copy(pUW);
  if(is_void(rge) || allof(rge == 0)) rge=[-300,300];
  
  if (_DVV_RESAMPLE_DISPLAY) pUW= img_resample(pUW, _DVV_RESAMPLE_DISPLAY);

  vel= dvv_Vel(pUW, scale= scale);
  if (nfringes) vel= dvv_ShiftShock(vel, nfringes, [ctr, radius], bkbox);
  if (!is_void(box)) vel= img_extract(vel, box);
  if (!nfringes) (*vel.data) -= (*vel.data)(*)(avg);

  if(is_void(nlevs)) nlevs= 5;
  if (numberof(rge) == 2) pls, *vel.data, *vel.yscale, *vel.xscale, cmin= rge(1), cmax= rge(2), 
            levs= span(rge(1), rge(2), nlevs), marks= marks, height= 18, color= color,type= type;
  else pls, *vel.data, *vel.yscale, *vel.xscale, nlevs= nlevs, marks= marks,
             height= 18, color= color, type= type;

  if(is_void(ncolors)) ncolors= 1000;
  color_bar1, span(rge(1), rge(2), ncolors+1), span(1,200,ncolors+2),
    vert= 1, labs= [1, ncolors/(nlevs-1)], edges= 0;
  plt, vel.z_label+" ("+vel.z_unit+")", 0.70, 0.64, tosys= 0,
    justify= "CC", orient= 1, height= 18;
  
  xytitles, vel.x_label+" ("+vel.x_unit+")", vel.y_label+" ("+vel.y_unit+")";
  pltitle, shot+": "+suffix+" velocity";
  
  if(is_void(pal)) pal= "earth.gp";
  palette, pal;
  //wprt, w, dir+"/"+shot+"-"+suffix+"-velContour.ps";
  wpdf, w, dir+"/"+shot+"-"+suffix+"-velContour";
}

func dvv_ShowSpectV(pUW, &pV, &pVaz, az0=, w=, suffix=, lims=, pal=, box=, bxctr=,
                    rge=, rmslims=, lgscl=, prlo=, norm=, pl1D=, pl1Drge=, nodc=, opfilt=)
/* DOCUMENT dvv_ShowSpectV, pUW, &pV, &pVaz, az0=, w=, suffix=, lims=, pal=, box=, bxctr=,
                  rge=, rmslims=, lgscl=, prlo=, norm=, pl1D=, pl1Drge=, nodc=, opfilt=


  SEE ALSO:
  dvv_ShowUWPhase, dvv_ShowUWVelocity
  
*/
{
  extern _DVV_RESAMPLE_DISPLAY;
  extern dir, reg, shot;
  extern ctr;
  extern _VPF, _TAU;

  if(is_void(suffix)) suffix= "";
  if(is_void(shot)) shot= "";
  if(is_void(w)) w= 0;
  if(is_void(nodc)) nodc= 2;
  window, w, legends= 0;

  if (!is_void(box)) box1= box;  else box1= reg;

  if (!is_void(bxctr)) {
    box1.x1+=bxctr(1);
    box1.x2+=bxctr(1);
    box1.y1+=bxctr(2);
    box1.y2+=bxctr(2);
  }
  
  print, "SpectV: box = ", box1;

  phs2D= img_extract(pUW, box1);
  if (!is_void(opfilt))
    phs2D= dvv_Filter(phs2D, dvv_GenFilterMask(phs2D, opfilt));
  V= dvv_Vel(phs2D, scale= 1000.0);
  Lx= (*V.xscale)([1, V.nx])(dif)(1);
  Ly= (*V.yscale)([1, V.ny])(dif)(1);
  if(!is_void(rmslims)) fm= dvv_GenFilterMask(phs2D, rmslims, type="band_pass");

  pV= dvv_AMP(V, wndw= "hanning", norm= norm, nodc= nodc);
  pVaz= pVazRMS= dvv_AMP_1D(pV, qavg= 1);
  phs= dvv_AMP_1D(dvv_AMP(phs2D, wndw= "hanning", norm= 1, nodc= nodc), qavg= 1);
  if(is_void(norm) || norm == 1) {
    pV.z_label="rms velocity";
    pV.z_unit="m/s";
    pVazRMS(,2)*=sqrt(Lx*Ly);
  } else if (norm == 2) {
    pV.z_label="velocity density";
    pV.z_unit="m/s-!mm";
  }
  
  if(numberof(rge) == 2) dvv_DisplayImg, pV, rge(1), rge(2), lgscl= lgscl, rmsunit="m/s",
                           w= w, pltrms= "amplitude", rmslims= rmslims, norm= norm,
      pltitl=shot+": "+suffix+" velocity spectrum", prlo= prlo;
  else dvv_DisplayImg, pV, w= w, lgscl= lgscl, pltrms= "amplitude", rmslims= rmslims,
         norm= norm, pltitl=shot+": "+suffix+" velocity spectrum",
      prlo= prlo, rmsunit="m/s";

  lofact= 1.0;
  plsys, 2;
  plg, lofact*pVaz(,2), pVaz(,1), color= "red", width= 3, marks= 0;
  plg, lofact*pVaz(,2), -pVaz(,1), color= "red", width= 3, marks= 0;
  plsys, 3; 
  plg, pVaz(,1), lofact*pVaz(,2), color= "red", width= 3, marks= 0;
  plg, -pVaz(,1), lofact*pVaz(,2), color= "red", width= 3, marks= 0;
  
  if(numberof(az0) > 0) {
    plsys, 2;
    plg, lofact*az0(,2), az0(,1), color= "blue", width= 3, marks= 0;
    plg, lofact*az0(,2), -az0(,1), color= "blue", width= 3, marks= 0;
    plsys, 3; 
    plg, az0(,1), lofact*az0(,2), color= "blue", width= 3, marks= 0;
    plg, -az0(,1), lofact*az0(,2), color= "blue", width= 3, marks= 0;
  }
  
  if(!is_void(prlo)) {
    fname= swrite(format="%s-AzAvg.txt", prlo);
    fh= open(fname, "w");
    write, fh, format="# %s", prlo+" - Velocity spectrum azimuthal average\n";
    write, fh, format="# Nx = %d, Ny = %d\n", pV.nx, pV.ny;
    write, fh, format="# Lx = %7.2f, Ly = %7.2f\n", Lx, Ly;
    write, fh, format="# Box limits: [%6.2f, %6.2f, %6.2f, %6.2f]\n",
      box1.x1, box1.x2, box1.y1, box1.y2;
    write, fh, format="# %d points in 1D array\n", dimsof(pVaz)(2);
    if(is_void(norm)) write, fh, format="# Normalization is %d\n", 1;
    else write, fh, format="# Normalization is %d\n", norm;
    if(!is_void(rmslims))
    write, fh, format="# Frequency range is [%6.4f, %6.4f]\n", rmslims(1), rmslims(2);
    write, fh, format="# RMS phase is %12.5e [unfiltered %12.5e]\n",
      dvv_rms(phs2D, mask= fm), dvv_rms(phs2D);
    write, fh, format="# RMS phase is %12.5e [unfiltered %12.5e]\n",
      dvv_AMP_1D_rms(phs, rmslims, qavg= 1, norm= 1),
      dvv_AMP_1D_rms(phs, qavg= 1, norm= 1);
    write, fh, format="# RMS velocity is %12.5e [unfiltered %12.5e]\n",
      dvv_rms(V, mask= fm), dvv_rms(V);
    write, fh, format="# RMS velocity is %12.5e [unfiltered %12.5e]\n",
      dvv_AMP_1D_rms(pVaz, rmslims, qavg= 1, norm= norm),
      dvv_AMP_1D_rms(pVaz, qavg= 1, norm= norm);
    write, fh, format="# VPF = %6.4f km/s/fringe, TAU = %7.3f ps\n", _VPF, _TAU*1000.0;
    write, fh, format="# %s(%s) %s(%s) %s rms_phase(mrad)\n", pV.x_label, pV.x_unit,
      pV.z_label, pV.z_unit, "Normalized_velocity(m/s-um)";
    for(j= 1; j<= dimsof(pVaz)(2); j++)
      write, fh, pVaz(j,1), pVaz(j,2), pVazRMS(j,2), phs(j,2)*1000.0;
    close, fh;
  }
  plsys, 1;
  if(is_void(lims)) limits, -0.5, 0.5, -0.5, 0.5;
  else limits, lims(1), lims(2), lims(3), lims(4);
  limsReset;
  if (numberof(rge)==2) {
    plsys, 2;
    range, rge(1), rge(2);
    plsys, 3;
    limits, rge(1), rge(2);
  }
  if(is_void(pal) || pal == "") pal= "earth.gp";
  palette, pal;
  
  //wprt, w, dir+"/"+shot+"-"+suffix+"-2D-spectrum.ps";
  if (!is_void(dir) && !is_void(suffix) && !is_void(shot))
    wpdf, w, dir+"/"+shot+"-"+suffix+"-2D-spectrum";
  if(!is_void(pl1D)) {
    if(is_void(pl1Drge)) pl1Drge=[10.0, 1000.0];
    window, w+20, legends= 0;
    fma;
    plg, pVazRMS(,2), pVaz(,1), color= "red"; width= 3;
    if(!is_void(az0)) {
      if(norm == 1) plg, az0(,2)*sqrt(Lx*Ly), az0(,1), type= 2, width= 3;
      else plg, az0(,2), az0(,1), type= 2, width= 3;
    }
    logxy, 1, 1;
    limits, 0.002, 0.5, pl1Drge(1), pl1Drge(2);
    gridxy, 1, 1;
    xytitles, "Spatial Frequency (!mm^-1^)", "Velocity density (m/s-!mm)";
    pltitle, shot+suffix+"- velocity spectrum";
    //wprt, w+20, dir+"/"+shot+"-"+suffix+"-1D-spectrum.ps";
    if (!is_void(dir) && !is_void(suffix) && !is_void(shot))
      wpdf, w+20, dir+"/"+shot+"-"+suffix+"-1D-spectrum";
    logxy, 1, 0;
    limits, 0.01, 0.5, 0.0, pl1Drge(2);
    if (!is_void(dir) && !is_void(suffix) && !is_void(shot))
      wpdf, w+20, dir+"/"+shot+"-"+suffix+"-1D-spectrum-B";
  }
}

func dvv_ShowIContour(fs, &sig, w=, suffix=, pal=, rge=, box=, marks=,
                          nlevs=, ncolors=, color=, type=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _DVV_RESAMPLE_DISPLAY;
  extern dir, shot, reg;
  extern ctr, radius;
  
  if(is_void(w)) w= 0;
  if(is_void(suffix)) suffix= "";
  window, w, legends= 0;
  if(is_void(pal)) pal= "earth.gp";
  palette, pal;
  if(!is_void(radius)) rr= radius;
  if(is_void(rge) || allof(rge == 0)) rge=[0,10000];

  if (_DVV_RESAMPLE_DISPLAY) fs= img_resample(fs, _DVV_RESAMPLE_DISPLAY);

  if(!is_void(box)) sig= img_extract(dvv_Ampl(fs), box);
  else sig= dvv_Ampl(fs);

  if(is_void(nlevs)) nlevs= 5;
  if (numberof(rge) == 2) pls, *sig.data, *sig.yscale, *sig.xscale, cmin= rge(1), cmax= rge(2), 
            levs= span(rge(1), rge(2), nlevs), marks= marks, height= 18, color= color,type= type;
  else pls, *sig.data, *sig.yscale, *sig.xscale, nlevs= nlevs, marks= marks,
             height= 18, color= color, type= type;

  if(is_void(ncolors)) ncolors= 1000;
  color_bar, span(rge(1), rge(2), ncolors+1), span(1,200,ncolors+2),
    vert= 1, labs= [1, ncolors/(nlevs-1)], edges= 0;
  plt, sig.z_label+" ("+sig.z_unit+")", 0.70, 0.64, tosys= 0,
    justify= "CC", orient= 1, height= 18;
  
  xytitles, sig.x_label+" ("+sig.x_unit+")", sig.y_label+" ("+sig.y_unit+")";
  pltitle, shot+": "+suffix+" intensity";
  
  if(is_void(pal)) pal= "earth.gp";
  palette, pal;
  //wprt, w, dir+"/"+shot+"-"+suffix+"-IContour.ps";
  wpdf, w, dir+"/"+shot+"-"+suffix+"-IContour";
}

func dvv_CheckNorm(V, wndw=, norm=, type=, w=, rmsl=)
/* DOCUMENT dvv_CheckNorm, V, wndw=, norm=, type=, w=, rmsl=

   Generates plots to demonstrate correct normalization of a spectral
   decomposition of the input data set, V.  The set is resampled in various
   ways, as well as a subregion is extracted, and various versions
   of the power spectra are calculated. All of the curves should be
   in good agreement.

   The normalized spectrum gives rms amplitude per unit frequency, through
   3 different measures:
     (1) azimuthal sum over ...
     (2) azimuthal average over the mode power
     (3) azimuthal average over the mode amplitude
     
   SEE ALSO:
 */
{
  if(is_void(norm)) norm= 1;
  if(is_void(w)) w= 0;

  //Generate resampled & subsampled versions of
  //this data set
  xlims= (*V.xscale)([1,V.nx]);
  ylims= (*V.yscale)([1,V.ny]);
  //Resampled at 1/2 resolution (blue)
  V1= img_resample(V, V.nx/2, V.ny/2);
  //Resampled at double resolution (green)
  V2= img_resample(V, V.nx*2, V.ny*2);
  //Extract 0.7 x 0.7 of central area (magenta)
  Vx= img_extract(V, REGION(x1= xlims(1)*0.7, x2= xlims(2)*0.7,
                        y1= ylims(1)*0.7, y2= ylims(2)*0.7));
  //Resample extracted image to initial resolution (red)
  Vx1= img_resample(Vx, V.nx, V.ny);

  //Generate filter masks if needed
  if (!is_void(rmsl)) {
    m10= dvv_GenFilterMask(V, rmsl, type="band_pass");
    m11= dvv_GenFilterMask(V1, rmsl, type="band_pass");
    m12= dvv_GenFilterMask(V2, rmsl, type="band_pass");
    m1x= dvv_GenFilterMask(Vx, rmsl, type="band_pass");
    m1x1= dvv_GenFilterMask(Vx1, rmsl, type="band_pass");
  }

  write, format="%s\n", " ";
  write, format="%s%d\n", " Normalization: ", norm;
  write, format="%s\n",   "-------------------";
  write, format="%s\n", " ";
  write, format="%s\n", "RMS within frequency range computed using 2D RMS ... ";
  if (!is_void(rmsl))
  write, format="%s[%6.4f,%6.4f]\n", " Frequency ramge: ", rmsl(1), rmsl(2);
  write, format="%s\n", "Case        V      V1     V2      Vx      Vx1";
  write, format="%s\n", "=============================================";

  //Compute RMS measures for each case
  rms0= dvv_rms(V, mask= m10);
  rms1= dvv_rms(V1, mask= m11);
  rms2= dvv_rms(V2, mask= m12);
  rmsx= dvv_rms(Vx, mask= m1x);
  rmsx1= dvv_rms(Vx1, mask= m1x1);
  write, format="%s %6.4f %6.4f %6.4f %6.4f %6.4f\n", "Direct    ", rms0, rms1, rms2, rmsx, rmsx1;

  //Compute 2D psd for each case
  psd= dvv_PSD(V, wndw= wndw, norm= norm);
  psd1= dvv_PSD(V1, wndw= wndw, norm= norm);
  psd2= dvv_PSD(V2, wndw= wndw, norm= norm);
  psdx= dvv_PSD(Vx, wndw= wndw, norm= norm);
  psdx1= dvv_PSD(Vx1, wndw= wndw, norm= norm);

  //Compute RMS measures for 2D psd
  rms0= dvv_PSD_rms(psd, norm= norm, mask= m10);
  rms1= dvv_PSD_rms(psd1, norm= norm, mask= m11);
  rms2= dvv_PSD_rms(psd2, norm= norm, mask= m12);
  rmsx= dvv_PSD_rms(psdx, norm= norm, mask= m1x);
  rmsx1= dvv_PSD_rms(psdx1, norm= norm, mask= m1x1);
  write, format="%s %6.4f %6.4f %6.4f %6.4f %6.4f\n", "2D PSD    ", rms0, rms1, rms2, rmsx, rmsx1;

  //Power-averaged rms per mode
  p1d= dvv_PSD_1D(psd, qavg= 1);
  p1d1= dvv_PSD_1D(psd1, qavg= 1);
  p1d2= dvv_PSD_1D(psd2,qavg= 1);
  p1dx= dvv_PSD_1D(psdx, qavg= 1);
  p1dx1= dvv_PSD_1D(psdx1, qavg= 1);
  window, w, legends= 0;
  plg, p1d(,2), p1d(,1), color= "black", width= 3, type= type;
  plg, p1d1(,2), p1d1(,1), color= "blue", width= 3, type= type;
  plg, p1d2(,2), p1d2(,1), color= "green", width= 3, type= type;
  plg, p1dx(,2), p1dx(,1), color= "magenta", width= 3, type= type;
  plg, p1dx1(,2), p1dx1(,1), color= "red", width= 3, type= type;
  pltitle, swrite(format="Power-averaged rms per mode\n norm = %d", norm);
  logxy, 1, 1;

  //Compute RMS measures for each case
  rms0= dvv_PSD_1D_rms(p1d, rmsl, qavg= 1, norm= norm);
  rms1= dvv_PSD_1D_rms(p1d1, rmsl, qavg= 1, norm= norm);
  rms2= dvv_PSD_1D_rms(p1d2, rmsl, qavg= 1, norm= norm);
  rmsx= dvv_PSD_1D_rms(p1dx, rmsl, qavg= 1, norm= norm);
  rmsx1= dvv_PSD_1D_rms(p1dx1, rmsl, qavg= 1, norm= norm);
  write, format="%s %6.4f %6.4f %6.4f %6.4f %6.4f\n", "1D PSD AVG", rms0, rms1, rms2, rmsx, rmsx1;

  //limits, 0.02,2,0.1,100;

  //Power-summed rms per mode
  q1d= dvv_PSD_1D(psd);
  q1d1= dvv_PSD_1D(psd1);
  q1d2= dvv_PSD_1D(psd2);
  q1dx= dvv_PSD_1D(psdx);
  q1dx1= dvv_PSD_1D(psdx1);
  window, w+1, legends= 0; 
  plg, q1d(,2), q1d(,1), color= "red", width= 3, type= type;
  plg, q1d1(,2), q1d1(,1), color= "blue", width= 3, type= type;
  plg, q1d2(,2), q1d2(,1), color= "green", width= 3, type= type;
  plg, q1dx(,2), q1dx(,1), color= "magenta", width= 3, type= type;
  plg, q1dx1(,2), q1dx1(,1), color= "black", width= 3, type= type;
  pltitle, swrite(format="Power-summed rms per mode\n norm = %d", norm);
  logxy, 1, 1;
  //limits, 0.02,2,0.1,100;
  
  //Compute RMS measures for each case
  rms0= dvv_PSD_1D_rms(q1d, rmsl, norm= norm);
  rms1= dvv_PSD_1D_rms(q1d1, rmsl, norm= norm);
  rms2= dvv_PSD_1D_rms(q1d2, rmsl, norm= norm);
  rmsx= dvv_PSD_1D_rms(q1dx, rmsl, norm= norm);
  rmsx1= dvv_PSD_1D_rms(q1dx1, rmsl, norm= norm);
  write, format="%s %6.4f %6.4f %6.4f %6.4f %6.4f\n", "1D PSD SUM", rms0, rms1, rms2, rmsx, rmsx1;

  //Amplitude-averaged rms per mode
  p1d= dvv_AMP_1D(dvv_AMP(V, norm= norm), qavg= 1);
  p1d1= dvv_AMP_1D(dvv_AMP(V1, norm= norm), qavg= 1);
  p1d2= dvv_AMP_1D(dvv_AMP(V2, norm= norm), qavg= 1);
  p1dx= dvv_AMP_1D(dvv_AMP(Vx, norm= norm), qavg= 1);
  p1dx1= dvv_AMP_1D(dvv_AMP(Vx1, norm= norm), qavg= 1);
  window, w+2, legends= 0
  plg, p1d(,2), p1d(,1), color= "red", width= 3, type= type;
  plg, p1d1(,2), p1d1(,1), color= "blue", width= 3, type= type;
  plg, p1d2(,2), p1d2(,1), color= "green", width= 3, type= type;
  plg, p1dx(,2), p1dx(,1), color= "magenta", width= 3, type= type;
  plg, p1dx1(,2), p1dx1(,1), color= "black", width= 3, type= type;
  pltitle, swrite(format="Amplitude-averaged rms per mode\n norm = %d", norm);
  logxy, 1, 1;  //limits, 0.02,2,0.1,100;

  //Compute RMS measures for each case
  rms0= dvv_AMP_1D_rms(p1d, rmsl, qavg= 1, norm= norm);
  rms1= dvv_AMP_1D_rms(p1d1, rmsl, qavg= 1, norm= norm);
  rms2= dvv_AMP_1D_rms(p1d2, rmsl, qavg= 1, norm= norm);
  rmsx= dvv_AMP_1D_rms(p1dx, rmsl, qavg= 1, norm= norm);
  rmsx1= dvv_AMP_1D_rms(p1dx1, rmsl, qavg= 1, norm= norm);
  write, format="%s %6.4f %6.4f %6.4f %6.4f %6.4f\n", "1D AMP AVG", rms0, rms1, rms2, rmsx, rmsx1;

  //Amplitude of mode sum
  p1d= dvv_AMP_1D(dvv_AMP(V, norm= norm));
  p1d1= dvv_AMP_1D(dvv_AMP(V1, norm= norm));
  p1d2= dvv_AMP_1D(dvv_AMP(V2, norm= norm));
  p1dx= dvv_AMP_1D(dvv_AMP(Vx, norm= norm));
  p1dx1= dvv_AMP_1D(dvv_AMP(Vx1, norm= norm));
  window, w+3, legends= 0
  plg, p1d(,2), p1d(,1), color= "red", width= 3, type= type;
  plg, p1d1(,2), p1d1(,1), color= "blue", width= 3, type= type;
  plg, p1d2(,2), p1d2(,1), color= "green", width= 3, type= type;
  plg, p1dx(,2), p1dx(,1), color= "magenta", width= 3, type= type;
  plg, p1dx1(,2), p1dx1(,1), color= "black", width= 3, type= type;
  pltitle, swrite(format="Amplitude-summed rms per mode\n norm = %d", norm);
  logxy, 1, 1;
  //limits, 0.02,2,0.1,100;

  //Compute RMS measures for each case
  rms0= dvv_AMP_1D_rms(p1d, rmsl, norm= norm);
  rms1= dvv_AMP_1D_rms(p1d1, rmsl, norm= norm);
  rms2= dvv_AMP_1D_rms(p1d2, rmsl, norm= norm);
  rmsx= dvv_AMP_1D_rms(p1dx, rmsl, norm= norm);
  rmsx1= dvv_AMP_1D_rms(p1dx1, rmsl, norm= norm);
  write, format="%s %6.4f %6.4f %6.4f %6.4f %6.4f\n", "1D AMP SUM", rms0, rms1, rms2, rmsx, rmsx1;
}

func dvv_ShowChannels(fs, shotid=, filtered=, w=)
/* DOCUMENT 
     
   SEE ALSO:
 */
{
  extern _DVV_CH1S_FILTER, _DVV_CH2S_FILTER;
  extern _DVV_CH1P_FILTER, _DVV_CH2P_FILTER;
  extern _DVV_CHS_FILTER, _DVV_CHP_FILTER;
  extern _CSPHI, _SNPHI, _CSPHI0, _SNPHI0, _CSPHI1, _SNPHI1;
  extern _ms, _mp;
  extern smask, pmask, sdiff, pdiff, sdiff0, pdiff0;
  
  //if(is_void(shotid)) shotid= fs.raw.shotid;
  if(is_void(shotid)) shotid= fs.ch1S.shotid;
  if(is_void(w)) w= 0;

  //Filter out the secondary interference mode for each channel
  if(filtered == 1) {
    ftxt=",filtered";
    ftxt1=", filtered";
    ch1S= dvv_Filter(fs.ch1S, dvv_GenFilterMaskSet(fs.ch1S, _DVV_CH1S_FILTER));
    ch2S= dvv_Filter(fs.ch2S, dvv_GenFilterMaskSet(fs.ch2S, _DVV_CH2S_FILTER));
    ch1P= dvv_Filter(fs.ch1P, dvv_GenFilterMaskSet(fs.ch1P, _DVV_CH1P_FILTER));
    ch2P= dvv_Filter(fs.ch2P, dvv_GenFilterMaskSet(fs.ch2P, _DVV_CH2P_FILTER));
  } else {
    ftxt="";
    ftxt1="";
    ch1S= img_copy(fs.ch1S);
    ch2S= img_copy(fs.ch2S);
    ch1P= img_copy(fs.ch1P);
    ch2P= img_copy(fs.ch2P);
  }

  //Spectra of each of the 4 channels
  dvv_DisplayImg, dvv_PSD(ch1S, wndw="hanning"), 10, 1e5, lgscl= 1,
    w= w, pltitl= shotid+", Ch1S spectrum"+ftxt1;
  plsys, 1; limits, -0.2, 0.2, -0.2, 0.2;
  limsReset;
  //wprt, w, shotid+"-spectrum,ch1S"+ftxt+".ps";
  wpdf, w, shotid+"-spectrum,ch1S"+ftxt;

  w+=1;
  dvv_DisplayImg, dvv_PSD(ch2S, wndw="hanning"), 10, 1e5, lgscl= 1,
    w= w, pltitl= shotid+", Ch2S spectrum"+ftxt1;
  plsys, 1; limits, -0.2, 0.2, -0.2, 0.2;
  limsReset;
  //wprt, w, shotid+"-spectrum,ch2S"+ftxt+".ps";
  wpdf, w, shotid+"-spectrum,ch2S"+ftxt;

  w+=1;
  dvv_DisplayImg, dvv_PSD(ch1P, wndw="hanning"), 10, 1e5, lgscl= 1,
    w= w, pltitl= shotid+", Ch1P spectrum"+ftxt1;
  plsys, 1; limits, -0.2, 0.2, -0.2, 0.2;
  limsReset;
  //wprt, w, shotid+"-spectrum,ch1P"+ftxt+".ps";
  wpdf, w, shotid+"-spectrum,ch1P"+ftxt;

  w+= 1;
  dvv_DisplayImg, dvv_PSD(ch2P, wndw="hanning"), 10, 1e5, lgscl= 1,
    w= w, pltitl= shotid+", Ch2P spectrum"+ftxt1;
  plsys, 1; limits, -0.2, 0.2, -0.2, 0.2;
  limsReset;
  //wprt, w, shotid+"-spectrum,ch2P"+ftxt+".ps";
  wpdf, w, shotid+"-spectrum,ch2P"+ftxt;

  //Sum over s-channels
  w+= 1;
  stot= img_add(ch1S, ch2S);
  if(filtered == 2) {
    stot= dvv_Filter(stot, dvv_GenFilterMaskSet(stot, _(_DVV_CH1S_FILTER, _DVV_CH2S_FILTER, _DVV_CHS_FILTER)));
    ftxt=",filtered-2";
    ftxt1=", filtered-type-2";
  }
  dvv_DisplayImg, stot, w= w, pltitl=shotid+", s-channel sum"+ftxt1;
  dvv_DisplayImg, dvv_PSD(stot, wndw= "hanning"), lgscl= 1, 100, 1e6,
    w= w+1, pltitl=shotid+", spectrum, s-channel sum"+ftxt1;
  plsys, 1; limits, -0.2, 0.2, -0.2, 0.2;
  limsReset;
  //wprt, w, shotid+"-s-channel-sum"+ftxt+".ps";
  //wprt, w+1, shotid+"-spectrum,s-channel-sum"+ftxt+".ps";
  wpdf, w, shotid+"-s-channel-sum"+ftxt;
  wpdf, w+1, shotid+"-spectrum,s-channel-sum"+ftxt;

  //Sum over p-channels
  w+= 2;
  ptot= img_add(ch1P, ch2P);
  if(filtered == 2) {
    ptot= dvv_Filter(ptot, dvv_GenFilterMaskSet(ptot, _(_DVV_CH1P_FILTER, _DVV_CH2P_FILTER, _DVV_CHP_FILTER)));
    ftxt=",filtered-2";
    ftxt1=", filtered-type-2";
  }
  dvv_DisplayImg, ptot, w= w, pltitl=shotid+", p-channel sum"+ftxt1;
  dvv_DisplayImg, dvv_PSD(ptot, wndw= "hanning"), lgscl= 1, 100, 1e6,
    w= w+1, pltitl=shotid+", spectrum, p-channel sum"+ftxt1;
  plsys, 1; limits, -0.2, 0.2, -0.2, 0.2;
  limsReset;
  //wprt, w, shotid+"-p-channel-sum"+ftxt+".ps";
  //wprt, w+1, shotid+"-spectrum,p-channel-sum"+ftxt+".ps";
  wpdf, w, shotid+"-p-channel-sum"+ftxt;
  wpdf, w+1, shotid+"-spectrum,p-channel-sum"+ftxt;

  //Difference of stot and ptot
  w+= 2;
  stot= img_add(ch1S, ch2S);
  ptot= img_add(ch1P, ch2P);
  spdiff= img_div(img_sub(ptot, stot), img_add(stot, ptot, average= 1));
  if(filtered == 2) {
    m= dvv_GenFilterMaskSet(spdiff, _(_DVV_CH1S_FILTER, _DVV_CH2S_FILTER,
                                     _DVV_CH1P_FILTER, _DVV_CH2P_FILTER,
                                      _DVV_CHS_FILTER, _DVV_CHP_FILTER));
    spdiff= dvv_Filter(spdiff, m, wndw= "hanning");
    ftxt=",filtered-2";
    ftxt1=", filtered-type-2";
  }
  dvv_DisplayImg, spdiff, w= w, -0.5, 0.5, pltitl=shotid+", stot,ptot-difference"+ftxt1;
  dvv_DisplayImg, dvv_PSD(spdiff, wndw= "hanning"), lgscl= 1, 1e-8, 1e-5, 
    w= w+1, pltitl=shotid+", spectrum,stot,ptot difference"+ftxt1;
  plsys, 1; limits, -0.2, 0.2, -0.2, 0.2;
  limsReset;
  //wprt, w, shotid+"-stot,ptot-difference"+ftxt+".ps";
  //wprt, w+1, shotid+"-spectrum,stot,ptot-difference"+ftxt+".ps";
  wpdf, w, shotid+"-stot,ptot-difference"+ftxt;
  wpdf, w+1, shotid+"-spectrum,stot,ptot-difference"+ftxt;

  //s-channel difference images
  w+= 2;
  sdiff= img_div(img_sub(ch2S, ch1S), img_add(ch1S, ch2S, average= 1));
  //  print, "Comparing sdiff0 with _CSPHI0: ", allof(*sdiff0.data == *_CSPHI0.data);
  if(filtered == 2) {
    smask= dvv_GenFilterMaskSet(sdiff, _(_DVV_CH1S_FILTER, _DVV_CH2S_FILTER, _DVV_CHS_FILTER));
    //    print, "Comparing smask with _ms: ", allof(smask == _ms);
    //    print, "Comparing sdiff0 with _CSPHI0: ", allof(*sdiff0.data == *_CSPHI0.data);
    sdiff= dvv_Filter(sdiff, smask, wndw= "hanning");
    //    print, "Comparing sdiff with _CSPHI: ", allof(*sdiff.data == *_CSPHI.data);
    ftxt=",filtered-2";
    ftxt1=", filtered-type-2";
  }
  dvv_DisplayImg, sdiff, w= w, -1.5, 1.5, pltitl=shotid+", s-channel difference, normalized"+ftxt1;
  dvv_DisplayImg, dvv_PSD(sdiff, wndw= "hanning"), lgscl= 1, 1e-8, 1e-5,
    w= w+1, pltitl=shotid+", spectrum, s-channel difference, normalized"+ftxt1;
  plsys, 1; limits, -0.2, 0.2, -0.2, 0.2;
  limsReset;
  //wprt, w, shotid+"-s-diff,normalized"+ftxt+".ps";
  //wprt, w+1, shotid+"-spectrum,s-diff,normalized"+ftxt+".ps";
  wpdf, w, shotid+"-s-diff,normalized"+ftxt;
  wpdf, w+1, shotid+"-spectrum,s-diff,normalized"+ftxt;

  //p-channel difference images
  w+= 2;
  pdiff= img_div(img_sub(ch2P, ch1P), img_add(ch1P, ch2P, average= 1));
  //  print, "Comparing pdiff0 with _SNPHI0: ", allof(*pdiff0.data == *_SNPHI0.data);
  if(filtered == 2) {
    pmask= dvv_GenFilterMaskSet(pdiff, _(_DVV_CH1P_FILTER, _DVV_CH2P_FILTER, _DVV_CHP_FILTER));
    //    print, "Comparing pmask with _mp: ", allof(pmask == _mp);
    //    print, "Comparing pdiff0 with _SNPHI0: ", allof(*pdiff0.data == *_SNPHI0.data);
    pdiff= dvv_Filter(pdiff, pmask, wndw= "hanning");
    //   print, "Comparing pdiff with _SNPHI: ", allof(*pdiff.data == *_SNPHI.data);
    ftxt=",filtered-2";
    ftxt1=", filtered-type-2";
  }
  dvv_DisplayImg, pdiff, w= w, -1.5, 1.5, pltitl=shotid+", p-channel difference, normalized"+ftxt1;
  dvv_DisplayImg, dvv_PSD(pdiff, wndw= "hanning"), lgscl= 1, 1e-8, 1e-5,
    w= w+1, pltitl=shotid+", spectrum, p-channel difference, normalized"+ftxt1;
  plsys, 1; limits, -0.2, 0.2, -0.2, 0.2;
  limsReset;
  //wprt, w, shotid+"-p-diff,normalized"+ftxt+".ps";
  //wprt, w+1, shotid+"-spectrum,p-diff,normalized"+ftxt+".ps";
  wpdf, w, shotid+"-p-diff,normalized"+ftxt;
  wpdf, w+1, shotid+"-spectrum,p-diff,normalized"+ftxt;
}

func dvv_ImprintPlots(cmin, cmax)
{
  vv= dvv_Vel(dtpUWf);
  pldefault, width= 3, color= "black";
  dvv_DisplayImg, vv, cmin, cmax, w= next_window();
  palette, "yarg.gp";
  pldefault, width= 3, color= "red";
  plsys, 1; plbox, [-1,1,-1,1]*50;
  pldefault, width= 3, color= "black";
  dvv_DisplayImg, img_extract(vv, [-1,1,-1,1]*50), cmin, cmax, w= next_window(), noresample= 1;
  palette, "yarg.gp";
  pldefault, width= 6, color= "red";
  plsys, 1; plbox, [-1,1,-1,1]*50;
  limits, -50, 50, -50, 50;
  pldefault, width= 3, color= "black";
}

func dvv_window(n)
{
  extern L;

  L--;
  cm= current_mouse();
  print, cm;
  if (L) after, 0.5, dvv_window;
}
