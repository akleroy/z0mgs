pro plot_gswlc_sfruv $
   , show=show

  thresh = 10
  ston = 3.0
  @constants.bat
  lsun_3p4 = 1.83d18
  restore, '../gswlc/gswlc_data.idl'

  nu_fuv = c/(154.d-9*1d2)
  nu_nuv = c/(231.d-9*1d2)
  nu_w1 = c/(3.4d-6*1d2)
  nu_w2 = c/(4.5d-6*1d2)
  nu_w3 = c/(12.d-6*1d2)
  nu_w4 = c/(22.d-6*1d2)

  plot, findgen(10), xtitle='!6'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE SUBSET OF DATA TO USE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  ind = where(gws_logmstar gt 0 $
              and w1_lum gt 0 $
              and gws_w1 gt ston*gws_ew1 $
              and gws_w3 gt ston*gws_ew3 $
              and gws_w4 gt ston*gws_ew4 $
              and gws_nuv gt ston*gws_enuv $
              and gws_fuv gt ston*gws_efuv $
              and gws_flagsed eq 0 $
             )
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ESTIMATE THE COEFFICIENTS TO CONVERT BROAD BAND TO STAR FORMATION RATE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  coef_fuv = ((10.d^(gws_logsfrsed))) / (nu_fuv*fuv_lum*10.^(gws_afuv/2.5))
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; HISTOGRAMS OF THE FAR UV CALIBRATION TERM
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  c_salim = 10.d^(-43.45)
  c_ke12 = 10.d^(-43.35)

  xmin = -44.0
  xmax = -43.0
  binsize = 0.025

  histpsfile = '../plots/gswlc_cfuv.eps'
  histpnfile = '../plots/gswlc_cfuv.png'
  
  coef = coef_fuv
  tag = '!6FUV+A!dFUV!n'

  ssfr = gws_ssfr[ind]
  vec = alog10(coef_fuv[ind])

  vec = vec[where(ssfr lt -11.)]
  vec = vec[sort(vec)]
  n = n_elements(vec)

  red_med = vec[long(0.5*n)]
  red_lo = vec[long(0.16*n)]
  red_hi = vec[long(0.84*n)]

  red_bins = bin_data(vec, vec*0.0+1.0 $
                      , xmin=xmin, xmax=xmax, binsize=binsize, /nan)
  
  ssfr = gws_ssfr[ind]
  vec = alog10(coef_fuv[ind])

  vec = vec[where(ssfr gt -11.)]
  vec = vec[sort(vec)]
  n = n_elements(vec)

  blue_med = vec[long(0.5*n)]
  blue_lo = vec[long(0.16*n)]
  blue_hi = vec[long(0.84*n)]

  blue_bins = bin_data(vec, vec*0.0+1.0 $
                       , xmin=xmin, xmax=xmax, binsize=binsize, /nan)

  c_z0mgs = blue_med
  
; REPORT OUR BEST ESTIMATES
  print, 'FUV SFR calibration for star forming galaxies: '
  print, blue_med
  print, blue_lo, ' to ', blue_hi
  
  print, 'FUV SFR calibration for RED galaxies: '
  print, red_med
  print, red_lo, ' to ', red_hi
  
; HISTOGRAM
  ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
      , file=histpsfile
  
  loadct, 0
  ylo = 0.0
  yhi = 0.25
  plot $
     , [0], [0], /nodata $`
     , xtitle='!6log!d10!n C [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
     , ytitle='!6Fraction of Galaxies' $
     , xthick=5, ythick=5, charthick=3, charsize=1.25 $
     , xrange=[xmin, xmax], yrange=[ylo, yhi] $
     , /xstyle
  
  lo = blue_lo
  hi = blue_hi
  polyfill, [lo, hi, hi, lo, lo], [ylo, ylo, yhi, yhi, ylo], /fill, /clip $
            , color=cgcolor('lightgray')
  oplot, (c_z0mgs)*[1,1], [-100, 100], lines=0, thick=10, color=cgcolor('charcoal')   
  oplot, alog10(c_salim)*[1,1], [-100, 100], lines=1, thick=10, color=cgcolor('firebrick')
  oplot, alog10(c_ke12)*[1,1], [-100, 100], lines=2, thick=10, color=cgcolor('darkgreen') 

  for kk = -100, 100 do $
     oplot, kk*binsize*10.*[1,1]+xmin, [-10, 10], lines=1, color=cgcolor('charcoal')
  
  for kk = -100, 100 do $
     oplot, [-1e6, 1e6], kk*0.05*[1,1], lines=1, color=cgcolor('charcoal')
  
  histplot $
     , red_bins.xmid, ((red_bins.counts/total(red_bins.counts,/nan) > (0.))) $
     , /overplot $
     , lthick=3, /fill $
     , fcolor=cgcolor('salmon') $
     , lcolor=cgcolor('firebrick')

  histplot $
     , blue_bins.xmid, ((blue_bins.counts/total(blue_bins.counts,/nan) > (0.))) $
     , /overplot $
     , lthick=3, /fill, /fline, forient=45, /nobars $
     , fcolor=cgcolor('cornflowerblue') $
     , lcolor=cgcolor('navy')
  
  al_legend, /top, /left, lines=-99, charsize=1.25, charthick=3 $
             , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
             , textcolor=[cgcolor('black'),cgcolor('royalblue'),cgcolor('firebrick')] $
             , [tag,'high sSFR', 'low sSFR']

  al_legend, /bottom, /left, lines=[1,2,0], charsize=1.25, charthick=3 $
             , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
             , textcolor=[cgcolor('firebrick'),cgcolor('darkgreen'), cgcolor('charcoal')] $
             , color=[cgcolor('firebrick'),cgcolor('darkgreen'), cgcolor('charcoal')] $
             , ['S07','KE12','GSWLC'], pspacing=1.0, thick=7

  plot $
     , [0], [0], /nodata $`
     , xtitle='!6log!d10!n C [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
     , ytitle='!6Fraction of Galaxies' $
     , xthick=5, ythick=5, charthick=3, charsize=1.25 $
     , xrange=[xmin, xmax], yrange=[ylo, yhi] $
     , /xstyle, /noerase
  
  ps, /xw
  if keyword_set(show) then $
     spawn, 'evince '+histpsfile+' &'
  spawn, 'convert -density 300x300 '+histpsfile+' '+histpnfile
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT THE COEFFICIENT ON THE EXTINCTION CORRECTED CASE VS MASS AND SSFR
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  mstar_psfile = '../plots/gswlc_cfuv_vs_mass.eps'
  mstar_pnfile = '../plots/gswlc_cfuv_vs_mass.png'
  
  ssfr_psfile = '../plots/gswlc_cfuv_vs_ssfr.eps'
  ssfr_pnfile = '../plots/gswlc_cfuv_vs_ssfr.png'

  coef = coef_fuv
  tag = '!6FUV+A!dFUV!n'
  
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; STELLAR MASS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  ps, /def, /ps, xs=5, ys=3.5, /color, /encaps, file=mstar_psfile
  
  loadct, 0
  reversect
  x = gws_logmstar[ind]
  y = alog10(coef[ind])
  ymin = -44.0
  ymax = -43.0
  grid = $
     grid_data(x, y, /nan $
               , xaxis_out = x_axis, yaxis_out = y_axis $
               , xmin=9., xmax=12., binsize_x=0.05 $
               , ymin=ymin, ymax=ymax, binsize_y=0.025 $
              )
  disp, alog10(grid), x_axis, y_axis $
        , ps=3 $
        , xrange=[9., 12.], yrange=[ymin, ymax] $
        , xthick=5, ythick=5, charsize=1.25, charthick=3 $
        , xtitle='!6log!d10!n Stellar Mass from CIGALE [M!d!9n!6!n]' $
        , ytitle='!6log!d10!n C [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
        , /xstyle, /ystyle, reserve=50, color=cgcolor('black', 255) $
        , max=3.0, position=[0.2,0.2,0.95, 0.95]

  oplot, [-100., 100], (c_z0mgs)*[1,1],  lines=0, thick=10, color=cgcolor('charcoal')   
  oplot,[-100., 100],  alog10(c_salim)*[1,1], lines=1, thick=10, color=cgcolor('firebrick')
  oplot,[-100., 100],  alog10(c_ke12)*[1,1], lines=2, thick=10, color=cgcolor('darkgreen') 

  for jj = -100, 100 do $
     oplot, jj*0.5*[1,1], [-1d6, 1d6], lines=1, color=cgcolor('black')
  for jj = -100, 100 do $
     oplot, [-1d6, 1d6], jj*0.1*[1,1]-44., lines=1, color=cgcolor('black')

  contour, /overplot, alog10(grid), x_axis, y_axis $
           , lev=[1., 1.5, 2.0, 2.5, 3.0], c_color=cgcolor('black')
  
  xfid = findgen(101)/100.*100.
  oplot, xfid, xfid*0.0, thick=13, color=cgcolor('charcoal')
  oplot, xfid, xfid*0.0, thick=5, color=cgcolor('gray')
  
  bins = bin_data(x, y, /nan, xmin=9.5, xmax=11.5, binsize=0.05)
  oploterror, bins.xmid $
              , bins.ymed, bins.ymad $
              , color=cgcolor('red') $
              , psym=cgsymcat('filledsquare') $
              , errthick=5, /nohat   
  
  al_legend, /bottom, /right $
             , box=1, clear=1, lines=-99, background=['lightgray'] $
             , charthick=3, charsize=1.25 $
             , tag
  
  ps, /xw
  if keyword_set(show) then $
     spawn, 'evince '+mstar_psfile+' &'
  spawn, 'convert -density 300x300 '+mstar_psfile+' '+mstar_pnfile

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; SFR/MSTAR
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  ps, /def, /ps, xs=5, ys=3.5, /color, /encaps, file=ssfr_psfile
  
  loadct, 0
  reversect
  x = gws_ssfr[ind]
  xlo = -12.0
  xhi = -8.0
  y = alog10(coef[ind])
  ymin = -44.0
  ymax = -43.0
  grid = $
     grid_data(x, y, /nan $
               , xaxis_out = x_axis, yaxis_out = y_axis $
               , xmin=xlo, xmax=xhi, binsize_x=0.05 $
               , ymin=ymin, ymax=ymax, binsize_y=0.025 $
              )
  disp, alog10(grid), x_axis, y_axis $
        , ps=3 $
        , xrange=[xlo,xhi], yrange=[ymin, ymax] $
        , xthick=5, ythick=5, charsize=1.25, charthick=3 $
        , xtitle='!6log!d10!n SFR/M!d*!n from CIGALE [yr!u-1!n]' $
        , ytitle='!6log!d10!n C [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
        , /xstyle, /ystyle, reserve=50, color=cgcolor('black', 255) $
        , max=3.0, position=[0.2,0.2,0.95, 0.95]

  oplot, [-100., 100], (c_z0mgs)*[1,1],  lines=0, thick=10, color=cgcolor('charcoal')   
  oplot,[-100., 100],  alog10(c_salim)*[1,1], lines=1, thick=10, color=cgcolor('firebrick')
  oplot,[-100., 100],  alog10(c_ke12)*[1,1], lines=2, thick=10, color=cgcolor('darkgreen') 

  for jj = -100, 100 do $
     oplot, jj*0.5*[1,1], [-1d6, 1d6], lines=1, color=cgcolor('black')
  for jj = -100, 100 do $
     oplot, [-1d6, 1d6], jj*0.1*[1,1]-44., lines=1, color=cgcolor('black')

  contour, /overplot, alog10(grid), x_axis, y_axis $
           , lev=[1., 1.5, 2.0, 2.5, 3.0], c_color=cgcolor('black')
  
  xfid = findgen(101)/100.*100.
  oplot, xfid, xfid*0.0, thick=13, color=cgcolor('charcoal')
  oplot, xfid, xfid*0.0, thick=5, color=cgcolor('gray')
  
  bins = bin_data(x, y, /nan, xmin=-10.75, xmax=-9.0, binsize=0.05)
  oploterror, bins.xmid $
              , bins.ymed, bins.ymad $
              , color=cgcolor('red') $
              , psym=cgsymcat('filledsquare') $
              , errthick=5, /nohat   
  
  al_legend, /bottom, /right $
             , box=1, clear=1, lines=-99, background=['lightgray'] $
             , charthick=3, charsize=1.25 $
             , tag
  
  ps, /xw
  if keyword_set(show) then $
     spawn, 'evince '+ssfr_psfile+' &'
  spawn, 'convert -density 300x300 '+ssfr_psfile+' '+ssfr_pnfile

  stop

end
