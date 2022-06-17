pro plot_mainseq
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  index15 = mrdfits('../measurements/delivery_index_gauss15.fits',1,h15)
  z0mgs_x = index15.logmass
  z0mgs_y = index15.logsfr - index15.logmass
; super restrictive:
;  ind = where(index15.has_wise4 and index15.has_fuv and
;  abs(index15.gb_deg gt 30.) and index15.dist_mpc le 50.)

  ind = where(index15.has_wise4 and index15.has_fuv)
  z0mgs_x = z0mgs_x[ind]
  z0mgs_y = z0mgs_y[ind] 
  ind = sort(z0mgs_x)
  z0mgs_x = z0mgs_x[ind]
  z0mgs_y = z0mgs_y[ind]
  
  @constants.bat
  lsun_3p4 = 1.83d18
  restore, '../gswlc/gswlc_data.idl'  

  gws_ind = $
     where(gws_logmstar gt 0 $
           ;and w1_lum gt 0 $
           ;and gws_w1 gt 5.*gws_ew1 $
           ;and gws_w3 gt 5.*gws_ew3 $
           ;and gws_w4 gt 5.*gws_ew4 $
           ;and gws_nuv gt 5.*gws_enuv $
           ;and gws_fuv gt 5.*gws_efuv $
           and gws_flagsed eq 0 $
           and gws_z le 0.05 $
          )

  gws_x = gws_logmstar[gws_ind]
  gws_y = (gws_logsfrsed-gws_logmstar)[gws_ind]
  gws_x = gws_x[gws_ind]
  gws_y = gws_y[gws_ind]   
  ind = sort(gws_x)
  gws_x = gws_x[ind]
  gws_y = gws_y[ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CALCULATIONS  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  ythresh = -11

  med_ind = where(z0mgs_y gt ythresh)
  z0mgs_x_for_med = z0mgs_x[med_ind]
  z0mgs_y_for_med = z0mgs_y[med_ind]
  z0mgs_med_x = z0mgs_x_for_med
  z0mgs_med_y = z0mgs_x_for_med*!values.f_nan
  z0mgs_lo_y = z0mgs_x_for_med*!values.f_nan
  z0mgs_hi_y = z0mgs_x_for_med*!values.f_nan
  n_z0mgs = n_elements(z0mgs_med_x)
  win = 250
  for kk = win, n_z0mgs-win-1 do begin
     lo = kk-win
     hi = kk+win
     vec = z0mgs_y_for_med[lo:hi]
     vec = vec[sort(vec)]
     z0mgs_med_y[kk] = median(vec)
     z0mgs_lo_y[kk] = vec[(win*2+1)*0.16]
     z0mgs_hi_y[kk] = vec[(win*2+1)*0.84]
  endfor

  med_ind = where(gws_y gt ythresh)
  gws_x_for_med = gws_x[med_ind]
  gws_y_for_med = gws_y[med_ind]
  gws_med_x = gws_x[med_ind]
  gws_med_y = gws_x_for_med*!values.f_nan
  gws_lo_y = gws_x_for_med*!values.f_nan
  gws_hi_y = gws_x_for_med*!values.f_nan
  n_gws = n_elements(gws_med_x)
  win = 250
  for kk = win, n_gws-win-1 do begin
     lo = kk-win
     hi = kk+win
     vec = gws_y_for_med[lo:hi]
     vec = vec[sort(vec)]
     gws_med_y[kk] = median(vec)
     gws_lo_y[kk] = vec[(win*2+1)*0.16]
     gws_hi_y[kk] = vec[(win*2+1)*0.84]
  endfor

  xnorm = 10.0

  fit_ind = where(z0mgs_y gt ythresh and finite(z0mgs_x) and $
                  z0mgs_x gt 9.5 and z0mgs_x le 11.0)
  sixlin, (z0mgs_x[fit_ind]-xnorm), z0mgs_y[fit_ind], z0mgs_a, z0mgs_sa, z0mgs_b, z0mgs_sb

  z0mgs_pred = (z0mgs_x-xnorm)*z0mgs_b[0]+z0mgs_a[0]
  z0mgs_resid = z0mgs_y-z0mgs_pred
  z0mgs_active = where(z0mgs_y gt -11 and z0mgs_x gt 9.5 and z0mgs_x le 11.0)

  print, '--------------------------------'
  print, 'z0mgs star forming main sequence'
  print, 'Norm at 1e10 Msun: ', z0mgs_a[0]
  print, 'Slope: ', z0mgs_b[0]
  print, 'Scatter: ', mad(z0mgs_resid[z0mgs_active])
  print, '--------------------------------'

  gws_fit_ind = where(gws_y gt -11.0 and finite(gws_x) and $
                      gws_x gt 9.5 and gws_x le 11.0)
  sixlin, (gws_x[gws_fit_ind]-xnorm), gws_y[gws_fit_ind], gws_a, gws_sa, gws_b, gws_sb
    

  gws_pred = (gws_x-xnorm)*gws_b[0]+gws_a[0]
  gws_resid = gws_y-gws_pred
  gws_active = where(gws_y gt -11 and gws_x gt 9.5 and gws_x le 11.0)

  print, '--------------------------------'
  print, 'gwslc star forming main sequence'
  print, 'Norm at 1e10 Msun: ', gws_a[0]
  print, 'Slope: ', gws_b[0]
  print, 'Scatter: ', mad(gws_resid[gws_active])
  print, '--------------------------------'

  ms_fid = findgen(101)/100.*3.0+8.0

; speagle et al. 2014
  t = 13.7
  s14_fid = (0.84 -0.026*t)*ms_fid-(6.51-0.11*t) - ms_fid

; z0mgs fit
  z0mgs_fid = (ms_fid-xnorm)*z0mgs_b[0]+z0mgs_a[0]

; gwslc fit
  gws_fid = (ms_fid-xnorm)*gws_b[0]+gws_a[0]
 
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GRID
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  gws_grid = grid_data(gws_x, gws_y $
                       , xaxis_out = xaxis_gws, yaxis_out = yaxis_gws $
                       , ymin=-12.0, ymax=-8.5, binsize_x=0.025 $
                       , xmin=9.0, xmax=11.5, binsize_y=0.025)
  
  psf = psf_gaussian(fwhm=[7], npix=21, /norm)
  gws_smooth = convol(gws_grid*1.0, psf, /edge_zero)
  gws_smooth /= total(gws_smooth, /nan)*(0.025*0.025)

  gws_all = grid_data(gws_logmstar $
                      , (gws_logsfrsed-gws_logmstar) $
                      , xaxis_out = xaxis_gws, yaxis_out = yaxis_gws $
                      , ymin=-12.0, ymax=-8.5, binsize_x=0.025 $
                      , xmin=9.0, xmax=11.5, binsize_y=0.025)

  z0mgs_grid = grid_data(z0mgs_x, z0mgs_y $
                         , xaxis_out = xaxis_z0mgs, yaxis_out = yaxis_z0mgs $
                         , ymin=-12.0, ymax=-8.5, binsize_x=0.025 $
                         , xmin=9.0, xmax=11.5, binsize_y=0.025 $
                         , /nan)
  
  psf = psf_gaussian(fwhm=[7], npix=21, /norm)
  z0mgs_smooth = convol(z0mgs_grid*1.0, psf, /edge_zero)
  z0mgs_smooth /= total(z0mgs_smooth, /nan)*(0.025*0.025)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; z0mgs DENSITY ONLY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  xmin = 9.0
  xmax = 11.5
  ymin = -12.0
  ymax = -8.5

  levs = 0.05*2.0^(0.5*findgen(14))

  psfile = '../plots/z0mgs_mainseq_contours.eps'
  pnfile = '../plots/z0mgs_mainseq_contours.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile

  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Stellar Mass [M!d!9n!6!n]' $
     , ytitle='!6log!d10!n SFR / M!d*!n [yr!u-1!n]' $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.25 $
     , xrange=[xmin, xmax], yrange=[ymin, ymax], /xstyle, /ystyle $
     , color=cgcolor('black') $
     , position=[0.2, 0.2, 0.95, 0.95]

  loadct, 3
  reversect
  contour $
     , z0mgs_smooth, xaxis_z0mgs, yaxis_z0mgs, /overplot $
     , lev=levs, /fill

  for jj = -100, 100 do $
     oplot, [-100, 100], jj*0.5*[1,1], lines=1, color=cgcolor('black')
  for jj = -100, 100 do $
     oplot, jj*0.5*[1,1], [-100, 100], lines=1, color=cgcolor('black')

  loadct, 0
  contour $
     , z0mgs_smooth, xaxis_z0mgs, yaxis_z0mgs, /overplot $
     , lev=levs, color=cgcolor('black'), c_thick=2

  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Stellar Mass [M!d!9n!6!n]' $
     , ytitle='!6log!d10!n SFR / M!d*!n [yr!u-1!n]' $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.25 $
     , xrange=[xmin, xmax], yrange=[ymin, ymax], /xstyle, /ystyle $
     , color=cgcolor('black') $
     , position=[0.2, 0.2, 0.95, 0.95] $
     , /noerase

  al_legend, /top, /right $
             , lines=-99, charsize=1.25, charthick=3 $
             , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
             , ['!6Z0MGS density']

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile   

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GSWLC + z0MGS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  psfile = '../plots/z0mgs+gswlc_mainseq_contours.eps'
  pnfile = '../plots/z0mgs+gswlc_mainseq_contours.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile

  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Stellar Mass [M!d!9n!6!n]' $
     , ytitle='!6log!d10!n SFR / M!d*!n [yr!u-1!n]' $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.25 $
     , xrange=[xmin, xmax], yrange=[ymin, ymax], /xstyle, /ystyle $
     , color=cgcolor('black') $
     , position=[0.2, 0.2, 0.95, 0.95]

  loadct, 1
  reversect
  contour $
     , gws_smooth, xaxis_gws, yaxis_gws, /overplot $
     , lev=levs, /fill

  for jj = -100, 100 do $
     oplot, [-100, 100], jj*0.5*[1,1], lines=1, color=cgcolor('black')
  for jj = -100, 100 do $
     oplot, jj*0.5*[1,1], [-100, 100], lines=1, color=cgcolor('black')

  loadct, 0
  contour $
     , z0mgs_smooth, xaxis_z0mgs, yaxis_z0mgs, /overplot $
     , lev=levs, color=cgcolor('firebrick'), c_thick=2

  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Stellar Mass [M!d!9n!6!n]' $
     , ytitle='!6log!d10!n SFR / M!d*!n [yr!u-1!n]' $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.25 $
     , xrange=[xmin, xmax], yrange=[ymin, ymax], /xstyle, /ystyle $
     , color=cgcolor('black') $
     , position=[0.2, 0.2, 0.95, 0.95] $
     , /noerase

  al_legend, /top, /right $
             , lines=-99, charsize=1.25, charthick=3 $
             , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
             , ['!6GSWLC density, Z0MGS contours']

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile   

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MARKED UP VERSION WITH POINTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  xmin = 9.0
  xmax = 11.0
  ymin = -11.0
  ymax = -9.0

  levs = 0.05*2.0^(0.5*findgen(10))

  psfile = '../plots/z0mgs_mainseq.eps'
  pnfile = '../plots/z0mgs_mainseq.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile

  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Stellar Mass [M!d!9n!6!n]' $
     , ytitle='!6log!d10!n SFR / M!d*!n [yr!u-1!n]' $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.25 $
     , xrange=[xmin, xmax] $
     , yrange=[ymin, ymax] $
     , /xstyle, /ystyle $
     , color=cgcolor('black') $
     , position=[0.2, 0.2, 0.95, 0.95]

;  polyx = [xmin, xmin, xmax, xmax, xmin]
;  polyy = [ymin, ythresh, ythresh, ymin, ymin]  
;  polyfill, polyx, polyy, noclip=0, color=cgcolor('lightgray')

  oplot, z0mgs_x, z0mgs_y, ps=cgsymcat('filledcircle') $
         , color=cgcolor('gray'), symsize=0.50

  for jj = -100, 100 do $
     oplot, [-100, 100], jj*0.5*[1,1], lines=1, color=cgcolor('black')
  for jj = -100, 100 do $
     oplot, jj*0.5*[1,1], [-100, 100], lines=1, color=cgcolor('black')

  oplot, z0mgs_med_x, z0mgs_med_y $
         , lines=0, color=cgcolor('salmon'), thick=7
  oplot, z0mgs_med_x, z0mgs_lo_y $
         , lines=0, color=cgcolor('salmon'), thick=7
  oplot, z0mgs_med_x, z0mgs_hi_y $
         , lines=0, color=cgcolor('salmon'), thick=7

; z0mgs
  oplot, ms_fid, z0mgs_fid, color=cgcolor('firebrick'), thick=15, lines=0

; Legend
  al_legend, /top, /right $
             , lines=[-99,0,0], charsize=1.25, charthick=3 $
             , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
             , ['!6Z0MGS points','Sliding percentile','Fit'] $
             , color=[cgcolor('black'), cgcolor('salmon'), cgcolor('firebrick')] $
             , textcolor=cgcolor('black') $
             , thick=7, pspacing=1.0

  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Stellar Mass [M!d!9n!6!n]' $
     , ytitle='!6log!d10!n SFR / M!d*!n [yr!u-1!n]' $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.25 $
     , xrange=[xmin, xmax] $
     , yrange=[ymin, ymax] $
     , /xstyle, /ystyle $
     , color=cgcolor('black') $
     , position=[0.2, 0.2, 0.95, 0.95] $
     , /noerase

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile   

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SCHEMATIC
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  xmin = 9.0
  xmax = 11.0
  ymin = -11.0
  ymax = -9.0

  levs = 0.05*2.0^(0.5*findgen(10))

  psfile = '../plots/markup_mainseq.eps'
  pnfile = '../plots/markup_mainseq.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile

  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Stellar Mass [M!d!9n!6!n]' $
     , ytitle='!6log!d10!n SFR / M!d*!n [yr!u-1!n]' $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.25 $
     , xrange=[xmin, xmax] $
     , yrange=[ymin, ymax] $
     , /xstyle, /ystyle $
     , color=cgcolor('black') $
     , position=[0.2, 0.2, 0.95, 0.95]

  ind_for_poly = where(finite(z0mgs_med_y))
  polyx = [z0mgs_med_x[ind_for_poly] $
           , reverse(z0mgs_med_x[ind_for_poly]) $
           , z0mgs_med_x[ind_for_poly[0]]]
  y_for_poly = reverse(z0mgs_lo_y[ind_for_poly])
  nan_ind = where(finite(y_for_poly) eq 0, nan_ct)
  if nan_ct gt 0 then y_for_poly[nan_ind] = -11.0
  polyy = [z0mgs_hi_y[ind_for_poly], y_for_poly, z0mgs_hi_y[0]]
  polyfill, polyx, polyy, noclip=0, color=cgcolor('salmon')

  for jj = -100, 100 do $
     oplot, [-100, 100], jj*0.5*[1,1], lines=1, color=cgcolor('black')
  for jj = -100, 100 do $
     oplot, jj*0.5*[1,1], [-100, 100], lines=1, color=cgcolor('black')

  oplot, gws_med_x, gws_med_y $
         , lines=0, color=cgcolor('navy'), thick=5
  oplot, gws_med_x, gws_lo_y $
         , lines=0, color=cgcolor('navy'), thick=5
  oplot, gws_med_x, gws_hi_y $
         , lines=0, color=cgcolor('navy'), thick=5

; z0mgs
  oplot, ms_fid, z0mgs_fid, color=cgcolor('black'), thick=12, lines=0
  oplot, ms_fid, z0mgs_fid, color=cgcolor('firebrick'), thick=7, lines=0

; gws
  oplot, ms_fid, gws_fid, color=cgcolor('black'), thick=12, lines=2
  oplot, ms_fid, gws_fid, color=cgcolor('royalblue'), thick=7, lines=2

; speagle et al. 2014
  oplot, ms_fid, s14_fid, color=cgcolor('black'), thick=12, lines=3
  oplot, ms_fid, s14_fid, color=cgcolor('darkgreen'), thick=7, lines=3

; Legend
  al_legend, /top, /right $
             , lines=-99, charsize=1.25, charthick=3 $
             , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
             , ['!6Z0MGS data','GSWLC (Salim et al. 2018)', 'Speagle et al. (2014)'] $
             , textcolor=[cgcolor('salmon'),cgcolor('royalblue'),cgcolor('darkgreen')]

  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Stellar Mass [M!d!9n!6!n]' $
     , ytitle='!6log!d10!n SFR / M!d*!n [yr!u-1!n]' $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.25 $
     , xrange=[xmin, xmax] $
     , yrange=[ymin, ymax] $
     , /xstyle, /ystyle $
     , color=cgcolor('black') $
     , position=[0.2, 0.2, 0.95, 0.95] $
     , /noerase

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile   

  stop


end
