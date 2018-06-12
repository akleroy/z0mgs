pro compare_photometry_to_s4g

  s4g = read_ipac_table('~/idl/galbase/gal_data/s4g_ipac_table.txt')
  s4g_pgc = get_pgc(s4g.object)

; 4.11.3
  corrfac_irac1 = 0.91*1.01
  corrfac_irac2 = 0.94*1.01

  s4g_irac1_flux = 3631.*10.^(-1.*s4g.mabs1/2.5)* $
                   (10./(s4g.dmean*10.^6))^2.*corrfac_irac1
  s4g_irac2_flux = 3631.*10.^(-1.*s4g.mabs2/2.5)* $
                   (10./(s4g.dmean*10.^6))^2.*corrfac_irac2
  
  us = mrdfits('../measurements/z0mgs_photometry.fits', 1, h)
  band = strcompress(us.band, /rem)

  n_s4g = n_elements(s4g_pgc)
  us_wise1_ind = -1 + lindgen(n_s4g)*0L
  for ii = 0, n_s4g-1 do $
     us_wise1_ind[ii] = where(us.pgc eq s4g_pgc[ii] and band eq 'WISE1')
  use_w1 = where(us_wise1_ind ne -1)

  us_wise2_ind = -1 + lindgen(n_s4g)*0L
  for ii = 0, n_s4g-1 do $
     us_wise2_ind[ii] = where(us.pgc eq s4g_pgc[ii] and band eq 'WISE2')
  use_w2 = where(us_wise2_ind ne -1)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PRINT COMPARISONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  print, ""
  
  us_to_s4g_wise1 = $
     us[us_wise1_ind[use_w1]].mean[1]/s4g_irac1_flux[use_w1]
    
  vec = us_to_s4g_wise1
  vec = vec[where(finite(vec))]
  vec = vec[sort(vec)]
  n = n_elements(vec)
  print, "WISE 1 - using only integral inside r25:"
  print, "... 5-16-50-84-95: "
  print, vec[0.05*n], vec[0.16*n], vec[0.5*n], vec[0.84*n], vec[0.95*n]

  us_to_s4g_wise1 = $
     ((us.mean[1]+us.med[5]-us.med[1])[us_wise1_ind[use_w1]])/s4g_irac1_flux[use_w1]
  
  vec = us_to_s4g_wise1
  vec = vec[where(finite(vec))]
  vec = vec[sort(vec)]
  n = n_elements(vec)
  print, "WISE 1 - adding the median outside r25:"
  print, "... 5-16-50-84-95: "
  print, vec[0.05*n], vec[0.16*n], vec[0.5*n], vec[0.84*n], vec[0.95*n]

  us_to_s4g_wise1_chisq = $
     (us[us_wise1_ind[use_w1]].mean[1] - s4g_irac1_flux[use_w1])/us[us_wise1_ind[use_w1]].unc_conf[1]
  
  vec = us_to_s4g_wise1_chisq
  vec = vec[where(finite(vec))]
  vec = vec[sort(vec)]
  n = n_elements(vec)
  print, "WISE 1 - difference, normalized by uncertainty:"
  print, "... 5-16-50-84-95: "
  print, vec[0.05*n], vec[0.16*n], vec[0.5*n], vec[0.84*n], vec[0.95*n]

  us_to_s4g_wise2 = $
     us[us_wise2_ind[use_w2]].mean[1]/s4g_irac2_flux[use_w2]

  vec = us_to_s4g_wise2
  vec = vec[where(finite(vec))]
  vec = vec[sort(vec)]
  n = n_elements(vec)
  print, "WISE 2 - using only integral inside r25:"
  print, "... 5-16-50-84-95: "
  print, vec[0.05*n], vec[0.16*n], vec[0.5*n], vec[0.84*n], vec[0.95*n]

  us_to_s4g_wise2 = $
     ((us.mean[1]+us.med[5]-us.med[1])[us_wise2_ind[use_w2]])/s4g_irac2_flux[use_w2]

  vec = us_to_s4g_wise2
  vec = vec[where(finite(vec))]
  vec = vec[sort(vec)]
  n = n_elements(vec)
  print, "WISE 2 - Adding the median outside r25:"
  print, "... 5-16-50-84-95: "
  print, vec[0.05*n], vec[0.16*n], vec[0.5*n], vec[0.84*n], vec[0.95*n]

  us_to_s4g_wise2_chisq = $
     (us[us_wise2_ind[use_w2]].mean[1] - s4g_irac2_flux[use_w2])/us[us_wise2_ind[use_w2]].unc_conf[1]
  
  vec = us_to_s4g_wise2_chisq
  vec = vec[where(finite(vec))]
  vec = vec[sort(vec)]
  n = n_elements(vec)
  print, "WISE 2 - difference, normalized by uncertainty:"
  print, "... 5-16-50-84-95: "
  print, vec[0.05*n], vec[0.16*n], vec[0.5*n], vec[0.84*n], vec[0.95*n]


; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT SCALINGS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  plot, findgen(10), xtitle='!6'

  psfile = '../plots/s4g_z0mgs_wise1.eps'
  pnfile = '../plots/s4g_z0mgs_wise1.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
    
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n S4G Flux [Jy]' $
     , ytitle='!6log!d10!n Z0MGS "Hybrid" Flux [Jy]' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[-4., 2.], yrange=[-4., 2.0] $
     , /ystyle, /xstyle

  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')

  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.5*[1,1], lines=1, color=cgcolor('charcoal')

  xfid = findgen(1001)/1000.*20 - 10.
  oplot, xfid, xfid, thick=5, color=cgcolor('gray')

  oplot, alog10(s4g_irac1_flux[use_w1]) $
         , alog10(us[us_wise1_ind[use_w1]].mean[1]) $
         , ps=cgsymcat('opencircle'), color=cgcolor('charcoal')

  al_legend, /top, /left $
             , box=1, clear=1 $
             , lines=-99 $
             , background=cgcolor('lightgray') $
             , textcolor=cgcolor('black') $
             , ['!6WISE1 vs. IRAC1'] $
             , charthick=3.0, charsize=1.25

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile+' &'

  plot, findgen(10), xtitle='!6'

  psfile = '../plots/s4g_z0mgs_wise2.eps'
  pnfile = '../plots/s4g_z0mgs_wise2.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
    
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n S4G Flux [Jy]' $
     , ytitle='!6log!d10!n Z0MGS "Hybrid" Flux [Jy]' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[-4., 2.], yrange=[-4., 2.0] $
     , /ystyle, /xstyle

  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')

  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.5*[1,1], lines=1, color=cgcolor('charcoal')

  xfid = findgen(1001)/1000.*20 - 10.
  oplot, xfid, xfid, thick=5, color=cgcolor('gray')

  oplot, alog10(s4g_irac2_flux[use_w2]) $
         , alog10(us[us_wise2_ind[use_w2]].mean[1]) $
         , ps=cgsymcat('opencircle'), color=cgcolor('charcoal')

  al_legend, /top, /left $
             , box=1, clear=1 $
             , lines=-99 $
             , background=cgcolor('lightgray') $
             , textcolor=cgcolor('black') $
             , ['!6WISE2 vs. IRAC2'] $
             , charthick=3.0, charsize=1.25

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile+' &'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT RATIO HISTOGRAMS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  plot, findgen(10), xtitle='!6'

  psfile = '../plots/s4g_z0mgs_wise1_hist.eps'
  pnfile = '../plots/s4g_z0mgs_wise1_hist.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
    
  wise1_bins = $
     bin_data(alog10(us[us_wise1_ind[use_w1]].mean[1]/s4g_irac1_flux[use_w1]) $
                     , s4g_irac1_flux[use_w1]*0.0+1.0, /nan $
                     , xmin=-0.5, xmax=0.5, binsize=0.025)

  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Z0MGS/S4G Flux' $
     , ytitle='!6Number of Galaxies' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[-0.5, 0.5], yrange=[0., 500.0] $
     , /ystyle, /xstyle

  for ii = -100, 100 do $
     oplot, ii*0.1*[1,1], [-1d6, 1d6], lines=1, color=cgcolor('charcoal')

  for ii = -100, 100 do $
     oplot, [-1d6, 1d6], ii*50.*[1,1], lines=1, color=cgcolor('charcoal')
  
  histplot $
     , wise1_bins.xmid, ((wise1_bins.counts)) $
     , /overplot $
     , lthick=3, /fill $
     , fcolor=cgcolor('lightgray') $
     , lcolor=cgcolor('charcoal')

  al_legend, /top, /left $
             , box=1, clear=1 $
             , lines=-99 $
             , background=cgcolor('lightgray') $
             , textcolor=cgcolor('black') $
             , ['!6WISE1 vs. IRAC1'] $
             , charthick=3.0, charsize=1.25

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile+' &'

  plot, findgen(10), xtitle='!6'

  psfile = '../plots/s4g_z0mgs_wise2_hist.eps'
  pnfile = '../plots/s4g_z0mgs_wise2_hist.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
    
  wise2_bins = $
     bin_data(alog10(us[us_wise2_ind[use_w2]].mean[1]/s4g_irac2_flux[use_w2]) $
                     , s4g_irac2_flux[use_w2]*0.0+1.0, /nan $
                     , xmin=-0.5, xmax=0.5, binsize=0.025)

  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Z0MGS/S4G Flux' $
     , ytitle='!6Number of Galaxies' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[-0.5, 0.5], yrange=[0., 500.0] $
     , /ystyle, /xstyle

  for ii = -100, 100 do $
     oplot, ii*0.1*[1,1], [-1d6, 1d6], lines=1, color=cgcolor('charcoal')

  for ii = -100, 100 do $
     oplot, [-1d6, 1d6], ii*50.*[1,1], lines=1, color=cgcolor('charcoal')
  
  histplot $
     , wise2_bins.xmid, ((wise2_bins.counts)) $
     , /overplot $
     , lthick=3, /fill $
     , fcolor=cgcolor('lightgray') $
     , lcolor=cgcolor('charcoal')

;  xfid = findgen(1001)/1000.-0.5
;  yfid = alog10(exp(-1.*(xfid-0.0)^2/2./0.1^2.)*max(wise2_bins.counts))*0.9
;  oplot, xfid, yfid $
;         , thick=3, color=cgcolor('black'), lines=2

  al_legend, /top, /left $
             , box=1, clear=1 $
             , lines=-99 $
             , background=cgcolor('lightgray') $
             , textcolor=cgcolor('black') $
             , ['!6WISE2 vs. IRAC2'] $
             , charthick=3.0, charsize=1.25

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile+' &'

  stop

end
