pro plot_distance_unc

  tab = mrdfits('~/idl/galbase/gal_data/dist_base.fits',1,h)
  
  gold = where(tab.edd_code eq 1)
  silver = where(tab.edd_code eq 2)
  copper = where(tab.edd_code gt 2)
  lead =  where(tab.edd_code eq -1)

  plot, findgen(101), title='!6Hello'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; JUST PLOT DISTANCE VS DISTANCE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  rat = alog10(tab[copper].leda_bestdist_mpc / $
               tab[copper].edd_dist_mpc)
  print, "For best LEDA distance: "
  print, "... median and mean ratio with EDD: ", median(rat), mean(rat,/nan)
  print, "... m.a.d. and stddev ratio with EDD: ", mad(rat), stddev(rat,/nan)

  psfile = '../plots/edd_vs_leda.eps'
  pnfile = '../plots/edd_vs_leda.png'
  ps, /ps, /def, file=psfile, xs=8, ys=8, /color, /encaps

  plot $
     , [0], [0] $
     , ps=1, symsize=0.5 $
     , xtitle='!6log!d10!n LEDA Best Distance [Mpc]' $
     , ytitle='!6log!d10!n EDD Best Distance [Mpc]' $
     , xrange=[0, 2] $
     , yrange=[0, 2] $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.5

  for ii = -10, 25 do $
     oplot, ii*0.25*[1,1], [-1d6, 1d6], lines=1
  for ii = -10, 25 do $
     oplot, [-1d6, 1d6], ii*0.25*[1,1], lines=1

  xfid = findgen(101)/100.*2.
  
  delta = 0.1 ; dex typical unc.
  poly_x = [xfid, reverse(xfid), xfid[0]]
  poly_y = [xfid+delta, reverse(xfid)-delta, xfid[0]+delta]
  polyfill, poly_x, poly_y, color=cgcolor('gray'), noclip=0B

  oplot, xfid, xfid, lines=0 $
         , color=cgcolor('charcoal') $
         , thick=10

  circle, /fill
  oplot, alog10(tab[copper].leda_dist0_mpc) $
         , alog10(tab[copper].edd_dist_mpc) $
         , ps=8, symsize=0.8 $
         , color=cgcolor('firebrick')
  oplot, alog10(tab[copper].leda_dist0_mpc) $
         , alog10(tab[copper].edd_dist_mpc) $
         , ps=8, symsize=0.45 $
         , color=cgcolor('salmon')

  oplot, alog10(tab[silver].leda_dist0_mpc) $
         , alog10(tab[silver].edd_dist_mpc) $
         , ps=8, symsize=1.1 $
         , color=cgcolor('navy')
  oplot, alog10(tab[silver].leda_dist0_mpc) $
         , alog10(tab[silver].edd_dist_mpc) $
         , ps=8, symsize=0.70 $
         , color=cgcolor('dodgerblue')

  oplot, alog10(tab[gold].leda_dist0_mpc) $
         , alog10(tab[gold].edd_dist_mpc) $
         , ps=8, symsize=1.1 $
         , color=cgcolor('black')
  oplot, alog10(tab[gold].leda_dist0_mpc) $
         , alog10(tab[gold].edd_dist_mpc) $
         , ps=8, symsize=0.70 $
         , color=cgcolor('goldenrod')

  al_legend, /top, /left $
             , box=1, clear=1, background=cgcolor('lightgray') $
             , ['!6Has EDD TRGB Distance', '!6Has EDD Quality Distance', '!6Lacks Either'] $
             , lines=0, psym=[8,8,8] $
             , color=[cgcolor('goldenrod'), cgcolor('dodgerblue'), cgcolor('salmon')] $
             , charthick=3, charsize=1.5 $
             , symsize=[1.5,1.5,1.5]

  al_legend, /bottom, /right $
             , box=1, clear=1, background=cgcolor('lightgray') $
             , ['!6Line Shows Equality !9+!60.1 dex'] $
             , lines=-99 $
             , charthick=3, charsize=1.5

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 150x150 '+psfile+' '+pnfile

  rat = alog10(tab[copper].cf_dist_mpc / $
               tab[copper].edd_dist_mpc)
  print, "For Cosmic Flows distance: "
  print, "... median and mean ratio with EDD: ", median(rat), mean(rat,/nan)
  print, "... m.a.d. and stddev ratio with EDD: ", mad(rat), stddev(rat,/nan)

  psfile = '../plots/edd_vs_cf3.eps'
  pnfile = '../plots/edd_vs_cf3.png'
  ps, /ps, /def, file=psfile, xs=8, ys=8, /color, /encaps

  plot $
     , [0], [0] $
     , ps=1, symsize=0.5 $
     , xtitle='!6log!d10!n Cosmic Flows 3 Distance [Mpc]' $
     , ytitle='!6log!d10!n EDD Best Distance [Mpc]' $
     , xrange=[0, 2] $
     , yrange=[0, 2] $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.5

  for ii = -10, 25 do $
     oplot, ii*0.25*[1,1], [-1d6, 1d6], lines=1
  for ii = -10, 25 do $
     oplot, [-1d6, 1d6], ii*0.25*[1,1], lines=1

  xfid = findgen(101)/100.*2.
  
  delta = 0.1 ; dex typical unc.
  poly_x = [xfid, reverse(xfid), xfid[0]]
  poly_y = [xfid+delta, reverse(xfid)-delta, xfid[0]+delta]
  polyfill, poly_x, poly_y, color=cgcolor('gray'), noclip=0B

  oplot, xfid, xfid, lines=0 $
         , color=cgcolor('charcoal') $
         , thick=10

  circle, /fill
  oplot, alog10(tab[copper].cf_dist_mpc) $
         , alog10(tab[copper].edd_dist_mpc) $
         , ps=8, symsize=0.8 $
         , color=cgcolor('firebrick')
  oplot, alog10(tab[copper].cf_dist_mpc) $
         , alog10(tab[copper].edd_dist_mpc) $
         , ps=8, symsize=0.45 $
         , color=cgcolor('salmon')

  oplot, alog10(tab[silver].cf_dist_mpc) $
         , alog10(tab[silver].edd_dist_mpc) $
         , ps=8, symsize=1.1 $
         , color=cgcolor('navy')
  oplot, alog10(tab[silver].cf_dist_mpc) $
         , alog10(tab[silver].edd_dist_mpc) $
         , ps=8, symsize=0.70 $
         , color=cgcolor('dodgerblue')

  oplot, alog10(tab[gold].cf_dist_mpc) $
         , alog10(tab[gold].edd_dist_mpc) $
         , ps=8, symsize=1.1 $
         , color=cgcolor('black')
  oplot, alog10(tab[gold].cf_dist_mpc) $
         , alog10(tab[gold].edd_dist_mpc) $
         , ps=8, symsize=0.70 $
         , color=cgcolor('goldenrod')

  al_legend, /top, /left $
             , box=1, clear=1, background=cgcolor('lightgray') $
             , ['!6Has EDD TRGB Distance', '!6Has EDD Quality Distance', '!6Lacks Either'] $
             , lines=0, psym=[8,8,8] $
             , color=[cgcolor('goldenrod'), cgcolor('dodgerblue'), cgcolor('salmon')] $
             , charthick=3, charsize=1.5 $
             , symsize=[1.5,1.5,1.5]

  al_legend, /bottom, /right $
             , box=1, clear=1, background=cgcolor('lightgray') $
             , ['!6Line Shows Equality !9+!60.1 dex'] $
             , lines=-99 $
             , charthick=3, charsize=1.5

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 150x150 '+psfile+' '+pnfile

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; NOW CALCULATE THE SCATTER AS HISTOGRAMS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; LEAD

  leda_cf3_rat_lead = $
     alog10(tab[lead].leda_bestdist_mpc / $
            tab[lead].cf_dist_mpc)

  rat_lead = [leda_cf3_rat_lead]

; COPPER

  leda_edd_rat_copper = $
     alog10(tab[copper].leda_bestdist_mpc / $
            tab[copper].edd_dist_mpc)

  leda_cf3_rat_copper = $
     alog10(tab[copper].leda_bestdist_mpc / $
            tab[copper].cf_dist_mpc)

  edd_cf3_copper = $
     alog10(tab[copper].edd_dist_mpc / $
            tab[copper].cf_dist_mpc)

  rat_copper = [leda_edd_rat_copper, leda_cf3_rat_copper, edd_cf3_copper]

; SILVER

  leda_edd_rat_silver = $
     alog10(tab[silver].leda_bestdist_mpc / $
            tab[silver].edd_dist_mpc)

  leda_cf3_rat_silver = $
     alog10(tab[silver].leda_bestdist_mpc / $
            tab[silver].cf_dist_mpc)

  edd_cf3_silver = $
     alog10(tab[silver].edd_dist_mpc / $
            tab[silver].cf_dist_mpc)

  rat_silver = [leda_edd_rat_silver, leda_cf3_rat_silver, edd_cf3_silver]


; GOLD

  leda_edd_rat_gold = $
     alog10(tab[gold].leda_bestdist_mpc / $
            tab[gold].edd_dist_mpc)

  leda_cf3_rat_gold = $
     alog10(tab[gold].leda_bestdist_mpc / $
            tab[gold].cf_dist_mpc)

  edd_cf3_gold = $
     alog10(tab[gold].edd_dist_mpc / $
            tab[gold].cf_dist_mpc)

  rat_gold = [leda_edd_rat_gold, leda_cf3_rat_gold, edd_cf3_gold]

; CALCULATE SOME NUMBERS

  rms_gold = stddev(rat_gold, /nan)
  rms_silver = stddev(rat_silver, /nan)
  rms_copper = stddev(rat_copper, /nan)

  psfile = '../plots/dist_hist.eps'
  pnfile = '../plots/dist_hist.png'
  ps, /ps, /def, file=psfile, xs=8, ys=8, /color, /encaps
 
  !p.multi=[0,1,1]

  plot $
     , [0], [0] $
     , ps=1, symsize=0.5 $
     , xtitle='!6log!d10!n Difference in Distance Estimates [dex]' $
     , ytitle='!6log!d10!n Probability [Galaxies / dex]' $
     , xrange=[-0.5, 0.5] $
     , yrange=[2, 5.5] $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.5

  for ii = -10, 25 do $
     oplot, ii*0.1*[1,1], [-1d6, 1d6], lines=1
  for ii = -10, 25 do $
     oplot, [-1d6, 1d6], ii*0.25*[1,1], lines=1

  binsize = 0.01
  bins = bin_data(rat_copper, rat_copper*0.0+1.0, xmin=-2., xmax=2., binsize=binsize,/nan)
  histplot, bins.xmid, alog10(bins.ysum/(binsize)) $
            , /fill, fcolor=cgcolor('salmon') $
            , lcolor=cgcolor('black') $
            , lthick=3 $
            , /overplot, /nobar

  binsize = 0.01
  bins = bin_data(rat_silver, rat_silver*0.0+1.0, xmin=-2., xmax=2., binsize=binsize,/nan)
  histplot, bins.xmid, alog10(bins.ysum/(binsize)) $
            , /fill, fcolor=cgcolor('dodgerblue') $
            , lcolor=cgcolor('black') $
            , lthick=3 $
            , /overplot, /nobar

  binsize = 0.005
  bins = bin_data(rat_gold, rat_gold*0.0+1.0, xmin=-2., xmax=2., binsize=binsize,/nan)
  histplot, bins.xmid, alog10(bins.ysum/(binsize)) $
            , /fill, fcolor=cgcolor('goldenrod') $
            , lcolor=cgcolor('black') $
            , lthick=3 $
            , /overplot, /nobar
    
  xyouts, [0.], [5.1], '!9+!61!7r!6 range', align=0.5, charthick=3, charsize=1.5
  oplot, [-1.,1.]*rms_copper, 5.0*[1,1], thick=10, color=cgcolor('salmon')
  oplot, [-1.,1.]*rms_silver, 4.9*[1,1], thick=10, color=cgcolor('dodgerblue')
  oplot, [-1.,1.]*rms_gold, 4.8*[1,1], thick=10, color=cgcolor('goldenrod')

  al_legend, /top, /left $
             , box=1, clear=1, background=cgcolor('lightgray') $
             , ['!6Has EDD TRGB Distance', '!6Has EDD Quality Distance', '!6Lacks Either'] $
             , lines=0, psym=[8,8,8] $
             , color=[cgcolor('goldenrod'), cgcolor('dodgerblue'), cgcolor('salmon')] $
             , charthick=3, charsize=1.5 $
             , symsize=[1.5,1.5,1.5]

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 150x150 '+psfile+' '+pnfile  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT A SKETCH OF OUR STRATEGY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  psfile = '../plots/dist_sketch.eps'
  pnfile = '../plots/dist_sketch.png'
  ps, /ps, /def, file=psfile, xs=8, ys=8, /color, /encaps
  
  !p.multi=[0,1,1]

  plot $
     , [0], [0] $
     , ps=1, symsize=0.5 $
     , xtitle='!6log!d10!n Virgocentric-Flow Corrected Velocity' $
     , ytitle='!6log!d10!n Distance Uncertainty [dex]' $
     , xrange=[2., 4.5], /xstyle $
     , yrange=[0.0, 0.5] $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.5
  
  xlo = 0.
  xhi = alog10(3500.)
  ylo = 0.
  yhi = 10.
  polyfill, [xlo, xhi, xhi, xlo, xlo], [ylo, ylo, yhi, yhi, ylo] $
            , noclip=0B, color=cgcolor('lightgray')
  oplot, alog10(3500.)*[1., 1.], [-10., 10.], thick=10, lines=0, color=cgcolor('charcoal')
  xyouts, [2.65], [0.25], '!6Use redshift!cindependent distance', align=0.5, charsize=1.5, charthick=3

  for ii = -10, 25 do $
     oplot, ii*0.25*[1,1], [-1d6, 1d6], lines=1
  for ii = -10, 25 do $
     oplot, [-1d6, 1d6], ii*0.1*[1,1], lines=1

  xyouts, [4.0], [0.25], '!6Use Hubble flow', align=0.5, charsize=1.5, charthick=3

  oplot, alog10(tab.vvir_kms), tab.e_leda_vdist, ps=1

  oplot, [0., 10], 0.125*[1., 1.], thick=10, lines=2, color=cgcolor('salmon')
  xyouts, [2.65], [0.125]-0.004, '!6RMS Other', align=0.5, charsize=1.5, charthick=5, color=cgcolor('black')
;  xyouts, [2.65], [0.125]-0.004, '!6RMS Other', align=0.5, charsize=1.5, charthick=1, color=cgcolor('salmon')

  oplot, [0., 10], 0.06*[1, 1], thick=10, lines=2, color=cgcolor('dodgerblue') 
  xyouts, [2.65], [0.06]-0.004, '!6RMS Quality', align=0.5, charsize=1.5, charthick=5, color=cgcolor('black')
;  xyouts, [2.65], [0.06]-0.004, '!6RMS Quality', align=0.5, charsize=1.5, charthick=1, color=cgcolor('dodgerblue')

  oplot, [0., 10], 0.035*[1, 1], thick=10, lines=2, color=cgcolor('goldenrod') 
  xyouts, [2.65], [0.035]-0.004, '!6RMS TRGB', align=0.5, charsize=1.5, charthick=5, color=cgcolor('black')
;  xyouts, [2.65], [0.035]-0.004, '!6RMS TRGB', align=0.5, charsize=1.5, charthick=1, color=cgcolor('goldenrod')
  
  xyouts, [3.1], [0.40], '!6LEDA Hubble Flow', charsize=1.5, charthick=3, orient=-75.

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 150x150 '+psfile+' '+pnfile  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PRINT SOME UNCERTAINTY RELATED NUMBERS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  print, "Mean uncertainty on low quality values from LEDA: ", mean(tab[copper].e_leda_dist0, /nan)
  print, "Mean uncertainty on quality values from LEDA: ", mean(tab[silver].e_leda_dist0, /nan)
  print, "Mean uncertainty on TRGB values from LEDA: ", mean(tab[gold].e_leda_dist0, /nan)

  stop

end
