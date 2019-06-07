pro motivate_selection $
   , write=do_write

  selection_thresh = 0.1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE FULL GALAXY DATABASE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  tab = gal_data(/all, /full)
  tab = tab[where(tab.dist_mpc le 150)]
  absb = tab.btc-tab.distmod

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MONTE CARLO SELECTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  ngal = n_elements(tab)
  unc_btc = 0.5
  niter = 1000L
  selected = lonarr(ngal)
  absb_cutoff = -18.
  dist_cutoff = 50.
  for ii = 0, niter-1 do begin
     counter, ii, niter, 'Monte ...'
     rand_btc = tab.btc+unc_btc*randomn(seed,ngal)
     rand_dist = tab.dist_mpc*10.^(tab.e_dist_dex*randomn(seed,ngal))
     rand_distmod = 5.*alog10(rand_dist*1d6/10.)
     rand_absb = rand_btc-rand_distmod
     select = rand_absb le absb_cutoff and rand_dist le dist_cutoff
     selected += select
  endfor
  prob = selected*1.0/niter

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT THE SELECTION FUNCTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  plot, findgen(10), xtitle='!6'

  psfile = '../plots/select_wedge.eps'
  pnfile = '../plots/select_wedge.png'
  ps, /ps, /def, file=psfile, xs=5, ys=5, /color, /encaps
  
  !p.multi=[0,1,1]
 
  loadct, 0

  xmin = 1.0
  xmax = 2.0
  ymin = 11.
  ymax = 16.

  plot $
     , [0], [0] $
     , ps=1, symsize=0.5 $
     , xtitle='!6log!d10!n Distance [Mpc]' $
     , ytitle='!6LEDA Corrected Apparent B Magnitude [mag]' $
     , xrange=[xmin, xmax], /xstyle $
;     , yrange=[-4., -1.] $
     , yrange=[ymin, ymax], /ystyle $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.0 $
     , position=[0.15,0.1,0.95,0.85]
 
  xfid = findgen(1000)+1.
  dmfid = 5.*alog10(xfid*1d6/10.)
  cutoff = -18.+dmfid
  cutoff = cutoff*(xfid le 50.)

  polyfill, [alog10(xfid), reverse(alog10(xfid)), (alog10(xfid))[0]] $
            , [cutoff, ymax+xfid*0.0, cutoff[0]] $
            , color=cgcolor('lightgray'), noclip=0
  
  noise_x = randomn(seed, n_elements(tab))*0.01
  noise_y = randomn(seed, n_elements(tab))*0.05
  ;oplot, alog10(tab.dist_mpc)+noise_x, tab.btc+noise_y, ps=1 $
  ;       , symsize=0.25
  colorvec = round(prob*255)
  loadct, 33
  plots, alog10(tab.dist_mpc)+noise_x, tab.btc+noise_y, ps=cgsymcat('filledcircle') $
         , color=colorvec, symsize=0.50, noclip=0

  loadct, 0
  for ii = -10, 25 do $
     oplot, ii*0.2*[1,1], [-1d6, 1d6], lines=1, thick=3, color=cgcolor('charcoal')
  for ii = -40, 40 do $
     oplot, [0, 1e3], ii*1.0*[1,1], lines=1, thick=3, color=cgcolor('charcoal')

  oplot, alog10(xfid), cutoff, thick=10, color=cgcolor('gray')

  !p.charthick=3
  loadct, 33
  cgColorbar, Range=[0., 1.] $
              , Position=[0.15, 0.95, 0.95, 0.975] $
              , xtickformat='(F5.2)' $
              , title="!6Probability of Selection (!8p!6)" $
              , annotatecolor=cgcolor("black", 255) $
              , charthick = 3, charsize=1.0 $
              , textthick=3
  !p.charthick=1
  loadct, 0

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CUMULATIVE PROBABILITY DISTRIBUTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  plot, findgen(10), xtitle='!6'

  psfile = '../plots/select_cdf.eps'
  pnfile = '../plots/select_cdf.png'
  ps, /ps, /def, file=psfile, xs=5, ys=5, /color, /encaps
  
  !p.multi=[0,1,1]
 
  loadct, 0

  xmin = 0.0
  xmax = 1.0
  ymin = 1.0
  ymax = 5.0

  plot $
     , [0], [0] $
     , ps=1, symsize=0.5 $
     , xtitle='!6Probability of Selection !8p!6' $
     , ytitle='!6log!d10!n Number of Galaxies With Probability > !8p!6' $
     , xrange=[xmin, xmax], /xstyle $
     , yrange=[ymin, ymax], /ystyle $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.0 $
     , position=[0.15,0.1,0.95,0.95]

  pvec = prob[where(finite(prob))]
  pvec = pvec[sort(pvec)]
  cumdist = (n_elements(pvec) - total(pvec*0.0+1L, /cumul))
   
  loadct, 0
  for ii = -10, 25 do $
     oplot, ii*0.2*[1,1], [-1d6, 1d6], lines=1, thick=3, color=cgcolor('charcoal')
  for ii = -40, 40 do $
     oplot, [0, 1e3], ii*1.0*[1,1], lines=1, thick=3, color=cgcolor('charcoal')

  oplot, pvec, alog10(cumdist), lines=0, thick=7

  oplot, 0.1*[1,1], [-10., 10.], thick=10 $
         , color=cgcolor('navy')
  xyouts,  0.1*[1]+0.05, [2.0], orient=90, align=0.5, '!6!8p!6 = 0.1' $
           , color=cgcolor('navy'), charsize=1.0, charthick=3

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IDENTIFY SELECTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  select_ind = where(prob gt selection_thresh, select_ct)
  stab = tab[select_ind]
  selected_absb = absb[select_ind]

; COMPLETENESS ONLY FOR M_B < -18

  dvec = stab.dist_mpc
  tvec = stab.t
  bvec = selected_absb
  
  bind = where(bvec lt -18., comp_ct)
  dvec = dvec[bind]
  tvec = tvec[bind]
  bvec = bvec[bind]
  
  dind = sort(dvec)
  dvec = dvec[dind]
  tvec = tvec[dind]
  bvec = bvec[dind]

; COMPLETENESS IN B

  bthresh_lo = -18
  bthresh = -19
  bthresh_mid = -20
  frac_lowb = fltarr(select_ct)*!values.f_nan
  frac_midb = fltarr(select_ct)*!values.f_nan
  frac_highb = fltarr(select_ct)*!values.f_nan
  
  window = 500
  for ii = window/2, comp_ct-1-window/2-1 do begin
     lo = ii-window/2
     hi = ii+window/2-1
     frac_lowb[ii] = total(bvec[lo:hi] gt bthresh and bvec[lo:hi] lt bthresh_lo,/nan)/total(bvec[lo:hi]*0.+1.,/nan)
     frac_midb[ii] = total(bvec[lo:hi] gt bthresh_mid and bvec[lo:hi] lt bthresh,/nan)/total(bvec[lo:hi]*0.+1.,/nan)
     frac_highb[ii] = total(bvec[lo:hi] lt bthresh_mid,/nan)/total(bvec[lo:hi]*0.+1.,/nan)
  endfor

; COMPLETENESS IN T

  tthresh_early = 2
  tthresh_late = 6  
  frac_latet = fltarr(select_ct)*!values.f_nan
  frac_midt = fltarr(select_ct)*!values.f_nan
  frac_earlyt = fltarr(select_ct)*!values.f_nan

  window = 500
  for ii = window/2, comp_ct-1-window/2-1 do begin
     lo = ii-window/2
     hi = ii+window/2-1
     frac_latet[ii] = $
        total(tvec[lo:hi] gt tthresh_late,/nan)/total(tvec[lo:hi]*0.+1.,/nan) & $
        frac_midt[ii] = total(tvec[lo:hi] le tthresh_late $
                              and tvec[lo:hi] ge tthresh_early,/nan)/total(tvec[lo:hi]*0.+1.,/nan)
     frac_earlyt[ii] = total(tvec[lo:hi] lt 2,/nan)/total(tvec[lo:hi]*0.+1.,/nan)
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT FRACTION BY MORPHOLOGY AND B MAGNITUDE VS DISTANCE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  plot, findgen(10), xtitle='!6'

  psfile = '../plots/morph_vs_dist.eps'
  pnfile = '../plots/morph_vs_dist.png'
  ps, /ps, /def, file=psfile, xs=5, ys=5, /color, /encaps
  
  !p.multi=[0,1,1]
 
  loadct, 0

  xmin = 10.0
  xmax = 50.
  ymin = 0.0
  ymax = 0.8

  plot $
     , [0], [0] $
     , ps=1, symsize=0.5 $
     , xtitle='!6Estimated Distance [Mpc]' $
     , ytitle='!6Fraction of Galaxies With Given Morphology' $
     , xrange=[xmin, xmax], /xstyle $
     , yrange=[ymin, ymax], /ystyle $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.0 $
     , position=[0.15,0.1,0.975,0.95]
   
  loadct, 0
  for ii = -10, 25 do $
     oplot, ii*1e1*[1,1], [-1d6, 1d6], lines=1, thick=3, color=cgcolor('charcoal')
  for ii = -40, 40 do $
     oplot, [0, 1e6], ii*0.1*[1,1], lines=1, thick=3, color=cgcolor('charcoal')

  oplot, dvec, frac_earlyt, thick=5, color=cgcolor('firebrick')
  oplot, dvec, frac_midt, thick=5, color=cgcolor('darkgreen')
  oplot, dvec, frac_latet, thick=5, color=cgcolor('royalblue')

  al_legend, /top, /right, clear=1, outline=cgcolor('black') $
             , background=cgcolor('lightgray') $
             , lines=[-99,0,0,0] $
             , color=[0, cgcolor('firebrick'), cgcolor('darkgreen'), cgcolor('royalblue')] $
             , textcolor=[cgcolor('black'), cgcolor('firebrick'), cgcolor('darkgreen'), cgcolor('royalblue')] $
             , ['!6Fraction per 500-galaxy window ...' $
               , '!6... with !8t!6 < 2' $
               , '!6... with !8t!6 = 2 to 6' $
               , '!6... with !8t!6 > 6' $
               ] $
             , pspacing=1.0, thick=5 $
             , charthick=2, charsize=1.0

  al_legend, ['!8M!dB!n!6 < -18 mag only'], lines=-99, charthick=2, charsize=1.0 $
             , /bottom, /left, box=0, clear=0

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

  psfile = '../plots/absb_vs_dist.eps'
  pnfile = '../plots/absb_vs_dist.png'
  ps, /ps, /def, file=psfile, xs=5, ys=5, /color, /encaps
  
  !p.multi=[0,1,1]
 
  loadct, 0

  xmin = 10.0
  xmax = 50.
  ymin = 0.0
  ymax = 0.7

  plot $
     , [0], [0] $
     , ps=1, symsize=0.5 $
     , xtitle='!6Estimated Distance [Mpc]' $
     , ytitle='!6Fraction of Galaxies With Given Luminosity' $
     , xrange=[xmin, xmax], /xstyle $
     , yrange=[ymin, ymax], /ystyle $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.0 $
     , position=[0.15,0.1,0.975,0.95]
   
  loadct, 0
  for ii = -10, 25 do $
     oplot, ii*1e1*[1,1], [-1d6, 1d6], lines=1, thick=3, color=cgcolor('charcoal')
  for ii = -40, 40 do $
     oplot, [0, 1e6], ii*0.1*[1,1], lines=1, thick=3, color=cgcolor('charcoal')

  oplot, dvec, frac_highb, thick=5, color=cgcolor('firebrick')
  oplot, dvec, frac_midb, thick=5, color=cgcolor('darkgreen')
  oplot, dvec, frac_lowb, thick=5, color=cgcolor('royalblue')

  al_legend, /top, /left, clear=1, outline=cgcolor('black') $
             , background=cgcolor('lightgray') $
             , lines=[-99,0,0,0] $
             , color=[0, cgcolor('firebrick'), cgcolor('darkgreen'), cgcolor('royalblue')] $
             , textcolor=[cgcolor('black'), cgcolor('firebrick'), cgcolor('darkgreen'), cgcolor('royalblue')] $
             , ['!6Fraction per 500-galaxy window ...' $
               , '!6... with !8M!dB!n!6 < -20' $
               , '!6... with !8M!dB!n!6 = -19 to -20' $
               , '!6... with !8M!dB!n!6 = -18 to -19' $
               ] $
             , pspacing=1.0, thick=5 $
             , charthick=2, charsize=1.0
  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

  stop

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT SIZE VS DISTANCE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  psfile = '../plots/size_vs_dist.eps'
  pnfile = '../plots/size_vs_dist.png'
  ps, /ps, /def, file=psfile, xs=5, ys=5, /color, /encaps
    
  !p.multi=[0,1,1]
 
  loadct, 0

  xmin = 0.75
  xmax = 2.0
  ymin = 1.0
  ymax = 3.0

  plot $
     , [0], [0] $
     , ps=1, symsize=0.5 $
     , xtitle='!6log!d10!n Estimated Distance [Mpc]' $
     , ytitle='!6log!d10!n Optical diameter !8d!d25!n!6 [arcseconds]' $
     , xrange=[xmin, xmax], /xstyle $
     , yrange=[ymin, ymax], /ystyle $
     , xthick=5, ythick=5 $
     , charthick=3, charsize=1.0 $
     , position=[0.15,0.1,0.975,0.95]

  pvec = prob[where(finite(prob))]
  pvec = pvec[sort(pvec)]
  cumdist = (n_elements(pvec) - total(pvec*0.0+1L, /cumul))

  for ii = -10, 25 do $
     oplot, ii*0.2*[1,1], [-1d6, 1d6], lines=1, thick=3, color=cgcolor('charcoal')
  for ii = -30, 30 do $
     oplot, [-1d6, 1d6], ii*0.2*[1,1], lines=1
  
  ;xnoise = randomn(seed,select_ct)*1.0
  ;ynoise = 0.0
  ;oplot, (tab[select_ind].dist_mpc)+xnoise $
  ;       , (2.0*tab[select_ind].r25_deg*3600.)+ynoise $
  ;       , ps=cgsymcat('filledcircle') $
  ;       , color=cgcolor('gray'), symsize=0.50      

  ind = where(absb le -18.)
  x = alog10(tab[ind].dist_mpc)
  y = alog10(2.0*tab[ind].r25_deg*3600.)
  oplot, x, y $
         , ps=1 $
         , color=cgcolor('gray')
  
  bins = bin_data(x, y, binsize=0.075, xmin=1.0, xmax=2.0, /nan)
  oploterror, bins.xmid, bins.ymed, bins.ymad $
              , psym=cgsymcat('filledcircle'), /nohat $
              , errthick=5, symsize=1.5, color=cgcolor('firebrick')

  oplot, [0, 1e6], alog10(15.)*[1,1], lines=2, thick=5, color=cgcolor('black')
  xyouts, [1.0], alog10(15.)*[1]+0.05, '!61 beam' $
          , align=0.5, charsize=1.0, charthick=3, color=cgcolor('black')

  oplot, [0, 1e6], alog10(5.*15.)*[1,1], lines=2, thick=5, color=cgcolor('black')
  xyouts, [1.0], alog10(5.*15.)*[1]+0.05, '!65 beams' $
          , align=0.5, charsize=1.0, charthick=3, color=cgcolor('black')

  oplot, [0, 1e6], alog10(10.*15.)*[1,1], lines=2, thick=5, color=cgcolor('black')
  xyouts, [1.0], alog10(10.*15.)*[1]+0.05, '!610 beams' $
          , align=0.5, charsize=1.0, charthick=3, color=cgcolor('black')

  oplot, alog10(50.)*[1,1], [-10., 10.], thick=10 $
         , color=cgcolor('navy')
  xyouts,  alog10(50.)*[1]-0.025, [2.75], orient=90, align=0.5, '!6!8d!6 = 50 Mpc' $
           , color=cgcolor('navy'), charsize=1.0, charthick=3

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PRINT OUR TARGET LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_write) then begin
     openw, 1, '../measurements/selection.txt'
     for ii = 0L, select_ct-1 do $
        printf, 1, (stab[ii].pgc)  
     close, 1
  endif

  stop

end
