pro plot_s4g_comp

  restore, '../measurements/s4g_z0mgs_comp.idl', /v
  wise1 = where(band_comp eq 'WISE1')
  wise2 = where(band_comp eq 'WISE2')

; WISE1
  
  high_sn = where(s4g_comp[wise1] gt 10.*noise_comp[wise1] and $
                  z0mgs_comp[wise1] gt 10.*noise_comp[wise1])
  x = alog10(z0mgs_comp[wise1[high_sn]])
  y = alog10(s4g_comp[wise1[high_sn]]/z0mgs_comp[wise1[high_sn]])
  
  wise1_grid = grid_data(x, y $
                       , xaxis_out = xaxis, yaxis_out = yaxis $
                       , ymin=-0.5, ymax=0.5, binsize_x=0.025 $
                       , xmin=-1.5, xmax=2.5, binsize_y=0.025)
  
  psfile = '../plots/s4g_comp_wise1.eps'
  pnfile = '../plots/s4g_comp_wise1.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile

  loadct, 0
  reversect      
  disp, alog10(wise1_grid), xaxis, yaxis $
        , xtitle='!6log!d10!n This Atlas Intensity [MJy sr!u-1!n]' $
        , ytitle='!6log!d10!n S4G / This Atlas' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xstyle=1, ystyle=1, reserve=5, color=cgcolor('black',255) $
        , position=[0.2, 0.2, 0.95, 0.95]
  
  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')
  
  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.1*[1,1], lines=1, color=cgcolor('charcoal')

  oplot, [-1d6, 1d6], [0,0], thick=10, color=cgcolor('salmon')

  bins = bin_data(x, y $
                  , xmin=-1.0, xmax=1.5, binsize=0.125, /nan)
  oploterror, (bins.xmid), (bins.ymed), bins.ymad $
              , color=cgcolor('dodgerblue'), psym=cgsymcat('filledcircle') $
              , /nohat, errthick=5
  
  al_legend, /top, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Typical offset: '+string(median(abs(bins.ymed)), format='(F5.3)'), $
              'Typical scatter: '+string(median(bins.ymad), format='(F5.3)')] $
             , lines=-99, charthick=3, charsize=1.25

  al_legend, /bottom, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , lines=-99, ['!6WISE1/IRAC1'], charthick=3, charsize=1.5

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

; WISE2
  
  high_sn = where(s4g_comp[wise2] gt 10.*noise_comp[wise2] and $
                  z0mgs_comp[wise2] gt 10.*noise_comp[wise2])
  x = alog10(z0mgs_comp[wise2[high_sn]])
  y = alog10(s4g_comp[wise2[high_sn]]/z0mgs_comp[wise2[high_sn]])
  
  wise2_grid = grid_data(x, y $
                       , xaxis_out = xaxis, yaxis_out = yaxis $
                       , ymin=-0.5, ymax=0.5, binsize_x=0.025 $
                       , xmin=-1.5, xmax=2.5, binsize_y=0.025)
  
  psfile = '../plots/s4g_comp_wise2.eps'
  pnfile = '../plots/s4g_comp_wise2.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile

  loadct, 0
  reversect      
  disp, alog10(wise2_grid), xaxis, yaxis $
        , xtitle='!6log!d10!n This Atlas Intensity [MJy sr!u-1!n]' $
        , ytitle='!6log!d10!n S4G / This Atlas' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xstyle=1, ystyle=1, reserve=5, color=cgcolor('black',255) $
        , position=[0.2, 0.2, 0.95, 0.95]
  
  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')
  
  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.1*[1,1], lines=1, color=cgcolor('charcoal')

  oplot, [-1d6, 1d6], [0,0], thick=10, color=cgcolor('salmon')

  bins = bin_data(x, y $
                  , xmin=-1.0, xmax=1.25, binsize=0.125, /nan)
  oploterror, (bins.xmid), (bins.ymed), bins.ymad $
              , color=cgcolor('dodgerblue'), psym=cgsymcat('filledcircle') $
              , /nohat, errthick=5
  
  al_legend, /top, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Typical offset: '+string(median(abs(bins.ymed)), format='(F5.3)'), $
              'Typical scatter: '+string(median(bins.ymad), format='(F5.3)')] $
             , lines=-99, charthick=3, charsize=1.25

  al_legend, /bottom, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , lines=-99, ['!6WISE2/IRAC2'], charthick=3, charsize=1.5

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

; OFFSET

  xmin = -1.0
  xmax = 1.0
  binsize = 0.025
  uniq_ind = uniq(gal_comp[wise1])
  vec = (offset2_comp)[wise1[uniq_ind]]*1d2
  bins_wise1 = $
     bin_data((vec), vec*0.0+1.0 $
              , xmin=xmin, xmax=xmax, binsize=binsize, /nan)
  
  uniq_ind = uniq(gal_comp[wise2])
  vec = (offset2_comp)[wise2[uniq_ind]]*1d2
  bins_wise2 = $
     bin_data((vec), vec*0.0+1.0 $
              , xmin=xmin, xmax=xmax, binsize=binsize, /nan)
  
  psfile = '../plots/s4g_hist.eps'
  pnfile = '../plots/s4g_hist.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6This Atlas - S4G [10!u-2!n MJy sr!u-1!n]' $
     , ytitle='!6Number of Galaxies' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[xmin, xmax], yrange=[0., 2d2], ystyle=1
      
  for ii = -100, 100 do $
     oplot, ii*binsize*5.*[1,1], [-1d10, 1d10], lines=1, color=cgcolor('charcoal')

  for ii = -100, 100 do $
     oplot, [-10, 10], ii*50*[1,1], lines=1, color=cgcolor('charcoal')
 

  histplot $
     , bins_wise2.xmid, ((bins_wise2.counts) > (0.)) $
     , /overplot $
     , lthick=3, /nobar, /fill $
     , fcolor=cgcolor('salmon') $
     , lcolor=cgcolor('firebrick')

  histplot $
     , bins_wise1.xmid, ((bins_wise1.counts) > (0.)) $
     , /overplot $
     , lthick=3, /nobar, /fline, forient=45 $
     , fcolor=cgcolor('royalblue') $
     , lcolor=cgcolor('royalblue')
  
  histplot $
     , bins_wise1.xmid, ((bins_wise1.counts) > (0.)) $
     , /overplot $
     , lthick=3, /nobar, /fline, forient=-45 $
     , fcolor=cgcolor('royalblue') $
     , lcolor=cgcolor('royalblue')
  
  oplot, [0., 0.], [-1d10, 1d10], thick=5

  al_legend, /top, /left, box=1, clear=1, background=cgcolor('lightgray') $
             , lines=-99, ['!6Zero','Point'], charthick=3, charsize=1.5

  al_legend $
     , /top, /right $
     , box=1, clear=1 $
     , background=cgcolor('lightgray') $
     , charsize=1.5, charthick=3 $
     , lines=-99 $
     , ['WISE2/IRAC2','WISE1/IRAC1'] $
     , textcolor=[cgcolor('salmon'), cgcolor('royalblue')]

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile

  stop

end
