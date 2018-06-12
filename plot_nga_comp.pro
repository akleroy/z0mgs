pro plot_nga_comp

  restore, '../measurements/nga_z0mgs_comp.idl', /v
  fuv = where(band_comp eq 'FUV' and gal_comp ne 'PGC2557' and gal_comp ne 'PGC5818')
  nuv = where(band_comp eq 'NUV' and gal_comp ne 'PGC2557' and gal_comp ne 'PGC5818')

; FUV
  
  high_sn = where(nga_comp[fuv] gt 10.*noise_comp[fuv] and $
                  z0mgs_comp[fuv] gt 10.*noise_comp[fuv])
  x = alog10(z0mgs_comp[fuv[high_sn]])
  y = alog10(nga_comp[fuv[high_sn]]/z0mgs_comp[fuv[high_sn]])
  
  fuv_grid = grid_data(x, y $
                       , xaxis_out = xaxis, yaxis_out = yaxis $
                       , ymin=-0.5, ymax=0.5, binsize_x=0.025 $
                       , xmin=-3.5, xmax=0.5, binsize_y=0.025)
  
  psfile = '../plots/nga_comp_fuv.eps'
  pnfile = '../plots/nga_comp_fuv.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile

  loadct, 0
  reversect      
  disp, alog10(fuv_grid), xaxis, yaxis $
        , xtitle='!6log!d10!n This Atlas Intensity [MJy sr!u-1!n]' $
        , ytitle='!6log!d10!n Nearby Galaxy Atlas / This Atlas' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xstyle=1, ystyle=1, reserve=5, color=cgcolor('black',255) $
        , position=[0.2, 0.2, 0.95, 0.95]
  
  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')
  
  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.1*[1,1], lines=1, color=cgcolor('charcoal')

  oplot, [-1d6, 1d6], [0,0], thick=10, color=cgcolor('salmon')

  bins = bin_data(x, y $
                  , xmin=-3.5, xmax=0., binsize=0.25, /nan)
  oploterror, (bins.xmid), (bins.ymed), bins.ymad $
              , color=cgcolor('dodgerblue'), psym=cgsymcat('filledcircle') $
              , /nohat, errthick=5
  
  al_legend, /bottom, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Typical offset: '+string(median(abs(bins.ymed)), format='(F5.3)'), $
              'Typical scatter: '+string(median(bins.ymad), format='(F5.3)')] $
             , lines=-99, charthick=3, charsize=1.25

  al_legend, /top, /left, box=1, clear=1, background=cgcolor('lightgray') $
             , lines=-99, ['!6FUV'], charthick=3, charsize=1.5

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

; NUV
  
  high_sn = where(nga_comp[nuv] gt 10.*noise_comp[nuv] and $
                  z0mgs_comp[nuv] gt 10.*noise_comp[nuv])
  x = alog10(z0mgs_comp[nuv[high_sn]])
  y = alog10(nga_comp[nuv[high_sn]]/z0mgs_comp[nuv[high_sn]])
  
  nuv_grid = grid_data(x, y $
                       , xaxis_out = xaxis, yaxis_out = yaxis $
                       , ymin=-0.5, ymax=0.5, binsize_x=0.025 $
                       , xmin=-3.5, xmax=0.5, binsize_y=0.025)
  
  psfile = '../plots/nga_comp_nuv.eps'
  pnfile = '../plots/nga_comp_nuv.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile

  loadct, 0
  reversect      
  disp, alog10(nuv_grid), xaxis, yaxis $
        , xtitle='!6log!d10!n This Atlas Intensity [MJy sr!u-1!n]' $
        , ytitle='!6log!d10!n Nearby Galaxy Atlas / This Atlas' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xstyle=1, ystyle=1, reserve=5, color=cgcolor('black',255) $
        , position=[0.2, 0.2, 0.95, 0.95]
  
  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')
  
  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.1*[1,1], lines=1, color=cgcolor('charcoal')

  oplot, [-1d6, 1d6], [0,0], thick=10, color=cgcolor('salmon')

  bins = bin_data(x, y $
                  , xmin=-3., xmax=0.5, binsize=0.25, /nan)
  oploterror, (bins.xmid), (bins.ymed), bins.ymad $
              , color=cgcolor('dodgerblue'), psym=cgsymcat('filledcircle') $
              , /nohat, errthick=5

  al_legend, /top, /left, box=1, clear=1, background=cgcolor('lightgray') $
             , lines=-99, ['!6NUV'], charthick=3, charsize=1.5

  al_legend, /bottom, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Typical offset: '+string(median(abs(bins.ymed)), format='(F5.3)'), $
              'Typical scatter: '+string(median(bins.ymad), format='(F5.3)')] $
             , lines=-99, charthick=3, charsize=1.25

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

; OFFSET
  xmin = -1.0
  xmax = 1.0
  binsize = 0.025
  uniq_ind = uniq(gal_comp[fuv])
  vec = (offset2_comp)[fuv[uniq_ind]]*1d3
  print, "FUV background scatter: ", mad(vec)
  bins_fuv = $
     bin_data((vec), vec*0.0+1.0 $
              , xmin=xmin, xmax=xmax, binsize=binsize, /nan)
  
  uniq_ind = uniq(gal_comp[nuv])
  vec = (offset2_comp)[nuv[uniq_ind]]*1d3
  print, "NUV background scatter: ", mad(vec)
  bins_nuv = $
     bin_data((vec), vec*0.0+1.0 $
              , xmin=xmin, xmax=xmax, binsize=binsize, /nan)
  
  psfile = '../plots/nga_hist.eps'
  pnfile = '../plots/nga_hist.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6This Atlas - NGA [10!u-3!n MJy sr!u-1!n]' $
     , ytitle='!6Number of Galaxies' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[xmin, xmax], yrange=[0., 1d2], ystyle=1
      
  for ii = -100, 100 do $
     oplot, ii*binsize*5.*[1,1], [-1d10, 1d10], lines=1, color=cgcolor('charcoal')

  for ii = -100, 100 do $
     oplot, [-10, 10], ii*10*[1,1], lines=1, color=cgcolor('charcoal')
 

  histplot $
     , bins_nuv.xmid, ((bins_nuv.counts) > (0.)) $
     , /overplot $
     , lthick=3, /nobar, /fill $
     , fcolor=cgcolor('salmon') $
     , lcolor=cgcolor('firebrick')

  histplot $
     , bins_fuv.xmid, ((bins_fuv.counts) > (0.)) $
     , /overplot $
     , lthick=3, /nobar, /fline, forient=45 $
     , fcolor=cgcolor('royalblue') $
     , lcolor=cgcolor('royalblue')
  
  histplot $
     , bins_fuv.xmid, ((bins_fuv.counts) > (0.)) $
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
     , charsize=1.75, charthick=3 $
     , lines=-99 $
     , ['NUV','FUV'] $
     , textcolor=[cgcolor('salmon'), cgcolor('royalblue')]

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile

  stop

end
