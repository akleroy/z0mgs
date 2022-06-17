pro plot_nga_comp $
   , res_str=res_str

  if n_elements(res_str) eq 0 then $
     res_str = 'gauss15'

  if res_str eq 'gauss7p5' then begin
     label = '7.5" Beam'
  endif

  if res_str eq 'gauss15' then begin
     label = '15" Beam'
  endif

  restore, '../measurements/nga_z0mgs_comp_'+res_str+'.idl', /v
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
  
  psfile = '../plots/nga_comp_fuv_'+res_str+'.eps'
  pnfile = '../plots/nga_comp_fuv_'+res_str+'.png'
  ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
      , file=psfile

  loadct, 0
  reversect      
  disp, alog10(fuv_grid), xaxis, yaxis $
        , xtitle='!6log!d10!n This Atlas Intensity [MJy sr!u-1!n]' $
        , ytitle='!6log!d10!n NGA / This Atlas' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xstyle=1, ystyle=1, reserve=50, color=cgcolor('black',255) $
        , position=[0.2, 0.2, 0.95, 0.95]
  
  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')
  
  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.1*[1,1], lines=1, color=cgcolor('charcoal')

  oplot, [-1d6, 1d6], [0,0], thick=10, color=cgcolor('cornflowerblue')

  bins = bin_data(x, y $
                  , xmin=-3.5, xmax=0., binsize=0.25, /nan)
  oploterror, (bins.xmid), (bins.ymed), bins.ymad $
              , color=cgcolor('salmon'), psym=cgsymcat('filledcircle') $
              , /nohat, errthick=7, symsize=1.5
  oplot, (bins.xmid), (bins.ymed) $
         , color=cgcolor('firebrick'), psym=cgsymcat('filledcircle')
         
  al_legend, /bottom, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Offset: '+string(median(abs(bins.ymed)), format='(F5.3)'), $
              'Scatter: '+string(median(bins.ymad), format='(F5.3)')] $
             , lines=-99, charthick=3, charsize=1.25

  al_legend, /top, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , lines=-99, ['!6FUV',label], charthick=3, charsize=1.25

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
  
  psfile = '../plots/nga_comp_nuv_'+res_str+'.eps'
  pnfile = '../plots/nga_comp_nuv_'+res_str+'.png'
  ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
      , file=psfile

  loadct, 0
  reversect      
  disp, alog10(nuv_grid), xaxis, yaxis $
        , xtitle='!6log!d10!n This Atlas Intensity [MJy sr!u-1!n]' $
        , ytitle='!6log!d10!n NGA / This Atlas' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xstyle=1, ystyle=1, reserve=50, color=cgcolor('black',255) $
        , position=[0.2, 0.2, 0.95, 0.95]
  
  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')
  
  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.1*[1,1], lines=1, color=cgcolor('charcoal')

  oplot, [-1d6, 1d6], [0,0], thick=10, color=cgcolor('cornflowerblue')

  bins = bin_data(x, y $
                  , xmin=-3., xmax=0.5, binsize=0.25, /nan)
  oploterror, (bins.xmid), (bins.ymed), bins.ymad $
              , color=cgcolor('salmon'), psym=cgsymcat('filledcircle') $
              , /nohat, errthick=5,symsize=1.5
  oplot, (bins.xmid), (bins.ymed) $
         , color=cgcolor('firebrick'), psym=cgsymcat('filledcircle')

  al_legend, /top, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , lines=-99, ['!6NUV',label], charthick=3, charsize=1.25

  al_legend, /bottom, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Offset: '+string(median(abs(bins.ymed)), format='(F5.3)'), $
              'Scatter: '+string(median(bins.ymad), format='(F5.3)')] $
             , lines=-99, charthick=3, charsize=1.25

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

  stop

end
