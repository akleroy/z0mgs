pro plot_s4g_comp $
   , res_str=res_str

  if n_elements(res_str) eq 0 then $
     res_str = 'gauss15'

  if res_str eq 'gauss7p5' then begin
     label = '7.5" Beam'
  endif

  if res_str eq 'gauss15' then begin
     label = '15" Beam'
  endif

  restore, '../measurements/s4g_z0mgs_comp_'+res_str+'.idl', /v
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
  
  psfile = '../plots/s4g_comp_wise1_'+res_str+'.eps'
  pnfile = '../plots/s4g_comp_wise1_'+res_str+'.png'
  ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
      , file=psfile

  loadct, 0
  reversect      
  disp, alog10(wise1_grid), xaxis, yaxis $
        , xtitle='!6log!d10!n This Atlas Intensity [MJy sr!u-1!n]' $
        , ytitle='!6log!d10!n S4G / This Atlas' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xstyle=1, ystyle=1, reserve=50, color=cgcolor('black',255) $
        , position=[0.2, 0.2, 0.95, 0.95]
  
  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')
  
  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.1*[1,1], lines=1, color=cgcolor('charcoal')

  oplot, [-1d6, 1d6], [0,0], thick=10, color=cgcolor('cornflowerblue')

  bins = bin_data(x, y $
                  , xmin=-1.0, xmax=1.5, binsize=0.125, /nan)
  oploterror, (bins.xmid), (bins.ymed), bins.ymad $
              , color=cgcolor('salmon'), psym=cgsymcat('filledcircle') $
              , /nohat, errthick=7,symsize=1.5
  oplot, (bins.xmid), (bins.ymed) $
         , color=cgcolor('firebrick'), psym=cgsymcat('filledcircle')
         
  al_legend, /top, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Offset: '+string(median(abs(bins.ymed)), format='(F5.3)'), $
              'Scatter: '+string(median(bins.ymad), format='(F5.3)')] $
             , lines=-99, charthick=3, charsize=1.25

  al_legend, /bottom, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , lines=-99, ['!6WISE1/IRAC1',label], charthick=3, charsize=1.25

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
  
  psfile = '../plots/s4g_comp_wise2_'+res_str+'.eps'
  pnfile = '../plots/s4g_comp_wise2_'+res_str+'.png'
  ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
      , file=psfile

  loadct, 0
  reversect      
  disp, alog10(wise2_grid), xaxis, yaxis $
        , xtitle='!6log!d10!n This Atlas Intensity [MJy sr!u-1!n]' $
        , ytitle='!6log!d10!n S4G / This Atlas' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xstyle=1, ystyle=1, reserve=50, color=cgcolor('black',255) $
        , position=[0.2, 0.2, 0.95, 0.95]
  
  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')
  
  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.1*[1,1], lines=1, color=cgcolor('charcoal')


  oplot, [-1d6, 1d6], [0,0], thick=10, color=cgcolor('cornflowerblue')

  bins = bin_data(x, y $
                  , xmin=-1.0, xmax=1.25, binsize=0.125, /nan)

  oploterror, (bins.xmid), (bins.ymed), bins.ymad $
              , color=cgcolor('salmon'), psym=cgsymcat('filledcircle') $
              , /nohat, errthick=7,symsize=1.5
  oplot, (bins.xmid), (bins.ymed) $
         , color=cgcolor('firebrick'), psym=cgsymcat('filledcircle')

  al_legend, /top, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Offset: '+string(median(abs(bins.ymed)), format='(F5.3)'), $
              'Scatter: '+string(median(bins.ymad), format='(F5.3)')] $
             , lines=-99, charthick=3, charsize=1.25

  al_legend, /bottom, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , lines=-99, ['!6WISE2/IRAC2', label], charthick=3, charsize=1.25

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

  stop

end
