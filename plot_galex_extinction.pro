pro plot_galex_extinction

  tab = mrdfits('../measurements/delivery_index.fits',1,h)
  b = tab.gb_deg

  plot, findgen(10), title='!6Test'

  bmin = -87.5
  bmax = +87.5
  bbinsize = 5.0
  afuv_bins = $
     bin_data(b, tab.afuv, /nan $
              , xmin=bmin, xmax=bmax, binsize=bbinsize)
  anuv_bins = $
     bin_data(b, tab.anuv, /nan $
              , xmin=bmin, xmax=bmax, binsize=bbinsize)
     
  psfile = '../plots/galex_extinction.eps'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  ymin = 0.0
  ymax = 2.0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6Galactic Latitude' $
     , ytitle='!6A!dFUV!n and A!dNUV!n [mag]' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[-100., 100.], yrange=[ymin, ymax]
  
  for ii = -100, 100 do $
     oplot, ii*15.*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')
  
  for ii = -100, 100 do $
     oplot, [-100, 100], ii*0.25*[1,1], lines=1, color=cgcolor('charcoal')
  
  circle, /fill
  
  oploterror, anuv_bins.xmid-1.0, (anuv_bins.ymed), anuv_bins.ymad_log $
              , ps=cgsymcat('filledsquare') $
              , symsize=2.0, color=cgcolor('black'), errthick=3
  oplot, anuv_bins.xmid-1.0, (anuv_bins.ymed), ps=cgsymcat('filledsquare') $
         , symsize=1.0 $
         , color=cgcolor('salmon')
  
  oploterror, afuv_bins.xmid+1.0, (afuv_bins.ymed), afuv_bins.ymad_log $
              , ps=cgsymcat('filledcircle'), symsize=2.0, color=cgcolor('black'), errthick=3
  oplot, afuv_bins.xmid+1.0, (afuv_bins.ymed), ps=cgsymcat('filledcircle') $
         , symsize=1.0, color=cgcolor('dodgerblue')
  
  al_legend $
     , /top, /left $
     , box=0, clear=0, charsize=1.75, charthick=3 $
     , lines=-99 $
     , ['!6NUV', 'FUV'] $
     , textcolor=[cgcolor('salmon'), cgcolor('dodgerblue')]
  
  ps, /xw
  spawn, 'evince '+psfile+' &'

  stop

end
