pro plot_galex_extinction

  tab = mrdfits('../measurements/delivery_index_gauss15.fits',1,h)
  b = tab.gb_deg

  plot, findgen(10), title='!6Test'

  bmin = -2.5
  bmax = +87.5
  bbinsize = 5.0
  afuv_bins = $
     bin_data(abs(b), tab.afuv, /nan $
              , xmin=bmin, xmax=bmax, binsize=bbinsize)

  bmin = 0.0
  bmax = +85.0
  bbinsize = 5.0
  anuv_bins = $
     bin_data(abs(b), tab.anuv, /nan $
              , xmin=bmin, xmax=bmax, binsize=bbinsize)
     
  psfile = '../plots/galex_extinction.eps'
  ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
      , file=psfile
  
  ymin = 0.0
  ymax = 2.0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6Galactic Latitude' $
     , ytitle='!6A!dFUV!n and A!dNUV!n [mag]' $
     , xthick=5, ythick=5, charthick=3, charsize=1.25 $
     , xrange=[0., 100.], yrange=[ymin, ymax]
  
  for ii = -100, 100 do $
     oplot, ii*15.*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')
  
  for ii = -100, 100 do $
     oplot, [-100, 100], ii*0.25*[1,1], lines=1, color=cgcolor('charcoal')
  
  circle, /fill
  
  oploterror, anuv_bins.xmid, (anuv_bins.ymed), anuv_bins.ymad_log $
              , ps=cgsymcat('filledsquare') $
              , symsize=1.0, color=cgcolor('firebrick'), errthick=3 $
              , /nohat
;  oplot, anuv_bins.xmid-1.0, (anuv_bins.ymed), ps=cgsymcat('filledsquare') $
;         , symsize=1.0 $
;         , color=cgcolor('salmon')
  
  oploterror, afuv_bins.xmid, (afuv_bins.ymed), afuv_bins.ymad_log $
              , ps=cgsymcat('filledcircle'), symsize=1.0 $
              , color=cgcolor('navy'), errthick=3, /nohat
;  oplot, afuv_bins.xmid+1.0, (afuv_bins.ymed), ps=cgsymcat('filledcircle') $
;         , symsize=1.0, color=cgcolor('navy')
  
  al_legend $
     , /top, /right $
     , box=1, clear=1, background=cgcolor('lightgray') $
     , charsize=1.25, charthick=3 $
     , lines=-99 $
     , ['!6NUV', 'FUV'] $
     , textcolor=[cgcolor('firebrick'), cgcolor('navy')]
  
  ps, /xw
  spawn, 'evince '+psfile+' &'

  stop

end
