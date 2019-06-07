pro plot_detection_vs_absb

  tab = mrdfits('../measurements/delivery_index_gauss15.fits',1,h)
  dat = gal_data(pgc=tab.pgc)
  absb = dat.btc-dat.distmod

  ;ind = where(absb gt -19)
  ;x = dat[ind].btc

  plot, findgen(10), xtitle='!6'

  psfile = '../plots/ston_vs_absb.eps'
  pnfile = '../plots/ston_vs_absb.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile

  ind = where(finite(absb) and (dat.dist_mpc gt 25.))
  x = absb[ind]
  y = (alog10(tab.flux_wise4/tab.std_flux_wise4))[ind]
  plot, x, y $
        , ps=1 $
        ;, xrange=[12, 16] $
        , xrange=[-22,-17] $
        , yrange=[0.0, 3.0] $
        , xtitle='!6Absolute B Magnitude [mag]' $
        , ytitle='!6log!d10!n Signal-to-Noise at WISE4' $
        , xthick=5, ythick=5 $
        , charthick=3, charsize=1.0 $
        , /nodata
  oplot, x, y, ps=1, symsize=0.5, color=cgcolor('gray')
  for ii = -40, 40 do $
     oplot, ii*1.0*[1,1], [-10,10], lines=1, color=cgcolor('charcoal')
  for ii = -40, 40 do $
     oplot, [-100,100], ii*0.5*[1,1], lines=1, color=cgcolor('charcoal')

  bins = bin_data(x, y, /nan $
                  ;, xmin=12.5, xmax=15.5, binsize=0.25)
                  , xmin=-21.5, xmax=-18, binsize=0.25)
  oploterror, bins.xmid, bins.ymed, bins.ymad $
              , color=cgcolor('firebrick'), errcolor=cgcolor('firebrick') $
              , /nohat, errthick=5, ps=cgsymcat('filledcircle') $
              , symsize=1.5
  
  oplot, -18.*[1,1], [-100, 100], color=cgcolor('navy'), thick=10
  xyouts, -18.15, [2.5], align=0.5, orient=90, '!8M!dB!n!6 = -18' $
          , charthick=3, charsize=1.25, color=cgcolor('navy')
  
  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

  stop

end
