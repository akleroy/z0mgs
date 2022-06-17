pro plot_galex_noise_scaling 

  tab = mrdfits('../measurements/delivery_index_gauss15.fits',1,h)
  b = tab.gb_deg

  plot, findgen(10), title='!6Test'

  ind = where(abs(b) gt 40. and tab.has_fuv)
  fuv_t = alog10((tab.time_fuv)[ind])
  fuv_rms = alog10((tab.rms_fuv)[ind])

  ind = where(abs(b) gt 40. and tab.has_nuv)
  nuv_t = alog10((tab.time_nuv)[ind])
  nuv_rms = alog10((tab.rms_nuv)[ind])

  nuv_bins = $
     bin_data(nuv_t, nuv_rms, /nan $
              , xmin=1.75, xmax=4.5, binsize=0.25)

  fuv_bins = $
     bin_data(fuv_t, fuv_rms, /nan $
              , xmin=1.70, xmax=4.5, binsize=0.25)

  psfile = '../plots/galex_noise_vs_t.eps'
  ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
      , file=psfile

  plot $
     , [0], [0], /nodata $
     , ytitle='!6log!d10!n Robust Noise [MJy sr!u-1!n]' $
     , xtitle='!6log!d10!n Effective Integration [s]' $
     , xthick=5, ythick=5, charthick=3, charsize=1.25 $
     , xrange=[1., 5.], yrange=[-4.75, -3.25], /ystyle $
     , /xstyle

  for ii = 0, 20 do $
     oplot, ii*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')

  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.25*[1,1], lines=1, color=cgcolor('charcoal')

  xfid = findgen(101)/10.
  oplot, xfid, -0.5*(xfid-2.0)-3.6, color=cgcolor('gray'), thick=10

  oplot, nuv_t, nuv_rms, ps=3, color=cgcolor('salmon')
  oplot, fuv_t, fuv_rms, ps=3, color=cgcolor('dodgerblue')

  oploterror, nuv_bins.xmid, nuv_bins.ymed, nuv_bins.ymad $
              , ps=cgsymcat('filledcircle'), color=cgcolor('firebrick') $
              , errthick=5, /nohat
  oploterror, fuv_bins.xmid, fuv_bins.ymed, fuv_bins.ymad $
              , ps=cgsymcat('filledcircle'), color=cgcolor('navy') $
              , errthick=5, /nohat

  al_legend $
     , /top, /right $
     , box=1, clear=1, background=cgcolor('lightgray'), outline=cgcolor('black') $
     , charsize=1.25, charthick=3 $
     , lines=-99 $
     , ['!6NUV', 'FUV'] $
     , textcolor=[cgcolor('firebrick'), cgcolor('navy')]

  ps, /xw
  spawn, 'evince '+psfile+' &'

  stop

end
