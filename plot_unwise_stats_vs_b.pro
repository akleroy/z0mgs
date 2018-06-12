pro plot_unwise_stats_vs_b

  tab = mrdfits('../measurements/delivery_index.fits',1,h)
  b = tab.gb_deg

  plot, findgen(10), title='!6Test'

  for this_band = 0, 3 do begin

     if this_band eq 0 then begin
        this_rms = tab.rms_wise1
        this_std = tab.std_wise1
     endif
     if this_band eq 1 then begin
        this_rms = tab.rms_wise2
        this_std = tab.std_wise2
     endif
     if this_band eq 2 then begin
        this_rms = tab.rms_wise3
        this_std = tab.std_wise3
     endif
     if this_band eq 3 then begin
        this_rms = tab.rms_wise4
        this_std = tab.std_wise4
     endif

     bmin = -87.5
     bmax = +87.5
     bbinsize = 5.0

     rms_bins = $
        bin_data(b, this_rms, /nan $
                 , xmin=bmin, xmax=bmax, binsize=bbinsize)
     std_bins = $
        bin_data(b, this_std, /nan $
                 , xmin=bmin, xmax=bmax, binsize=bbinsize)
     
     psfile = '../plots/unwise_noise_lat_band'+str(this_band+1)+'.eps'
     ps, /def, /ps, xs=8, ys=4, /color, /encaps $
         , file=psfile
     
     ymin = -3.
     ymax = 1.0
     plot $
        , [0], [0], /nodata $
        , xtitle='!6Galactic Latitude' $
        , ytitle='!6log!d10!n Value [MJy sr!u-1!n]' $
        , xthick=5, ythick=5, charthick=3, charsize=1.5 $
        , xrange=[-100., 100.], yrange=[ymin, ymax]

     for ii = -100, 100 do $
        oplot, ii*15.*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')

     for ii = -100, 100 do $
        oplot, [-100, 100], ii*0.25*[1,1], lines=1, color=cgcolor('charcoal')

     circle, /fill

     oploterror, std_bins.xmid-1.0, alog10(std_bins.ymed), std_bins.ymad_log $
                 , ps=cgsymcat('filledsquare') $
                 , symsize=2.0, color=cgcolor('black'), errthick=3
     oplot, std_bins.xmid-1.0, alog10(std_bins.ymed), ps=cgsymcat('filledsquare') $
            , symsize=1.0 $
            , color=cgcolor('salmon')

     oploterror, rms_bins.xmid+1.0, alog10(rms_bins.ymed), rms_bins.ymad_log $
            , ps=cgsymcat('filledcircle'), symsize=2.0, color=cgcolor('black'), errthick=3
     oplot, rms_bins.xmid+1.0, alog10(rms_bins.ymed), ps=cgsymcat('filledcircle') $
            , symsize=1.0, color=cgcolor('dodgerblue')

     al_legend $
        , /top, /right $
        , box=0, clear=0, charsize=1.75, charthick=3 $
        , lines=-99 $
        , ['!6WISE '+str(this_band+1)]

     al_legend $
        , /top, /left $
        , box=0, clear=0, charsize=1.75, charthick=3 $
        , lines=-99 $
        , ['!6Std. Dev.', 'Robust Noise'] $
        , textcolor=[cgcolor('salmon'), cgcolor('dodgerblue')]

     ps, /xw
     spawn, 'evince '+psfile+' &'

  endfor

  stop

end
