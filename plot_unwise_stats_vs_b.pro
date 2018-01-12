pro plot_unwise_stats_vs_b $
   , wise=wise $
   , galex=galex

  restore, '../measurements/unwise_stats_with_dat.idl', /v

  plot, findgen(10), title='!6Test'

  for this_band = 0, 3 do begin

     if this_band eq 0 then begin
        ymin = -3.0
        ymax = 1.0
        binsize = 0.01
     endif
     if this_band eq 1 then begin
        ymin = -3.0
        ymax = 1.0
        binsize = 0.01
     endif
     if this_band eq 2 then begin
        ymin = -2.0
        ymax = 0.5
        binsize = 0.01
     endif
     if this_band eq 3 then begin
        ymin = -1.0
        ymax = 0.5
        binsize = 0.01
     endif

     this_med = (stats[*,2,this_band,0].med)
     this_mean = (stats[*,2,this_band,0].mean)
     this_rms =  (stats[*,2,this_band,0].rms)
     this_std =  (stats[*,2,this_band,0].std)
     this_rej =  (stats[*,2,this_band,0].rejfrac)

     bmin = -87.5
     bmax = +87.5
     bbinsize = 5.0
     mean_bins = $
        bin_data(b, this_mean, /nan $
                 , xmin=bmin, xmax=bmax, binsize=bbinsize)
     rms_bins = $
        bin_data(b, this_rms, /nan $
                 , xmin=bmin, xmax=bmax, binsize=bbinsize)
     std_bins = $
        bin_data(b, this_std, /nan $
                 , xmin=bmin, xmax=bmax, binsize=bbinsize)
     rej_bins = $
        bin_data(b, this_rej, /nan $
                 , xmin=bmin, xmax=bmax, binsize=bbinsize)
     
     psfile = '../plots/unwise_noise_lat_band'+str(this_band+1)+'.eps'
     ps, /def, /ps, xs=8, ys=8, /color, /encaps $
         , file=psfile
     
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

;     oplot, b, alog10(this_rms), ps=1, color=cgcolor('olive')
;     oplot, b, alog10(this_std), ps=1, color=cgcolor('orchid')
;     oplot, b, alog10(this_mean), ps=3, color=cgcolor('salmon')
;     oplot, b, alog10(this_rej), ps=3, color=cgcolor('goldenrod')

     circle, /fill
;     oplot, mean_bins.xmid, alog10(mean_bins.ymed), ps=8, symsize=1.2, color=cgcolor('black')
;     oplot, mean_bins.xmid, alog10(mean_bins.ymed), ps=8, symsize=0.8, color=cgcolor('firebrick')

     oploterror, std_bins.xmid-1.0, alog10(std_bins.ymed), std_bins.ymad_log $
                 , ps=8, symsize=2.0, color=cgcolor('black'), errthick=3
     oplot, std_bins.xmid-1.0, alog10(std_bins.ymed), ps=8, symsize=1.0, color=cgcolor('purple')

     oploterror, rms_bins.xmid+1.0, alog10(rms_bins.ymed), rms_bins.ymad_log $
            , ps=8, symsize=2.0, color=cgcolor('black'), errthick=3
     oplot, rms_bins.xmid+1.0, alog10(rms_bins.ymed), ps=8, symsize=1.0, color=cgcolor('darkgreen')

;     oplot, rms_bins.xmid, alog10(rej_bins.ymed), ps=8, symsize=1.2, color=cgcolor('black')
;     oplot, rms_bins.xmid, alog10(rej_bins.ymed), ps=8, symsize=0.8, color=cgcolor('goldenrod')

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
        , textcolor=[cgcolor('purple'), cgcolor('darkgreen')]

     ps, /xw
     spawn, 'evince '+psfile+' &'

  endfor

  stop

end
