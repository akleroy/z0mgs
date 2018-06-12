pro plot_galex_stats_vs_b

  tab = mrdfits('../measurements/delivery_index.fits',1,h)
  b = tab.gb_deg

;  restore, '../measurements/unwise_stats_with_dat.idl', /v

  plot, findgen(10), title='!6Test'

;  restore, file='../measurements/unwise_stats_with_dat.idl', /v
;  restore, file='../measurements/galex_stats.idl'

  plot, findgen(10), title='!6Test'

  for this_band = 0, 1 do begin

     if this_band eq 0 then begin
        band_str = 'NUV'     
        this_rms = tab.rms_nuv
        this_std = tab.std_nuv
     endif
     if this_band eq 1 then begin
        band_str = 'FUV'
        this_rms = tab.rms_fuv
        this_std = tab.std_fuv
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
     
     psfile = '../plots/galex_noise_lat_'+band_str+'.eps'
     ps, /def, /ps, xs=8, ys=4, /color, /encaps $
         , file=psfile
     
     ymin = -4.0
     ymax = 0.0
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
        , ['!6GALEX '+band_str]

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
