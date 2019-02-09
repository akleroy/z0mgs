pro plot_galex_stats_vs_b

  tab15 = mrdfits('../measurements/delivery_index_gauss15.fits',1,h)
  tab7p5 = mrdfits('../measurements/delivery_index_gauss7p5.fits',1,h)
  if total(tab15.pgc ne tab7p5.pgc) then begin
     print, "Some mismatch in tables."
     stop
  endif
  b = tab15.gb_deg

  plot, findgen(10), title='!6Test'

  for this_band = 0, 1 do begin

     if this_band eq 0 then begin
        band_str = 'NUV'     
        this_rms15 = tab15.rms_nuv
        this_std15 = tab15.std_nuv
        this_rms7p5 = tab7p5.rms_nuv
        this_std7p5 = tab7p5.std_nuv
        ymin = -4.5
        ymax = -1.75
     endif
     if this_band eq 1 then begin
        band_str = 'FUV'
        this_rms15 = tab15.rms_fuv
        this_std15 = tab15.std_fuv
        this_rms7p5 = tab7p5.rms_fuv
        this_std7p5 = tab7p5.std_fuv
        ymin = -4.5
        ymax = -2.0     
     endif

     bmin = 2.5
     bmax = 87.5
     bbinsize = 5.0

     rms_bins15 = $
        bin_data(b, this_rms15, /nan $
                 , xmin=bmin, xmax=bmax, binsize=bbinsize)
     std_bins15 = $
        bin_data(b, this_std15, /nan $
                 , xmin=bmin, xmax=bmax, binsize=bbinsize)
     rms_bins7p5 = $
        bin_data(b, this_rms7p5, /nan $
                 , xmin=bmin, xmax=bmax, binsize=bbinsize)
     std_bins7p5 = $
        bin_data(b, this_std7p5, /nan $
                 , xmin=bmin, xmax=bmax, binsize=bbinsize)
     
     psfile = '../plots/galex_noise_lat_'+band_str+'.eps'
     ps, /def, /ps, xs=8, ys=4, /color, /encaps $
         , file=psfile
     
     plot $
        , [0], [0], /nodata $
        , xtitle='!3|!6 Galactic Latitude !3|!6' $
        , ytitle='!6log!d10!n Value [MJy sr!u-1!n]' $
        , xthick=5, ythick=5, charthick=3, charsize=1.5 $
        , xrange=[0., 100.], yrange=[ymin, ymax] $
        , ystyle=1

     for ii = -100, 100 do $
        oplot, ii*15.*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')

     for ii = -100, 100 do $
        oplot, [-100, 100], ii*0.25*[1,1], lines=1, color=cgcolor('charcoal')

     circle, /fill

     oplot,  std_bins7p5.xmid, alog10(std_bins7p5.ymed), lines = 2 $
             , thick=10, color=cgcolor('firebrick')

     oplot,  std_bins15.xmid, alog10(std_bins15.ymed), lines = 0 $
             , thick=10, color=cgcolor('firebrick')

     oplot,  rms_bins7p5.xmid, alog10(rms_bins7p5.ymed), lines = 2 $
             , thick=10, color=cgcolor('navy')

     oplot,  rms_bins15.xmid, alog10(rms_bins15.ymed), lines = 0 $
             , thick=10, color=cgcolor('navy')

     al_legend $
        , /top, /right $
        , box=1, clear=1, background=cgcolor('lightgray'), charsize=1.25, charthick=3 $
        , lines=-99 $
        , ['!6GALEX '+band_str, '!6Std. Dev.', 'Robust Noise'] $
        , textcolor=[cgcolor('black'), cgcolor('firebrick'), cgcolor('navy')]

     al_legend $
        , /bottom, /left $
        , box=1, clear=1, background=cgcolor('lightgray'), charsize=1.25, charthick=3 $
        , ['!67.5"','15"'] $
        , color=cgcolor('black'), thick=5, pspacing=1.0, lines=[2,0]

     ps, /xw
     spawn, 'evince '+psfile+' &'

  endfor

  stop

end
