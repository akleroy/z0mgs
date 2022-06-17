pro plot_unwise_noise

  perc_lo = 2.0/1000.
  perc_hi = 1.0 - perc_lo

  tab15 = mrdfits('../measurements/delivery_index_gauss15.fits',1,h)
  tab7p5 = mrdfits('../measurements/delivery_index_gauss7p5.fits',1,h)
  if total(tab15.pgc ne tab7p5.pgc) then begin
     print, "Some mismatch in tables."
     stop
  endif

  b = tab15.gb_deg

  plot, findgen(10), title='!6Test'

  for this_band = 0, 3 do begin

     if this_band eq 0 then begin
        xmin = -3.0
        xmax = -1.0
        binsize = 0.01
        stat15 = tab15.rms_wise1
        stat7p5 = tab7p5.rms_wise1
     endif
     if this_band eq 1 then begin
        xmin = -3.0
        xmax = -1.0
        binsize = 0.01
        stat15 = tab15.rms_wise2
        stat7p5 = tab7p5.rms_wise2
     endif
     if this_band eq 2 then begin
        xmin = -3.0
        xmax = 0.5
        binsize = 0.01
        stat15 = tab15.rms_wise3
        stat7p5 = tab7p5.rms_wise3
     endif
     if this_band eq 3 then begin
        xmin = -1.5
        xmax = 0.0
        binsize = 0.01
        stat15 = tab15.rms_wise4
        stat7p5 = tab15.rms_wise4*!values.f_nan
     endif

     high_b = where(abs(b) gt 40.)
     low_b =  where(abs(b) le 40.)

     vec15 = (alog10(stat15))[high_b]
     bins15 = bin_data(vec15, vec15*0.0+1.0 $
                       , xmin=xmin, xmax=xmax, binsize=binsize, /nan)

     vec7p5 = (alog10(stat7p5))[high_b]
     bins7p5 = bin_data(vec7p5, vec7p5*0.0+1.0 $
                        , xmin=xmin, xmax=xmax, binsize=binsize, /nan)
     
     vec15 = vec15[sort(vec15)]
     n = n_elements(vec15)
     print, 'WISE BAND '+str(this_band+1)+' at 15"'
     print, '... 16: ', 10.^(vec15[round(n*0.16)])
     print, '... 50: ', 10.^(vec15[round(n*0.50)])
     print, '... 84: ', 10.^(vec15[round(n*0.84)])

     vec7p5 = vec7p5[sort(vec7p5)]
     n = n_elements(vec7p5)
     print, 'WISE BAND '+str(this_band+1)+' at 7p5"'
     print, '... 16: ', 10.^(vec7p5[round(n*0.16)])
     print, '... 50: ', 10.^(vec7p5[round(n*0.50)])
     print, '... 84: ', 10.^(vec7p5[round(n*0.84)])

     ;vec = (alog10(stat))[low_b]
     ;bins15_lowb = bin_data(vec, vec*0.0+1.0 $
     ;                        , xmin=xmin, xmax=xmax, binsize=binsize, /nan)

     psfile = '../plots/unwise_noise_band'+str(this_band+1)+'.eps'
     ps, /def, /ps, xs=8, ys=4, /color, /encaps $
         , file=psfile

     plot $
        , [0], [0], /nodata $
        , xtitle='!6log!d10!n Robust Noise Estimate [MJy sr!u-1!n]' $
        , ytitle='!6log!d10!n Number of Images' $
        , xthick=5, ythick=5, charthick=3, charsize=1.5 $
        , xrange=[xmin, xmax], yrange=[0., 4.]

     for ii = -100, 100 do $
        oplot, ii*binsize*10.*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')

     for ii = -100, 100 do $
        oplot, [-10, 10], ii*0.25*[1,1], lines=1, color=cgcolor('charcoal')

     histplot $
        , bins15.xmid, (alog10(bins15.counts) > (0.)) $
        , /overplot $
        , lthick=3, /nobar, /fill $
        , fcolor=cgcolor('salmon') $
        , lcolor=cgcolor('firebrick')

     histplot $
        , bins7p5.xmid, (alog10(bins7p5.counts) > (0.)) $
        , /overplot $
        , lthick=3, /nobar, /fline, forient=45 $
        , fcolor=cgcolor('royalblue') $
        , lcolor=cgcolor('royalblue')

     histplot $
        , bins7p5.xmid, (alog10(bins7p5.counts) > (0.)) $
        , /overplot $
        , lthick=3, /nobar, /fline, forient=-45 $
        , fcolor=cgcolor('royalblue') $
        , lcolor=cgcolor('royalblue')

     al_legend $
        , /top, /right $
        , box=1, clear=1 $
        , background=cgcolor('lightgray') $
        , charsize=1.75, charthick=3 $
        , lines=-99 $
        , ['!6WISE '+str(this_band+1)+', !3|!8b!6!9!3|!6 > 40'+textoidl('\circ') $
           , '... at 15"' $
           , '... at 7.5"'] $
        , textcolor=[cgcolor('black'), cgcolor('salmon') $
                     , cgcolor('royalblue')]

     ps, /xw
     spawn, 'evince '+psfile+' &'
     
  endfor

  stop

end
