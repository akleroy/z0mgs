pro plot_unwise_noise

  perc_lo = 2.0/1000.
  perc_hi = 1.0 - perc_lo

;  restore, '../measurements/unwise_stats_with_dat.idl', /v
  tab = mrdfits('../measurements/delivery_index.fits',1,h)
  b = tab.gb_deg

  plot, findgen(10), title='!6Test'

  for this_band = 0, 3 do begin

     if this_band eq 0 then begin
        xmin = -3.0
        xmax = -1.0
        binsize = 0.01
        stat = tab.rms_wise1
     endif
     if this_band eq 1 then begin
        xmin = -3.0
        xmax = -1.0
        binsize = 0.01
        stat = tab.rms_wise2
     endif
     if this_band eq 2 then begin
        xmin = -2.5
        xmax = -0.5
        binsize = 0.01
        stat = tab.rms_wise3
     endif
     if this_band eq 3 then begin
        xmin = -1.5
        xmax = 0.0
        binsize = 0.01
        stat = tab.rms_wise4
     endif

     high_b = where(abs(b) gt 40.)
     low_b = where(abs(b) le 40.)

     vec = (alog10(stat))[high_b]
     bins_15 = bin_data(vec, vec*0.0+1.0 $
                        , xmin=xmin, xmax=xmax, binsize=binsize, /nan)
     vec = vec[sort(vec)]
     n = n_elements(vec)
     print, 'WISE BAND '+str(this_band+1)
     print, '... 16: ', 10.^(vec[round(n*0.16)])
     print, '... 50: ', 10.^(vec[round(n*0.50)])
     print, '... 84: ', 10.^(vec[round(n*0.84)])

     vec = (alog10(stat))[low_b]
     bins_15_lowb = bin_data(vec, vec*0.0+1.0 $
                             , xmin=xmin, xmax=xmax, binsize=binsize, /nan)

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
        , bins_15.xmid, (alog10(bins_15.counts) > (0.)) $
        , /overplot $
        , lthick=3, /nobar, /fill $
        , fcolor=cgcolor('salmon') $
        , lcolor=cgcolor('firebrick')

     histplot $
        , bins_15_lowb.xmid, (alog10(bins_15_lowb.counts) > (0.)) $
        , /overplot $
        , lthick=3, /nobar, /fline, forient=45 $
        , fcolor=cgcolor('royalblue') $
        , lcolor=cgcolor('royalblue')

     histplot $
        , bins_15_lowb.xmid, (alog10(bins_15_lowb.counts) > (0.)) $
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
        , ['!6WISE '+str(this_band+1)+' at 15"' $
           , '!3|!8b!6!9!3|!6 > 40'+textoidl('\circ'), '!3|!8b!3|!6 < 40'+textoidl('\circ')] $
        , textcolor=[cgcolor('black'), cgcolor('salmon') $
                     , cgcolor('royalblue')]

     ps, /xw
     spawn, 'evince '+psfile+' &'
     
  endfor

  stop

end
