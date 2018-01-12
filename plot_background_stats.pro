pro plot_background_stats $
   , wise=wise $
   , galex=galex

  if keyword_set(wise) then begin  

     spawn, 'rm -rf ../plots/outliers/wise_bkgrd/*.png'

     perc_lo = 2.0/1000.
     perc_hi = 1.0 - perc_lo

     readcol, '../measurements/unwise_bkgrd.txt', format='A,I,I,I,F,F,F' $
              , gal, band, x, y, c0, c1, c2

     for this_band = 1, 4 do begin
        ind = where(band eq this_band)
        
        if this_band eq 1 then begin
           xmin = -0.1
           xmax = 0.1
           binsize = 1d-3
        endif
        if this_band eq 2 then begin
           xmin = -0.1
           xmax = 0.1
           binsize = 1d-3
        endif
        if this_band eq 3 then begin
           xmin = -1.0
           xmax = 1.0
           binsize = 1d-2
        endif
        if this_band eq 4 then begin
           xmin = -1.5
           xmax = 1.5
           binsize = 1d-2
        endif

        vec = [(c1*x)[ind], (c2*y)[ind]]
        vec = vec[sort(vec)]
        bins_c1c2 = bin_data(vec, vec*0.0+1.0 $
                             , xmin=xmin, xmax=xmax, binsize=binsize)
        n = n_elements(vec)
        lo_c1c2 = vec[round(n*perc_lo)]
        hi_c1c2 = vec[round(n*perc_hi)]

        vec = [(c0)[ind]]
        vec = vec[sort(vec)]
        bins_c0 = bin_data(vec, vec*0.0+1.0 $
                           , xmin=xmin, xmax=xmax, binsize=binsize)
        n = n_elements(vec)
        lo_c0 = vec[round(n*perc_lo)]
        hi_c0 = vec[round(n*perc_hi)]

        psfile = '../plots/unwise_background_band'+str(this_band)+'.eps'
        ps, /def, /ps, xs=8, ys=8, /color, /encaps $
            , file=psfile

        plot $
           , [0], [0], /nodata $
           , xtitle='!6Background Magnitude [MJy sr!u-1!n]' $
           , ytitle='!6log!d10!n Number of Images' $
           , xthick=5, ythick=5, charthick=3, charsize=1.5 $
           , xrange=[xmin, xmax], yrange=[0., 4.]

        for ii = -100, 100 do $
           oplot, ii*binsize*10.*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')

        for ii = -100, 100 do $
           oplot, [-10, 10], ii*0.25*[1,1], lines=1, color=cgcolor('charcoal')

        oplot, lo_c1c2*[1,1], [-10, 10], lines=0, color=cgcolor('black'), thick=10
        oplot, hi_c1c2*[1,1], [-10, 10], lines=0, color=cgcolor('black'), thick=10
        oplot, lo_c1c2*[1,1], [-10, 10], lines=0, color=cgcolor('salmon'), thick=3
        oplot, hi_c1c2*[1,1], [-10, 10], lines=0, color=cgcolor('salmon'), thick=3

        oplot, lo_c0*[1,1], [-10, 10], lines=0, color=cgcolor('black'), thick=10
        oplot, hi_c0*[1,1], [-10, 10], lines=0, color=cgcolor('black'), thick=10
        oplot, lo_c0*[1,1], [-10, 10], lines=0, color=cgcolor('royalblue'), thick=3
        oplot, hi_c0*[1,1], [-10, 10], lines=0, color=cgcolor('royalblue'), thick=3

        histplot $
           , bins_c1c2.xmid, (alog10(bins_c1c2.counts) > (0.)) $
           , /overplot $
           , lthick=3, /nobar, /fill, fcolor=cgcolor('salmon') $
           , lcolor=cgcolor('firebrick')

        histplot $
           , bins_c0.xmid, (alog10(bins_c0.counts) > (0.)) $
           , /overplot $
           , lthick=3, /nobar, /fline, forient=45. $
           , fcolor=cgcolor('navy'), lcolor=cgcolor('navy')

        histplot $
           , bins_c0.xmid, (alog10(bins_c0.counts) > (0.)) $
           , /overplot $
           , lthick=3, /nobar, /fline, forient=-45. $
           , fcolor=cgcolor('navy'), lcolor=cgcolor('navy')
        
        al_legend $
           , /top, /left $
           , box=0, clear=0, charsize=1.75, charthick=3 $
           , lines=-99 $
           , ['!6WISE '+str(this_band) $
              , 'Tilt', 'Zero point'] $
           , textcolor=[cgcolor('black'), cgcolor('salmon') $
                        , cgcolor('royalblue')]

        ps, /xw
        spawn, 'evince '+psfile+' &'

        outliers = $
           where(band eq this_band and $
                 (c0 lt lo_c0 or $
                 x*c1 lt lo_c1c2 or $
                 y*c2 lt lo_c1c2 or $
                 c0 gt hi_c0 or $
                 x*c1 gt hi_c1c2 or $
                 y*c2 gt hi_c1c2) $
                 , outlier_ct)
        
        for jj = 0, outlier_ct-1 do begin

           im = $
              readfits('../unwise/atlas/'+gal[outliers[jj]]+ $
                       '_w'+str(this_band)+'_bksub.fits',hdr)
           
           loadct, 33
           disp, im, max=mad(im)*3., min=-3.*mad(im), /sq
           im = tvrd(true=1)
           write_png, '../plots/outliers/wise_bkgrd/'+ $
                      gal[outliers[jj]]+'_w'+str(this_band)+'_bksub.png' $
                      , im
        endfor
        
     endfor

  endif

end
  

end
