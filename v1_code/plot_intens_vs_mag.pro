pro plot_intens_vs_mag

  for jj = 0, 1 do begin
     
     if jj eq 0 then begin
        res_str = 'gauss7p5'
        legend_string = '7.5" Resolution'
     endif
     if jj eq 1 then begin
        res_str = 'gauss15'
        legend_string = '15" Resolution'
     endif

     restore, '../measurements/star_stacks/star_stack_w1_'+res_str+'.idl'
     w1_val = val
     w1_mag = mag
     w1_bins = bin_data(alog10(w1_val), w1_mag, /nan $
                        , xmin=-1.0, xmax=2.5, binsize=0.1)
     w1_coef = median(w1_val/(10.^(-1.*w1_mag/2.5)))

     restore, '../measurements/star_stacks/star_stack_w2_'+res_str+'.idl'
     w2_val = val
     w2_mag = mag
     w2_bins = bin_data(alog10(w2_val), w2_mag, /nan $
                        , xmin=-1.0, xmax=2.5, binsize=0.1)
     w2_coef = median(w2_val/(10.^(-1.*w2_mag/2.5)))

     restore, '../measurements/star_stacks/star_stack_w3_'+res_str+'.idl'
     w3_val = val
     w3_mag = mag
     w3_bins = bin_data(alog10(w3_val), w3_mag, /nan $
                        , xmin=-1.0, xmax=2.5, binsize=0.1)
     w3_coef = median(w3_val/(10.^(-1.*w3_mag/2.5)))

     if res_str eq 'gauss15' then begin
        restore, '../measurements/star_stacks/star_stack_w4_'+res_str+'.idl'
        w4_val = val
        w4_mag = mag
        w4_bins = bin_data(alog10(w4_val), w4_mag, /nan $
                           , xmin=-1.0, xmax=2.5, binsize=0.1)
        w4_coef = median(w4_val/(10.^(-1.*w4_mag/2.5)))
     endif

     plot, findgen(10)
     psfile = '../plots/mag_vs_intens_wise_'+res_str+'.eps'
     pnfile = '../plots/mag_vs_intens_wise_'+res_str+'.png'
     ps, /def, /ps, xs=5, ys=5, /color, /encaps $
         , file=psfile
     
     plot $
        , [0], [0], /nodata $
        , ytitle='!6log!d10!n Intensity [MJy sr!u-1!n]' $
        , xtitle='!6log!d10!n Ks (2MASS) Magnitude' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , yrange=[-1.0, 2.5], xrange=[11., 0.] $
        , /xstyle
     
     for ii = -100, 100 do $
        oplot, ii*1.0*[1,1], [-100,100], lines=1
     for ii = -100, 100 do $
        oplot, [-100,100], ii*0.5*[1,1], lines=1

     if res_str eq 'gauss15' then begin
        oplot, w4_mag, alog10(w4_val), ps=1, symsize=0.5, color=cgcolor('lightseagreen')
     endif

     oplot, w3_mag, alog10(w3_val), ps=1, symsize=0.5, color=cgcolor('goldenrod')
     oplot, w2_mag, alog10(w2_val), ps=1, symsize=0.5, color=cgcolor('salmon')
     oplot, w1_mag, alog10(w1_val), ps=1, symsize=0.5, color=cgcolor('firebrick')
     
     fid = 10^(findgen(101)/100.*5.-2.)
     
     oplot, -1.*alog10(1./w1_coef*fid)*2.5, alog10(fid) $
            , thick=3, lines=2, color=cgcolor('black')
     oplot, -1.*alog10(1./w2_coef*fid)*2.5, alog10(fid) $
            , thick=3, lines=2, color=cgcolor('black')
     oplot, -1.*alog10(1./w3_coef*fid)*2.5, alog10(fid) $
            , thick=3, lines=2, color=cgcolor('black')
     
     if res_str eq 'gauss15' then begin
        oplot, -1.*alog10(1./w4_coef*fid)*2.5, alog10(fid) $
               , thick=3, lines=2, color=cgcolor('black')
     endif     
     
     al_legend, /bottom, /right, [legend_string,'WISE1','WISE2','WISE3','WISE4'] $
                , box=1, clear=1, lines=-99 $
                , background=cgcolor('lightgray') $
                , textcolor=[cgcolor('black'),cgcolor('firebrick'),cgcolor('salmon'),cgcolor('goldenrod'),cgcolor('lightseagreen')] $
                , charsize=1.25, charthick=3

     ps, /xw
     spawn, 'evince '+psfile+' &'     
     spawn, 'convert -density 300x300 '+psfile+' '+pnfile     

     restore, '../measurements/star_stacks/star_stack_nuv_'+res_str+'.idl'
     nuv_val = val
     nuv_mag = mag
     nuv_bins = bin_data(alog10(nuv_val), nuv_mag, /nan $
                        , xmin=-1.0, xmax=2.5, binsize=0.1)
     nuv_coef = median(nuv_val/(10.^(-1.*nuv_mag/2.5)))

     restore, '../measurements/star_stacks/star_stack_fuv_'+res_str+'.idl'
     fuv_val = val
     fuv_mag = mag
     fuv_bins = bin_data(alog10(fuv_val), fuv_mag, /nan $
                        , xmin=-1.0, xmax=2.5, binsize=0.1)
     fuv_coef = median(fuv_val/(10.^(-1.*fuv_mag/2.5)))
     
     plot, findgen(10)
     psfile = '../plots/mag_vs_intens_galex_'+res_str+'.eps'
     pnfile = '../plots/mag_vs_intens_galex_'+res_str+'.png'
     ps, /def, /ps, xs=5, ys=5, /color, /encaps $
         , file=psfile
     
     plot $
        , [0], [0], /nodata $
        , ytitle='!6log!d10!n Intensity [MJy sr!u-1!n]' $
        , xtitle='!6log!d10!n Ks (2MASS) Magnitude' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , yrange=[-4.0, 1.0], xrange=[11., 0.] $
        , /xstyle
     
     for ii = -100, 100 do $
        oplot, ii*1.0*[1,1], [-100,100],  lines=1
     for ii = -100, 100 do $
        oplot, [-100,100], ii*0.5*[1,1], lines=1

     oplot, nuv_mag, alog10(nuv_val), ps=1, symsize=0.5, color=cgcolor('dodgerblue')
     oplot, fuv_mag, alog10(fuv_val), ps=1, symsize=0.5, color=cgcolor('orchid')
     
     fid = 10^(findgen(101)/100.*10.-5.)
     
     oplot $
        , -1.*alog10(1./nuv_coef*fid)*2.5, alog10(fid), thick=3, lines=2, color=cgcolor('black')
     oplot $
        , -1.*alog10(1./fuv_coef*fid)*2.5, alog10(fid), thick=3, lines=2, color=cgcolor('black')
     
     al_legend, /bottom, /right, [legend_string,'FUV','NUV'] $
                , box=1, clear=1, lines=-99 $
                , background=[cgcolor('lightgray')] $
                , textcolor=[cgcolor('black'),cgcolor('orchid'),cgcolor('dodgerblue')] $
                , charsize=1.25, charthick=3

     ps, /xw
     spawn, 'evince '+psfile+' &'     
     spawn, 'convert -density 300x300 '+psfile+' '+pnfile     

     print, "Best Fit Coefficients (Note Scatter for UV):"
     print, 'WISE 1: ', w1_coef
     print, 'WISE 2: ', w2_coef
     print, 'WISE 3: ', w3_coef
     if res_str eq 'gauss15' then begin
        print, 'WISE 4: ', w4_coef
     endif
     print, 'NUV: ', nuv_coef
     print, 'FUV: ', fuv_coef

  endfor


  stop

end
