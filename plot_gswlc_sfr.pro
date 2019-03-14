pro plot_gswlc_sfr

  thresh = 10
  @constants.bat
  lsun_3p4 = 1.83d18
  restore, '../gswlc/gswlc_data.idl'

  nu_fuv = c/(154.d-9*1d2)
  nu_nuv = c/(231.d-9*1d2)
  nu_w1 = c/(3.4d-6*1d2)
  nu_w2 = c/(4.5d-6*1d2)
  nu_w3 = c/(12.d-6*1d2)
  nu_w4 = c/(22.d-6*1d2)

  plot, findgen(10), xtitle='!6'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE SUBSET OF DATA TO USE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sfr_ind = where(gws_logmstar gt 0 $
                  and w1_lum gt 0 $
                  and gws_w1 gt 3.*gws_ew1 $
                  and gws_w3 gt 3.*gws_ew3 $
                  and gws_w4 gt 3.*gws_ew4 $
                  and gws_nuv gt 3.*gws_enuv $
                  and gws_fuv gt 3.*gws_efuv $
                  and mtol_w1 gt 0.02 and mtol_w1 lt 1.0 $
                  and gws_ssfr gt -10.5 $
                 )

  red_ind = where(gws_logmstar gt 0 $
                  and w1_lum gt 0 $
                  and gws_w1 gt 3.*gws_ew1 $
                  and gws_w3 gt 3.*gws_ew3 $
                  and gws_w4 gt 3.*gws_ew4 $
                  and gws_nuv gt 3.*gws_enuv $
                  and gws_fuv gt 3.*gws_efuv $
                  and mtol_w1 gt 0.02 and mtol_w1 lt 1.0 $
                  and gws_ssfr lt -10.5 $
                 )

  nowise_ind = where(gws_logmstar gt 0 $
                     and w1_lum gt 0 $
                     and gws_w1 gt 3.*gws_ew1 $
                     ;and gws_w3 gt 3.*gws_ew3 $
                     ;and gws_w4 gt 3.*gws_ew4 $
                     and gws_nuv gt 3.*gws_enuv $
                     and gws_fuv gt 3.*gws_efuv $
                     and mtol_w1 gt 0.02 and mtol_w1 lt 1.0 $
                     ;and gws_ssfr lt -10.75 $
                    )

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ESTIMATE THE WISE COEFFICIENTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  coef_fuv = ((10.d^(gws_logsfrsed))) / (nu_fuv*fuv_lum*10.^(gws_afuv/2.5))
  coef_justfuv = ((10.d^(gws_logsfrsed))) / (nu_fuv*fuv_lum)
  coef_justnuv = ((10.d^(gws_logsfrsed))) / (nu_nuv*nuv_lum)

  coef_justw3 = ((10.d^(gws_logsfrsed))) / (nu_w3*w3_lum)
  coef_justw4 = ((10.d^(gws_logsfrsed))) / (nu_w4*w4_lum)
  
  coef_w3nuv = ((10.d^(gws_logsfrsed)) - sfr_nuv_z19) /  (nu_w3*w3_lum)
  coef_w4nuv = ((10.d^(gws_logsfrsed)) - sfr_nuv_z19) /  (nu_w4*w4_lum)

  coef_w3fuv = ((10.d^(gws_logsfrsed)) - sfr_fuv_z19) /  (nu_w3*w3_lum)
  coef_w4fuv = ((10.d^(gws_logsfrsed)) - sfr_fuv_z19) /  (nu_w4*w4_lum)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER BAND COMBINATIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  loadct, 0

  for ii = 0, 5 do begin
     
     x = gws_logmstar[sfr_ind]

     if ii eq 0 then begin
        psfile = '../plots/gswlc_justw3_coef.eps'
        pnfile = '../plots/gswlc_justw3_coef.png'

        histpsfile = '../plots/gswlc_justw3_hist.eps'
        histpnfile = '../plots/gswlc_justw3_hist.png'

        coef = coef_justw3
        tag = ['!6WISE3 only']

        ref = 10.d^(-42.9)

        xmin = -43.5
        xmax = -41.5
     endif

     if ii eq 1 then begin
        psfile = '../plots/gswlc_justw4_coef.eps'
        pnfile = '../plots/gswlc_justw4_coef.png'

        histpsfile = '../plots/gswlc_justw4_hist.eps'
        histpnfile = '../plots/gswlc_justw4_hist.png'

        coef = coef_justw4
        tag = ['!6WISE4 only']

        ref = 10.d^(-42.7)*22./24. ; MIPS -> WISE

        xmin = -43.5
        xmax = -41.5
     endif

     if ii eq 2 then begin
        psfile = '../plots/gswlc_nuvw3_coef.eps'
        pnfile = '../plots/gswlc_nuvw3_coef.png'

        histpsfile = '../plots/gswlc_nuvw3_hist.eps'
        histpnfile = '../plots/gswlc_nuvw3_hist.png'

        coef = coef_w3nuv
        tag = ['!6WISE3','in NUV+WISE3']

        ref = 10.d^(-42.9)

        xmin = -43.5
        xmax = -41.5
     endif

     if ii eq 3 then begin
        psfile = '../plots/gswlc_nuvw4_coef.eps'
        pnfile = '../plots/gswlc_nuvw4_coef.png'

        histpsfile = '../plots/gswlc_nuvw4_hist.eps'
        histpnfile = '../plots/gswlc_nuvw4_hist.png'

        coef = coef_w4nuv
        tag = ['!6WISE4','in NUV+WISE4']

        ref = 10.d^(-42.81)*22./25. ; KENNICUTT/HAO -> WISE4 from IRAS2

        xmin = -43.5
        xmax = -41.5
     endif

     if ii eq 4 then begin
        psfile = '../plots/gswlc_fuvw3_coef.eps'
        pnfile = '../plots/gswlc_fuvw3_coef.png'

        histpsfile = '../plots/gswlc_fuvw3_hist.eps'
        histpnfile = '../plots/gswlc_fuvw3_hist.png'

        coef = coef_w3fuv
        tag = ['!6WISE3','in FUV+WISE3']

        ref = 10.d^(-42.9)

        xmin = -43.5
        xmax = -41.5
     endif

     if ii eq 5 then begin
        psfile = '../plots/gswlc_fuvw4_coef.eps'
        pnfile = '../plots/gswlc_fuvw4_coef.png'

        histpsfile = '../plots/gswlc_fuvw4_hist.eps'
        histpnfile = '../plots/gswlc_fuvw4_hist.png'

        coef = coef_w4fuv
        tag = ['!6WISE4','in FUV+WISE4']

        ref = 10.d^(-42.76)*22./25. ; KENNICUTT IRAS -> WISE

        xmin = -43.5
        xmax = -41.5
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT HISTOGRAMS OF THE SOLVED FOR COEFFICIENTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
     binsize = 0.05

     vec = alog10(coef[sfr_ind])
     vec = vec[sort(vec)]
     n = n_elements(vec)
     lo = vec[long(0.16*n)]
     hi = vec[long(0.84*n)]
     c_z0mgs = median(vec)
     
     print, tag
     print, c_z0mgs
     print, '+/- ', mad(vec), ' range ', lo, hi

     bins = bin_data(vec, vec*0.0+1.0 $
                     , xmin=xmin, xmax=xmax, binsize=binsize, /nan)
     normhist = ((bins.counts/total(bins.counts,/nan) > (0.)))
     
     ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
         , file=histpsfile
     
     loadct, 0
     plot $
        , [0], [0], /nodata $`
        , xtitle='!6log!d10!n C [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
        , ytitle='!6Fraction of Galaxies' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xrange=[xmin, xmax], yrange=[0., 0.25] $
        , /xstyle
     
     for kk = -100, 100 do $
        oplot, kk*binsize*10.*[1,1]+xmin, [-10, 10], lines=1, color=cgcolor('charcoal')
     
     for kk = -100, 100 do $
        oplot, [-1e6, 1e6], kk*0.05*[1,1], lines=1, color=cgcolor('charcoal')

     polyfill, [lo, hi, hi, lo, lo], [ylo, ylo, yhi, yhi, ylo], /fill, /clip $
               , color=cgcolor('lightgray')
     oplot, (c_z0mgs)*[1,1], [-100, 100], lines=0, thick=10, color=cgcolor('charcoal')   

     oplot, alog10(ref)*[1,1],  [-1d6, 1d6], lines=2, thick=10, color=cgcolor('darkgreen')
     
     histplot $
        , bins.xmid, normhist $
        , /overplot $
        , lthick=3, /fill $
        , fcolor=cgcolor('cornflowerblue') $
        , lcolor=cgcolor('royalblue')
     
     al_legend, /top, /left, lines=-99, charsize=1.25, charthick=3 $
                , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
                , textcolor=[cgcolor('black')] $
                , tag
     
;  al_legend, /top, /right, lines=-99, charsize=1.25, charthick=3 $
;             , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
;             , textcolor=[cgcolor('black')] $
;             , ['!6med: '+string(median(vec),format='(F5.3)'),'!6(!9+!6'+string(mad(vec),format='(F5.3)')+')']

     al_legend, /top, /right, lines=[2], charsize=1.25, charthick=3 $
                , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
                , textcolor=[cgcolor('darkgreen')] $
                , color=[cgcolor('darkgreen')] $
                , ['J13/KE12'], pspacing=1.0, thick=7

     loadct, 0
     plot $
        , [0], [0], /nodata $`
        , xtitle='!6log!d10!n C [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
        , ytitle='!6Fraction of Galaxies' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xrange=[xmin, xmax], yrange=[0., 0.25] $
        , /xstyle, /noerase     

     ps, /xw
     spawn, 'evince '+histpsfile+' &'
     spawn, 'convert -density 300x300 '+histpsfile+' '+histpnfile

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT THE SOLVED-FOR COEFFICIENTS VS STELLAR MASS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
         , file=psfile
     
     loadct, 0
     reversect
     y = alog10(coef[sfr_ind])
     grid = $
        grid_data(x, y, /nan $
                  , xaxis_out = x_axis, yaxis_out = y_axis $
                  , xmin=9., xmax=12., binsize_x=0.05 $
                  , ymin=xmin, ymax=xmax, binsize_y=0.025 $
                 )

     disp, alog10(grid), x_axis, y_axis $
           , ps=3 $
           , xrange=[9., 12.], yrange=[xmin, xmax] $
           , xthick=5, ythick=5, charsize=1.25, charthick=3 $
           , xtitle='!6log!d10!n Stellar Mass from CIGALE [M!d!9n!6!n]' $
           , ytitle='!6log!d10!n C [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
           , /xstyle, /ystyle, reserve=50, color=cgcolor('black', 255) $
           , max=3.0, position=[0.2,0.2,0.95, 0.95]

     for jj = -100, 100 do $
        oplot, jj*0.5*[1,1], [-1d6, 1d6], lines=1, color=cgcolor('black')
     for jj = -100, 100 do $
        oplot, [-1d6, 1d6], jj*0.25*[1,1]-43., lines=1, color=cgcolor('black')

     contour, /overplot, alog10(grid), x_axis, y_axis $
              , lev=[1., 1.5, 2.0, 2.5, 3.0], c_color=cgcolor('black')
     
     oplot, [-1d6, 1d6], (c_z0mgs)*[1,1], lines=0, thick=10, color=cgcolor('charcoal')   
     oplot, [-1d6, 1d6], alog10(ref)*[1,1], lines=2, thick=10, color=cgcolor('darkgreen')

     ;xfid = findgen(101)/100.*100.
     ;oplot, xfid, xfid*0.0, thick=13, color=cgcolor('charcoal')
     ;oplot, xfid, xfid*0.0, thick=5, color=cgcolor('gray')
     
     bins = bin_data(x, y, /nan, xmin=9.5, xmax=11.5, binsize=0.05)
     oploterror, bins.xmid $
                 , bins.ymed, bins.ymad $
                 , color=cgcolor('red') $
                 , psym=cgsymcat('filledsquare') $
                 , errthick=5, /nohat   
     
     al_legend, /top, /right $
                , box=1, clear=1, lines=-99, background=['lightgray'] $
                , charthick=3, charsize=1.25 $
                , tag
     
     ps, /xw
;     spawn, 'evince '+psfile+' &'
     spawn, 'convert -density 300x300 '+psfile+' '+pnfile
     
  endfor

  stop

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT PAIRS OF SFR MEASUREMENTS AGAINST EACH OTHER
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  loadct, 0

  for ii = 0, 6 do begin
     
     x = gws_logmstar[sfr_ind]

     if ii eq 0 then begin
        psfile = '../plots/gswlc_sfr_fuvw4_estimate.eps'
        pnfile = '../plots/gswlc_sfr_fuvw4_estimate.png'
        rat = alog10((sfr_fuvw4_z19/10.d^gws_logsfrsed)[sfr_ind])
        tag = '!6FUV+WISE4 / CIGALE'
     endif

     if ii eq 1 then begin
        psfile = '../plots/gswlc_sfr_nuvw3_estimate.eps'
        pnfile = '../plots/gswlc_sfr_nuvw3_estimate.png'
        rat = alog10((sfr_nuvw3_z19/10.d^gws_logsfrsed)[sfr_ind])
        tag = '!6NUV+WISE3 / CIGALE'
     endif

     if ii eq 2 then begin
        psfile = '../plots/gswlc_sfr_fuvw4nuvw3_estimate.eps'
        pnfile = '../plots/gswlc_sfr_fuvw4nuvw3_estimate.png'
        rat = alog10((sfr_nuvw3_z19/sfr_fuvw4_z19)[sfr_ind])
        tag = '!6NUV+WISE3 / FUV+WISE4'
     endif

     if ii eq 3 then begin
        psfile = '../plots/gswlc_sfr_w4fuvw4_estimate.eps'
        pnfile = '../plots/gswlc_sfr_w4fuvw4_estimate.png'
        rat = alog10((sfr_w4_z19/sfr_fuvw4_z19)[sfr_ind])
        tag = '!6Only WISE4 / FUV+WISE4'
     endif

     if ii eq 4 then begin
        psfile = '../plots/gswlc_sfr_w3nuvw3_estimate.eps'
        pnfile = '../plots/gswlc_sfr_w3nuvw3_estimate.png'
        rat = alog10((sfr_w3_z19/sfr_nuvw3_z19)[sfr_ind])
        tag = '!6Only WISE3 / NUV+WISE3'
     endif

     if ii eq 5 then begin
        psfile = '../plots/gswlc_sfr_w3w4_estimate.eps'
        pnfile = '../plots/gswlc_sfr_w3w4_estimate.png'
        rat = alog10((sfr_w3_z19/sfr_w4_z19)[sfr_ind])
        tag = '!6Only WISE3 / Only WISE4'
     endif

     if ii eq 1 then begin
        psfile = '../plots/gswlc_sfr_fuvw3_estimate.eps'
        pnfile = '../plots/gswlc_sfr_fuvw3_estimate.png'
        rat = alog10((sfr_fuvw3_z19/10.d^gws_logsfrsed)[sfr_ind])
        tag = '!6FUV+WISE3 / CIGALE'
     endif

     ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
         , file=psfile
     
     loadct, 0
     reversect
     y = rat
     grid = $
        grid_data(x, y, /nan $
                  , xaxis_out = x_axis, yaxis_out = y_axis $
                  , xmin=9., xmax=12., binsize_x=0.05 $
                  , ymin=-0.5, ymax=+0.5, binsize_y=0.01 $
                 )
     disp, alog10(grid), x_axis, y_axis $
           , ps=3 $
           , xrange=[9., 12.], yrange=[-0.5, 0.5] $
           , xthick=5, ythick=5, charsize=1.25, charthick=3 $
           , xtitle='!6log!d10!n Stellar Mass from CIGALE [M!d!9n!6!n]' $
           , ytitle='!6log!d10!n Residual SFR [M!d!9n!6!n yr!u-1!n]' $
           , /xstyle, /ystyle, reserve=50, color=cgcolor('black', 255) $
           , max=3.0, position=[0.2,0.2,0.95, 0.95]

     for jj = -100, 100 do $
        oplot, jj*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('black')
     for jj = -100, 100 do $
        oplot, [-100, 100], jj*0.1*[1,1], lines=1, color=cgcolor('black')

     contour, /overplot, alog10(grid), x_axis, y_axis $
              , lev=[1., 1.5, 2.0, 2.5, 3.0], c_color=cgcolor('black')
     
     xfid = findgen(101)/100.*100.
     oplot, xfid, xfid*0.0, thick=13, color=cgcolor('charcoal')
     oplot, xfid, xfid*0.0, thick=5, color=cgcolor('gray')
     
     bins = bin_data(x, rat, /nan, xmin=9.5, xmax=11.5, binsize=0.05)
     oploterror, bins.xmid $
                 , bins.ymed, bins.ymad $
                 , color=cgcolor('red') $
                 , psym=cgsymcat('filledsquare') $
                 , errthick=5, /nohat   
     
     al_legend, /top, /right $
                , box=1, clear=1, lines=-99, background=['lightgray'] $
                , charthick=3, charsize=1.25 $
                , tag
     
     ps, /xw
     spawn, 'evince '+psfile+' &'
     spawn, 'convert -density 300x300 '+psfile+' '+pnfile
     
  endfor

  stop

end
