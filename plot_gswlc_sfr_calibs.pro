pro plot_gswlc_sfr_calibs $
   , show=show

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

  all_ind = where(gws_logmstar gt 0 $
                  and w1_lum gt 0 $
                  and gws_w1 gt 3.*gws_ew1 $
                  and gws_w3 gt 3.*gws_ew3 $
                  and gws_w4 gt 3.*gws_ew4 $
                  and gws_nuv gt 3.*gws_enuv $
                  and gws_fuv gt 3.*gws_efuv $
                  and mtol_w1 gt 0.02 and mtol_w1 lt 1.0 $
                 )

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
     
     if ii eq 0 then begin
        mstar_psfile = '../plots/gswlc_justw3_coef.eps'
        mstar_pnfile = '../plots/gswlc_justw3_coef.png'

        ssfr_psfile = '../plots/gswlc_justw3_ssfr.eps'
        ssfr_pnfile = '../plots/gswlc_justw3_ssfr.png'

        histpsfile = '../plots/gswlc_justw3_hist.eps'
        histpnfile = '../plots/gswlc_justw3_hist.png'

        coef = coef_justw3
        tag = ['!6WISE3 only']

        ref = 10.d^(-42.9)

        xmin = -43.5
        xmax = -41.5
     endif

     if ii eq 1 then begin
        mstar_psfile = '../plots/gswlc_justw4_coef.eps'
        mstar_pnfile = '../plots/gswlc_justw4_coef.png'

        ssfr_psfile = '../plots/gswlc_justw4_ssfr.eps'
        ssfr_pnfile = '../plots/gswlc_justw4_ssfr.png'

        histpsfile = '../plots/gswlc_justw4_hist.eps'
        histpnfile = '../plots/gswlc_justw4_hist.png'

        coef = coef_justw4
        tag = ['!6WISE4 only']

        ref = 10.d^(-42.7)*22./24. ; MIPS -> WISE

        xmin = -43.5
        xmax = -41.5
     endif

     if ii eq 2 then begin
        mstar_psfile = '../plots/gswlc_nuvw3_coef.eps'
        mstar_pnfile = '../plots/gswlc_nuvw3_coef.png'

        ssfr_psfile = '../plots/gswlc_nuvw3_ssfr.eps'
        ssfr_pnfile = '../plots/gswlc_nuvw3_ssfr.png'

        histpsfile = '../plots/gswlc_nuvw3_hist.eps'
        histpnfile = '../plots/gswlc_nuvw3_hist.png'

        coef = coef_w3nuv
        tag = ['!6WISE3','in NUV+WISE3']

        ref = 10.d^(-42.9)

        xmin = -43.5
        xmax = -41.5
     endif

     if ii eq 3 then begin
        mstar_psfile = '../plots/gswlc_nuvw4_coef.eps'
        mstar_pnfile = '../plots/gswlc_nuvw4_coef.png'

        ssfr_psfile = '../plots/gswlc_nuvw4_ssfr.eps'
        ssfr_pnfile = '../plots/gswlc_nuvw4_ssfr.png'

        histpsfile = '../plots/gswlc_nuvw4_hist.eps'
        histpnfile = '../plots/gswlc_nuvw4_hist.png'

        coef = coef_w4nuv
        tag = ['!6WISE4','in NUV+WISE4']

        ref = 10.d^(-42.81)*22./25. ; KENNICUTT/HAO -> WISE4 from IRAS2

        xmin = -43.5
        xmax = -41.5
     endif

     if ii eq 4 then begin
        mstar_psfile = '../plots/gswlc_fuvw3_coef.eps'
        mstar_pnfile = '../plots/gswlc_fuvw3_coef.png'

        ssfr_psfile = '../plots/gswlc_fuvw3_ssfr.eps'
        ssfr_pnfile = '../plots/gswlc_fuvw3_ssfr.png'

        histpsfile = '../plots/gswlc_fuvw3_hist.eps'
        histpnfile = '../plots/gswlc_fuvw3_hist.png'

        coef = coef_w3fuv
        tag = ['!6WISE3','in FUV+WISE3']

        ref = 10.d^(-42.9)

        xmin = -43.5
        xmax = -41.5
     endif

     if ii eq 5 then begin
        mstar_psfile = '../plots/gswlc_fuvw4_coef.eps'
        mstar_pnfile = '../plots/gswlc_fuvw4_coef.png'

        ssfr_psfile = '../plots/gswlc_fuvw4_ssfr.eps'
        ssfr_pnfile = '../plots/gswlc_fuvw4_ssfr.png'

        histpsfile = '../plots/gswlc_fuvw4_hist.eps'
        histpnfile = '../plots/gswlc_fuvw4_hist.png'

        coef = coef_w4fuv
        tag = ['!6WISE4','in FUV+WISE4']

        ref = 10.d^(-42.76)*22./25. ; KENNICUTT IRAS -> WISE

        xmin = -43.5
        xmax = -41.5
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT HISTOGRAMS OF THE SOLVED-FOR COEFFICIENTS
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
     
     ylo = 0.0
     yhi = 0.25

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
     if keyword_set(show) then $
        spawn, 'evince '+histpsfile+' &'
     spawn, 'convert -density 300x300 '+histpsfile+' '+histpnfile

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT THE SOLVED-FOR COEFFICIENTS VS STELLAR MASS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
         , file=mstar_psfile
     
     loadct, 0
     reversect
     x = gws_logmstar[all_ind] ; [sfr_ind]
     y = alog10(coef[all_ind]) ; [sfr_ind]
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
     if keyword_set(show) then $
        spawn, 'evince '+mstar_psfile+' &'
     spawn, 'convert -density 300x300 '+mstar_psfile+' '+mstar_pnfile

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT THE SOLVED-FOR COEFFICIENT VS SPECIFIC STAR FORMATION RATE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
         , file=ssfr_psfile
     
     loadct, 0
     reversect
     x = gws_ssfr[all_ind] ; [sfr_ind]
     y = alog10(coef[all_ind]) ; [sfr_ind]
     grid = $
        grid_data(x, y, /nan $
                  , xaxis_out = x_axis, yaxis_out = y_axis $
                  , xmin=-12., xmax=-8., binsize_x=0.05 $
                  , ymin=xmin, ymax=xmax, binsize_y=0.025 $
                 )

     disp, alog10(grid), x_axis, y_axis $
           , ps=3 $
           , xrange=[-12., -8.], yrange=[xmin, xmax] $
           , xthick=5, ythick=5, charsize=1.25, charthick=3 $
           , xtitle='!6log!d10!n SFR/M!d*!n from CIGALE [M!d!9n!6!n]' $
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
     
     bins = bin_data(x, y, /nan, xmin=-11.25, xmax=-9.0, binsize=0.05)
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
     if keyword_set(show) then $
        spawn, 'evince '+ssfr_psfile+' &'
     spawn, 'convert -density 300x300 '+ssfr_psfile+' '+ssfr_pnfile
     
  endfor

  stop

end
