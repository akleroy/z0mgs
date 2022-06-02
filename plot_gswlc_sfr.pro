pro plot_gswlc_sfr $
  , show=show

  thresh = 10
  ston = 3.0
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

  ind = where(gws_logmstar gt 0 $
              and w1_lum gt 0 $
              and gws_w1 gt ston*gws_ew1 $
              and gws_w3 gt ston*gws_ew3 $
              and gws_w4 gt ston*gws_ew4 $
              and gws_nuv gt ston*gws_enuv $
              and gws_fuv gt ston*gws_efuv $
              and gws_flagsed eq 0 $
             )

  fuvw1 = (alog10(fuv_lum/w1_lum))[ind]
  nuvw1 = (alog10(fuv_lum/w1_lum))[ind]
  w4w1 = (alog10(w4_lum/w1_lum))[ind]
  w3w1 = (alog10(w4_lum/w1_lum))[ind]

  mstar = ((gws_logmstar))[ind]
  ssfr = ((gws_ssfr))[ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ESTIMATE THE WISE COEFFICIENTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

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

  y = alog10(coef_w4fuv[ind])
  ymin = -43.5
  ymax = -41.75
  ytitle = '!6log!d10!n C [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]'

  fid = -42.68

  y_justw4 = alog10(coef_justw4[ind])

  for ii = 0, 3 do begin
     
     if ii eq 0 then begin
        x = (ssfr)
        xtitle = '!6log!d10!n SFR/M!d*!n from CIGALE [yr!u-1!n]'
        tag = 'ssfr'
        xmin = -12.0
        xmax = -9.0
        binmin = -11.75
        binmax = -8.5
        binsize = 0.1
     endif

     if ii eq 1 then begin
        x = alog10(ssfr_like_fuvw4[ind])
        xtitle = '!6log!d10!n SFR(FUV+W4)/W1 Lum. [M!d!9n!6!n yr!u-1!n L!d!9n!6!n!u-1!n]'
        tag = 'ssfrlike'
        xmin = -12.5
        xmax = -9.5
        binmin = -12.5
        binmax = -10.0
        binsize = 0.1
     endif

     if ii eq 2 then begin
        x = w4w1
        xtitle = '!6log!d10!n W4-to-W1'
        tag = 'w4w1'
        xmin = -1.5
        xmax = 1.25
        binmin = -1.25
        binmax = 1.0
        binsize = 0.1
     endif

     if ii eq 3 then begin
        x = nuvw1
        xtitle = '!6log!d10!n NUV-to-W1'
        tag = 'nuvw1'
        xmin = -3.5
        xmax = 0.0
        binmin = -3.0
        binmax = -0.25
        binsize = 0.1
     endif

     psfile = '../plots/gswlc_cfuvw4_'+tag+'.eps'
     pnfile = '../plots/gswlc_cfuvw4_'+tag+'.png'

     txtfile = '../measurements/gswlc_cfuvw4_'+tag+'.txt'
     txtfile_justw4 = '../measurements/gswlc_justw4_'+tag+'.txt'

     ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
         , file=psfile
     
     loadct, 0
     reversect

     cfuvw4_grid = $
        grid_data(x, y, /nan $
                  , xaxis_out = xaxis_grid, yaxis_out = yaxis_grid $
                  , ymin=ymin, ymax=ymax, binsize_y=0.025 $
                  , xmin=xmin, xmax=xmax, binsize_x=(xmax-xmin)/50.0)

     disp, alog10(cfuvw4_grid), x_axis, y_axis $
           , ps=3 $
           , xrange=[xmin, xmax], yrange=[ymin, ymax] $
           , xthick=5, ythick=5, charsize=1.25, charthick=3 $
           , xtitle=xtitle $
           , ytitle=ytitle $
           , /xstyle, /ystyle, reserve=50, color=cgcolor('black', 255) $
           , max=3.0, position=[0.2,0.2,0.95, 0.95]

     for jj = -100, 100 do $
        oplot, jj*0.5*[1,1], [-1d6, 1d6], lines=1, color=cgcolor('black')
     for jj = -100, 100 do $
        oplot, [-1d6, 1d6], jj*0.25*[1,1]-43., lines=1, color=cgcolor('black')

     oplot, [-1d6, 1d6], fid*[1,1], thick=10, color=cgcolor('dodgerblue')
     oplot, [-1d6, 1d6], fid*[1,1], thick=3, color=cgcolor('royalblue')

     contour, /overplot, alog10(cfuvw4_grid), xaxis_grid, yaxis_grid $
              , lev=[1., 1.5, 2.0, 2.5, 3.0], c_color=cgcolor('black')
     
     bins = $
        bin_data(x, y $
                 , xmin=binmin, xmax=binmax, binsize=binsize, /nan)

     bins_justw4 = $
        bin_data(x, y_justw4 $
                 , xmin=binmin, xmax=binmax, binsize=binsize, /nan)

     oploterror, bins.xmid $
                 , bins.ymed, bins.ymad $
                 , color=cgcolor('red') $
                 , psym=cgsymcat('filledsquare') $
                 , errthick=5, /nohat            

     al_legend, /top, /left $
                , box=1, clear=1, lines=-99, background=['lightgray'] $
                , charthick=3, charsize=1.25 $
                , ['!6C in WISE4+FUV']
     
     ps, /xw
     if keyword_set(show) then begin
        spawn, 'evince '+psfile+' &'
     endif
     spawn, 'convert -density 300x300 '+psfile+' '+pnfile
     
     get_lun, lun
     openw, lun, txtfile
     printf, lun, '# coeff on w4 in FUV+W4 SFR'
     printf, lun, '# binned by '+tag
     printf, lun, '# format x, coeff, scatter'
     printf, lun, 'x coeff scatter'
     for yy = 0, n_elements(bins)-1 do $
        printf, lun, bins[yy].xmid, bins[yy].ymed, bins[yy].ymad
     close, lun

     get_lun, lun
     openw, lun, txtfile_justw4
     printf, lun, '# coeff on w4 in just W4 SFR'
     printf, lun, '# binned by '+tag
     printf, lun, '# format x, coeff, scatter'
     printf, lun, 'x coeff scatter'
     for yy = 0, n_elements(bins)-1 do $
        printf, lun, bins_justw4[yy].xmid, bins_justw4[yy].ymed, bins[yy].ymad
     close, lun
     
  endfor

  stop

end
