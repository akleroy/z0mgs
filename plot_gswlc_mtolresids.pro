pro plot_gswlc_mtolresids

  ston = 3.0
  @constants.bat
  lsun_3p4 = 1.83d18
  restore, '../gswlc/gswlc_data.idl'
    
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE SUBSET OF DATA TO USE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  mtol_ind = where(gws_logmstar gt 0 $
                   and w1_lum gt 0 $
                   and gws_w1 gt ston*gws_ew1 $
                   and gws_w3 gt ston*gws_ew3 $
                   and gws_w4 gt ston*gws_ew4 $
                   and gws_nuv gt ston*gws_enuv $
                   and gws_fuv gt ston*gws_efuv $
                   and gws_flagsed eq 0 $
                  )

  mass_bins_min = [9., 9.5, 10., 10.5, 11.]
  mass_bins_max = [9.5, 10., 10.5, 11., 11.5]
  n_mass_bins = n_elements(mass_bins_min)
  mass_bins_color = reverse(['red','salmon','goldenrod','lightseagreen','dodgerblue'])

  ssfrlike_fuvw4 = alog10(ssfr_like_fuvw4[mtol_ind])
  ssfrlike_nuvw4 = alog10(ssfr_like_nuvw4[mtol_ind])
  ssfr = gws_ssfr[mtol_ind]
  nuvw1 = nuvw1[mtol_ind]
  fuvw1 = fuvw1[mtol_ind]
  w3w1 = w3w1[mtol_ind]
  w4w1 = w4w1[mtol_ind]
  mtol_w1 = mtol_w1[mtol_ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PREDICT MTOL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  loadct, 0

  mtol_fix = $
     lookup_mtol()

  mtol_ssfr = $
     lookup_mtol(truessfr=ssfr)

  mtol_ssfrlike_fuvw4 = $
     lookup_mtol(ssfrlike=ssfrlike_fuvw4)

  mtol_ssfrlike_nuvw4 = $
     lookup_mtol(ssfrlike=ssfrlike_nuvw4)

  mtol_nuvw3w1 = $
     lookup_mtol(nuvw1=nuvw1,w3w1=w3w1)

  mtol_fuvw4w1 = $
     lookup_mtol(fuvw1=fuvw1,w4w1=w4w1)
  
  mtol_w4w1 = $
     lookup_mtol(w4w1=w4w1)

  mtol_w3w1 = $
     lookup_mtol(w3w1=w3w1)

  mtol_fuvw1 = $
     lookup_mtol(fuvw1=fuvw1)

  mtol_nuvw1 = $
     lookup_mtol(nuvw1=nuvw1)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COMPARE TO GSWLC CALCULATION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  x = ssfr

  xtitle = '!6log!d10!n SFR/M!d*!n from CIGALE [yr!u-1!n]'
  tag = 'ssfr'
  xmin = -12.0
  xmax = -9.0
  binmin = -11.75
  binmax = -8.5
  binsize = 0.1

  for ii = 0, 10 do begin

     if ii eq 0 then begin
        psfile = '../plots/mtol_fix_estimate.eps'
        pnfile = '../plots/mtol_fix_estimate.png'
        rat = alog10(mtol_fix/mtol_w1)
        tag = '!6Fixed !7T!6!d*!n'
     endif
  
     if ii eq 1 then begin
        psfile = '../plots/mtol_ssfr_estimate.eps'
        pnfile = '../plots/mtol_ssfr_estimate.png'
        rat = alog10(mtol_ssfr/mtol_w1)
        tag = '!6SFR/M!d*!n'
     endif

     if ii eq 2 then begin
        psfile = '../plots/mtol_ssfrlikefuvw4_estimate.eps'
        pnfile = '../plots/mtol_ssfrlikefuvw4_estimate.png'
        rat = alog10(mtol_ssfrlike_fuvw4/mtol_w1)
        tag = '!6SFR(FUV+W4)/W1'
     endif

     if ii eq 3 then begin
        psfile = '../plots/mtol_ssfrlikenuvw4_estimate.eps'
        pnfile = '../plots/mtol_ssfrlikenuvw4_estimate.png'
        rat = alog10(mtol_ssfrlike_nuvw4/mtol_w1)
        tag = '!6SFR(NUV+W4)/W1'
     endif

     if ii eq 4 then begin
        print, "GRIDS TURNED OFF AS TOO UNSTABLE."
        continue
        psfile = '../plots/mtol_nuvw3w1_estimate.eps'
        pnfile = '../plots/mtol_nuvw3w1_estimate.png'
        rat = alog10(mtol_nuvw3w1/mtol_w1)
        tag = '!6NUV+W3+W1'
     endif

     if ii eq 5 then begin
        print, "GRIDS TURNED OFF AS TOO UNSTABLE."
        continue
        psfile = '../plots/mtol_fuvw4w1_estimate.eps'
        pnfile = '../plots/mtol_fuvw4w1_estimate.png'
        rat = alog10(mtol_fuvw4w1/mtol_w1)
        tag = '!6FUV+W4+W1'
     endif

     if ii eq 6 then begin
        psfile = '../plots/mtol_w4w1_estimate.eps'
        pnfile = '../plots/mtol_w4w1_estimate.png'
        rat = alog10(mtol_w4w1/mtol_w1)
        tag = '!6WISE4/WISE1'
     endif

     if ii eq 7 then begin
        psfile = '../plots/mtol_w3w1_estimate.eps'
        pnfile = '../plots/mtol_w3w1_estimate.png'
        rat = alog10(mtol_w3w1/mtol_w1)
        tag = '!6WISE3/WISE1'
     endif

     if ii eq 8 then begin
        psfile = '../plots/mtol_nuvw1_estimate.eps'
        pnfile = '../plots/mtol_nuvw1_estimate.png'
        rat = alog10(mtol_nuvw1/mtol_w1)
        tag = '!6NUV/WISE1'
     endif

     if ii eq 9 then begin
        psfile = '../plots/mtol_fuvw1_estimate.eps'
        pnfile = '../plots/mtol_fuvw1_estimate.png'
        rat = alog10(mtol_fuvw1/mtol_w1)
        tag = '!6FUV/WISE1'
     endif

     ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
         , file=psfile
     
     loadct, 0
     reversect
     y = rat

     grid = $
        grid_data(x, y, /nan $
                  , xaxis_out = x_axis, yaxis_out = y_axis $
                  , xmin=xmin, xmax=xmax, binsize_x=(xmax-xmin)/100.0 $
                  , ymin=-0.5, ymax=0.5, binsize_y=0.02 $
                 )

     disp, alog10(grid), x_axis, y_axis $
           , ps=3 $
           , xrange=[xmin, xmax], yrange=[-0.5, 0.5] $
           , xthick=5, ythick=5, charsize=1.25, charthick=3 $
           , xtitle=xtitle $
           , ytitle='!6log!d10!n Residual !7T!6!d*!n' $
           , /xstyle, /ystyle, reserve=50, color=cgcolor('black', 255) $
           , max=3.0

     for jj = -100, 100 do $
        oplot, jj*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('black')
     for jj = -100, 100 do $
        oplot, [-100, 100], jj*0.1*[1,1], lines=1, color=cgcolor('black')

     contour, /overplot, alog10(grid), x_axis, y_axis $
              , lev=[1., 1.5, 2.0, 2.5, 3.0], c_color=cgcolor('black')
     
     xfid = findgen(101)/100.*100.-50.
     oplot, xfid, xfid*0.0, thick=13, color=cgcolor('salmon')
     oplot, xfid, xfid*0.0, thick=5, color=cgcolor('rose')
     
     bins = bin_data(x, y, /nan, xmin=binmin, xmax=binmax, binsize=binsize)
     oploterror, bins.xmid $
                 , bins.ymed, bins.ymad $
                 , color=cgcolor('red') $
                 , psym=cgsymcat('filledsquare') $
                 , errthick=5, /nohat   
     
     al_legend, /top, /left $
                , box=1, clear=1, lines=-99, background=['lightgray'] $
                , charthick=3, charsize=1.25 $
                , tag
     
     ps, /xw
     spawn, 'evince '+psfile+' &'
     spawn, 'convert -density 300x300 '+psfile+' '+pnfile
     
  endfor

  stop

end
