pro plot_gswlc_resids

  @constants.bat
  lsun_3p4 = 1.83d18
  restore, '../gswlc/gswlc_data.idl'
    
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE SUBSET OF DATA TO USE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  mtol_ind = where(gws_logmstar gt 0 $ ; sanity check
                   and w1_lum gt 0 $ ; sanity check
                   and gws_w1 gt 3.*gws_ew1 $ ; detected
                   and gws_w3 gt 3.*gws_ew3 $ ; detected
                   and gws_nuv gt 3.*gws_enuv $ ; detected
                   and mtol_w1 gt 0.02 and mtol_w1 lt 1.0 $ ; reasonable m/l
                   and gws_flagsed eq 0 $ ; good fit
                  )

  nuvw1 = (alog10(nuv_lum/w1_lum))[mtol_ind]
  w3w1 = (alog10(w3_lum/w1_lum))[mtol_ind]
  lumw1 = (alog10(w1_lum/lsun_3p4))[mtol_ind]
  mtol = mtol_w1[mtol_ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT PREDICTION VS MEASUREMENT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  loadct, 0

  mtol_fix = $
     lumw1*0.0+0.45

  mtol_nuvw1 = $
     lookup_mtol(nuvw1=nuvw1 $
                 , unc=unc)

  mtol_w3w1 = $
     lookup_mtol(w3w1=w3w1 $
                 , unc=unc)

  mtol_justcolor = $
     lookup_mtol(nuvw1=nuvw1 $
                 , w3w1=w3w1 $
                 , unc=unc)

  mtol_justwise = $
     lookup_mtol(w3w1=w3w1 $
                 , lumw1=lumw1 $
                 , unc=unc)

  mtol_3d = $
     lookup_mtol(nuvw1=nuvw1 $
                 , w3w1=w3w1 $
                 , lumw1=lumw1 $
                 , unc=unc)

  for ii = 0, 5 do begin
  
     if ii eq 0 then begin
        psfile = '../plots/mtol_nuvw1_estimate.eps'
        pnfile = '../plots/mtol_nuvw1_estimate.png'
        rat = alog10(mtol_nuvw1/mtol)
        tag = '!6Just NUV/W1'
     endif

     if ii eq 1 then begin
        psfile = '../plots/mtol_w3w1_estimate.eps'
        pnfile = '../plots/mtol_w3w1_estimate.png'
        rat = alog10(mtol_w3w1/mtol)
        tag = '!6Just W3/W1'
     endif
  
     if ii eq 2 then begin
        psfile = '../plots/mtol_justcolor_estimate.eps'
        pnfile = '../plots/mtol_justcolor_estimate.png'
        rat = alog10(mtol_justcolor/mtol)
        tag = '!6NUV/W1 + W3/W1'
     endif

     if ii eq 3 then begin
        psfile = '../plots/mtol_justwise_estimate.eps'
        pnfile = '../plots/mtol_justwise_estimate.png'
        rat = alog10(mtol_justwise/mtol)
        tag = '!6W3/W1+W1'
     endif

     if ii eq 4 then begin
        psfile = '../plots/mtol_cube_estimate.eps'
        pnfile = '../plots/mtol_cube_estimate.png'
        rat = alog10(mtol_3d/mtol)
        tag = '!6NUV/W1 + W3/W1 + W1'
     endif

     if ii eq 5 then begin
        psfile = '../plots/mtol_fix_estimate.eps'
        pnfile = '../plots/mtol_fix_estimate.png'
        rat = alog10(mtol_fix/mtol)
        tag = '!6Fixed !7T!6!d*!n'
     endif

     ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
         , file=psfile
     
     loadct, 0
     reversect
     x = gws_logmstar[mtol_ind]
     y = rat
     grid = $
        grid_data(x, y, /nan $
                  , xaxis_out = x_axis, yaxis_out = y_axis $
                  , xmin=9., xmax=12., binsize_x=0.05 $
                  , ymin=-0.3, ymax=0.3, binsize_y=0.01 $
                 )
     disp, alog10(grid), x_axis, y_axis $
           , ps=3 $
           , xrange=[9., 12.], yrange=[-0.3, 0.3] $
           , xthick=5, ythick=5, charsize=1.25, charthick=3 $
           , xtitle='!6log!d10!n Stellar Mass from CIGALE [M!d!9n!6!n]' $
           , ytitle='!6log!d10!n Residual !7T!6!d*!n' $
           , /xstyle, /ystyle, reserve=50, color=cgcolor('black', 255) $
           , max=3.0

     for jj = -100, 100 do $
        oplot, jj*0.5*[1,1], [-10, 10], lines=1, color=cgcolor('black')
     for jj = -100, 100 do $
        oplot, [-100, 100], jj*0.1*[1,1], lines=1, color=cgcolor('black')

     contour, /overplot, alog10(grid), x_axis, y_axis $
              , lev=[1., 1.5, 2.0, 2.5, 3.0], c_color=cgcolor('black')
     
     xfid = findgen(101)/100.*100.
     oplot, xfid, xfid*0.0, thick=13, color=cgcolor('charcoal')
     oplot, xfid, xfid*0.0, thick=5, color=cgcolor('gray')
     
     bins = bin_data(gws_logmstar[mtol_ind] $
                     , rat, /nan, xmin=9.5, xmax=11.5, binsize=0.05)
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
