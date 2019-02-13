pro plot_gswlc_grid

  @constants.bat
  lsun_3p4 = 1.83d18
  thresh = 10

  plot, findgen(10), xtitle='!6'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PREPARE THE GWSLC MEASUREMENTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  gws = mrdfits('~/idl/galbase/gal_data/'+$
                'hlsp_gswlc_galex-sdss-wise_multi_x1_multi_v1_cat.fits' $
                , 1, h_gws)
  readcol, '/data/kant/0/leroy.42/allsky/gswlc/galex_unwise_fluxes_GSWLC-X.dat' $
           , id, gws_fuv, gws_efuv, gws_nuv, gws_enuv, gws_w1_nm, gws_ew1 $
           , gws_w2_nm, gws_ew2, gws_w3_nm, gws_ew3, gws_w4_nm, gws_ew4 $
           , format='L,F,F,F,F,F,F,F,F,F,F,F,F'

  gws_w1 = 3631.*10^(-0.4*(22.5+2.683))*gws_w1_nm
  gws_w2 = 3631*10^(-0.4*(22.5+3.319))*gws_w2_nm
  gws_w3 = 3631*10^(-0.4*(22.5+5.242))*gws_w3_nm
  gws_w4 = 3631*10^(-0.4*(22.5+6.604))*gws_w4_nm

  nu_fuv = c/(154.d-9*1d2)
  nu_nuv = c/(231.d-9*1d2)
  nu_w1 = c/(3.4d-6*1d2)
  nu_w2 = c/(4.5d-6*1d2)
  nu_w3 = c/(12.d-6*1d2)
  nu_w4 = c/(22.d-6*1d2)

  fuv_lum = gws_fuv/1d3*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2
  nuv_lum = gws_nuv/1d3*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2
  w1_lum = gws_w1*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2
  w2_lum = gws_w2*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2
  w3_lum = gws_w3*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2
  w4_lum = gws_w4*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2

  sfr_fuv_ke12 = $
     lum_to_sfr(band='FUV', cal='KE12', lum=fuv_lum*nu_fuv)
  sfr_nuv_ke12 = $
     lum_to_sfr(band='NUV', cal='KE12', lum=nuv_lum*nu_nuv)
  sfr_w3_j13 = $
     lum_to_sfr(band='WISE3', cal='J13', lum=w3_lum*nu_w3)
  sfr_w4_j13 = $
     lum_to_sfr(band='WISE4', cal='J13', lum=w4_lum*nu_w4)

  sfr_fuvw4_ke12 = $
     lum_to_sfr(band='FUV', cal='KE12' $
                , lum=(fuv_lum*nu_fuv + 3.89*w4_lum*nu_w4))  
  sfr_nuvw3 = $
     sfr_nuv_ke12+sfr_w3_j13

  mtol_w1 = 10.^(gws.logmstar - alog10(w1_lum/lsun_3p4))
  gws_ssfr = gws.logsfrsed-gws.logmstar
  ssfr_like = sfr_nuvw3 / (w1_lum / lsun_3p4)
  w2w1 = alog10(w2_lum/w1_lum)
  w3w1 = alog10(w3_lum/w1_lum)
  nuvw1 = alog10(nuv_lum/w1_lum)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE SUBSET OF DATA TO USE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sane_ind = where(gws.logmstar gt 0 $
                   and w1_lum gt 0 $
                   and gws_w1_nm gt 3.*gws_ew1 $
                   and gws_w3_nm gt 3.*gws_ew3 $
                   and gws_nuv gt 3.*gws_enuv $
                   and mtol_w1 gt 0.02 and mtol_w1 lt 1.0 $
                  )

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT PREDICTION VS MEASUREMENT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  nuvw1 = (alog10(nuv_lum/w1_lum))[sane_ind]
  w3w1 = (alog10(w3_lum/w1_lum))[sane_ind]
  lumw1 = (alog10(w1_lum/lsun_3p4))[sane_ind]
  mtol = mtol_w1[sane_ind]

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
     x = gws[sane_ind].logmstar
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

     contour, /overplot, alog10(grid), x_axis, y_axis $
              , lev=[1., 1.5, 2.0, 2.5, 3.0], c_color=cgcolor('black')
     
     xfid = findgen(101)/100.*100.
     oplot, xfid, xfid*0.0, thick=13, color=cgcolor('charcoal')
     oplot, xfid, xfid*0.0, thick=5, color=cgcolor('gray')
     
     bins = bin_data(gws[sane_ind].logmstar $
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
