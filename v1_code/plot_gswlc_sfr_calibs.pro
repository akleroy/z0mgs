pro plot_gswlc_sfr_calibs $
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
; LOAD THE COEFFICIENT GRIDS IN MSTAR-SFR/MSTAR SPACE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  counts = readfits('../measurements/ssfr_mstar_count_grid.fits')

  grid_justw3 = readfits('../measurements/ssfr_mstar_cjustw3_grid.fits')
  grid_justw4 = readfits('../measurements/ssfr_mstar_cjustw4_grid.fits')

  grid_w3nuv = readfits('../measurements/ssfr_mstar_cnuvw3_grid.fits')
  grid_w3fuv = readfits('../measurements/ssfr_mstar_cfuvw3_grid.fits')

  grid_w4nuv = readfits('../measurements/ssfr_mstar_cnuvw4_grid.fits')
  grid_w4fuv = readfits('../measurements/ssfr_mstar_cfuvw4_grid.fits')

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CALCULATE AND PRINT THE MEDIAN C AND THE SCATTER BY GAL AND BY CELL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  loadct, 0

  print, 'Median Coefficients (all log10 C):'

  psfile = '../plots/gswlc_cwise.eps'
  pnfile = '../plots/gswlc_cwise.png'

  ps, /def, /ps, xs=7, ys=4, /color, /encaps $
      , file=psfile
  
  loadct, 0

  plot $
     , [0], [0] $
     , xtitle='!6Tracer', xrange=[-1,6], xstyle=5 $
     , yrange=[-43.5, -41.75], ystyle=1 $
     , ytitle='!6log!d10!n C [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
     , xthick=5, ythick=5, charthick=3, charsize=1.25 $
     , /nodata

  for ii = -200, 0 do $
     oplot, [-100,100], ii*0.25*[1,1], lines=1
  for ii = 0, 10 do $
     oplot, ii*[1,1], [-100,100], lines=1
  oplot, [-100, 100], -43.5*[1,1], lines=0, thick=5
  oplot, [-100, 100], -41.75*[1,1], lines=0, thick=5

  for ii = 0, 5 do begin

     if ii eq 0 then begin
        tag = "WISE4+FUV"
        rat = alog10(coef_w4fuv)
        grid = grid_w4fuv
     endif

     if ii eq 1 then begin
        tag = "WISE4+NUV"
        rat = alog10(coef_w4nuv)
        grid = grid_w4nuv
     endif

     if ii eq 2 then begin
        tag = "Just WISE4"
        rat = alog10(coef_justw4)
        grid = grid_justw4
     endif

     if ii eq 3 then begin
        tag = "WISE3+FUV"
        rat = alog10(coef_w3fuv)
        grid = grid_w3fuv
     endif

     if ii eq 4 then begin
        tag = "WISE3+NUV"
        rat = alog10(coef_w3nuv)
        grid = grid_w3nuv
     endif

     if ii eq 5 then begin
        tag = "Just WISE3"
        rat = alog10(coef_justw3)
        grid = grid_justw3
     endif
     print, tag
     rat = rat[ind]
     grid = alog10(grid)

     medrat = median(rat)
     madrat = mad(rat)
     medgrid = median(grid)
     madgrid = mad(grid)

     print, "... median and scatter BY GALAXY: ", medrat, " +/- ", madrat
     print, "... median and scatter IN GRID: ", medgrid, " +/- ", madgrid

     oploterror, (ii-0.1)*[1], medrat*[1], [madrat] $
                 , psym=cgsymcat('filledcircle'), symsize=2, errthick=10 $
                 , color=cgcolor('firebrick')
     oploterror, (ii+0.1)*[1], medgrid*[1], [madgrid] $
                 , psym=cgsymcat('filledsquare'), symsize=2, errthick=10 $
                 , color=cgcolor('royalblue')
     xyouts, ii*1.0, -42.05-(ii mod 2)*0.10 $
             , align=0.5, tag, charthick=3, charsize=1.25

  endfor

  al_legend, /bottom, /left $
             , box=1, clear=1, lines=-99, background=['lightgray'] $
             ;, box=0, clear=0, lines=-99 $
             , charthick=3, charsize=1.25 $
             , ['!6Individual galaxies','!6Grid cells in sSFR-M!d*!n'] $
             , textcolor=[cgcolor('firebrick'), cgcolor('royalblue')]
   
  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile

  stop

end
