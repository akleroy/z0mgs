pro build_gswlc_grid

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
; DEFINE SUBSET OF DATA TO USE FOR MASS TO LIGHT RATIOS
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
; CALCULATE THE M/L AND WISE-TO-SFR COEFFICIENTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  mtol = mtol_w1[ind]

  coef_fuv = (((10.d^(gws_logsfrsed))) / (nu_fuv*fuv_lum*10.^(gws_afuv/2.5)))[ind]
  coef_justfuv = (((10.d^(gws_logsfrsed))) / (nu_fuv*fuv_lum))[ind]
  coef_justnuv = (((10.d^(gws_logsfrsed))) / (nu_nuv*nuv_lum))[ind]

  coef_justw3 = (((10.d^(gws_logsfrsed))) / (nu_w3*w3_lum))[ind]
  coef_justw4 = (((10.d^(gws_logsfrsed))) / (nu_w4*w4_lum))[ind]
  
  coef_w3nuv = (((10.d^(gws_logsfrsed)) - sfr_nuv_z19) /  (nu_w3*w3_lum))[ind]
  coef_w4nuv = (((10.d^(gws_logsfrsed)) - sfr_nuv_z19) /  (nu_w4*w4_lum))[ind]

  coef_w3fuv = (((10.d^(gws_logsfrsed)) - sfr_fuv_z19) /  (nu_w3*w3_lum))[ind]
  coef_w4fuv = (((10.d^(gws_logsfrsed)) - sfr_fuv_z19) /  (nu_w4*w4_lum))[ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE AXES AND BUILD GRIDS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  min_fuvw1 = -3.5
  max_fuvw1 = 0.5
  bin_fuvw1 = 0.05
  tol_fuvw1 = 0.05;0.1

  min_nuvw1 = -3.5
  max_nuvw1 = 0.5
  bin_nuvw1 = 0.05
  tol_nuvw1 = 0.05;0.1

  min_w4w1 = -2.0
  max_w4w1 = 2.0
  bin_w4w1 = 0.05
  tol_w4w1 = 0.05;0.1

  min_w3w1 = -2.0
  max_w3w1 = 2.0
  bin_w3w1 = 0.05
  tol_w3w1 = 0.05;0.1

  min_ssfr = -12.
  max_ssfr = -8.
  bin_ssfr = 0.125
  tol_ssfr = 0.125 ;0.25

  min_mstar = 9.0
  max_mstar = 11.0
  bin_mstar = 0.125
  tol_mstar = 0.125 ;0.25  

  count_fuvw1_w4w1 = $
     grid_data(fuvw1, w4w1, /nan $
               , xaxis_out = fuvw1_axis, yaxis_out = w4w1_axis $
               , xmin=min_fuvw1, xmax=max_fuvw1, binsize_x=binsize_fuvw1 $
               , ymin=min_w4w1, ymax=max_w4w1, binsize_y=binsize_w4w1 $
              )
  n_x = n_elements(fuvw1_axis)
  n_y = n_elements(w4w1_axis)

  count_nuvw1_w3w1 = $
     grid_data(nuvw1, w3w1, /nan $
               , xaxis_out = nuvw1_axis, yaxis_out = w3w1_axis $
               , xmin=min_nuvw1, xmax=max_nuvw1, binsize_x=binsize_nuvw1 $
               , ymin=min_w3w1, ymax=max_w3w1, binsize_y=binsize_w3w1 $
              )
  n_u = n_elements(nuvw1_axis)
  n_v = n_elements(w3w1_axis)

  count_mstar_ssfr = $
     grid_data(mstar, ssfr, /nan $
               , xaxis_out = mstar_axis, yaxis_out = ssfr_axis $
               , xmin=min_mstar, xmax=max_mstar, binsize_x=binsize_mstar $
               , ymin=min_ssfr, ymax=max_ssfr, binsize_y=binsize_ssfr $
              )
  n_p = n_elements(mstar_axis)
  n_q = n_elements(ssfr_axis)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE GRIDS FOR FUV/W1+W4/W1
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  count_xy_grid = fltarr(n_x, n_y)*!values.f_nan

  mtol_xy_grid = fltarr(n_x, n_y)*!values.f_nan
  mtol_xy_mad_grid = mtol_xy_grid*!values.f_nan

  cfuvw4_xy_grid = fltarr(n_x, n_y)*!values.f_nan
  cfuvw4_xy_mad_grid = cfuvw4_xy_grid*!values.f_nan

  cnuvw4_xy_grid = fltarr(n_x, n_y)*!values.f_nan
  cnuvw4_xy_mad_grid = cnuvw4_xy_grid*!values.f_nan

  cfuvw3_xy_grid = fltarr(n_x, n_y)*!values.f_nan
  cfuvw3_xy_mad_grid = cfuvw3_xy_grid*!values.f_nan

  cnuvw3_xy_grid = fltarr(n_x, n_y)*!values.f_nan
  cnuvw3_xy_mad_grid = cnuvw3_xy_grid*!values.f_nan

  cjustw4_xy_grid = fltarr(n_x, n_y)*!values.f_nan
  cjustw4_xy_mad_grid = cjustw4_xy_grid*!values.f_nan

  cjustw3_xy_grid = fltarr(n_x, n_y)*!values.f_nan
  cjustw3_xy_mad_grid = cjustw3_xy_grid*!values.f_nan

  for ii = 0, n_x-1 do begin
     for jj = 0, n_y-1 do begin

        use_this = abs(fuvw1 - fuvw1_axis[ii]) lt tol_fuvw1
        use_this *= abs(w4w1 - w4w1_axis[jj]) lt tol_w4w1

        this_ind = where(use_this, this_ct)

        if this_ct gt thresh then begin

           count_xy_grid[ii,jj] = this_ct*1.0

           mtol_xy_grid[ii,jj] = median(mtol[this_ind])
           mtol_xy_mad_grid[ii,jj] = mad(mtol[this_ind])

           cfuvw4_xy_grid[ii,jj] = median(coef_w4fuv[this_ind])
           cfuvw4_xy_mad_grid[ii,jj] = mad(coef_w4fuv[this_ind])

           cnuvw4_xy_grid[ii,jj] = median(coef_w4nuv[this_ind])
           cnuvw4_xy_mad_grid[ii,jj] = mad(coef_w4nuv[this_ind])

           cfuvw3_xy_grid[ii,jj] = median(coef_w3fuv[this_ind])
           cfuvw3_xy_mad_grid[ii,jj] = mad(coef_w3fuv[this_ind])

           cnuvw3_xy_grid[ii,jj] = median(coef_w3nuv[this_ind])
           cnuvw3_xy_mad_grid[ii,jj] = mad(coef_w3nuv[this_ind])

           cjustw4_xy_grid[ii,jj] = median(coef_justw4[this_ind])
           cjustw4_xy_mad_grid[ii,jj] = mad(coef_justw4[this_ind])

           cjustw3_xy_grid[ii,jj] = median(coef_justw3[this_ind])
           cjustw3_xy_mad_grid[ii,jj] = mad(coef_justw3[this_ind])

        endif

     endfor
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE GRIDS FOR NUV/W1+W3/W1
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  count_uv_grid = fltarr(n_u, n_v)*!values.f_nan

  mtol_uv_grid = fltarr(n_u, n_v)*!values.f_nan
  mtol_uv_mad_grid = mtol_uv_grid*!values.f_nan

  cfuvw4_uv_grid = fltarr(n_u, n_v)*!values.f_nan
  cfuvw4_uv_mad_grid = cfuvw4_uv_grid*!values.f_nan

  cnuvw4_uv_grid = fltarr(n_u, n_v)*!values.f_nan
  cnuvw4_uv_mad_grid = cnuvw4_uv_grid*!values.f_nan

  cfuvw3_uv_grid = fltarr(n_u, n_v)*!values.f_nan
  cfuvw3_uv_mad_grid = cfuvw3_uv_grid*!values.f_nan

  cnuvw3_uv_grid = fltarr(n_u, n_v)*!values.f_nan
  cnuvw3_uv_mad_grid = cnuvw3_uv_grid*!values.f_nan

  cjustw4_uv_grid = fltarr(n_u, n_v)*!values.f_nan
  cjustw4_uv_mad_grid = cjustw4_uv_grid*!values.f_nan

  cjustw3_uv_grid = fltarr(n_u, n_v)*!values.f_nan
  cjustw3_uv_mad_grid = cjustw3_uv_grid*!values.f_nan

  for ii = 0, n_u-1 do begin
     for jj = 0, n_v-1 do begin

        use_this = abs(nuvw1 - nuvw1_axis[ii]) lt tol_nuvw1
        use_this *= abs(w3w1 - w3w1_axis[jj]) lt tol_w3w1
        
        this_ind = where(use_this, this_ct)
        if this_ct gt thresh then begin

           count_uv_grid[ii,jj] = this_ct*1.0

           mtol_uv_grid[ii,jj] = median(mtol[this_ind])
           mtol_uv_mad_grid[ii,jj] = mad(mtol[this_ind])

           cfuvw4_uv_grid[ii,jj] = median(coef_w4fuv[this_ind])
           cfuvw4_uv_mad_grid[ii,jj] = mad(coef_w4fuv[this_ind])

           cnuvw4_uv_grid[ii,jj] = median(coef_w4nuv[this_ind])
           cnuvw4_uv_mad_grid[ii,jj] = mad(coef_w4nuv[this_ind])

           cfuvw3_uv_grid[ii,jj] = median(coef_w3fuv[this_ind])
           cfuvw3_uv_mad_grid[ii,jj] = mad(coef_w3fuv[this_ind])

           cnuvw3_uv_grid[ii,jj] = median(coef_w3nuv[this_ind])
           cnuvw3_uv_mad_grid[ii,jj] = mad(coef_w3nuv[this_ind])

           cjustw4_uv_grid[ii,jj] = median(coef_justw4[this_ind])
           cjustw4_uv_mad_grid[ii,jj] = mad(coef_justw4[this_ind])

           cjustw3_uv_grid[ii,jj] = median(coef_justw3[this_ind])
           cjustw3_uv_mad_grid[ii,jj] = mad(coef_justw3[this_ind])
        endif

        
     endfor
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE TWOD GRIDS FOR SSFR VS MSTAR
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  count_pq_grid = fltarr(n_p, n_q)*!values.f_nan

  mtol_pq_grid = fltarr(n_p, n_q)*!values.f_nan
  mtol_pq_mad_grid = mtol_pq_grid*!values.f_nan

  cfuvw4_pq_grid = fltarr(n_p, n_q)*!values.f_nan
  cfuvw4_pq_mad_grid = cfuvw4_pq_grid*!values.f_nan

  cnuvw4_pq_grid = fltarr(n_p, n_q)*!values.f_nan
  cnuvw4_pq_mad_grid = cnuvw4_pq_grid*!values.f_nan

  cfuvw3_pq_grid = fltarr(n_p, n_q)*!values.f_nan
  cfuvw3_pq_mad_grid = cfuvw3_pq_grid*!values.f_nan

  cnuvw3_pq_grid = fltarr(n_p, n_q)*!values.f_nan
  cnuvw3_pq_mad_grid = cnuvw3_pq_grid*!values.f_nan

  cjustw4_pq_grid = fltarr(n_p, n_q)*!values.f_nan
  cjustw4_pq_mad_grid = cjustw4_pq_grid*!values.f_nan

  cjustw3_pq_grid = fltarr(n_p, n_q)*!values.f_nan
  cjustw3_pq_mad_grid = cjustw3_pq_grid*!values.f_nan

  cfuv_pq_grid = fltarr(n_p, n_q)*!values.f_nan
  cfuv_pq_mad_grid = cfuv_pq_grid*!values.f_nan

  for ii = 0, n_p-1 do begin
     for jj = 0, n_q-1 do begin

        use_this = abs(ssfr - ssfr_axis[jj]) lt tol_ssfr

        if ii eq 0 then $
           use_this *= (mstar lt tol_mstar+mstar_axis[ii]) $
        else if ii eq n_q-1 then $
           use_this *= (mstar gt mstar_axis[ii]-tol_mstar) $
        else $
           use_this *= (mstar gt mstar_axis[ii]-tol_mstar and $
                        mstar lt mstar_axis[ii]+tol_mstar)
        
        this_ind = where(use_this, this_ct)
        if this_ct gt thresh then begin

           mtol_pq_grid[ii,jj] = median(mtol[this_ind])
           mtol_pq_mad_grid[ii,jj] = mad(mtol[this_ind])
           count_pq_grid[ii,jj] = this_ct*1.0

           cfuv_pq_grid[ii,jj] = median(coef_fuv[this_ind])
           cfuv_pq_mad_grid[ii,jj] = mad(coef_fuv[this_ind])

           cfuvw4_pq_grid[ii,jj] = median(coef_w4fuv[this_ind])
           cfuvw4_pq_mad_grid[ii,jj] = mad(coef_w4fuv[this_ind])

           cnuvw4_pq_grid[ii,jj] = median(coef_w4nuv[this_ind])
           cnuvw4_pq_mad_grid[ii,jj] = mad(coef_w4nuv[this_ind])

           cfuvw3_pq_grid[ii,jj] = median(coef_w3fuv[this_ind])
           cfuvw3_pq_mad_grid[ii,jj] = mad(coef_w3fuv[this_ind])

           cnuvw3_pq_grid[ii,jj] = median(coef_w3nuv[this_ind])
           cnuvw3_pq_mad_grid[ii,jj] = mad(coef_w3nuv[this_ind])

           cjustw4_pq_grid[ii,jj] = median(coef_justw4[this_ind])
           cjustw4_pq_mad_grid[ii,jj] = mad(coef_justw4[this_ind])

           cjustw3_pq_grid[ii,jj] = median(coef_justw3[this_ind])
           cjustw3_pq_mad_grid[ii,jj] = mad(coef_justw3[this_ind])

        endif
        
     endfor
  endfor
 
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SAVE THE GRIDS AS FITS FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; SAVE THE XY GRIDS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; MAKE HEADERS

  mkhdr, xy_hdr, mtol_xy_grid
  
  sxaddpar, xy_hdr, 'CTYPE1', 'LOGFUVW1'
  sxaddpar, xy_hdr, 'CRPIX1', 1.0
  sxaddpar, xy_hdr, 'CRVAL1', fuvw1_axis[0]
  sxaddpar, xy_hdr, 'CDELT1', fuvw1_axis[1]-fuvw1_axis[0]

  sxaddpar, xy_hdr, 'CTYPE2', 'LOGW4W1'
  sxaddpar, xy_hdr, 'CRPIX2', 1.0
  sxaddpar, xy_hdr, 'CRVAL2', w4w1_axis[0]
  sxaddpar, xy_hdr, 'CDELT2', w4w1_axis[1]-w4w1_axis[0]

  mtol_xy_hdr = xy_hdr
  sxaddpar, mtol_xy_hdr, 'BUNIT', 'MSUN/LSUN', 'at 3.4um'

  sfr_xy_hdr = xy_hdr
  sxaddpar, sfr_xy_hdr, 'BUNIT', 'COEF', 'WISE->SFR'

  ct_xy_hdr = xy_hdr
  sxaddpar, sfr_xy_hdr, 'BUNIT', 'COUNTS', 'GSWLC GALAXIES'

; WRITE GRIDS

; ... COUNTS

  writefits $
     , '../measurements/fuvw1_w4w1_count_grid.fits' $
     , count_xy_grid, ct_xy_hdr

; ... M/L

  writefits $
     , '../measurements/fuvw1_w4w1_mtol_grid.fits' $
     , mtol_xy_grid, mtol_xy_hdr

  writefits $
     , '../measurements/fuvw1_w4w1_madmtol_grid.fits' $
     , mtol_xy_mad_grid, mtol_xy_hdr

; ... FUV

  writefits $
     , '../measurements/fuvw1_w4w1_cfuvw4_grid.fits' $
     , cfuvw4_xy_grid, sfr_xy_hdr

  writefits $
     , '../measurements/fuvw1_w4w1_madcfuvw4_grid.fits' $
     , cfuvw4_xy_mad_grid, sfr_xy_hdr

  writefits $
     , '../measurements/fuvw1_w4w1_cfuvw3_grid.fits' $
     , cfuvw3_xy_grid, sfr_xy_hdr

  writefits $
     , '../measurements/fuvw1_w4w1_madcfuvw3_grid.fits' $
     , cfuvw3_xy_mad_grid, sfr_xy_hdr

  writefits $
     , '../measurements/fuvw1_w4w1_cjustw4_grid.fits' $
     , cjustw4_xy_grid, sfr_xy_hdr

  writefits $
     , '../measurements/fuvw1_w4w1_madcjustw4_grid.fits' $
     , cjustw4_xy_mad_grid, sfr_xy_hdr

; ... NUV

  writefits $
     , '../measurements/fuvw1_w4w1_cnuvw4_grid.fits' $
     , cnuvw4_xy_grid, sfr_xy_hdr

  writefits $
     , '../measurements/fuvw1_w4w1_madcnuvw4_grid.fits' $
     , cnuvw4_xy_mad_grid, sfr_xy_hdr

  writefits $
     , '../measurements/fuvw1_w4w1_cnuvw3_grid.fits' $
     , cnuvw3_xy_grid, sfr_xy_hdr

  writefits $
     , '../measurements/fuvw1_w4w1_madcnuvw3_grid.fits' $
     , cnuvw3_xy_mad_grid, sfr_xy_hdr

  writefits $
     , '../measurements/fuvw1_w4w1_cjustw3_grid.fits' $
     , cjustw3_xy_grid, sfr_xy_hdr

  writefits $
     , '../measurements/fuvw1_w4w1_madcjustw3_grid.fits' $
     , cjustw3_xy_mad_grid, sfr_xy_hdr

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; SAVE THE UV GRIDS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; MAKE HEADERS

  mkhdr, uv_hdr, mtol_uv_grid

  sxaddpar, uv_hdr, 'CTYPE1', 'LOGNUVW1'
  sxaddpar, uv_hdr, 'CRPIX1', 1.0
  sxaddpar, uv_hdr, 'CRVAL1', nuvw1_axis[0]
  sxaddpar, uv_hdr, 'CDELT1', nuvw1_axis[1]-nuvw1_axis[0]

  sxaddpar, uv_hdr, 'CTYPE2', 'LOGW3W1'
  sxaddpar, uv_hdr, 'CRPIX2', 1.0
  sxaddpar, uv_hdr, 'CRVAL2', w3w1_axis[0]
  sxaddpar, uv_hdr, 'CDELT2', w3w1_axis[1]-w3w1_axis[0]

  mtol_uv_hdr = uv_hdr
  sxaddpar, mtol_uv_hdr, 'BUNIT', 'MSUN/LSUN', 'at 3.4um'

  sfr_uv_hdr = uv_hdr
  sxaddpar, sfr_uv_hdr, 'BUNIT', 'COEF', 'WISE->SFR'

  ct_uv_hdr = uv_hdr
  sxaddpar, sfr_uv_hdr, 'BUNIT', 'COUNTS', 'GSWLC GALAXIES'

; WRITE GRIDS

; ... COUNTS

  writefits $
     , '../measurements/nuvw1_w3w1_count_grid.fits' $
     , count_uv_grid, ct_uv_hdr

; ... M/L

  writefits $
     , '../measurements/nuvw1_w3w1_mtol_grid.fits' $
     , mtol_uv_grid, mtol_uv_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_madmtol_grid.fits' $
     , mtol_uv_mad_grid, mtol_uv_hdr

; ... FUV

  writefits $
     , '../measurements/nuvw1_w3w1_cfuvw4_grid.fits' $
     , cfuvw4_uv_grid, sfr_uv_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_madcfuvw4_grid.fits' $
     , cfuvw4_uv_mad_grid, sfr_uv_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_cfuvw3_grid.fits' $
     , cfuvw3_uv_grid, sfr_uv_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_madcfuvw3_grid.fits' $
     , cfuvw3_uv_mad_grid, sfr_uv_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_cjustw4_grid.fits' $
     , cjustw4_uv_grid, sfr_uv_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_madcjustw4_grid.fits' $
     , cjustw4_uv_mad_grid, sfr_uv_hdr

; ... NUV

  writefits $
     , '../measurements/nuvw1_w3w1_cnuvw4_grid.fits' $
     , cnuvw4_uv_grid, sfr_uv_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_madcnuvw4_grid.fits' $
     , cnuvw4_uv_mad_grid, sfr_uv_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_cnuvw3_grid.fits' $
     , cnuvw3_uv_grid, sfr_uv_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_madcnuvw3_grid.fits' $
     , cnuvw3_uv_mad_grid, sfr_uv_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_cjustw3_grid.fits' $
     , cjustw3_uv_grid, sfr_uv_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_madcjustw3_grid.fits' $
     , cjustw3_uv_mad_grid, sfr_uv_hdr

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; SAVE THE SSFR-MSTAR GRIDS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; MAKE HEADERS

  mkhdr, pq_hdr, mtol_pq_grid

  sxaddpar, pq_hdr, 'CTYPE1', 'LOGMSTAR'
  sxaddpar, pq_hdr, 'CRPIX1', 1.0
  sxaddpar, pq_hdr, 'CRVAL1', mstar_axis[0]
  sxaddpar, pq_hdr, 'CDELT1', mstar_axis[1]-mstar_axis[0]

  sxaddpar, pq_hdr, 'CTYPE2', 'LOGSSFR'
  sxaddpar, pq_hdr, 'CRPIX2', 1.0
  sxaddpar, pq_hdr, 'CRVAL2', ssfr_axis[0]
  sxaddpar, pq_hdr, 'CDELT2', ssfr_axis[1]-ssfr_axis[0]

  mtol_pq_hdr = pq_hdr
  sxaddpar, mtol_pq_hdr, 'BUNIT', 'MSUN/LSUN', 'at 3.4um'

  sfr_pq_hdr = pq_hdr
  sxaddpar, sfr_pq_hdr, 'BUNIT', 'COEF', 'WISE->SFR'

  ct_pq_hdr = pq_hdr
  sxaddpar, sfr_pq_hdr, 'BUNIT', 'COUNTS', 'GSWLC GALAXIES'

; WRITE GRIDS

; ... COUNTS

  writefits $
     , '../measurements/ssfr_mstar_count_grid.fits' $
     , count_pq_grid, ct_pq_hdr

; ... M/L

  writefits $
     , '../measurements/ssfr_mstar_mtol_grid.fits' $
     , mtol_pq_grid, mtol_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_madmtol_grid.fits' $
     , mtol_pq_mad_grid, mtol_pq_hdr

; ... FUV

  writefits $
     , '../measurements/ssfr_mstar_cfuv_grid.fits' $
     , cfuv_pq_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_madcfuv_grid.fits' $
     , cfuv_pq_mad_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_cfuvw4_grid.fits' $
     , cfuvw4_pq_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_madcfuvw4_grid.fits' $
     , cfuvw4_pq_mad_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_cfuvw3_grid.fits' $
     , cfuvw3_pq_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_madcfuvw3_grid.fits' $
     , cfuvw3_pq_mad_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_cjustw4_grid.fits' $
     , cjustw4_pq_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_madcjustw4_grid.fits' $
     , cjustw4_pq_mad_grid, sfr_pq_hdr

; ... NUV

  writefits $
     , '../measurements/ssfr_mstar_cnuvw4_grid.fits' $
     , cnuvw4_pq_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_madcnuvw4_grid.fits' $
     , cnuvw4_pq_mad_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_cnuvw3_grid.fits' $
     , cnuvw3_pq_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_madcnuvw3_grid.fits' $
     , cnuvw3_pq_mad_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_cjustw3_grid.fits' $
     , cjustw3_pq_grid, sfr_pq_hdr

  writefits $
     , '../measurements/ssfr_mstar_madcjustw3_grid.fits' $
     , cjustw3_pq_mad_grid, sfr_pq_hdr

  stop

end
