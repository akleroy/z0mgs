pro build_gswlc_grid

  @constants.bat
  lsun_3p4 = 1.83d18
  thresh = 10

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

  nuvw1 = (alog10(nuv_lum/w1_lum))[sane_ind]
  w3w1 = (alog10(w3_lum/w1_lum))[sane_ind]
  lumw1 = (alog10(w1_lum/lsun_3p4))[sane_ind]
  mtol = mtol_w1[sane_ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE AXES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  min_nuvw1 = -3.5
  max_nuvw1 = 0.5
  bin_nuvw1 = 0.05
  tol_nuvw1 = 0.1

  min_w3w1 = -2.0
  max_w3w1 = 1.5
  bin_w3w1 = 0.05
  tol_w3w1 = 0.1

  min_lumw1 = 9.7
  max_lumw1 = 11.95
  binsize_lumw1 = 0.15
  tol_lumw1 = 0.3

  count_nuvw1_w3w1 = $
     grid_data(nuvw1, w3w1, /nan $
               , xaxis_out = nuvw1_axis, yaxis_out = w3w1_axis $
               , xmin=min_nuvw1, xmax=max_nuvw1, binsize_x=binsize_nuvw1 $
               , ymin=min_w3w1, ymax=max_w3w1, binsize_y=binsize_w3w1 $
              )
  n_x = n_elements(nuvw1_axis)
  n_y = n_elements(w3w1_axis)

  count_w3w1_lumw1 = $
     grid_data(w3w1, w1_lum, /nan $
               , xaxis_out = w3w1_axis, yaxis_out = lumw1_axis $
               , xmin=min_w3w1, xmax=max_w3w1, binsize_x=binsize_w3w1 $
               , ymin=min_lumw1, ymax=max_lumw1, binsize_y=binsize_lumw1 $
              )
  n_z = n_elements(lumw1_axis)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE TWOD GRIDS FOR NUV/W1+W3/W1 and W3/W1+LUM_W1
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  xy_grid = fltarr(n_x, n_y)*!values.f_nan
  xy_mad_grid = xy_grid*!values.f_nan

  for ii = 0, n_x-1 do begin
     for jj = 0, n_y-1 do begin

        use = abs(nuvw1 - nuvw1_axis[ii]) lt tol_nuvw1
        use *= abs(w3w1 - w3w1_axis[jj]) lt tol_w3w1
        ind = where(use, ct)
        if ct gt thresh then begin
           xy_grid[ii,jj] = median(mtol[ind])
           xy_mad_grid[ii,jj] = mad(mtol[ind])
        endif

     endfor
  endfor

  yz_grid = fltarr(n_y, n_z)*!values.f_nan
  yz_mad_grid = yz_grid*!values.f_nan

  for ii = 0, n_y-1 do begin
     for jj = 0, n_z-1 do begin

        use = abs(w3w1 - w3w1_axis[ii]) lt tol_w3w1

        if jj eq 0 then $
           use *= (lumw1 lt tol_lumw1+lumw1_axis[jj]) $
        else if jj eq n_z-1 then $
           use *= (lumw1 gt lumw1_axis[jj]-tol_lumw1) $
        else $
           use *= (lumw1 gt lumw1_axis[jj]-tol_lumw1 and $
                   lumw1 lt lumw1_axis[jj]+tol_lumw1)
        
        ind = where(use, ct)
        if ct gt thresh then begin
           yz_grid[ii,jj] = median(mtol[ind])
           yz_mad_grid[ii,jj] = mad(mtol[ind])
        endif
        
     endfor
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE A THREE D CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  xyz_cube = fltarr(n_x, n_y, n_z)*!values.f_nan
  xyz_mad_cube = xyz_cube*!values.f_nan

  for ii = 0, n_x-1 do begin
     for jj = 0, n_y-1 do begin
        for kk = 0, n_z-1 do begin

           use = abs(nuvw1 - nuvw1_axis[ii]) lt tol_nuvw1
           use *= abs(w3w1 - w3w1_axis[jj]) lt tol_w3w1

           if kk eq 0 then $
              use *= (lumw1 lt tol_lumw1+lumw1_axis[kk]) $
           else if kk eq n_z-1 then $
              use *= (lumw1 gt lumw1_axis[kk]-tol_lumw1) $
           else $
              use *= (lumw1 gt lumw1_axis[kk]-tol_lumw1 and $
                      lumw1 lt lumw1_axis[kk]+tol_lumw1)

           ind = where(use, ct)           
           if ct gt thresh then begin
              xyz_cube[ii,jj,kk] = median(mtol[ind])
              xyz_mad_cube[ii,jj,kk] = mad(mtol[ind])
           endif
        endfor

     endfor
  endfor
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SAVE THE GRIDS AS FITS FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; SAVE THE XY GRID
  mkhdr, xy_hdr, xy_grid
  
  sxaddpar, xy_hdr, 'CTYPE1', 'LOGNUVW1'
  sxaddpar, xy_hdr, 'CRPIX1', 1.0
  sxaddpar, xy_hdr, 'CRVAL1', nuvw1_axis[0]
  sxaddpar, xy_hdr, 'CDELT1', nuvw1_axis[1]-nuvw1_axis[0]

  sxaddpar, xy_hdr, 'CTYPE2', 'LOGW3W1'
  sxaddpar, xy_hdr, 'CRPIX2', 1.0
  sxaddpar, xy_hdr, 'CRVAL2', w3w1_axis[0]
  sxaddpar, xy_hdr, 'CDELT2', w3w1_axis[1]-w3w1_axis[0]

  sxaddpar, xy_hdr, 'BUNIT', 'MSUN/LSUN', 'at 3.4um'

  writefits $
     , '../measurements/nuvw1_w3w1_mtol_grid.fits' $
     , xy_grid, xy_hdr

  writefits $
     , '../measurements/nuvw1_w3w1_madmtol_grid.fits' $
     , xy_mad_grid, xy__hdr

; SAVE THE YZ GRIDS
  mkhdr, yz_hdr, yz_grid

  sxaddpar, yz_hdr, 'CTYPE1', 'LOGW3W1'
  sxaddpar, yz_hdr, 'CRPIX1', 1.0
  sxaddpar, yz_hdr, 'CRVAL1', w3w1_axis[0]
  sxaddpar, yz_hdr, 'CDELT1', w3w1_axis[1]-w3w1_axis[0]

  sxaddpar, yz_hdr, 'CTYPE2', 'LOGW1LSUN'
  sxaddpar, yz_hdr, 'CRPIX2', 1.0
  sxaddpar, yz_hdr, 'CRVAL2', lumw1_axis[0]
  sxaddpar, yz_hdr, 'CDELT2', lumw1_axis[1]-lumw1_axis[0]

  sxaddpar, yz_hdr, 'BUNIT', 'MSUN/LSUN', 'at 3.4um'

  writefits $
     , '../measurements/w3w1_lumw1_mtol_grid.fits' $
     , yz_grid, yz_hdr

  writefits $
     , '../measurements/w3w1_lumw1_madmtol_grid.fits' $
     , yz_mad_grid, yz__hdr

; SAVE THE THREE-D GRID
  mkhdr, xyz_hdr, xyz_cube

  sxaddpar, xyz_hdr, 'CTYPE1', 'LOGNUVW1'
  sxaddpar, xyz_hdr, 'CRPIX1', 1.0
  sxaddpar, xyz_hdr, 'CRVAL1', nuvw1_axis[0]
  sxaddpar, xyz_hdr, 'CDELT1', nuvw1_axis[1]-nuvw1_axis[0]

  sxaddpar, xyz_hdr, 'CTYPE2', 'LOGW3W1'
  sxaddpar, xyz_hdr, 'CRPIX2', 1.0
  sxaddpar, xyz_hdr, 'CRVAL2', w3w1_axis[0]
  sxaddpar, xyz_hdr, 'CDELT2', w3w1_axis[1]-w3w1_axis[0]

  sxaddpar, xyz_hdr, 'CTYPE3', 'LOGW1LSUN'
  sxaddpar, xyz_hdr, 'CRPIX3', 1.0
  sxaddpar, xyz_hdr, 'CRVAL3', lumw1_axis[0]
  sxaddpar, xyz_hdr, 'CDELT3', lumw1_axis[1]-lumw1_axis[0]

  sxaddpar, xyz_hdr, 'BUNIT', 'MSUN/LSUN', 'at 3.4um'

  writefits $
     , '../measurements/mtol_cube.fits' $
     , xyz_cube, xyz_hdr

  writefits $
     , '../measurements/madmtol_cube.fits' $
     , xyz_mad_cube, xyz_hdr

  stop

end
