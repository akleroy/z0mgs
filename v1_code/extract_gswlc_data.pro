pro extract_gswlc_data

  @constants.bat
  lsun_3p4 = 1.83d18

  nu_fuv = c/(154.d-9*1d2)
  nu_nuv = c/(231.d-9*1d2)
  nu_w1 = c/(3.4d-6*1d2)
  nu_w2 = c/(4.5d-6*1d2)
  nu_w3 = c/(12.d-6*1d2)
  nu_w4 = c/(22.d-6*1d2)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GSWLC CATALOG
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;  gws = mrdfits('~/idl/galbase/gal_data/'+$
;                'hlsp_gswlc_galex-sdss-wise_multi_x1_multi_v1_cat.fits' $
;                , 1, h_gws)
; gws_logsfrsed = gws.logsfrsed
; gws_logmstar = gws.logmstar
; gws_z = gws.z

  readcol, '/data/kant/0/leroy.42/allsky/gswlc/GSWLC-X2.dat' $
           , gws_objid, gws_glxid, gws_plate, gws_mjd, gws_fiber $
           , gws_ra, gws_dec, gws_z, gws_chisq $
           , gws_logmstar, gws_elogmstar, gws_logsfrsed, gws_elogsfrsed $
           , gws_afuv, gws_eafuv, gws_ab, gws_eab, gws_av, gws_eav $
           , gws_flagsed, gws_uvsurvey, gws_uvflag, gws_midir, gws_mgs $
           , format='LL,LL,LL,L,L'+',F,F,F,F'+',F,F,F,F'+',F,F,F,F,F,F'+',L,L,L,L,L'
  n_gws = n_elements(gws_objid)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GWSLC FLUXES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, '/data/kant/0/leroy.42/allsky/gswlc/galex_unwise_fluxes_GSWLC-X.dat' $
           , flux_id $
           , flux_gws_fuv, flux_gws_efuv $
           , flux_gws_nuv, flux_gws_enuv $
           , flux_gws_w1_nm, flux_gws_w1_invar $
           , flux_gws_w2_nm, flux_gws_w2_invar $
           , flux_gws_w3_nm, flux_gws_w3_invar $
           , flux_gws_w4_nm, flux_gws_w4_invar $
           , format='LL,F,F,F,F,F,F,F,F,F,F,F,F'
  
  flux_gws_w1 = 3631.*10^(-0.4*(22.5+2.683))*flux_gws_w1_nm
  flux_gws_w2 = 3631*10^(-0.4*(22.5+3.319))*flux_gws_w2_nm
  flux_gws_w3 = 3631*10^(-0.4*(22.5+5.242))*flux_gws_w3_nm
  flux_gws_w4 = 3631*10^(-0.4*(22.5+6.604))*flux_gws_w4_nm

  flux_gws_ew1 = 3631.*10^(-0.4*(22.5+2.683))*(1./sqrt(flux_gws_w1_invar))
  flux_gws_ew2 = 3631*10^(-0.4*(22.5+3.319))*(1./sqrt(flux_gws_w2_invar))
  flux_gws_ew3 = 3631*10^(-0.4*(22.5+5.242))*(1./sqrt(flux_gws_w3_invar))
  flux_gws_ew4 = 3631*10^(-0.4*(22.5+6.604))*(1./sqrt(flux_gws_w4_invar))

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CROSS MATCH
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  gws_fuv = !values.f_nan*gws_logmstar
  gws_efuv = !values.f_nan*gws_logmstar
  gws_nuv = !values.f_nan*gws_logmstar
  gws_enuv = !values.f_nan*gws_logmstar
  gws_w1 = !values.f_nan*gws_logmstar
  gws_ew1 = !values.f_nan*gws_logmstar
  gws_w2 = !values.f_nan*gws_logmstar
  gws_ew2 = !values.f_nan*gws_logmstar
  gws_w3 = !values.f_nan*gws_logmstar
  gws_ew3 = !values.f_nan*gws_logmstar
  gws_w4 = !values.f_nan*gws_logmstar
  gws_ew4 = !values.f_nan*gws_logmstar

  flux_ind = lonarr(n_gws)
  for ii = 0, n_gws-1 do begin
     counter, ii, n_gws, 'Cross match '
     flux_ind[ii] = where(flux_id eq gws_objid[ii], flux_ct)
     if flux_ct gt 1 then begin
        print, 'Double match. I do not think this should happen.'
        stop
     endif
  endfor
  use = where(flux_ind ne -1, keep_ct)
  flux_ind = flux_ind[use]

  gws_fuv[use] = flux_gws_fuv[flux_ind]
  gws_efuv[use] = flux_gws_efuv[flux_ind]

  gws_nuv[use] = flux_gws_nuv[flux_ind]
  gws_enuv[use] = flux_gws_enuv[flux_ind]

  gws_w1[use] = flux_gws_w1[flux_ind]
  gws_ew1[use] = flux_gws_ew1[flux_ind]

  gws_w2[use] = flux_gws_w2[flux_ind]
  gws_ew2[use] = flux_gws_ew2[flux_ind]

  gws_w3[use] = flux_gws_w3[flux_ind]
  gws_ew3[use] = flux_gws_ew3[flux_ind]

  gws_w4[use] = flux_gws_w4[flux_ind]
  gws_ew4[use] = flux_gws_ew4[flux_ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVERT UNITS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  h0 = 70.0d
  dist_cm = lumdist(gws_z, h0=h0, omega_m = 0.27, lambda=0.73)*1d6*pc
  dumb_dist = gws_z*c/1d5/70.*1d6*pc

  fuv_lum = gws_fuv/1d3*1d-23*4.*!pi*(dist_cm)^2
  nuv_lum = gws_nuv/1d3*1d-23*4.*!pi*(dist_cm)^2
  w1_lum = gws_w1*1d-23*4.*!pi*(dist_cm)^2
  w2_lum = gws_w2*1d-23*4.*!pi*(dist_cm)^2
  w3_lum = gws_w3*1d-23*4.*!pi*(dist_cm)^2
  w4_lum = gws_w4*1d-23*4.*!pi*(dist_cm)^2

  sfr_fuv_ke12 = $
     lum_to_sfr(band='FUV', cal='KE12', lum=fuv_lum*nu_fuv)
  sfr_nuv_ke12 = $
     lum_to_sfr(band='NUV', cal='KE12', lum=nuv_lum*nu_nuv)

  sfr_fuv_s07 = $
     lum_to_sfr(band='FUV', cal='S07', lum=fuv_lum*nu_fuv)
  sfr_nuv_s07 = $
     lum_to_sfr(band='NUV', cal='S07', lum=nuv_lum*nu_nuv)

  sfr_fuvw4_ke12 = $
     lum_to_sfr(band='FUV', cal='KE12' $
                , lum=(fuv_lum*nu_fuv + 3.89*w4_lum*nu_w4))  

  sfr_fuv_z19 = $
     lum_to_sfr(band='FUV', cal='Z19', lum=fuv_lum*nu_fuv)
  sfr_nuv_z19 = $
     lum_to_sfr(band='NUV', cal='Z19', lum=nuv_lum*nu_nuv)

  sfr_w3_j13 = $
     lum_to_sfr(band='WISE3', cal='J13', lum=w3_lum*nu_w3)
  sfr_w4_j13 = $
     lum_to_sfr(band='WISE4', cal='J13', lum=w4_lum*nu_w4)

  sfr_w3_z19 = $
     lum_to_sfr(band='WISE3', cal='Z19', lum=w3_lum*nu_w3)
  sfr_w4_z19 = $
     lum_to_sfr(band='WISE4', cal='Z19', lum=w4_lum*nu_w4)

  sfr_fuvw4_z19 = $
     lum_to_sfr(band='FUV', cal='Z19', lum=fuv_lum*nu_fuv) + $
     lum_to_sfr(band='WISE4+FUV', cal='Z19', lum=w4_lum*nu_w4)
  
  sfr_fuvw3_z19 = $
     lum_to_sfr(band='FUV', cal='Z19', lum=fuv_lum*nu_fuv) + $
     lum_to_sfr(band='WISE3+FUV', cal='Z19', lum=w3_lum*nu_w3)

  sfr_nuvw4_z19 = $
     lum_to_sfr(band='NUV', cal='Z19', lum=nuv_lum*nu_nuv) + $
     lum_to_sfr(band='WISE4+NUV', cal='Z19', lum=w4_lum*nu_w4)
  
  sfr_nuvw3_z19 = $
     lum_to_sfr(band='NUV', cal='Z19', lum=nuv_lum*nu_nuv) + $
     lum_to_sfr(band='WISE3+NUV', cal='Z19', lum=w3_lum*nu_w3)

  mtol_w1 = 10.^(gws_logmstar - alog10(w1_lum/lsun_3p4))
  gws_ssfr = gws_logsfrsed-gws_logmstar
  w2w1 = alog10(w2_lum/w1_lum)
  w3w1 = alog10(w3_lum/w1_lum)
  w4w1 = alog10(w3_lum/w1_lum)
  nuvw1 = alog10(nuv_lum/w1_lum)
  fuvw1 = alog10(fuv_lum/w1_lum)
  
  ssfr_like_fuvw4 = sfr_fuvw4_z19 / (w1_lum/lsun_3p4)
  ssfr_like_nuvw4 = sfr_nuvw4_z19 / (w1_lum/lsun_3p4)
  ssfr_like_nuvw3 = sfr_nuvw3_z19 / (w1_lum/lsun_3p4)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SAVE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  save $
     , gws_objid, gws_glxid, gws_plate, gws_mjd, gws_fiber $
     , gws_ra, gws_dec, gws_z, gws_chisq $
     , gws_logmstar, gws_elogmstar, gws_logsfrsed, gws_elogsfrsed $
     , gws_afuv, gws_eafuv, gws_ab, gws_eab, gws_av, gws_eav $
     , gws_flagsed, gws_uvsurvey, gws_uvflag, gws_midir, gws_mgs $
     , gws_fuv, gws_efuv, gws_nuv, gws_enuv, gws_w1, gws_ew1 $
     , gws_w2, gws_ew2, gws_w3, gws_ew3, gws_w4, gws_ew4 $
     , fuv_lum, nuv_lum, w1_lum, w2_lum, w3_lum, w4_lum $
     , sfr_fuv_ke12, sfr_nuv_ke12, sfr_fuv_s07, sfr_nuv_s07, sfr_fuv_z19, sfr_nuv_z19 $
     , sfr_w3_j13, sfr_w4_j13, sfr_w3_z19, sfr_w4_z19 $
     , sfr_fuvw4_ke12, sfr_fuvw4_z19, sfr_fuvw3_z19, sfr_nuvw4_z19, sfr_nuvw3_z19 $
     , mtol_w1, gws_ssfr, ssfr_like_fuvw4, ssfr_like_nuvw4, ssfr_like_nuvw3 $
     , w2w1, w3w1, w4w1, nuvw1, fuvw1 $
     , file='../gswlc/gswlc_data.idl'  
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RECAST INTO A TABLE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  nan = !values.f_nan  
  empty = $
     { $
     gws_objid:-1L, $
     gws_glxid:-1L, $
     gws_plate:-1L, $
     gws_mjd:-1L, $
     gws_fiber:-1L, $
     gws_ra: nan, $
     gws_dec: nan, $
     gws_z: nan, $ 
     gws_chisq: nan, $
     gws_logmstar: nan, $     
     gws_elogmstar: nan, $    
     gws_logsfrsed: nan, $   
     gws_elogsfrsed: nan, $    
     gws_afuv: nan, $
     gws_eafuv: nan, $    
     gws_ab: nan, $
     gws_eab: nan, $    
     gws_av: nan, $  
     gws_eav: nan, $    
     gws_flagsed: -1L, $
     gws_uvsurvey: -1L, $
     gws_uvflag: -1L, $
     gws_midir: -1L, $
     gws_mgs: -1L, $
     gws_fuv: nan, $
     gws_efuv: nan, $
     gws_nuv: nan, $
     gws_enuv: nan, $
     gws_w1: nan, $
     gws_ew1: nan, $
     gws_w2: nan, $
     gws_ew2: nan, $
     gws_w3: nan, $
     gws_ew3: nan, $
     gws_w4: nan, $
     gws_ew4: nan, $
     dist_cm: nan, $
     fuv_lum: nan, $
     nuv_lum: nan, $
     w1_lum: nan, $
     w2_lum: nan, $
     w3_lum: nan, $
     w4_lum: nan $
     }

  n = n_elements(w1_lum)
  dat = replicate(empty, n)
  for ii = 0, n-1 do begin
     counter, ii, n, 'Filling entry: '

     dat[ii].gws_objid = gws_objid[ii]
     dat[ii].gws_glxid = gws_glxid[ii]
     dat[ii].gws_plate = gws_plate[ii]
     dat[ii].gws_mjd = gws_mjd[ii]
     dat[ii].gws_fiber = gws_fiber[ii]
     dat[ii].gws_ra = gws_ra[ii]
     dat[ii].gws_dec = gws_dec[ii]
     dat[ii].gws_z = gws_z[ii]
     dat[ii].gws_chisq = gws_chisq[ii]
     dat[ii].gws_logmstar = gws_logmstar[ii]
     dat[ii].gws_elogmstar = gws_elogmstar[ii]
     dat[ii].gws_logsfrsed = gws_logsfrsed[ii]
     dat[ii].gws_elogsfrsed = gws_elogsfrsed[ii]
     dat[ii].gws_afuv = gws_afuv[ii]
     dat[ii].gws_eafuv = gws_eafuv[ii]
     dat[ii].gws_ab = gws_ab[ii]
     dat[ii].gws_eab = gws_eab[ii]
     dat[ii].gws_av = gws_av[ii]
     dat[ii].gws_eav = gws_eav[ii]
     dat[ii].gws_flagsed = gws_flagsed[ii]
     dat[ii].gws_uvsurvey = gws_uvsurvey[ii]
     dat[ii].gws_uvflag = gws_uvflag[ii]
     dat[ii].gws_midir = gws_midir[ii]
     dat[ii].gws_mgs = gws_mgs[ii]
     dat[ii].gws_fuv = gws_fuv[ii]
     dat[ii].gws_efuv = gws_efuv[ii]
     dat[ii].gws_nuv = gws_nuv[ii]
     dat[ii].gws_enuv = gws_enuv[ii]
     dat[ii].gws_w1 = gws_w1[ii]
     dat[ii].gws_ew1 = gws_ew1[ii]
     dat[ii].gws_w2 = gws_w2[ii]
     dat[ii].gws_ew2 = gws_ew2[ii]
     dat[ii].gws_w3 = gws_w3[ii]
     dat[ii].gws_ew3 = gws_ew3[ii]
     dat[ii].gws_w4 = gws_w4[ii]
     dat[ii].gws_ew4 = gws_ew4[ii]
     dat[ii].dist_cm = dist_cm[ii]
     dat[ii].fuv_lum = fuv_lum[ii]
     dat[ii].nuv_lum = nuv_lum[ii]
     dat[ii].w1_lum = w1_lum[ii]
     dat[ii].w2_lum = w2_lum[ii]
     dat[ii].w3_lum = w3_lum[ii]
     dat[ii].w4_lum = w4_lum[ii]

  endfor

  mwrfits, dat, '../gswlc/gswlc_with_fluxes.fits', /create

  stop

end
