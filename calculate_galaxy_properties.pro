pro calculate_galaxy_properties

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; READ THE Z0MGS PHOTOMETRY FROM THE DELIVERY FILE
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  tab = mrdfits('../measurements/delivery_index_gauss15.fits', 1, h)

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; INITIALIZAE THE OUTPUT
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
  
  nan = !values.d_nan

  empty = $
     { $
     pgc:0L $
     , pgc_name:'' $
     , dist_mpc: nan $
     , e_dist_dex: nan $
     , fuv: nan $
     , e_fuv: nan $
     , nuv: nan $
     , e_nuv: nan $
     , w1: nan $
     , e_w1: nan $
     , w2: nan $
     , e_w2: nan $
     , w3: nan $
     , e_w3: nan $
     , w4: nan $
     , e_w4: nan $
     , ha: nan $
     , e_ha: nan $
     , sfr_fuv: nan $
     , e_sfr_fuv: nan $
     , sfr_nuv: nan $
     , e_sfr_nuv: nan $
     , sfr_w3: nan $
     , e_sfr_w3: nan $
     , sfr_w4: nan $
     , e_sfr_w4: nan $
     , sfr_ha: nan $
     , e_sfr_ha: nan $
     , mstar: nan $
     , e_mstar: nan $
     , mtol3p6: nan $
     , mtol_code: -1 $     
     , ssfr: nan $
     , e_ssfr: nan $
     , delta_ssfr: nan $
     , e_delta_ssfr: nan $
     }

;    OPTICAL
;    2MASS
;    SPITZER
;    HERSCHEL
;    IRAS?
;    AKARI?

  props = replicate(empty, n_elements(tab))

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; CALCULATE LUMINOSITIES
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  @constants.bat
  lsun_3p6 = 1.7d18
  mtol_3p6 = 0.5
  
  pgc = tab.pgc
  n_pgc = n_elements(pgc)
  pgcname = 'PGC'+str(pgc)
  s = gal_data(pgcname)
  dist_pc = s.dist_mpc*1d6

; FILL IN OUTPUT
  props.pgc = pgc
  props.pgc_name = pgcname
  props.dist_mpc = s.dist_mpc
  props.e_dist_dex = s.e_dist_dex

  fuv = tab.flux_fuv
  e_fuv = tab.rms_flux_fuv
  nu_fuv = c/(154.d-9*1d2)
  lum_fuv = 4.*!pi*(dist_pc*pc)^2*nu_fuv*1d-23*fuv 

  nuv = tab.flux_nuv
  e_nuv = tab.rms_flux_nuv
  nu_nuv = c/(231.d-9*1d2)
  lum_nuv = 4.*!pi*(dist_pc*pc)^2*nu_nuv*1d-23*nuv

  w1 = tab.flux_wise1
  e_w1 = tab.rms_flux_wise1
  nu_w1 = c/(3.6d-6*1d2)
  lum_w1 = 4.*!pi*(dist_pc*pc)^2*nu_w1*1d-23*w1

  w2 = tab.flux_wise2
  e_w2 = tab.rms_flux_wise2
  nu_w2 = c/(4.5d-6*1d2)
  lum_w2 = 4.*!pi*(dist_pc*pc)^2*nu_w2*1d-23*w2

  w3 = tab.flux_wise3
  e_w3 = tab.rms_flux_wise3
  nu_w3 = c/(12.d-6*1d2)
  lum_w3 = 4.*!pi*(dist_pc*pc)^2*nu_w3*1d-23*w3

  w4 = tab.flux_Wise4
  e_w4 = tab.rms_flux_wise4
  nu_w4 = c/(22.d-6*1d2)
  lum_w4 = 4.*!pi*(dist_pc*pc)^2*nu_w4*1d-23*w4

; FILL IN OUTPUT

  props.fuv = lum_fuv
  props.nuv = lum_nuv
  props.w1 = lum_w1
  props.w2 = lum_w2
  props.w3 = lum_w3
  props.w4 = lum_w4

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; CALCULATE STAR FORMATION RATES
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  sfr_fuv_ke12 = $
     lum_to_sfr(band='FUV', cal='KE12', lum=lum_fuv)
  props.sfr_fuv = sfr_fuv_ke12
  props.e_sfr_fuv = e_fuv/fuv * props.sfr_fuv

  sfr_nuv_ke12 = $
     lum_to_sfr(band='NUV', cal='KE12', lum=lum_nuv)
  props.sfr_nuv = sfr_nuv_ke12
  props.e_sfr_nuv = e_nuv/nuv * props.sfr_nuv

  sfr_w3_j13 = $
     lum_to_sfr(band='WISE3', cal='J13', lum=lum_w3)
  props.sfr_w3 = sfr_w3_j13
  props.e_sfr_w3 = e_w3/w3 * props.sfr_w3

  sfr_w4_j13 = $
     lum_to_sfr(band='WISE4', cal='J13', lum=lum_w4)
  props.sfr_w4 = sfr_w4_j13
  props.e_sfr_w4 = e_w4/w4 * props.sfr_w4

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; CALCULATE STELLAR MASSES
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  lumw1 = alog10(lum_w1/lsun_3p6/nu_w1)
  w3w1 = alog10(lum_w3/nu_w3 / (lum_w1/nu_w1))
  nuvw1 = alog10(lum_nuv/nu_nuv / (lum_w1/nu_w1))

  mtol_best = lookup_mtol(nuvw1=nuvw1, w3w1=w3w1, lumw1=lumw1)
  mtol_w3w1 = lookup_mtol(w3w1=w3w1, lumw1=lumw1)
  mtol_nuvw1 = lookup_mtol(nuvw1=nuvw1, lumw1=lumw1)

  mtol = lumw1*nan
  mtol_code = finite(lumw1)-1B

  ind = where(finite(mtol_best))
  mtol[ind] = mtol_best[ind]
  mtol_code[ind] = 1

  ind = where((finite(mtol) eq 0) and finite(mtol_w3w1))
  mtol[ind] = mtol_w3w1[ind]
  mtol_code[ind] = 2

  ind = where((finite(mtol) eq 0) and finite(mtol_nuvw1))
  mtol[ind] = mtol_nuvw1[ind]
  mtol_code[ind] = 3

  mstar = lum_w1/lsun_3p6/nu_w1*mtol

  props.mtol3p6 = mtol
  props.mtol_code = mtol_code
  props.mstar = mstar
  props.e_mstar = e_w1/w1*mstar

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; WRITE TO DISK
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  mwrfits, props, '../measurements/galaxy_properties.fits', /create

  stop

end
