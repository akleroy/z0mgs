pro add_galaxy_properties_to_index

;+
;
;-
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONSTANTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  @constants.bat
  lsun_3p4 = 1.83d18

  nu_fuv = c/(154.d-9*1d2)
  nu_nuv = c/(231.d-9*1d2)
  nu_w1 = c/(3.4d-6*1d2)
  nu_w2 = c/(4.5d-6*1d2)
  nu_w3 = c/(12.d-6*1d2)
  nu_w4 = c/(22.d-6*1d2)

  fixed_mtol = 0.5
 
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOAD INDEX FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  index15 = mrdfits('../measurements/delivery_index_gauss15.fits',1,h15)
  index7p5 = mrdfits('../measurements/delivery_index_gauss7p5.fits',1,h7p5)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for jj = 0, 1 do begin

;    ... select this index
     
     if jj eq 0 then this_index = index15
     if jj eq 1 then this_index = index7p5

     n_pgc = n_elements(this_index)

;    ... look up galaxy index

     dat = gal_data(pgc=this_index.pgc, /full)     
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Calculate luminosities, first SFR estimates, and colors
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    ... ... distance

     this_index.dist_mpc = dat.dist_mpc
     this_index.e_dist_dex = dat.e_dist_dex

;    ... ... b band luminosity
     this_b_apparent = dat.btc
     this_btc_abs = this_b_apparent - dat.distmod
     this_index.absbtc = this_btc_abs
     this_index.complete_sample = $
        (this_index.absbtc le -18.) and $
        (this_index.dist_mpc le 50.)

;    ... ... luminosity

     fuv_lum = this_index.flux_fuv*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2
     e_fuv_lum = this_index.std_flux_fuv*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2

     nuv_lum = this_index.flux_nuv*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2
     e_nuv_lum = this_index.std_flux_nuv*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2

     w1_lum = this_index.flux_wise1*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2
     e_w1_lum = this_index.std_flux_wise1*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2

     w2_lum = this_index.flux_wise2*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2
     e_w2_lum = this_index.std_flux_wise2*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2

     w3_lum = this_index.flux_wise3*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2    
     e_w3_lum = this_index.std_flux_wise3*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2    

     w4_lum = this_index.flux_wise4*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2    
     e_w4_lum = this_index.std_flux_wise4*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2    
     
;    ... ... sfr estimates

     sfr_nuvw4 = $
        lum_to_sfr(band='NUV', cal='Z19', lum=nuv_lum*nu_nuv) + $
        lum_to_sfr(band='WISE4+NUV', cal='Z19', lum=w4_lum*nu_w4)
     e_sfr_nuvw4 = $
        lum_to_sfr(band='NUV', cal='Z19', lum=e_nuv_lum*nu_nuv) + $
        lum_to_sfr(band='WISE4+NUV', cal='Z19', lum=e_w4_lum*nu_w4)
     isul_nuvw4 = $
        sfr_nuvw4 le 0.0
     
     sfr_fuvw4 = $
        lum_to_sfr(band='FUV', cal='Z19', lum=fuv_lum*nu_fuv) + $
        lum_to_sfr(band='WISE4+FUV', cal='Z19', lum=w4_lum*nu_w4)
     e_sfr_fuvw4 = $
        lum_to_sfr(band='FUV', cal='Z19', lum=e_fuv_lum*nu_fuv) + $
        lum_to_sfr(band='WISE4+FUV', cal='Z19', lum=e_w4_lum*nu_w4)     
     isul_fuvw4 = $
        sfr_fuvw4 le 0.0

     sfr_justw4 = $
        lum_to_sfr(band='WISE4', cal='Z19', lum=w4_lum*nu_w4)
     e_sfr_justw4 = $
        lum_to_sfr(band='WISE4', cal='Z19', lum=e_w4_lum*nu_w4)
     isul_justw4 = $
        sfr_justw4 le 0.0
     
;    ... ... colors

     w4w1 = alog10(w4_lum/w1_lum)
     ssfrlike_fuvw4 = alog10(sfr_fuvw4 / (w1_lum/lsun_3p4))
     ssfrlike_nuvw4 = alog10(sfr_nuvw4 / (w1_lum/lsun_3p4))

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Mass-to-light ratio
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    WE HAVE NO WISE4
     ind = where(this_index.has_wise4 eq 0.0, ct)
     if ct gt 0 then begin        
        this_index[ind].mtol = fixed_mtol
        this_index[ind].method_mtol = 'FIXEDUL'
        this_index[ind].logmass = $
           alog10((w1_lum/lsun_3p4*this_index.mtol)[ind])
     endif

;    WE ONLY HAVE W4/W1 

;    ... AND IT'S NOT AN UPPER LIMIT
     ind = where(this_index.has_wise1 and $
                 this_index.has_wise4 and $
                 this_index.has_fuv eq 0 and $
                 this_index.has_nuv eq 0 and $
                 this_index.flux_wise1 gt 0 and $
                 isul_justw4 eq 0 $
                 , ct)
     if ct gt 0 then begin        
        this_index[ind].mtol = $
           lookup_mtol(w4w1=w4w1[ind])
        this_index[ind].logmass = $
           alog10((w1_lum/lsun_3p4*this_index.mtol)[ind])
        this_index[ind].method_mtol = 'W4W1'
     endif

;    ... AND IT IS AN UPPER LIMIT
     ind = where(this_index.has_wise1 and $
                 this_index.has_wise4 and $
                 this_index.has_fuv eq 0 and $
                 this_index.has_nuv eq 0 and $
                 this_index.flux_wise1 gt 0 and $
                 isul_justw4 eq 1 $
                 , ct)
     if ct gt 0 then begin        
        this_index[ind].mtol = fixed_mtol
        this_index[ind].logmass = $
           alog10((w1_lum/lsun_3p4*this_index.mtol)[ind])
        this_index[ind].method_mtol = 'FIXEDUL'
     endif

;    WE HAVE FUV AND W4

;    ... AND IT'S NOT AN UPPER LIMIT
     ind = where(this_index.has_wise1 and $
                 this_index.has_fuv and $
                 this_index.has_wise4 and $
                 this_index.flux_wise1 gt 0 and $
                 isul_fuvw4 eq 0 $
                 , ct)
     if ct gt 0 then begin
        this_index[ind].mtol = $
           lookup_mtol(ssfrlike=ssfrlike_fuvw4[ind])           
        this_index[ind].logmass = $
           alog10(w1_lum[ind]/lsun_3p4*this_index[ind].mtol)
        this_index[ind].method_mtol = 'SSFRLIKE'
     endif
     
;    ... IT IS AN UPPER LIMIT
     ind = where(this_index.has_wise1 and $
                 this_index.has_fuv and $
                 this_index.has_wise4 and $
                 this_index.flux_wise1 gt 0 and $
                 isul_fuvw4 eq 1 $
                 , ct)
     if ct gt 0 then begin
        this_index[ind].mtol = fixed_mtol
        this_index[ind].logmass = $
           alog10((w1_lum/lsun_3p4*this_index.mtol)[ind])
        this_index[ind].method_mtol = 'FIXEDUL'
     endif     

;    WE HAVE NUV AND W4
     ind = where(this_index.has_wise1 and $
                 this_index.has_nuv and $
                 this_index.has_fuv eq 0 and $
                 this_index.has_wise4 and $
                 this_index.flux_wise1 gt 0 and $
                 isul_nuvw4 eq 0, ct)
     if ct gt 0 then begin
        this_index[ind].mtol = $
           lookup_mtol(ssfrlike=ssfrlike_nuvw4[ind])
        this_index[ind].logmass = $
           alog10(w1_lum[ind]/lsun_3p4*this_index[ind].mtol)        
        this_index[ind].method_mtol = 'SSFRLIKE'
     endif     

     ind = where(this_index.has_wise1 and $
                 this_index.has_nuv and $
                 this_index.has_fuv eq 0 and $
                 this_index.has_wise4 and $
                 this_index.flux_wise1 gt 0 and $
                 isul_nuvw4 eq 1, ct)
     if ct gt 0 then begin
        this_index[ind].mtol = fixed_mtol
        this_index[ind].logmass = $
           alog10((w1_lum/lsun_3p4*this_index.mtol)[ind])
        this_index[ind].method_mtol = 'FIXEDUL'
     endif     

;    UNCERTAINTY
     this_index.e_logmass = $
        sqrt(0.1^2 + (alog10(this_index.std_flux_wise1/this_index.flux_wise1+1.))^2)

;    WISE1 IS LESS THAN 0
     ind = where(this_index.flux_wise1 le 0.0, ct)
     if ct gt 0 then begin        
        this_index[ind].mtol = fixed_mtol
        this_index[ind].method_mtol = 'FIXEDUL'
        this_index[ind].logmass = !values.f_nan
        this_index[ind].e_logmass = alog10((3.*e_w1_lum/lsun_3p4*this_index.mtol)[ind])
     endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%     
; Star formation rates
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    has W4 and FUV

;    ... and is not an upper limit
     ind = where(this_index.has_fuv and $
                 this_index.has_wise4 and $
                 isul_fuvw4 eq 0, ct)    
     if ct gt 0 then begin
        this_index[ind].logsfr = alog10(sfr_fuvw4[ind])
        this_index[ind].e_logsfr = $
           sqrt(alog10((e_sfr_fuvw4/sfr_fuvw4+1.0)[ind])^2+0.2^2)
        this_index[ind].method_sfr = 'FUV+WISE4'
     endif

;    ... and is an upper limit
     ind = where(this_index.has_fuv and $
                 this_index.has_wise4 and $
                 isul_fuvw4 eq 1, ct)    
     if ct gt 0 then begin
        this_index[ind].logsfr = !values.f_nan
        this_index[ind].e_logsfr = alog10(3.0*e_sfr_fuvw4[ind])
        this_index[ind].method_sfr = 'FUV+WISE4'
     endif

;    has W4 and NUV

;    ... and is not an upper limit

     ind = where(this_index.has_nuv and $
                 this_index.has_fuv eq 0 and $
                 this_index.has_wise4 and $
                 isul_nuvw4 eq 0, ct)    
     if ct gt 0 then begin
        this_index[ind].logsfr = alog10(sfr_nuvw4[ind])
        this_index[ind].e_logsfr = $
           sqrt(alog10((e_sfr_nuvw4/sfr_nuvw4+1.0)[ind])^2+0.2^2)
        this_index[ind].method_sfr = 'NUV+WISE4'
     endif

;    ... and is an upper limit

     ind = where(this_index.has_nuv and $
                 this_index.has_fuv eq 0 and $
                 this_index.has_wise4 and $
                 isul_nuvw4 eq 1, ct)    
     if ct gt 0 then begin
        this_index[ind].logsfr = !values.f_nan
        this_index[ind].e_logsfr = alog10(3.0*e_sfr_nuvw4[ind])
        this_index[ind].method_sfr = 'NUV+WISE4'
     endif

;    has only W4

;    ... and is not an upper limit

     ind = where(this_index.has_nuv eq 0  and $
                 this_index.has_fuv eq 0 and $
                 this_index.has_wise4 and $
                 isul_justw4 eq 0, ct)    
     if ct gt 0 then begin
        this_index[ind].logsfr = alog10(sfr_justw4[ind])
        this_index[ind].e_logsfr = $
           sqrt(alog10((e_sfr_justw4/sfr_justw4+1.0)[ind])^2+0.2^2)
        this_index[ind].method_sfr = 'WISE4'
     endif

;    ... and is an upper limit

     ind = where(this_index.has_nuv eq 0  and $
                 this_index.has_fuv eq 0 and $
                 this_index.has_wise4 and $
                 isul_justw4 eq 1, ct)    
     if ct gt 0 then begin
        this_index[ind].logsfr = !values.f_nan
        this_index[ind].e_logsfr = alog10(3.0*e_sfr_justw4[ind])
        this_index[ind].method_sfr = 'WISE4'
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Offset from the main sequence
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
     xnorm = 10.0
     ssfr_pred = $
        (-0.34)*(this_index.logmass-xnorm) - $
        10.11
     ssfr_measure = $
        this_index.logsfr - this_index.logmass
     resid = $
        ssfr_measure - ssfr_pred
     ind = where((this_index.has_nuv ne 0 or $
                  this_index.has_fuv ne 0) and $
                 this_index.has_wise4 and $
                 (ssfr_measure ge -11) $
                 , ct)         
     if ct gt 0 then begin
        this_index.deltams = !values.f_nan
        this_index[ind].deltams = $
           resid[ind]
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Save the results
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if jj eq 0 then index15 = this_index
     if jj eq 1 then index7p5 = this_index

  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  mwrfits, index15 $
           , '../measurements/delivery_index_gauss15.fits' $
           , /create

  mwrfits, index7p5 $
           , '../measurements/delivery_index_gauss7p5.fits' $
           , /create  

  stop

end
