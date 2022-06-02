pro z0mgs_photometry $
   , start=start $
   , stop = stop $
   , just_pgc = just_pgc $
   , skip_pgc = skip_pgc $
   , pause = pause $
   , verbose = verbose

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  out_dir = '../measurements/'
  sample_dir = '../measurements/samples/'
  atlas_dir = '../delivery/'
  mask_dir = '../masks/'

  index7p5 = mrdfits('../measurements/delivery_index_gauss7p5.fits',1,h)
  index15 = mrdfits('../measurements/delivery_index_gauss15.fits',1,h)
  n_pgc = n_elements(index7p5)
  s = gal_data(pgc=index7p5.pgc, /full)

  if n_elements(start) eq 0 then start = 0
  if n_elements(stop) eq 0 then stop = n_pgc-1

  bands = ['FUV','NUV','WISE1','WISE2','WISE3','WISE4']
  n_bands = n_elements(bands)
  
  spacing_deg = 7.5/3600.
  spacing_rad = spacing_deg*!dtor  
  min_limit = 30./3600.

  oversamp = 4.0                ; would be 5.23 for hex

  rej_thresh_wise = 0.2
  rej_thresh_galex = 0.75

  wise1_corr = 10.^(-0.034/2.5)
  wise2_corr = 10.^(-0.041/2.5)
  wise3_corr = 10.^(+0.03/2.5)
  wise4_corr = 10.^(-0.0294/2.5)

  print, 'WISE1 correction', wise1_corr
  print, 'WISE2 correction', wise2_corr
  print, 'WISE3 correction', wise3_corr
  print, 'WISE4 correction', wise4_corr

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE OUTPUT STRUCTURE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  limit_ra = [0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 4.0]
  n_limit = n_elements(limit_ra)

  nan = !values.f_nan
  empty = $
     {pgc:0L $
      , gl_deg: nan $
      , gb_deg: nan $
      , band: '' $
      , present: 0B $
      , rms: nan $
      , std: nan $
      , true_posang: nan $
      , adopted_posang: nan $
      , true_incl: nan $
      , adopted_incl: nan $
      , true_r25: nan $
      , fiducial_limit: nan $
      , limit_used: limit_ra $
      , mean: limit_ra*nan $
      , med: limit_ra*nan $
      , unc_stat: limit_ra*nan $
      , unc_conf:limit_ra* nan $
      , rejected_flux: limit_ra*nan $
     }
  phot_gauss15 = replicate(empty, n_pgc, n_bands)
  phot_gauss7p5 = phot_gauss15

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  for jj = 0, 1 do begin

     if jj eq 0 then begin
        res_str = 'gauss15'
        this_index = index15
     endif else begin
        res_str = 'gauss7p5'
        this_index = index7p5
     endelse

     for ii = start, stop do begin

        counter, ii, n_pgc, 'Photometry for galaxy '
        
        this_pgc = s[ii].pgc
        if this_pgc eq 0 then begin
           print, this_index[ii].pgc, ' not in galbase. Skipping.'
           continue
        endif

        if this_index[ii].has_wise1 eq 0 then begin
           print, this_pgc, ' lacks WISE1. Should not.'
           continue
        endif

        if n_elements(just_pgc) gt 0 then begin
           if total(just_pgc eq this_pgc) eq 0 then $
              continue
        endif

        if n_elements(skip_pgc) gt 0 then begin
           if total(just_pgc eq skip_pgc) gt 0 then $
              continue
        endif

        pgc_name = strcompress('PGC'+str(this_pgc), /rem)

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    FILENAMES
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        w1_root = atlas_dir+pgc_name+'_w1_'+res_str
        w2_root = atlas_dir+pgc_name+'_w2_'+res_str
        w3_root = atlas_dir+pgc_name+'_w3_'+res_str
        w4_root = atlas_dir+pgc_name+'_w4_'+res_str
        nuv_root = atlas_dir+pgc_name+'_nuv_'+res_str
        fuv_root = atlas_dir+pgc_name+'_fuv_'+res_str
        rgrid_file = atlas_dir+pgc_name+'_'+res_str+'_rgrid.fits'
        galmask_file = atlas_dir+pgc_name+'_'+res_str+'_galaxies.fits'

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    READ FILES
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        w1 = readfits(w1_root+'.fits', w1_hdr, /silent)
        if this_index[ii].has_wise2 then w2 = readfits(w2_root+'.fits', w2_hdr, /silent)
        if this_index[ii].has_wise3 then w3 = readfits(w3_root+'.fits', w3_hdr, /silent)        
        if this_index[ii].has_wise4 then w4 = readfits(w4_root+'.fits', w4_hdr, /silent)
        if this_index[ii].has_nuv then nuv = readfits(nuv_root+'.fits', nuv_hdr, /silent)
        if this_index[ii].has_fuv then fuv = readfits(fuv_root+'.fits', fuv_hdr, /silent)

        rad_deg = readfits(rgrid_file, rgrid_hdr, /silent)
        gal_mask = readfits(galmask_file, galmask_hdr, /silent)

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    WORK OUT THE ADOPTED ORIENTATION
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        area_sr = (sxpar(w1_hdr, 'CD1_1')*!dtor)^2        
        true_r25 = s[ii].r25_deg
        this_fiducial_limit = sxpar(rgrid_hdr, 'FIDRAD')

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    LOOP OVER THE BANDS
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        for kk = 0, n_bands-1 do begin

           if keyword_set(verbose) then $
              print, "Band # "+str(kk)
           
           this_phot = empty

;          Metadata
           this_phot.pgc = s[ii].pgc
           this_phot.gl_deg = s[ii].gl_deg
           this_phot.gb_deg = s[ii].gb_deg
           this_phot.band = bands[kk]
           this_phot.present = 1B
           this_phot.true_r25 = true_r25
           this_phot.fiducial_limit = this_fiducial_limit
           
           if kk eq 0 then begin
              if this_index[ii].has_fuv eq 0 then begin
                 if res_str eq 'gauss15' then begin
                    phot_gauss15[ii,kk] = this_phot 
                 endif
                 
                 if res_str eq 'gauss7p5' then begin
                    phot_gauss7p5[ii,kk] = this_phot 
                 endif
                 
                 continue
              endif
              y = fuv
              h = fuv_hdr
              stars = readfits(fuv_root+'_stars.fits', star_hdr,/silent)
              corr = 1.0
           endif
           if kk eq 1 then begin
              if this_index[ii].has_nuv eq 0 then begin
                 if res_str eq 'gauss15' then begin
                    phot_gauss15[ii,kk] = this_phot 
                 endif
                 
                 if res_str eq 'gauss7p5' then begin
                    phot_gauss7p5[ii,kk] = this_phot 
                 endif
                 
                 continue
              endif
              y = nuv
              h = nuv_hdr
              stars = readfits(nuv_root+'_stars.fits', star_hdr,/silent)
              corr = 1.0
           endif
           if kk eq 2 then begin
              y = w1
              h = w1_hdr
              stars = readfits(w1_root+'_stars.fits', star_hdr,/silent)
              corr = wise1_corr
           endif
           if kk eq 3 then begin
              if this_index[ii].has_wise2 eq 0 then begin
                 if res_str eq 'gauss15' then begin
                    phot_gauss15[ii,kk] = this_phot 
                 endif
                 
                 if res_str eq 'gauss7p5' then begin
                    phot_gauss7p5[ii,kk] = this_phot 
                 endif

                 continue
              endif
              y = w2
              h = w2_hdr
              stars = readfits(w2_root+'_stars.fits', star_hdr,/silent)
              corr = wise2_corr
           endif
           if kk eq 4 then begin
              if this_index[ii].has_wise3 eq 0 then begin
                 if res_str eq 'gauss15' then begin
                    phot_gauss15[ii,kk] = this_phot 
                 endif
                 
                 if res_str eq 'gauss7p5' then begin
                    phot_gauss7p5[ii,kk] = this_phot 
                 endif
                 
                 continue
              endif
              y = w3
              h = w3_hdr
              stars = readfits(w3_root+'_stars.fits', star_hdr,/silent)
              corr = wise3_corr
           endif
           if kk eq 5 then begin
              if this_index[ii].has_wise4 eq 0 then begin

                 if res_str eq 'gauss15' then begin
                    phot_gauss15[ii,kk] = this_phot 
                 endif
                 
                 if res_str eq 'gauss7p5' then begin
                    phot_gauss7p5[ii,kk] = this_phot 
                 endif
                 
                 continue
              endif
              y = w4
              h = w4_hdr
              stars = readfits(w4_root+'_stars.fits', star_hdr,/silent)
              corr = wise4_corr
           endif

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    CONVERT UNITS AND STORE METADATA IN THE STRUCTURE
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

;          MJy -> Jy
           y *= 1d6

;          Apply correction factor
           y = y*corr

;          Noise
           this_rms = sxpar(h, 'RMS')*1d6 ; MJy -> Jy
           this_std = sxpar(h, 'STDDEV')*1d6 ; MJy -> Jy

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    CARRY OUT THE PHOTOMETRY
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
           
;          Identify indices to carry out photometry on           
           phot_ind = where(rad_deg le max(limit_ra)*this_fiducial_limit, phot_ct)
           if phot_ct eq 0 then begin           
              if s[ii].pgc eq 0 then stop
              print, "No valid data for PGC "+str(s[ii].pgc)

              if res_str eq 'gauss15' then begin
                 phot_gauss15[ii,kk] = this_phot 
              endif
              
              if res_str eq 'gauss7p5' then begin
                 phot_gauss7p5[ii,kk] = this_phot 
              endif

              continue
           endif                     
                      
           if total(finite(y[phot_ind])) eq 0 then begin
              print, "No valid data for PGC "+str(s[ii].pgc)

              if res_str eq 'gauss15' then begin
                 phot_gauss15[ii,kk] = this_phot 
              endif
              
              if res_str eq 'gauss7p5' then begin
                 phot_gauss7p5[ii,kk] = this_phot 
              endif

              continue
           endif

;          Two profiles: one to count pixels, the other to count flux

;          ... use the masks to set the image to NaN in places where
;          bright stars cause a problem.

;bad_ind = where(((bright_stars eq 1 or bright_stars eq 11) and (rad_deg ge 0.5*this_fiducial_limit)) or $
;                (gal_mask eq 1 and rad_deg ge this_fiducial_limit) or $
;                ((found_stars eq 1 or found_stars eq 11) and (rad_deg ge 0.5*this_fiducial_limit)) $
;                ,bad_ct)
           bad_ind = $
              where(((stars eq 1) and (rad_deg ge 0.5*this_fiducial_limit)) or $
                    (gal_mask eq 1 and rad_deg ge this_fiducial_limit) $
                    , bad_ct)
           if bad_ct gt 0 then y[bad_ind] = !values.f_nan

           ;loadct, 0
           ;disp, (y), /sq       
           ;contour, finite(y) eq 0, /overplot, lev=[1], color=cgcolor('red')
           ;ch = get_kbrd(1)

           bins_for_flux = $
              bin_data(rad_deg[phot_ind], (y)[phot_ind], /nan $
                       , xmin=-0.5*spacing_deg $
                       , xmax=max(limit_ra)*this_fiducial_limit $
                       , binsize=spacing_deg)

;          ... there will be no NaNs in this version, use it for the
;          area of each ring. Effectively, we interpolate.

           bins_for_area = $
              bin_data(rad_deg[phot_ind], rad_deg[phot_ind]*0.0+1.0, /nan $
                       , xmin=-0.5*spacing_deg $
                       , xmax=max(limit_ra)*this_fiducial_limit $
                       , binsize=spacing_deg)
           
           for mm = 0, n_limit-1 do begin

              this_limit = limit_ra[mm]*this_fiducial_limit
              bin_mask = bins_for_area.xmid le this_limit

              this_phot.mean[mm] = $
                 total(bins_for_flux.ymean*area_sr*bins_for_area.counts*bin_mask, /nan)
              this_phot.med[mm] = $
                 total(bins_for_flux.ymed*area_sr*bins_for_area.counts*bin_mask, /nan)  

              indep_meas = (total(bins_for_flux.counts*bin_mask)/oversamp)
              this_phot.unc_stat[mm] = this_rms * sqrt(indep_meas) * area_sr
              this_phot.unc_conf[mm] = this_std * sqrt(indep_meas) * area_sr

           endfor

           if res_str eq 'gauss15' then begin
              phot_gauss15[ii,kk] = this_phot 
           endif

           if res_str eq 'gauss7p5' then begin
              phot_gauss7p5[ii,kk] = this_phot 
           endif

        endfor

     endfor

  endfor

  if keyword_set(do_pause) then begin
     stop
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  save, file='../measurements/z0mgs_photometry.idl', phot_gauss15, phot_gauss7p5
  
  mwrfits, phot_gauss15 $
           , '../measurements/z0mgs_photometry_gauss15.fits' $
           , /create
  
  mwrfits, phot_gauss7p5 $
           , '../measurements/z0mgs_photometry_gauss7p5.fits' $
           , /create

  end
