pro z0mgs_photometry $
   , start=start $
   , stop = stop

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  out_dir = '../measurements/'
;  sample_dir = '../measurements/samples/'
  atlas_dir = '../delivery/'

  index = mrdfits('../measurements/delivery_index.fits',1,h)
  n_pgc = n_elements(index)
  s = gal_data(pgc=index.pgc)

  if n_elements(start) eq 0 then start = 0
  if n_elements(stop) eq 0 then stop = n_pgc-1

  bands = ['FUV','NUV','WISE1','WISE2','WISE3','WISE4']
  n_bands = n_elements(bands)
  
  spacing_deg = 7.5/3600.
  spacing_rad = spacing_deg*!dtor  
  min_limit = 30./3600.

  oversamp = 4.0                ;5.23

  rej_thresh_wise = 0.2
  rej_thresh_galex = 0.75

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
  phot = replicate(empty, n_pgc, n_bands)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = start, stop do begin

     counter, ii, n_pgc, 'Photometry for galaxy '
     
     this_pgc = s[ii].pgc
     if this_pgc eq 0 then begin
        print, index[ii].pgc, ' not in galbase. Skipping.'
        continue
     endif

     if index[ii].has_wise1 eq 0 then begin
        print, this_pgc, ' lacks WISE1. Should not.'
        continue
     endif

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    LOAD THE DATA
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     w1_file = atlas_dir+'PGC'+str(this_pgc)+'_w1.fits'
     w1_rejfile = atlas_dir+'PGC'+str(this_pgc)+'_w1_rejected.fits'

     w2_file = atlas_dir+'PGC'+str(this_pgc)+'_w2.fits'
     w2_rejfile = atlas_dir+'PGC'+str(this_pgc)+'_w2_rejected.fits'

     w3_file = atlas_dir+'PGC'+str(this_pgc)+'_w3.fits'
     w3_rejfile = atlas_dir+'PGC'+str(this_pgc)+'_w3_rejected.fits'

     w4_file = atlas_dir+'PGC'+str(this_pgc)+'_w4.fits'
     w4_rejfile = atlas_dir+'PGC'+str(this_pgc)+'_w4_rejected.fits'

     w1 = readfits(w1_file,w1hdr, /silent)
     w1_rej = readfits(w1_rejfile, /silent) ge rej_thresh_wise

     w2 = readfits(w2_file,w2hdr, /silent)
     w2_rej = readfits(w2_rejfile, /silent) ge rej_thresh_wise

     w3 = readfits(w3_file,w3hdr, /silent)
     w3_rej = readfits(w3_rejfile, /silent) ge rej_thresh_wise

     w4 = readfits(w4_file,w4hdr, /silent)
     w4_rej = readfits(w4_rejfile, /silent) ge rej_thresh_wise

     area_sr = (sxpar(w1hdr, 'CD1_1')*!dtor)^2

     if index[ii].has_nuv then begin
        nuv_file = atlas_dir+'PGC'+str(this_pgc)+'_nuv.fits'
        nuv_rejfile = atlas_dir+'PGC'+str(this_pgc)+'_nuv_rejected.fits'
        nuv = readfits(nuv_file, nuvhdr, /silent)
        nuv_rej = readfits(nuv_rejfile, /silent) ge rej_thresh_galex
     endif

     if index[ii].has_fuv then begin
        fuv_file = atlas_dir+'PGC'+str(this_pgc)+'_fuv.fits'
        fuv_rejfile = atlas_dir+'PGC'+str(this_pgc)+'_fuv_rejected.fits'
        fuv = readfits(fuv_file, fuvhdr, /silent)
        fuv_rej = readfits(fuv_rejfile, /silent) ge rej_thresh_galex
     endif

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    WORK OUT THE ADOPTED ORIENTATION
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     true_incl = s[ii].incl_deg
     this_incl = true_incl
     if finite(this_incl) eq 0 then $
        this_incl = 0.0
     if this_incl gt 60. then $
        this_incl = 60.

     true_posang = s[ii].posang_deg
     this_posang = s[ii].posang_deg
     if finite(this_posang) eq 0 then begin
        this_incl = 0.0
        this_posang = 0.0
     endif

     galpos = [this_posang, this_incl, s[ii].ra_deg, s[ii].dec_deg]

     make_axes, w1hdr, ri=ri, di=di
     deproject, ri, di, galpos, rgrid=rad_deg

     true_r25 = s[ii].r25_deg
     this_fiducial_limit = true_r25
     if finite(true_r25) eq 0 or true_r25 lt min_limit then $
        this_fiducial_limit = min_limit

     no_rej = rad_deg le true_r25

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    WORK OUT THE ADOPTED ORIENTATION
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     for jj = 0, n_bands-1 do begin
        
        this_phot = empty
        
        if jj eq 0 then begin
           if index[ii].has_fuv eq 0 then continue
           y = fuv
           r = fuv_rej
           h = fuvhdr
        endif
        if jj eq 1 then begin
           if index[ii].has_nuv eq 0 then continue
           y = nuv
           r = nuv_rej
           h = nuvhdr
        endif
        if jj eq 2 then begin
           y = w1
           r = w1_rej
           h = w1hdr
        endif
        if jj eq 3 then begin
           y = w2
           r = w2_rej
           h = w2hdr
        endif
        if jj eq 4 then begin
           y = w3
           r = w3_rej
           h = w3hdr
        endif
        if jj eq 5 then begin
           y = w4
           r = w4_rej
           h = w4hdr
        endif

        r = r and (no_rej eq 0)

;       MJy -> Jy
        y *= 1d6

        this_phot.pgc = s[ii].pgc
        this_phot.gl_deg = s[ii].gl_deg
        this_phot.gb_deg = s[ii].gb_deg
        this_phot.band = bands[jj]
        this_phot.present = 1B

        this_rms = sxpar(h, 'MADALL')*1d6 ; MJy -> Jy
        this_std = sxpar(h, 'STDALL')*1d6 ; MJy -> Jy

        this_phot.true_posang = true_posang
        this_phot.adopted_posang = this_posang

        this_phot.true_incl = true_incl
        this_phot.adopted_incl = this_incl

        this_phot.true_r25 = true_r25
        this_phot.fiducial_limit = this_fiducial_limit

;       Make the radial profiles once, for the maximum radius        
        phot_ind = where(rad_deg le max(limit_ra)*this_fiducial_limit, phot_ct)
        if phot_ct eq 0 then begin           
           if s[ii].pgc eq 0 then stop
           print, "No valid data for PGC "+str(s[ii].pgc)
           continue
        endif                     
           
        if total(finite(y[phot_ind])) eq 0 then $
           continue
           
        bins = $
           bin_data(rad_deg[phot_ind], y[phot_ind], /nan $
                    , xmin=-0.5*spacing_deg $
                    , xmax=max(limit_ra)*this_fiducial_limit $
                    , binsize=spacing_deg)

        rej_bins = $
           bin_data(rad_deg[phot_ind], (r*y)[phot_ind], /nan $
                    , xmin=-0.5*spacing_deg $
                    , xmax=max(limit_ra)*this_fiducial_limit $
                    , binsize=spacing_deg)
                
        for kk = 0, n_limit-1 do begin

           this_limit = limit_ra[kk]*this_fiducial_limit
           bin_mask = bins.xmid le this_limit

           this_phot.mean[kk] = $
              total(bins.ymean*area_sr*bins.counts*bin_mask, /nan)
           this_phot.med[kk] = $
              total(bins.ymed*area_sr*bins.counts*bin_mask, /nan)  
           this_phot.rejected_flux[kk] = $
              total(rej_bins.ymean*area_sr*rej_bins.counts*bin_mask, /nan)

           indep_meas = (total(bins.counts*bin_mask)/oversamp)
           this_phot.unc_stat[kk] = this_rms * sqrt(indep_meas) * area_sr
           this_phot.unc_conf[kk] = this_std * sqrt(indep_meas) * area_sr

        endfor

        phot[ii,jj] = this_phot

     endfor

  endfor

  stop
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  save, file='../measurements/z0mgs_photometry.idl', phot

  mwrfits, phot $
           , '../measurements/z0mgs_photometry.fits' $
           , /create

end
