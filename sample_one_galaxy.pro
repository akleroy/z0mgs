pro sample_one_galaxy $
   , pgc_name = pgc_name $
   , gal_data = gal_data $
   , index = index $
   , data_dir = data_dir $
   , outfile=outfile $
   , spacing=spacing $
   , beam_size=beam_size $
   , empty=empty

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINITIONS AND DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(spacing) eq 0 then $
     bin_size = 7.5

  if n_elements(beam_size) eq 0 then $
     beam_size = 15.

  max_bands = 6

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; EMPTY STRUCTURE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  nan = !values.f_nan
  sample = $
     { $
     pgc_name:'' $
     , ra_deg: nan $
     , dec_deg: nan $
     , beam_arcsec: nan $
     , fuv: nan $
     , fuv_rms: nan $
     , fuv_std: nan $
     , fuv_rej: nan $
     , fuv_wt: nan $
     , afuv: nan $
     , nuv: nan $
     , nuv_rms: nan $
     , nuv_std: nan $
     , nuv_rej: nan $
     , nuv_wt: nan $
     , anuv: nan $
     , wise1: nan $
     , wise1_rms: nan $
     , wise1_std: nan $
     , wise1_rej: nan $
     , wise2: nan $
     , wise2_rms: nan $
     , wise2_std: nan $
     , wise2_rej: nan $
     , wise3: nan $
     , wise3_rms: nan $
     , wise3_std: nan $
     , wise3_rej: nan $
     , wise4: nan $
     , wise4_rms: nan $
     , wise4_std: nan $
     , wise4_rej: nan $
     }

  empty = sample
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; USE THE INDEX FILE TO READ THE RELEVANT DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  found_ref = 0B
  
  if index.has_fuv eq 0 and $
     index.has_nuv eq 0 and $
     index.has_wise1 eq 0 and $
     index.has_wise2 eq 0 and $
     index.has_wise3 eq 0 and $
     index.has_wise4 eq 0 then begin
     print, "No files expected. Returning."
     return
  endif

  n_bands = 0

  if index.has_fuv then begin
     fuv_file = data_dir + pgc_name + '_fuv.fits'     
     test = file_search(fuv_file, count=fuv_ct)
     if fuv_ct eq 0 then begin
        print, "No FUV file found, but one expected. Stopping."
        stop
     endif
     fuv = readfits(fuv_file, fuv_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = fuv_hdr
        ref_map = fuv
     endif
     fuv_weight = readfits(data_dir + pgc_name + '_fuv_weight.fits')
     fuv_rej = readfits(data_dir + pgc_name + '_fuv_rejected.fits')     
     n_bands += 1
  endif

  if index.has_nuv then begin
     nuv_file = data_dir + pgc_name + '_nuv.fits'
     test = file_search(nuv_file, count=nuv_ct)
     if nuv_ct eq 0 then begin
        print, "No NUV file found, but one expected. Stopping."
        stop
     endif
     nuv = readfits(nuv_file, nuv_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = nuv_hdr
        ref_map = nuv
     endif
     nuv_weight = readfits(data_dir + pgc_name + '_nuv_weight.fits')
     nuv_rej = readfits(data_dir + pgc_name + '_nuv_rejected.fits')
     n_bands += 1
  endif

  if index.has_wise1 then begin
     wise1_file = data_dir + pgc_name + '_w1.fits'
     test = file_search(wise1_file, count=wise1_ct)
     if wise1_ct eq 0 then begin
        print, "No WISE1 file found, but one expected. Stopping."
        stop
     endif
     wise1 = readfits(wise1_file, wise1_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = wise1_hdr
        ref_map = wise1
     endif
     wise1_rej = readfits(data_dir + pgc_name + '_w1_rejected.fits')
     n_bands += 1
  endif

  if index.has_wise2 then begin
     wise2_file = data_dir + pgc_name + '_w2.fits'
     test = file_search(wise2_file, count=wise2_ct)
     if wise2_ct eq 0 then begin
        print, "No WISE2 file found, but one expected. Stopping."
        stop
     endif
     wise2 = readfits(wise2_file, wise2_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = wise2_hdr
        ref_map = wise2
     endif
     wise2_rej = readfits(data_dir + pgc_name + '_w2_rejected.fits')
     n_bands += 1
  endif

  if index.has_wise3 then begin
     wise3_file = data_dir + pgc_name + '_w3.fits'
     test = file_search(wise3_file, count=wise3_ct)
     if wise3_ct eq 0 then begin
        print, "No WISE3 file found, but one expected. Stopping."
        stop
     endif
     wise3 = readfits(wise3_file, wise3_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = wise3_hdr
        ref_map = wise3
     endif
     wise3_rej = readfits(data_dir + pgc_name + '_w3_rejected.fits')
     n_bands += 1
  endif

  if index.has_wise4 then begin
     wise4_file = data_dir + pgc_name + '_w4.fits'
     test = file_search(wise4_file, count=wise4_ct)
     if wise4_ct eq 0 then begin
        print, "No WISE4 file found, but one expected. Stopping."
        stop
     endif
     wise4 = readfits(wise4_file, wise4_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = wise4_hdr
        ref_map = wise4
     endif
     wise4_rej = readfits(data_dir + pgc_name + '_w4_rejected.fits')
     n_bands += 1
  endif

  if found_ref eq 0 then begin
     print, "I should have a reference header at this point. Stopping."
     stop
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD ASTROMETRY AND THE HEX GRID
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  pix_deg = get_pixel_scale(ref_hdr)
  sz = size(ref_map)
  max_rad = sz[1]*pix_deg*2.0

  hex_grid $
     , ctr_x = gal_data.ra_deg $
     , ctr_y = gal_data.dec_deg $
     , spacing = spacing/3600. $
     , /radec $
     , xout = samp_ra $
     , yout = samp_dec $
     , r_limit = max_rad $
     , /center

  adxy, ref_hdr, samp_ra, samp_dec, samp_x, samp_y

  keep = where(samp_x ge 0 and $
               samp_y ge 0 and $
               samp_x lt sz[1] and $
               samp_y lt sz[2], keep_ct)
  
  if keep_ct eq 0B then begin
     message, 'No sampling points survive inside mask. Returning.', /info
     return
  endif
  
  samp_ra = samp_ra[keep]
  samp_dec = samp_dec[keep]
  samp_x = samp_x[keep]
  samp_y = samp_y[keep]
  n_pts = keep_ct

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FILL IN HEADER
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sample.pgc_name = pgc_name
  sample.beam_arcsec = beam_size
  samples = replicate(sample, n_pts)
  samples.ra_deg = samp_ra
  samples.dec_deg = samp_dec
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; STEP THROUGH BANDS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
  if index.has_fuv then begin
     samples.fuv = interpolate(fuv, samp_x, samp_y $
                               , cubic=-0.5, missing=0.0, /double)
     samples.fuv_wt = interpolate(fuv_weight, samp_x, samp_y $
                                  , cubic=-0.5, missing=0.0, /double)
     samples.fuv_rej = interpolate(fuv_rej, samp_x, samp_y $
                                   , cubic=-0.5, missing=0.0, /double)
     samples.fuv_rms = index.flatrms_fuv * $
                       1.d/sqrt(samples.fuv_wt)
     samples.fuv_std = index.std_fuv
     samples.afuv = index.afuv
  endif

  if index.has_nuv then begin
     samples.nuv = interpolate(nuv, samp_x, samp_y $
                               , cubic=-0.5, missing=0.0, /double)
     samples.nuv_wt = interpolate(nuv_weight, samp_x, samp_y  $
                                  , cubic=-0.5, missing=0.0, /double)
     samples.nuv_rej = interpolate(nuv_rej, samp_x, samp_y  $
                                   , cubic=-0.5, missing=0.0, /double)
     samples.nuv_rms = index.flatrms_nuv * $
                       1.d/sqrt(samples.nuv_wt)
     samples.nuv_std = index.std_nuv
     samples.anuv = index.anuv
  endif

  if index.has_wise1 then begin
     samples.wise1 = interpolate(wise1, samp_x, samp_y $
                                 , cubic=-0.5, missing=0.0, /double)
     samples.wise1_rej = interpolate(wise1_rej, samp_x, samp_y $
                                     , cubic=-0.5, missing=0.0, /double)
     samples.wise1_rms = index.rms_wise1
     samples.wise1_std = index.std_wise1
  endif

  if index.has_wise2 then begin
     samples.wise2 = interpolate(wise2, samp_x, samp_y $
                                 , cubic=-0.5, missing=0.0, /double)
     samples.wise2_rej = interpolate(wise2_rej, samp_x, samp_y $
                                     , cubic=-0.5, missing=0.0, /double)
     samples.wise2_rms = index.rms_wise2
     samples.wise2_std = index.std_wise2
  endif

  if index.has_wise3 then begin
     samples.wise3 = interpolate(wise3, samp_x, samp_y $
                                 , cubic=-0.5, missing=0.0, /double)
     samples.wise3_rej = interpolate(wise3_rej, samp_x, samp_y $
                                     , cubic=-0.5, missing=0.0, /double)
     samples.wise3_rms = index.rms_wise3
     samples.wise3_std = index.std_wise3
  endif

  if index.has_wise4 then begin
     samples.wise4 = interpolate(wise4, samp_x, samp_y $
                                 , cubic=-0.5, missing=0.0, /double)
     samples.wise4_rej = interpolate(wise4_rej, samp_x, samp_y $
                                     , cubic=-0.5, missing=0.0, /double)
     samples.wise4_rms = index.rms_wise4
     samples.wise4_std = index.std_wise4
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SAVE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(outfile) gt 0 then begin    
     mwrfits, samples, outfile, /create
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RETURN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  return

end
