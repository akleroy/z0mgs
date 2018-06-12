function index_one_galaxy $
   , pgc_name=pgc_name $
   , gal_data=gal_data $
   , atlas_dir=atlas_dir $
   , empty=empty

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD THE STRUCTURE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  nan = !values.f_nan
  index = $
     { $
     pgc_name:'' $
     , pgc: 0L $
     , ra_deg: nan $
     , dec_deg: nan $
     , gl_deg: nan $
     , gb_deg: nan $
;    PRESENCE OF DATA
     , has_fuv: 0B $
     , has_nuv: 0B $
     , has_wise1: 0B $
     , has_wise2: 0B $
     , has_wise3: 0B $
     , has_wise4: 0B $     
;    DEPTH OF GALEX
     , time_fuv: nan $
     , time_nuv: nan $
     , afuv: nan $
     , anuv: nan $
;    NOISE IN IMAGE
     , rms_fuv: nan $
     , std_fuv: nan $
     , flatrms_fuv: nan $
     , rms_nuv: nan $
     , std_nuv: nan $
     , flatrms_nuv: nan $
     , rms_wise1: nan $
     , std_wise1: nan $
     , rms_wise2: nan $
     , std_wise2: nan $
     , rms_wise3: nan $
     , std_wise3: nan $
     , rms_wise4: nan $
     , std_wise4: nan $
;    FLAGS
     , psf_effects_fuv: 0B, 
     , psf_effects_nuv: 0B, 
     , psf_effects_wise1: 0B, 
     , psf_effects_wise2: 0B, 
     , psf_effects_wise3: 0B, 
     , psf_effects_wise4: 0B, 
     , overlap_star: 0B,
     , overlap_galaxy: 0B,
     , pathologies_fuv: '',
     , pathologies_nuv: '',
     , pathologies_wise1: '',
     , pathologies_wise2: '',
     , pathologies_wise3: '',
     , pathologies_wise4: '',
     , frac_blank_fuv: 0.0,
     , frac_blank_nuv: 0.0,
     , frac_blank_wise1: 0.0,
     , frac_blank_wise2: 0.0,
     , frac_blank_wise3: 0.0,
     , frac_blank_wise4: 0.0,
;    PHOTOMETRY
     , flux_fuv: nan $
     , rms_flux_fuv: nan $
     , std_flux_fuv: nan $
     , flux_nuv: nan $
     , rms_flux_nuv: nan $
     , std_flux_nuv: nan $
     , flux_wise1: nan $
     , rms_flux_wise1: nan $
     , std_flux_wise1: nan $
     , flux_wise2: nan $
     , rms_flux_wise2: nan $
     , std_flux_wise2: nan $
     , flux_wise3: nan $
     , rms_flux_wise3: nan $
     , std_flux_wise3: nan $
     , flux_wise4: nan $
     , rms_flux_wise4: nan $
     , std_flux_wise4: nan $
     }

  if keyword_set(empty) then begin
     return, index
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FILL IN THE METADATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  index.pgc_name = pgc_name
  index.pgc = gal_data.pgc
  index.ra_deg = gal_data.ra_deg
  index.dec_deg = gal_data.dec_deg
  index.gl_deg = gal_data.gl_deg
  index.gb_deg = gal_data.gb_deg

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GO FILE BY FILE AND CHECK EXISTENCE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  wise1_fname = atlas_dir+pgc_name+'_w1.fits'
  wise2_fname = atlas_dir+pgc_name+'_w2.fits'
  wise3_fname = atlas_dir+pgc_name+'_w3.fits'
  wise4_fname = atlas_dir+pgc_name+'_w4.fits'
  fuv_fname = atlas_dir+pgc_name+'_fuv.fits'
  nuv_fname = atlas_dir+pgc_name+'_nuv.fits'

  test = file_search(wise1_fname, count=w1_count)
  test = file_search(wise2_fname, count=w2_count)
  test = file_search(wise3_fname, count=w3_count)
  test = file_search(wise4_fname, count=w4_count)
  test = file_search(fuv_fname, count=fuv_count)
  test = file_search(nuv_fname, count=nuv_count)
  
  if w1_count eq 1 then index.has_wise1 = 1B
  if w2_count eq 1 then index.has_wise2 = 1B
  if w3_count eq 1 then index.has_wise3 = 1B
  if w4_count eq 1 then index.has_wise4 = 1B
  if fuv_count eq 1 then index.has_fuv = 1B
  if nuv_count eq 1 then index.has_nuv = 1B

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PARSE THE HEADERS FOR STATISTICS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if index.has_wise1 then begin
     hdr = headfits(wise1_fname)
     index.rms_wise1 = sxpar(hdr, 'MADALL')
     index.std_wise1 = sxpar(hdr, 'STDALL')
     index.rej_wise1 = sxpar(hdr, 'REJALL')
     index.max_wise1 = sxpar(hdr, 'MAXALL')
  endif

  if index.has_wise2 then begin
     hdr = headfits(wise2_fname)
     index.rms_wise2 = sxpar(hdr, 'MADALL')
     index.std_wise2 = sxpar(hdr, 'STDALL')
     index.rej_wise2 = sxpar(hdr, 'REJALL')
     index.max_wise2 = sxpar(hdr, 'MAXALL')
  endif

  if index.has_wise3 then begin
     hdr = headfits(wise3_fname)
     index.rms_wise3 = sxpar(hdr, 'MADALL')
     index.std_wise3 = sxpar(hdr, 'STDALL')
     index.rej_wise3 = sxpar(hdr, 'REJALL')
     index.max_wise3 = sxpar(hdr, 'MAXALL')
  endif

  if index.has_wise4 then begin
     hdr = headfits(wise4_fname)
     index.rms_wise4 = sxpar(hdr, 'MADALL')
     index.std_wise4 = sxpar(hdr, 'STDALL')
     index.rej_wise4 = sxpar(hdr, 'REJALL')
     index.max_wise4 = sxpar(hdr, 'MAXALL')
  endif

  if index.has_fuv then begin
     hdr = headfits(fuv_fname)
     index.time_fuv = sxpar(hdr, 'MEANINT')
     index.rms_fuv = sxpar(hdr, 'MADALL')
     index.std_fuv = sxpar(hdr, 'STDALL')
     index.rej_fuv = sxpar(hdr, 'REJALL')
     index.max_fuv = sxpar(hdr, 'MAXALL')
     index.flatrms_fuv = sxpar(hdr, 'FLATMADA')
     index.afuv = sxpar(hdr, 'AFUV')
  endif

  if index.has_nuv then begin
     hdr = headfits(nuv_fname)
     index.time_nuv = sxpar(hdr, 'MEANINT')
     index.rms_nuv = sxpar(hdr, 'MADALL')
     index.std_nuv = sxpar(hdr, 'STDALL')
     index.rej_nuv = sxpar(hdr, 'REJALL')
     index.max_nuv = sxpar(hdr, 'MAXALL')
     index.flatrms_nuv = sxpar(hdr, 'FLATMADA')
     index.anuv = sxpar(hdr, 'ANUV')
  endif

  return, index

end
