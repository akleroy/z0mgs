function index_one_galaxy $
   , pgc_name=pgc_name $
   , gal_data=gal_data $
   , atlas_dir=atlas_dir $
   , empty=empty $
   , res_str=res_str

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; TUNING PARAMETERS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(res_str) eq 0 then $
     res_str = 'gauss15'

  sat_thresh_wise1 = 1e2
  sat_thresh_wise2 = 1e2
  sat_thresh_wise3 = 3e2
  sat_thresh_wise4 = 3e2
  sat_thresh_nuv = 1e6
  sat_thresh_fuv = 1e6

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
;    BEAM
     , resolution: res_str $
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
;    PHOTOMETRY
     , flux_fuv: nan $
     , rms_flux_fuv: nan $
     , std_flux_fuv: nan $
     , outer_flux_fuv: nan $
     , flux_nuv: nan $
     , rms_flux_nuv: nan $
     , std_flux_nuv: nan $
     , outer_flux_nuv: nan $
     , flux_wise1: nan $
     , rms_flux_wise1: nan $
     , std_flux_wise1: nan $
     , outer_flux_wise1: nan $
     , flux_wise2: nan $
     , rms_flux_wise2: nan $
     , std_flux_wise2: nan $
     , outer_flux_wise2: nan $
     , flux_wise3: nan $
     , rms_flux_wise3: nan $
     , std_flux_wise3: nan $
     , outer_flux_wise3: nan $
     , flux_wise4: nan $
     , rms_flux_wise4: nan $
     , std_flux_wise4: nan $
     , outer_flux_wise4: nan $
;    MASS-TO-LIGHT RATIO
     , mtol_w1: nan $
     , mtol_unc: nan $
;    NOISE IN IMAGE
     , rms_fuv: nan $
     , std_fuv: nan $
     , maskfrac_fuv: nan $
     , rms_nuv: nan $
     , std_nuv: nan $
     , maskfrac_nuv: nan $
     , rms_wise1: nan $
     , std_wise1: nan $
     , maskfrac_wise1: nan $
     , rms_wise2: nan $
     , std_wise2: nan $
     , maskfrac_wise2: nan $
     , rms_wise3: nan $
     , std_wise3: nan $
     , maskfrac_wise3: nan $
     , rms_wise4: nan $
     , std_wise4: nan $
     , maskfrac_wise4: nan $
;    FLAGS
     , sat_effects_fuv: 0B $
     , sat_effects_nuv: 0B $
     , sat_effects_wise1: 0B $ 
     , sat_effects_wise2: 0B $
     , sat_effects_wise3: 0B $
     , sat_effects_wise4: 0B $
     , overlap_galaxy_mask: 0.0 $
     , overlap_bright_mask: 0.0 $
     , overlap_found_mask: 0.0 $
     , overlap_custom_mask: 0.0 $
     , pathologies_fuv: '' $
     , pathologies_nuv: '' $
     , pathologies_wise1: '' $
     , pathologies_wise2: '' $
     , pathologies_wise3: '' $
     , pathologies_wise4: '' $
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
  
  rgrid_fname = atlas_dir+pgc_name+'_'+res_str+'_rgrid.fits'
  rgrid = readfits(rgrid_fname, rgrid_hdr,/silent)

  wise1_fname = atlas_dir+pgc_name+'_w1_'+res_str+'.fits'
  wise2_fname = atlas_dir+pgc_name+'_w2_'+res_str+'.fits'
  wise3_fname = atlas_dir+pgc_name+'_w3_'+res_str+'.fits'
  wise4_fname = atlas_dir+pgc_name+'_w4_'+res_str+'.fits'
  fuv_fname = atlas_dir+pgc_name+'_fuv_'+res_str+'.fits'
  nuv_fname = atlas_dir+pgc_name+'_nuv_'+res_str+'.fits'

  if file_test(wise1_fname) then index.has_wise1 = 1B
  if file_test(wise2_fname) then index.has_wise2 = 1B
  if file_test(wise3_fname) then index.has_wise3 = 1B
  if file_test(wise4_fname) then index.has_wise4 = 1B
  if file_test(fuv_fname) then index.has_fuv = 1B
  if file_test(nuv_fname) then index.has_nuv = 1B
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PARSE THE HEADERS FOR STATISTICS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  footprint = rgrid le 1.25 * sxpar(rgrid_hdr, 'FIDRAD')

  if index.has_wise1 then begin
     map = readfits(wise1_fname, hdr,/silent)
     index.rms_wise1 = sxpar(hdr, 'RMS')
     index.std_wise1 = sxpar(hdr, 'STDDEV')
     index.maskfrac_wise1 = sxpar(hdr, 'MASKFRAC')
     index.sat_effects_wise1 = total((map gt sat_thresh_wise1)*footprint) ge 1
     sxaddpar, hdr, 'SATURATE', index.sat_effects_wise1, 'Saturation expected in image?'
     writefits, wise1_fname, map, hdr
  endif

  if index.has_wise2 then begin
     map = readfits(wise2_fname, hdr,/silent)
     index.rms_wise2 = sxpar(hdr, 'RMS')
     index.std_wise2 = sxpar(hdr, 'STDDEV')
     index.maskfrac_wise2 = sxpar(hdr, 'MASKFRAC')
     index.sat_effects_wise2 = total((map gt sat_thresh_wise2)*footprint) ge 1
     sxaddpar, hdr, 'SATURATE', index.sat_effects_wise2, 'Saturation expected in image?'
     writefits, wise2_fname, map, hdr
  endif

  if index.has_wise3 then begin
     map = readfits(wise3_fname, hdr,/silent)
     index.rms_wise3 = sxpar(hdr, 'RMS')
     index.std_wise3 = sxpar(hdr, 'STDDEV')
     index.maskfrac_wise3 = sxpar(hdr, 'MASKFRAC')
     index.sat_effects_wise3 = total((map gt sat_thresh_wise3)*footprint) ge 1
     sxaddpar, hdr, 'SATURATE', index.sat_effects_wise3, 'Saturation expected in image?'
     writefits, wise3_fname, map, hdr
  endif

  if index.has_wise4 then begin
     map = readfits(wise4_fname, hdr,/silent)
     index.rms_wise4 = sxpar(hdr, 'RMS')
     index.std_wise4 = sxpar(hdr, 'STDDEV')
     index.maskfrac_wise4 = sxpar(hdr, 'MASKFRAC')
     index.sat_effects_wise4 = total((map gt sat_thresh_wise4)*footprint) ge 1
     sxaddpar, hdr, 'SATURATE', index.sat_effects_wise4, 'Saturation expected in image?'
     writefits, wise4_fname, map, hdr
  endif

  if index.has_nuv then begin
     map = readfits(nuv_fname, hdr,/silent)
     index.time_nuv = sxpar(hdr, 'MEANINT')
     index.rms_nuv = sxpar(hdr, 'RMS')
     index.std_nuv = sxpar(hdr, 'STDDEV')
     index.maskfrac_nuv = sxpar(hdr, 'MASKFRAC')
     index.anuv = sxpar(hdr, 'MWEXT')
     index.sat_effects_nuv = total((map gt sat_thresh_nuv)*footprint) ge 1
     sxaddpar, hdr, 'SATURATE', index.sat_effects_nuv, 'Saturation expected in image?'
     writefits, nuv_fname, map, hdr
  endif

  if index.has_fuv then begin
     map = readfits(fuv_fname, hdr,/silent)
     index.time_fuv = sxpar(hdr, 'MEANINT')
     index.rms_fuv = sxpar(hdr, 'RMS')
     index.std_fuv = sxpar(hdr, 'STDDEV')
     index.maskfrac_fuv = sxpar(hdr, 'MASKFRAC')
     index.afuv = sxpar(hdr, 'MWEXT')
     index.sat_effects_fuv = total((map gt sat_thresh_fuv)*footprint) ge 1
     sxaddpar, hdr, 'SATURATE', index.sat_effects_fuv, 'Saturation expected in image?'
     writefits, fuv_fname, map, hdr
  endif

  return, index

end
