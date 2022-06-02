pro write_pared_down_gauss7p5_index

;+
;
;-
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONSTANTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  @constants.bat
  lsun_3p4 = 1.83d18

  infile = '../measurements/delivery_index_gauss7p5.fits'
  outfile = '../measurements/simple_index_gauss7p5.fits'
  index = mrdfits(infile,1,h7p5)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; REMOVE FIELDS THAT WE DON'T WANT TO DELIVER TWICE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  fid_entry = index[0]

; zap photometry 

  struct_delete_field, fid_entry, 'flux_fuv'
  struct_delete_field, fid_entry, 'rms_flux_fuv'
  struct_delete_field, fid_entry, 'std_flux_fuv'

  struct_delete_field, fid_entry, 'flux_nuv'
  struct_delete_field, fid_entry, 'rms_flux_nuv'
  struct_delete_field, fid_entry, 'std_flux_nuv'

  struct_delete_field, fid_entry, 'flux_wise1'
  struct_delete_field, fid_entry, 'rms_flux_wise1'
  struct_delete_field, fid_entry, 'std_flux_wise1'

  struct_delete_field, fid_entry, 'flux_wise2'
  struct_delete_field, fid_entry, 'rms_flux_wise2'
  struct_delete_field, fid_entry, 'std_flux_wise2'

  struct_delete_field, fid_entry, 'flux_wise3'
  struct_delete_field, fid_entry, 'rms_flux_wise3'
  struct_delete_field, fid_entry, 'std_flux_wise3'

  struct_delete_field, fid_entry, 'flux_wise4'
  struct_delete_field, fid_entry, 'rms_flux_wise4'
  struct_delete_field, fid_entry, 'std_flux_wise4'

; REMOVE MOST PHYSICAL PROPERTIES

  struct_delete_field, fid_entry, 'dist_mpc'
  struct_delete_field, fid_entry, 'e_dist_dex'

  struct_delete_field, fid_entry, 'mtol'
  struct_delete_field, fid_entry, 'method_mtol'
  struct_delete_field, fid_entry, 'logmass'
  struct_delete_field, fid_entry, 'e_logmass'

  struct_delete_field, fid_entry, 'logsfr'
  struct_delete_field, fid_entry, 'e_logsfr'  
  struct_delete_field, fid_entry, 'method_sfr'

  struct_delete_field, fid_entry, 'absbtc'
  struct_delete_field, fid_entry, 'complete_sample'  

  struct_delete_field, fid_entry, 'deltams'
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; NOW COPY AND ASSIGN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  new_index = replicate(fid_entry, n_elements(index))
  struct_assign, index, new_index, /verbose
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  mwrfits, new_index, outfile, /create  

  stop

end
