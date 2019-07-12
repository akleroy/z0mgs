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

; zap photometry 

  struct_delete_field, index, 'flux_fuv'
  struct_delete_field, index, 'rms_flux_fuv'
  struct_delete_field, index, 'std_flux_fuv'

  struct_delete_field, index, 'flux_nuv'
  struct_delete_field, index, 'rms_flux_nuv'
  struct_delete_field, index, 'std_flux_nuv'

  struct_delete_field, index, 'flux_wise1'
  struct_delete_field, index, 'rms_flux_wise1'
  struct_delete_field, index, 'std_flux_wise1'

  struct_delete_field, index, 'flux_wise2'
  struct_delete_field, index, 'rms_flux_wise2'
  struct_delete_field, index, 'std_flux_wise2'

  struct_delete_field, index, 'flux_wise3'
  struct_delete_field, index, 'rms_flux_wise3'
  struct_delete_field, index, 'std_flux_wise3'

  struct_delete_field, index, 'flux_wise4'
  struct_delete_field, index, 'rms_flux_wise4'
  struct_delete_field, index, 'std_flux_wise4'

; REMOVE MOST PHYSICAL PROPERTIES

  struct_delete_field, index, 'dist_mpc'
  struct_delete_field, index, 'e_dist_mpc'

  struct_delete_field, index, 'mtol'
  struct_delete_field, index, 'method_mtol'
  struct_delete_field, index, 'logmass'
  struct_delete_field, index, 'e_logmass'

  struct_delete_field, index, 'logsfr'
  struct_delete_field, index, 'e_logsfr'  
  struct_delete_field, index, 'method_sfr'

  struct_delete_field, index, 'absbtc'
  struct_delete_field, index, 'complete_sample'  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  mwrfits, index, outfile, /create  

  stop

end
