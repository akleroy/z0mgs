pro customize_masks

  readcol, 'custom_mask_specs.txt', format='A,F,F,F', comment='#' $
           , override_pgc, override_rad, override_pa, override_incl
  override_pgc = strcompress(override_pgc, /rem)
  n_pgc = n_elements(override_pgc)

;  start = 50
;  override_pgc = override_pgc[start:*]
;  override_rad = override_rad[start:*]
;  override_pa = override_pa[start:*]
;  override_incl = override_incl[start:*]

; LOOK OVER THE NEW MASKS
  compile_unwise_atlas $
     , just=override_pgc $
     , /show, /pause, /mask
  
  compile_galex_atlas $
     , just=override_pgc $
     , /show, /pause, /mask
  
; NOW RUN THE COMPILATION FOR THE NEW GALAXIES
  compile_galex_atlas $
     , just=override_pgc $
     , /show, /bksub, /special, /mask, /stat
  
  compile_unwise_atlas $
     , just=override_pgc $
     , /show, /bksub, /mask, /stat

  build_delivery $
     , just=override_pgc, /wise, /galex

end
