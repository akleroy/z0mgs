pro v2_convert_unwise_to_mjysr $
   , infile = infile $
   , outfile = outfile $
   , band = this_band

; CALIBRATION TO GO FROM VEGAS TO ABMAG
  w1_vtoab = 2.683
  w2_vtoab = 3.319
  w3_vtoab = 5.242
  w4_vtoab = 6.604

; NORMALIZATION OF UNITY IN VEGAS MAG
  norm_mag = 22.5
  pix_as = 2.75
  
; COUNTS -> JY CONVERSION
  w1_to_mjysr = 10.^((norm_mag+w1_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2
  w2_to_mjysr = 10.^((norm_mag+w2_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2
  w3_to_mjysr = 10.^((norm_mag+w3_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2
  w4_to_mjysr = 10.^((norm_mag+w4_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2

  im = readfits(infile, hdr, /silent)

  if typename(this_band) eq typename(0) then begin
     if this_band eq 1 then this_band = 'w1'
     if this_band eq 2 then this_band = 'w2'
     if this_band eq 3 then this_band = 'w3'
     if this_band eq 4 then this_band = 'w4'
  endif
  
  if this_band eq 'w1' then begin 
     im *= w1_to_mjysr
     sxaddpar, hdr, 'NMTOMJY', w1_to_mjysr
  endif

  if this_band eq 'w2' then begin
     im *= w2_to_mjysr
     sxaddpar, hdr, 'NMTOMJY', w2_to_mjysr
  endif

  if this_band eq 'w3' then begin
     im *= w3_to_mjysr
     sxaddpar, hdr, 'NMTOMJY', w3_to_mjysr
  endif

  if this_band eq 'w4' then begin
     im *= w4_to_mjysr
     sxaddpar, hdr, 'NMTOMJY', w4_to_mjysr
  endif
        
  sxaddpar, hdr, 'BUNIT', 'MJY/SR'
  writefits, outfile, im, hdr

end
