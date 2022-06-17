pro convert_unwise_to_mjysr $
   , infile = infile $
   , outfile = outfile $
   , band = band

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

  n_files = n_elements(infile)

  for ii = 0, n_elements(infile)-1 do begin
  
     im = readfits(infile, hdr, /silent)

     if n_elements(band) eq 1 then $
        this_band = band $
     else $
        this_band = band[ii]

     if this_band eq 1 then begin 
        im *= w1_to_mjysr
     endif
     if this_band eq 2 then begin
        im *= w2_to_mjysr
     endif
     if this_band eq 3 then begin
        im *= w3_to_mjysr
     endif
     if this_band eq 4 then begin
        im *= w4_to_mjysr
     endif
        
     sxaddpar, hdr, 'BUNIT', 'MJY/SR'
     writefits, outfile[ii], im, hdr

  endfor

end
