pro extract_unwise_stamp $
   , band=band $
   , ra_ctr=ra_ctr $
   , dec_ctr=dec_ctr $
   , size_deg=size_deg $
   , index=index $
   , image=outim $
   , hdr=target_hdr $
   , show=show $
   , pause=pause

  index_dir = '/data/tycho/0/leroy.42/allsky/z0mgs/'
  data_dir = '/data/tycho/0/leroy.42/allsky/unwise/sorted_tiles/'
  
 ; CALIBRATION TO GO FROM VEGAS TO ABMAG
  w1_vtoab = 2.683
  w2_vtoab = 3.319
  w3_vtoab = 5.242
  w4_vtoab = 6.604

; NORMALIZATION OF UNITY IN VEGAS MAG
  norm_mag = 22.5
  pix_as = 2.75
  
; COUNTS -> JY CONVERSION
  w1_to_mjysr = $
     10.^((norm_mag+w1_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2
  w2_to_mjysr = $
     10.^((norm_mag+w2_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2
  w3_to_mjysr = $
     10.^((norm_mag+w3_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2
  w4_to_mjysr = $
     10.^((norm_mag+w4_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2
 
; MAKE A HEADER
  pix_scale = 2.0/3600.
  pix_len = size_deg / pix_scale

  mkhdr, target_hdr, 2, [pix_len, pix_len]
  sxaddpar, target_hdr, 'NAXIS1', pix_len
  sxaddpar, target_hdr, 'NAXIS2', pix_len
  sxaddpar, target_hdr, 'CTYPE1', 'RA---TAN'
  sxaddpar, target_hdr, 'CRVAL1', double(ra_ctr)
  sxaddpar, target_hdr, 'CRPIX1', (pix_len/2.)*1.
  sxaddpar, target_hdr, 'CDELT1', -1.0*pix_scale
  sxaddpar, target_hdr, 'CTYPE2', 'DEC--TAN'
  sxaddpar, target_hdr, 'CRVAL2', double(dec_ctr)
  sxaddpar, target_hdr, 'CRPIX2', (pix_len/2.)*1.
  sxaddpar, target_hdr, 'CDELT2', pix_scale
  sxaddpar, target_hdr, 'EQUINOX', 2000.

; READ THE INDEX FILE (IF NOT PASSED IN)
  if n_elements(index) eq 0 then $
     index = mrdfits(index_dir+'unwise_index_file.fits',1,h)  

; CALCULATE TILE OVERLAP  
  tile_overlaps = $
     calc_tile_overlap( $
     ra_ctr=ra_ctr $
     , dec_ctr=dec_ctr $
     , pad=size_deg $
     , min_ra=index.min_ra $
     , max_ra=index.max_ra $
     , min_dec=index.min_dec $
     , max_dec=index.max_dec)

; FIND OVERLAPPING TILES WITH THE RIGHT BAND
  ind = where(index.band eq band and $
              tile_overlaps, ct_overlap)

; SET UP THE OUTPUT
  make_axes, target_hdr, ri=ri_targ, di=di_targ
  sz_out = size(ri_targ)
  outim = ri_targ*!values.f_nan

; LOOP OVER OVERLAPPING TILES AND STITCH ONTO TARGET HEADER

  for ii = 0, ct_overlap-1 do begin
     this_file = strcompress(data_dir+index[ind[ii]].fname, /rem)
     im = readfits(this_file, hdr)
     make_axes, hdr, ri=ri, di=di
     adxy, target_hdr, ri, di, x, y
     
     in_image = (x gt 0 and x lt (sz_out[1]-1)) and $
                (y gt 0 and y lt (sz_out[2]-1))
     if total(in_image) eq 0 then begin
        print, "No overlap. Proceeding."
        continue
     endif

     if band eq 1 then begin
        im *= w1_to_mjysr
     endif
     if band eq 2 then begin
        im *= w2_to_mjysr
     endif
     if band eq 3 then begin
        im *= w3_to_mjysr
     endif
     if band eq 4 then begin
        im *= w4_to_mjysr
     endif   
     
     hastrom, im, hdr, target_hdr $
              , interp=1, missing=!values.f_nan
     
     useful = where(finite(im), usect)
     outim[useful] = im[useful]

     sxaddpar, target_hdr, 'BUNIT', 'MJY/SR'
  
  endfor

; SHOW IF DESIRED
  
  if keyword_set(show) then begin
     disp, outim $
           , max=10*mad(outim) $
           , min=-5*mad(outim)
     if keyword_set(pause) then begin
        print, "Band "+str(band)+" hit a key to continue..."
        ch = get_kbrd(1)
     endif
  endif
  


  
end
