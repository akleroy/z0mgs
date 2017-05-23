pro extract_galex_stamp $
   , fuv=fuv $
   , ra_ctr=ra_ctr $
   , dec_ctr=dec_ctr $
   , size_deg=size_deg $
   , index=index $
   , image=combined $
   , weight=weight_image $
   , hdr=target_hdr $
   , show=show $
   , pause=pause $
   , bksub=bksub

  index_dir = '/data/tycho/0/leroy.42/allsky/code/'
  data_dir = '/data/tycho/0/leroy.42/allsky/galex/sorted_tiles/'
    
; MAKE A HEADER
  pix_scale = 1.5/3600.
  pix_len = size_deg / pix_scale

  if keyword_set(fuv) eq 0 then $
     fuv = 0B

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
     index = mrdfits(index_dir+'galex_index_file.fits',1,h)  

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
  ind = where(index.fuv eq fuv and $
              tile_overlaps, ct_overlap)

; SET UP THE OUTPUT
  make_axes, target_hdr, ri=ri_targ, di=di_targ
  sz_out = size(ri_targ)
  outim = ri_targ*!values.f_nan

; LOOP OVER OVERLAPPING TILES AND STITCH ONTO TARGET HEADER
  weight_image = ri_targ*0.0
  sum_image = ri_targ*0.0
  for ii = 0, ct_overlap-1 do begin
     this_file = strcompress(data_dir+index[ind[ii]].fname, /rem)
     this_rrhr = strcompress(data_dir+index[ind[ii]].rrhrfile, /rem)

     im = readfits(this_file, hdr)
     if total((im ne 0) and finite(im)) eq 0 then begin
        message, 'Nothing in image. Proceeding.', /nan
        continue
     endif

     make_axes, hdr, ri=ri, di=di
     adxy, target_hdr, ri, di, x, y
          
     in_image = (x gt 0 and x lt (sz_out[1]-1)) and $
                (y gt 0 and y lt (sz_out[2]-1))
     if total(in_image) eq 0 then begin
        print, "No overlap. Proceeding."
        continue
     endif

     rrhr = readfits(this_rrhr, rrhr_hdr)

     if keyword_set(bksub) then begin
        meanval = mean(im)
        im -= meanval
     endif

     xyad, hdr, 0, 0, ra1, de1
     xyad, hdr, 1, 0, ra2, de2
     pix_rad = sphdist(ra1, de1, ra2, de2, /deg)*!dtor
     pix_sr = pix_rad^2
     
     if keyword_set(fuv) then begin
        im = fuv_cps_to_jy(im)/1d6
        im /= pix_sr
     endif else begin
        im = nuv_cps_to_jy(im)/1d6
        im /= pix_sr
     endelse
     
     hastrom, im, hdr, target_hdr $
              , interp=1, missing=!values.f_nan

     hastrom, rrhr, rrhr_hdr, target_hdr $
              , interp=1, missing=!values.f_nan
     
     im_ind = where(finite(rrhr) and (rrhr gt 0), im_ct)
     if im_ct eq 0 then continue

     sum_image[im_ind] = sum_image[im_ind]+(rrhr[im_ind])*im[im_ind]
     weight_image[im_ind] = weight_image[im_ind]+(rrhr[im_ind])

     if keyword_set(show) then $
        disp, im, /sq, max=1d-3

     sxaddpar, target_hdr, 'BUNIT', 'MJY/SR'
  
  endfor

  combined = sum_image/weight_image

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
