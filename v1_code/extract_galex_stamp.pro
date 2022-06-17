pro extract_galex_stamp $
;  INPUTS
   , fuv=fuv $
   , ra_ctr=ra_ctr $
   , dec_ctr=dec_ctr $
   , size_deg=size_deg $
   , index=index $
   , useint=useint $
;  DISPLAY OPTIONS
   , show=show $
   , pause=pause $
   , bksub=bksub $
;  OUTPUTS
   , image=combined $
   , weight=weight_image $
   , noise=noise_image $
   , hdr=target_hdr

; FUV/NUV
  if keyword_set(fuv) eq 0 then $
     fuv = 0B

; DIRECTORIES   
  index_dir = '../z0mgs/'
  data_dir = '/data/fourier/leroy.42/allsky/all_tiles/'

; READ THE INDEX FILE (IF NOT PASSED IN)
  if n_elements(index) eq 0 then $
     index = mrdfits(index_dir+'galex_index_file.fits',1,h)  

; SPECIFY IF WE USE THE INT OR INTBGSUB CASES
  if keyword_set(useint) then $
     image_ext= 'int' $
  else $
     image_ext = 'intbgsub'

; MAKE A HEADER
  pix_scale = 1.5/3600.
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

; FIND CONTRIBUTING TILES
  ind = where(index.fuv eq fuv and $
              tile_overlaps, ct_overlap)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE THE OUTPUT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  make_axes, target_hdr, ri=ri_targ, di=di_targ

  sz_out = size(ri_targ)
  weight_image = ri_targ*0.0
  sum_image = ri_targ*0.0
;  sum_noise_image = ri_targ*0.0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER OVERLAPPING TILES AND STITCH ONTO TARGET HEADER
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = 0, ct_overlap-1 do begin

     counter, ii+1, ct_overlap, 'Processing image '

;    THE IMAGE
     this_image_fname = $
        data_dir+ $
        strcompress(str_replace((strsplit(index[ind[ii]].fname, '/', /extract))[1] $
                                , 'int', image_ext)+'.gz', /rem)
     
;    THE RELATIVE RESPONSE
     this_rrhr_fname = $
        str_replace(this_image_fname, image_ext, 'rrhr')

;    THE FLAG FILE
     this_flag_fname = $
        str_replace(this_image_fname, image_ext, 'flags')

;    READ THE IMAGE
     im = readfits(this_image_fname, hdr, /silent)
     if total((im ne 0) and finite(im)) eq 0 then begin
        message, 'Nothing in image. Proceeding.', /nan
        continue
     endif

     dummy = sxpar(hdr, 'CDELT1', count=cd_ct)
     if cd_ct ne 1 then begin
        message, 'No astrometry found. Proceeding.', /info
        continue
     endif

     rrhr = readfits(this_rrhr_fname, rrhr_hdr, /silent)     

     flag = readfits(this_flag_fname, flag_hdr, /silent)
     hastrom, flag, flag_hdr, hdr, interp=0, missing=1024

;    APPLY THE FLAGS
     edge = (((flag mod 256) mod 128) mod 64) ge 32
     im[where(edge)] = !values.f_nan
     rrhr[where(edge)] = !values.f_nan

;     make_axes, hdr, ri=ri, di=di
;     adxy, target_hdr, ri, di, x, y
          
;     in_image = (x gt 0 and x lt (sz_out[1]-1)) and $
;                (y gt 0 and y lt (sz_out[2]-1))
;     if total(in_image) eq 0 then begin
;        print, "No overlap. Proceeding."
;        continue
;     endif

     ;edge = where((flag mod 256 mod 128 mod 64) eq 32)
     ;im[edge] = !values.f_nan

;     if keyword_set(bksub) then begin
;        meanval = mean(im)
;        im -= meanval
;     endif

;    FIGURE THE PIXEL SIZE

     xyad, hdr, 0, 0, ra1, de1
     xyad, hdr, 1, 0, ra2, de2
     pix_rad = sphdist(ra1, de1, ra2, de2, /deg)*!dtor
     pix_sr = pix_rad^2

;    CONVERT TO MJY/SR
     
     if keyword_set(fuv) then begin
        im = fuv_cps_to_jy(im)/1d6
        im /= pix_sr
     endif else begin
        im = nuv_cps_to_jy(im)/1d6
        im /= pix_sr
     endelse

;    FIGURE THE RMS
     
;     im_for_noise = smooth(im, 5, /nan)
;     rms = mad(im_for_noise)
;     med = median(im_for_noise)
;     keep = where(abs(im_for_noise - med) lt 3.*rms, keep_ct)
;     if keep_ct eq 0 then stop
;     rms = mad(im_for_noise[keep])

;    ALIGN
     
     im = 1.0*im
     hastrom, im, hdr, target_hdr $
              , interp=1, missing=!values.f_nan $
              , errmsg=errmsg
     if strcompress(errmsg, /rem) ne '' then continue

     rrhr = 1.0*rrhr
     hastrom, rrhr, rrhr_hdr, target_hdr $
              , interp=1, missing=!values.f_nan
     if strcompress(errmsg, /rem) ne '' then continue
     
     im_ind = where(finite(rrhr) and (rrhr gt 0), im_ct)
     if im_ct eq 0 then continue

     sum_image[im_ind] = sum_image[im_ind]+(rrhr[im_ind])*im[im_ind]
     weight_image[im_ind] = weight_image[im_ind]+(rrhr[im_ind])
;     sum_noise_image[im_ind] = sum_noise_image[im_ind] + (rms*rrhr[im_ind])^2

     sxaddpar, target_hdr, 'BUNIT', 'MJY/SR'

  endfor

  combined = sum_image/weight_image
;  noise_image = sqrt(sum_noise_image)/weight_image

; RENORMALIZE THE NOISE
;  renorm = mad(combined)/median(noise_image)
;  noise_image = noise_image*renorm

; SHOW IF DESIRED
  
  if keyword_set(show) then begin
     disp, combined $
           , max=10*mad(combined) $
           , min=-5*mad(combined)
     if keyword_set(pause) then begin
        print, "hit a key to continue..."
        ch = get_kbrd(1)
     endif
  endif

end
