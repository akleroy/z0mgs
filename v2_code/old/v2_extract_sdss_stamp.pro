function v2_extract_sdss_stamp $
;  INPUTS
   , filter=filter $
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

; DIRECTORIES   
  index_dir = '../../working_data/sdss/index/'

; READ THE INDEX FILE (IF NOT PASSED IN)
  if n_elements(index) eq 0 then $
     index = mrdfits(index_dir+'sdss_frame_index.fits',1,h)  

; MAKE A HEADER
  pix_scale = 0.5/3600.
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

; Right now I don't record the corders in SDSS so we do this less
; precisely.
  
;  tile_overlaps = $
;     v2_calc_tile_overlap( $
;     ra_ctr=ra_ctr $
;     , dec_ctr=dec_ctr $
;     , pad=size_deg $
;     , max_ra=index.max_ra $
;     , min_ra=index.min_ra $
;     , max_dec=index.max_dec $
;     , min_dec=index.min_dec)

; FIND CONTRIBUTING TILES
  
  filter_matches = (strcompress(index.filter, /rem) eq filter)
  dist = sphdist(index.ra, index.dec, ra_ctr, dec_ctr, /deg)
  pad = 0.25
  dist_thresh = (size_deg*0.5+pad)
  
  ind = where(dist le dist_thresh and $
              filter_matches $
              , ct_overlap)

  if ct_overlap eq 0 then begin
     print, "No overlapping frames found. Returning."
     return, -1B
  endif
  
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

     this_image_fname = $
        strcompress(index[ind[ii]].fname, /rem)
     
;    READ THE IMAGE
     im = mrdfits(this_image_fname, 0, hdr)
     if total((im ne 0) and finite(im)) eq 0 then begin
        message, 'Nothing in image. Proceeding.', /nan
        continue
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVERT UNITS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%     

;    FIGURE THE PIXEL SIZE
     xyad, hdr, 0, 0, ra1, de1
     xyad, hdr, 1, 0, ra2, de2
     pix_rad = sphdist(ra1, de1, ra2, de2, /deg)*!dtor
     pix_sr = pix_rad^2

;    CONVERT UNITS TO MJY/SR or AB/AS2

; ... tbd     
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%          
; ALIGN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%          

     im = 1.0*im
     hastrom, im, hdr, target_hdr $
              , interp=1, missing=!values.f_nan $
              , errmsg=errmsg
     if strcompress(errmsg, /rem) ne '' then continue

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%               
; ACCUMULATE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%               
          
     im_fin_ind = where(finite(im), im_fin_ct)
     if im_fin_ct eq 0 then continue

     sum_image[im_fin_ind] = $
        sum_image[im_fin_ind] + im[im_fin_ind]
     weight_image[im_fin_ind] = $
        weight_image[im_fin_ind]+((finite(im)*1.0)[im_fin_ind])

     ; sxaddpar, target_hdr, 'BUNIT', 'MJY/SR'

  endfor

  combined = sum_image/weight_image

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

  return, 1B
  
end
