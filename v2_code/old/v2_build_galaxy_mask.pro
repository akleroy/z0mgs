pro v2_build_galaxy_mask $           
   , this_pgc=this_pgc $
   , pgc_to_skip=pgc_to_skip $
   , all_gal_data=all_gal_data $
   , infile=template_image_file $
   , outfile=outfile $
   , rad_to_blank=rad_to_blank $
   , min_rad=min_rad $
   , show=show $
   , pause=pause

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Initialize
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if n_elements(this_pgc) eq 0 then begin
     this_pgc = -1
  endif
  
  if n_elements(template_image_file) eq 0 then begin
     return
  endif
  
  if file_test(template_image_file) eq 0 then begin
     print, template_image_file, " not found."
     return
  endif
  
  if n_elements(all_gal_data) eq 0 then begin
     print, "... loading galaxy database."
     all_gal_data = gal_data(/all,/full)
  endif

  if n_elements(rad_to_blank) eq 0 then begin
     rad_to_blank = 1.0
  endif

  if n_elements(min_rad) eq 0 then begin
     min_rad = 7.5/3600.
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Find galaxies in the image  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
; Initialize the mask
  map = readfits(template_image_file, hdr)
  mask = finite(map)*0B

; Figure out the coordinates for the image

; ... pixel scale  
  xyad, hdr, [0,0], [0,1], ra, dec
  pix_deg = abs(sphdist(ra[0], dec[0], ra[1], dec[1], /deg))

; ... image edges  
  make_axes, hdr, ri=ri, di=di
  max_ra = max(ri, /nan)
  min_ra = min(ri, /nan)
  max_dec = max(di, /nan)
  min_dec = min(di, /nan)

; ... identify circular footprint of each galaxy
  cos_fac = cos(!dtor*(max_dec+min_dec)*0.5)
  tol = all_gal_data.r25_deg
  tol[where(finite(tol) eq 0)] = 30./3600.
  ra_tol = tol/cos_fac

; ... identify which galaxies are in the image  
  in_image = $
     (all_gal_data.ra_deg ge (min_ra - ra_tol)) and $
     (all_gal_data.ra_deg le (max_ra + ra_tol)) and $
     (all_gal_data.dec_deg ge (min_dec - tol)) and $
     (all_gal_data.dec_deg le (max_dec + tol))  

  gal_ind_to_mask = $
     where(in_image and $
           (all_gal_data.pgc ne this_pgc), gal_ct)

  print, "Found ", gal_ct, " galaxies in image ", template_image_file

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Loop over identified galaxies and add to mask
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  n_found = 0L
  for kk = 0, gal_ct-1 do begin

; Check if this is a galaxy we're skipping     
     this_other_pgc = all_gal_data[gal_ind_to_mask[kk]].pgc
     if n_elements(pgc_to_skip) gt 0 then $
        if total(this_other_pgc eq pgc_to_skip) gt 0 then $
           continue

; Get info for this galaxy     
     this_ra = all_gal_data[gal_ind_to_mask[kk]].ra_deg
     this_dec = all_gal_data[gal_ind_to_mask[kk]].dec_deg
     this_rad = all_gal_data[gal_ind_to_mask[kk]].r25_deg*rad_to_blank
     if finite(this_rad) eq 0 or (this_rad lt min_rad) then $
        this_rad = min_rad

     this_incl = all_gal_data[gal_ind_to_mask[kk]].incl_deg
     this_pa = all_gal_data[gal_ind_to_mask[kk]].posang_deg
     if finite(this_pa) eq 0 or finite(this_incl) eq 0 then begin
        this_incl = 0.0
        this_pa = 0.0              
     endif
     if this_incl gt 60. then this_incl = 60.

; Calculate the footprint of the galaxy in the image     
     this_pos_vec = $
        [this_pa, this_incl, this_ra, this_dec]
     deproject, ri, di, this_pos_vec, rgrid=this_rgrid

; Blank the mask image
     blank_ind = where(this_rgrid le (this_rad+pix_deg), blank_ct)
     if blank_ct gt 0 then begin
        mask[blank_ind] = 1B
     endif

; Record the found galaxy in the image     
     n_found += 1
     if n_elements(pgc_found) eq 0 then begin
        pgc_found = [this_other_pgc]
     endif else begin
        pgc_found = [pgc_found, this_other_pgc]
     endelse
           
  endfor   

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; If requested show the image and the mask
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(show) then begin

     !p.multi=0
     loadct, 0
     fin_ind = where(finite(map) and mask ne 1, fin_ct)
     if fin_ct lt 10 then begin
        rms = 0.0
     endif else begin
        rms = mad(map[fin_ind])
     endelse
     disp, map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=5, ystyle=5, reserve=5
     contour, mask, lev=[1], /overplot, color=cgcolor('blue', 255)
     if keyword_set(pause) then ch = get_kbrd(1)
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
; Write to disk
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if n_elements(outfile) ne 0 then begin
     sxaddpar, hdr, 'BUNIT', 'MASK'
     writefits, outfile, mask, hdr
  endif
  
end
