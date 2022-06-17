pro build_galaxy_mask $
   , infile = infile $
   , map = map $
   , hdr = hdr $
   , skip_pgc = pgc_to_skip $
   , outfile = outfile $
   , mask=mask $   
   , show=show $
   , pause=pause $
   , galdata = allgals $
   , rad_to_blank = rad_to_blank $
   , min_rad = min_rad $
   , n_found = n_found $
   , pgc_found = pgc_found
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if n_elements(map) eq 0 or n_elements(hdr) eq 0 then begin
     if file_test(infile) eq 0 then begin
        message, 'Target file not found and map and header not supplied.', /info
        return
     endif     
     map = readfits(infile, hdr, /silent)  
  endif
  xyad, hdr, [0,0], [0,1], ra, dec
  pix_deg = abs(sphdist(ra[0], dec[0], ra[1], dec[1], /deg))

  sz = size(map)
  mask = finite(map)*0B

  if n_elements(allgals) eq 0 then begin
     allgals = gal_data(/all)
  endif    

  if n_elements(rad_to_blank) eq 0 then begin
     rad_to_blank = 1.0
  endif

  if n_elements(min_rad) eq 0 then begin
     min_rad = 7.5/3600.
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FIND GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  make_axes, hdr, ri=ri, di=di
  max_ra = max(ri, /nan)
  min_ra = min(ri, /nan)
  max_dec = max(di, /nan)
  min_dec = min(di, /nan)
       
  cos_fac = cos(!dtor*mean(di,/nan))
  tol = allgals.r25_deg
  tol[where(finite(tol) eq 0)] = 30./3600.
  ra_tol = tol/cos_fac
  
  in_image = $
     (allgals.ra_deg ge (min_ra - ra_tol)) and $
     (allgals.ra_deg le (max_ra + ra_tol)) and $
     (allgals.dec_deg ge (min_dec - tol)) and $
     (allgals.dec_deg le (max_dec + tol))  
  
  gal_ind = where(in_image, gal_ct)
  print, "Found ", gal_ct, " galaxies in the image."
  
  n_found = 0L
  for kk = 0, gal_ct-1 do begin
     
     this_pgc = allgals[gal_ind[kk]].pgc
     if n_elements(pgc_to_skip) gt 0 then $
        if total(this_pgc eq pgc_to_skip) gt 0 then $
           continue

     this_ra = allgals[gal_ind[kk]].ra_deg
     this_dec = allgals[gal_ind[kk]].dec_deg
     this_rad = allgals[gal_ind[kk]].r25_deg*rad_to_blank
     if finite(this_rad) eq 0 or (this_rad lt min_rad) then $
        this_rad = min_rad

     this_incl = allgals[gal_ind[kk]].incl_deg
     this_pa = allgals[gal_ind[kk]].posang_deg
     if finite(this_pa) eq 0 or finite(this_incl) eq 0 then begin
        this_incl = 0.0
        this_pa = 0.0              
     endif
     if this_incl gt 60. then this_incl = 60.
           
     this_pos_vec = $
        [this_pa, this_incl, this_ra, this_dec]
     deproject, ri, di, this_pos_vec, rgrid=this_rgrid
           
     blank_ind = where(this_rgrid le (this_rad+pix_deg), blank_ct)
     if blank_ct gt 0 then begin
        mask[blank_ind] = 1B
     endif
     
     n_found += 1
     if n_elements(pgc_found) eq 0 then begin
        pgc_found = [this_pgc]
     endif else begin
        pgc_found = [pgc_found, this_pgc]
     endelse
     
  endfor   

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if keyword_set(show) then begin

     !p.multi=0
     loadct, 0
     rms = mad(map[where(finite(map) and mask ne 10)])
     disp, map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=5, ystyle=5, reserve=5
     contour, mask, lev=[1], /overplot, color=cgcolor('blue', 255)
     if keyword_set(pause) then ch = get_kbrd(1)

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if n_elements(outfile) gt 0 then begin
     sxaddpar, hdr, 'BUNIT', 'MASK'
     writefits, outfile, mask, hdr
  endif

end
