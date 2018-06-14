pro build_bright_star_mask $
   , infile = infile $
   , map = map $
   , hdr = hdr $
   , fwhm = fwhm $
   , outfile=outfile $
   , show=show $
   , pause=pause $
   , band=band $
   , minpeak=minpeak $
   , thresh=thresh $
   , star_tol = star_tol $
   , star_ra = star_ra $
   , star_dec = star_dec $
   , star_km = star_km $
   , n_found = n_found $
   , ra_found = ra_found $
   , dec_found = dec_found $
   , km_found = km_found

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

  if n_elements(fwhm) eq 0 then $
     fwhm = sxpar(hdr, 'BMAJ')

  sz = size(map)
  make_axes, hdr, ri=ri, di=di
  pix_deg = sphdist(ri[0,0], di[0,0], ri[1,0], di[1,0], /deg) 
  x = findgen(sz[1]) # (fltarr(sz[2])+1.0)
  y = (fltarr(sz[1])+1.0) # findgen(sz[2])
  
  if n_elements(star_ra) eq 0 then begin
     restore, '../measurements/2mass_stars.idl', /v
  endif
  
  if n_elements(star_tol) eq 0 then $
     star_tol = 0.0

  nan = !values.f_nan

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BLANK KNOWN STARS IN THE FIELD
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  star_mask = finite(map)*0B
  error_mask = finite(map)*0B
  
  ra_tol = star_tol/cos(!dtor*mean(di,/nan))

  max_ra = max(ri, /nan)
  min_ra = min(ri, /nan)
  max_dec = max(di, /nan)
  min_dec = min(di, /nan)
     
  in_image = $
     (star_ra ge (min_ra-ra_tol)) and $
     (star_ra le (max_ra+ra_tol)) and $
     (star_dec ge (min_dec - star_tol)) and $
     (star_dec le (max_dec + star_tol))
  
  star_ind = where(in_image, n_found)
  print, "Found ", n_found, " stars that I will blank."

  if n_found gt 0 then begin

     intens = fltarr(n_found)*nan
     adxy, hdr, star_ra[star_ind], star_dec[star_ind], star_x, star_y
     really_in_image = $
        ((star_x ge 0) and (star_y ge 0) and $
         (star_x le (sz[1]-1)) and (star_y le (sz[2]-1)))
     really_in_image_ind = where(really_in_image, really_ct)

     if really_ct gt 0 then begin
;       Note  that we subtract a background from the image before
;       calculating the intensity.
        bkval = median(map)
        intens[really_in_image_ind] = $
           map[star_x[really_in_image_ind], star_y[really_in_image_ind]] - bkval
     endif
     
     rad_to_blank = $
        get_rad_to_blank(mag=star_km[star_ind] $
                         , intens=intens $
                         , band=band $
                         , error_rad=error_rad $
                         , fwhm=fwhm $
                         , minpeak=minpeak $
                         , thresh=thresh)
     rad_to_blank_pix = rad_to_blank / pix_deg
     error_rad_pix = error_rad / pix_deg
     
     ra_found = star_ra[star_ind]
     dec_found = star_dec[star_ind]
     km_found = star_km[star_ind]
     
     for kk = 0, n_found-1 do begin
        counter, kk, n_found, 'Blanking star '

        this_x = round(star_x[kk])
        this_y = round(star_y[kk])
        this_rad = rad_to_blank_pix[kk]
        if finite(this_rad) eq 0 then continue
        this_error_rad = error_rad_pix[kk]

        footprint = this_rad + 5L
        if finite(this_error_rad) then $
           footprint = this_error_rad + 5L
        xlo = ((this_x - footprint) > 0)
        ylo = ((this_y - footprint) > 0)
        xhi = (this_x + footprint) < (sz[1]-1)
        yhi = (this_y + footprint) < (sz[2]-1)           
        if (xlo gt sz[1]-1) or (ylo gt sz[2]-1) or (xhi lt 0) or (yhi lt 0) then $
           continue

        dist = sqrt((x[xlo:xhi,ylo:yhi] - this_x)^2 + (y[xlo:xhi,ylo:yhi] - this_y)^2)        
        star_mask[xlo:xhi,ylo:yhi] = star_mask[xlo:xhi,ylo:yhi] or (dist le this_rad+1.0)
        
        if finite(this_error_rad) then begin
           error_mask[xlo:xhi,ylo:yhi] = (dist le this_error_rad) or error_mask[xlo:xhi,ylo:yhi]
        endif
        
     endfor
  
  endif

  mask = star_mask + 10.*error_mask

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if keyword_set(show) then begin

     !p.multi=[0,2,1]
     loadct, 0
     if total(finite(map)) lt 10 then $ 
        rms =0.0 $
     else $
        rms = mad(map[where(finite(map))])
     disp, map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, mask, lev=[1,11], /overplot, color=cgcolor('blue', 255)

     disp, map*(mask ne 1)*(mask ne 11) $
           , max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5

     if keyword_set(pause) then ch = get_kbrd(1)

     !p.multi=0

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if n_elements(outfile) gt 0 then begin
     sxaddpar, hdr, 'BUNIT', 'MASK'
     writefits, outfile, mask, hdr
  endif

end
