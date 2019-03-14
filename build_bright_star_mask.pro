pro build_bright_star_mask $
   , pgcname=pgcname $
   , galdata=this_data $
   , infile=infile $
   , band=band $
   , res=res $
   , outfile=outfile $
   , star_ra = star_ra $
   , star_dec = star_dec $
   , star_km = star_km $
   , n_found = n_found $
   , ra_found = ra_found $
   , dec_found = dec_found $
   , km_found = km_found $
   , nopad =nopad $
   , show=show $
   , pause=pause

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS AND READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  nan = !values.f_nan

  gaia_dir = '../stars/gaia/'

  if n_elements(res) eq 0 then  begin
     res = 'gauss7p5'
  endif
  
  valid_res = ['gauss7p5', 'gauss15']
  if total(res eq valid_res) eq 0 then begin
     print, "Defaulting to 7.5as Gaussian"
     res = 'gauss7p5'
  endif

  valid_bands = ['w1','w2','w3','w4','nuv','fuv']
  if total(band eq valid_bands) eq 0 then begin
     print, "Defaulting to WISE1"
     band = 'w1'
  endif

  if file_test(infile) eq 0 then begin
     message, 'Target file not found and map and header not supplied.', /info
     return
  endif     
  map = readfits(infile, hdr, /silent)  
  
  if n_elements(star_ra) eq 0 then begin
     print, "Loading 2MASS stars."
     restore, '../measurements/2mass_stars.idl', /v
  endif
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; TUNING PARAMETERS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if band eq 'w1' then begin
     if res eq 'gauss7p5' then begin
        min_peak = 1d-2
        thresh = 0.5d-2
     endif else begin
        min_peak = 5d-2
        thresh = 1d-2
     endelse
  endif

  if band eq 'w2' then begin
     if res eq 'gauss7p5' then begin
        min_peak = 1d-2
        thresh = 0.5d-2
     endif else begin
        min_peak = 5d-2
        thresh = 1d-2
     endelse
  endif

  if band eq 'w3' then begin
     if res eq 'gauss7p5' then begin
        min_peak = 5d-2
        thresh = 1d-2
     endif else begin
        min_peak = 5d-2
        thresh = 1d-2
     endelse
  endif

  if band eq 'w4' then begin
     if res eq 'gauss7p5' then begin
        min_peak = !values.f_nan
        thresh = !values.f_nan
     endif else begin
        min_peak = 1.0
        thresh = 5d-1
     endelse
  endif

  if band eq 'nuv' then begin
     if res eq 'gauss7p5' then begin
        min_peak = 5d-3
        thresh = 5d-4
     endif else begin
        min_peak = 5d-3
        thresh = 5d-4
     endelse
  endif

  if band eq 'fuv' then begin
     if res eq 'gauss7p5' then begin
        min_peak = 5d-3
        thresh = 5d-4
     endif else begin
        min_peak = 5d-3
        thresh = 5d-4
     endelse
  endif

  overlap_center_thresh = 5.0/3600.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sz = size(map)
  make_axes, hdr, ri=ri, di=di, ra=ra, da=da
  pix_deg = sphdist(ri[0,0], di[0,0], ri[1,0], di[1,0], /deg) 
  x = findgen(sz[1]) # (fltarr(sz[2])+1.0)
  y = (fltarr(sz[1])+1.0) # findgen(sz[2])
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IDENTIFY 2MASS STARS TO BLANK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; PLACE ALL STARS IN PIXEL COORDS
  adxy, hdr, star_ra, star_dec, star_x, star_y

; IDENTIFY STARS IN THE IMAGE
  in_image = $
     (star_x ge 0) and (star_x lt sz[1]) and $
     (star_y ge 0) and (star_y lt sz[2])
  
  star_ind = where(in_image, n_found)
;  print, "Found ", n_found, " stars that I will blank in ", infile

  if n_found gt 0 then begin

;    INITIALIZE THE VECTOR OF STARS TO BLANK     

     x_to_blank = star_x[star_ind]
     y_to_blank = star_y[star_ind]     
     intens = $
        pred_star_intensity(star_km[star_ind] $
                            , band=band $
                            , res=res)
     
;    USE THE MAP FOR UV (DISABLED FOR NOW)

     ;if band eq 'nuv' or band eq 'fuv' then begin
     ;   intens = map[star_x[star_ind], star_y[star_ind]]
     ;endif

;    GRAB THE FOUND STARS, TO PASS BACK OUT
     
     ra_found = star_ra[star_ind]
     dec_found = star_dec[star_ind]
     km_found = star_km[star_ind]

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IDENTIFY GAIA STARS TO BLANK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  gaia_file = gaia_dir+pgcname+'_gaia.txt'
  readcol $
     , gaia_file $
     , delim=',' $
     , gaia_id, gaia_ra, gaia_dec, gaia_gmag $
     , format='L,D,D,D,X,X,X,X,X,X' $
     , count=lines $
     , /silent
  
  if lines gt 0 then begin
     
     adxy, hdr, gaia_ra, gaia_dec, gaia_x, gaia_y
     
     gaia_intens = $
        pred_star_intensity(gaia_gmag $
                            , /gaia $
                            , band=band $
                            , res=res)
     
     if n_elements(x_to_blank) eq 0 then begin
        x_to_blank = gaia_x
        y_to_blank = gaia_y     
        intens = gaia_intens
     endif else begin
        x_to_blank = [x_to_blank, gaia_x]
        y_to_blank = [y_to_blank, gaia_y]
        intens = [intens, gaia_intens]
     endelse

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; EXCLUDE "STARS" THAT COULD BE THE GALAXY NUCLEUS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(this_data) gt 0 then begin

     if n_elements(x_to_blank) gt 0 then begin

        adxy, hdr, this_data[0].ra_deg, this_data[0].dec_deg, x_ctr, y_ctr
        
        dist_to_ctr = $
           sqrt((x_to_blank - x_ctr)^2 + (y_to_blank - y_ctr)^2) * $
           pix_deg
        
     endif

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CALCULATE RADII AND BLANK ALL FOUND STARS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  star_mask = finite(map)*0B

  if n_elements(x_to_blank) gt 0 then begin
     
     keep = where(intens ge min_peak, keep_ct) 
     
     if keep_ct gt 0 then begin
        
        print, "I will blank " + str(keep_ct) + " stars in this image at band "+band
        
        x_to_blank = x_to_blank[keep]
        y_to_blank = y_to_blank[keep]
        intens = intens[keep]
        
        rad_to_blank = $
           get_rad_to_blank( $
           intens=intens $
           , fwhm=fwhm $
           , thresh=thresh)
        
        rad_to_blank_pix = rad_to_blank / pix_deg
        
;    PADDING - TUNE THIS
        
        if keyword_set(nopad) eq 0 then begin
           rad_to_blank_pix += 1.0
        endif
        
;    LOOP OVER STARS TO BLANK
        n_to_blank = n_elements(x_to_blank)
        
        for kk = 0, n_to_blank-1 do begin
           
           counter, kk, n_to_blank, 'Blanking 2MASS star '
           
           this_x = round(x_to_blank[kk])
           this_y = round(y_to_blank[kk])
           this_rad = rad_to_blank_pix[kk]
           
           if finite(this_rad) eq 0 then continue

           footprint = this_rad + 5L
           xlo = ((this_x - footprint) > 0)
           ylo = ((this_y - footprint) > 0)
           xhi = (this_x + footprint) < (sz[1]-1)
           yhi = (this_y + footprint) < (sz[2]-1)           
           if (xlo gt sz[1]-1) or (ylo gt sz[2]-1) or (xhi lt 0) or (yhi lt 0) then $
              continue

           dist = sqrt((x[xlo:xhi,ylo:yhi] - this_x)^2 + (y[xlo:xhi,ylo:yhi] - this_y)^2)        

           star_mask[xlo:xhi,ylo:yhi] = $
              star_mask[xlo:xhi,ylo:yhi] or (dist le this_rad+1.0)
           
        endfor
        
     endif

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if keyword_set(show) then begin

     !p.multi=[0,2,2]

     if total(finite(map)) lt 10 then $ 
        rms =0.0 $
     else $
        rms = mad(map[where(finite(map))])

     loadct, 0
     disp, map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, star_mask, lev=[1,11], /overplot, color=cgcolor('blue', 255)

     loadct, 0
     disp, map, max=0.25, min=-0.01 $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, star_mask, lev=[1,11], /overplot, color=cgcolor('blue', 255)

     loadct, 0
     disp, map*(star_mask ne 1) $
           , max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5

     loadct, 0
     disp, map*(star_mask ne 1) $
           , max=0.25, min=-0.01 $
           , /sq, xstyle=1, ystyle=1, reserve=5

     if keyword_set(pause) then begin
        ch = get_kbrd(1)
     endif

     !p.multi=0

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if n_elements(outfile) gt 0 then begin
     sxaddpar, hdr, 'BUNIT', 'MASK'
     writefits, outfile, star_mask, hdr
  endif

end
