pro build_mask $
   , pgc=pgc_num $
   , mask=mask $
   , galdata=this_dat $
   , headerfile=headerfile $
   , map=map $
   , hdr=hdr $
   , outfile=outfile $
   , inner=inner $
   , outer=outer $
   , show=show $
   , pause=pause $
   , band=band $
   , do_flag_gals = do_flag_gals $
   , allgaldata = allgals $
   , galrad = rad_galflag $
   , do_flag_stars = do_flag_stars $
   , star_tol = star_tol $
   , star_ra = star_ra $
   , star_dec = star_dec $
   , star_km = star_km $
   , do_flag_sharp = do_flag_sharp

  if n_elements(this_dat) eq 0 then begin
     this_dat = gal_data(pgc=pgc_num)
  endif

  if n_elements(inner) eq 0 then begin
     inner = 2.0
  endif

  if n_elements(outer) eq 0 then begin
     outer = 4.0
  endif
  
  nan = !values.f_nan
  override_rad = nan
  override_pa = nan
  override_incl = nan
  if file_test('custom_mask_specs.txt') then begin
     override_file_found = 1B
     readcol, 'custom_mask_specs.txt', format='L,F,F,F', comment='#' $
              , override_pgc_list, override_rad_list, override_pa_list, override_incl_list
     ind = where(override_pgc_list eq pgc_num, ct)
     if ct eq 1 then begin
        override_rad = override_rad_list[ind]
        override_pa = override_pa_list[ind]
        override_incl = override_incl_list[ind]
     endif
  endif else begin
     override_file_found = 0B
  endelse
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  

  if n_elements(map) eq 0 or n_elements(hdr) eq 0 then begin
     if file_test(headerfile) eq 0 then begin
        message, 'Target astrometry file not found and map and header not supplied.', /info
        return
     endif
     
     map = readfits(headerfile, hdr, /silent)  
  endif

  sz = size(map)
  mask = finite(map)*0B
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IDENTIFY THE REGION NEAR THE GALAXY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

; Work out coordinates
  make_axes, hdr, ri=ri, di=di

; Set inclination and position angle, default to face on lacking this
; information. Even for edge on systems, assume an inclination of ~60
; to avoid issues with the thickness of the disk.

  pa = this_dat.posang_deg
  incl = this_dat.incl_deg        

  if n_elements(override_pa) gt 0 then begin
     if finite(override_pa) then $
        pa = override_pa
  endif

  if n_elements(override_incl) gt 0 then begin
     if finite(override_incl) then $
        incl = override_incl
  endif

  if finite(incl) eq 0 then begin
     incl = 0.0           
  endif

  if incl gt 60 then begin
     incl = 60.
  endif

  if finite(pa) eq 0 then begin
     pa = 0.0
     incl = 0.0
  endif
  
; Work out the galactocentric radius from the orientation and center

  xctr = this_dat.ra_deg
  yctr = this_dat.dec_deg

  gal_vec = [pa, incl, xctr, yctr]
  deproject, ri, di, gal_vec, rgrid=rgrid

; Pick the radius to define the area around the galaxy. This will be
; avoided in a background subtraction.

  fid_rad = this_dat.r25_deg        
  if n_elements(override_rad) gt 0 then begin
     if finite(override_rad) then $
        fid_rad = override_rad
  endif

  if fid_rad lt 30./3600. or finite(fid_rad) eq 0 then begin
     fid_rad = 30./3600.
  endif

; Label inner to outer r25 background, < inner galaxy

  gal_ind = where(rgrid lt fid_rad*inner, gal_ct)
  mask[gal_ind] = 10B
  bk_ind = where(rgrid ge fid_rad*inner and $
                 rgrid lt fid_rad*outer, bk_ct)
  mask[bk_ind] = 100B

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IF REQUESTED, BLANK ANY GALAXIES IN THE FIELD
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  gal_mask = mask*0B

  if keyword_set(do_flag_gals) then begin

;    IF WE DON'T HAVE THEM ALREADY LOAD THE GALAXIES FROM THE GALBASE
     if n_elements(allgals) eq 0 then begin
        allgals = gal_data(/all)
     endif

     max_ra = max(ri, /nan)
     min_ra = min(ri, /nan)
     max_dec = max(di, /nan)
     min_dec = min(di, /nan)

     if n_elements(rad_galflag) eq 0 then begin
        rad_galflag = 1.0
     endif
     
     tol = allgals.r25_deg
     ra_tol = tol/cos(!dtor*mean(di,/nan))

     in_image = $
        (allgals.ra_deg ge (min_ra-ra_tol)) and $
        (allgals.ra_deg le (max_ra+ra_tol)) and $
        (allgals.dec_deg ge (min_dec - ra_tol)) and $
        (allgals.dec_deg le (max_dec + ra_tol)) and $
        (allgals.pgc ne pgc_num)

     gal_ind = where(in_image, gal_ct)
     if gal_ct gt 0 then begin
        print, "Found ", gal_ct, " galaxies that I will blank."

        rad_to_blank = $
           allgals[gal_ind].r25_deg*rad_galflag
        
        for kk = 0, gal_ct-1 do begin
           
           this_ra = allgals[gal_ind[kk]].ra_deg
           this_dec = allgals[gal_ind[kk]].dec_deg
           this_rad = rad_to_blank[kk]
           
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
           
                                ;dist = sphdist(ri, di, this_ra, this_dec, /deg)
           blank_ind = where(this_rgrid le this_rad, blank_ct)
           if blank_ct gt 0 then begin
              gal_mask[blank_ind] = 1B
           endif
           
        endfor   
        
     endif

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IF REQUESTED, BLANK KNOWN STARS IN THE FIELD
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  star_mask = mask*0B

  if keyword_set(do_flag_stars) then begin

;    IF WE DON'T HAVE THEM ALREADY LOAD THE 2MASS STARS
     if n_elements(star_ra) eq 0 then begin
        restore, '../measurements/2mass_stars.idl', /v
     endif
     
     if n_elements(star_tol) eq 0 then $
        star_tol = 0.0
     
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

     star_ind = where(in_image, star_ct)
     if star_ct gt 0 then begin
        print, "Found ", star_ct, " stars that I will blank."

        intens = fltarr(star_ct)*nan
        adxy, hdr, star_ra[star_ind], star_dec[star_ind], star_x, star_y
        really_in_image = $
           ((star_x ge 0) and (star_y ge 0) and $
            (star_x le (sz[1]-1)) and (star_y le (sz[2]-1)))
        really_in_image_ind = where(really_in_image, really_ct)
        if really_ct gt 0 then $
           intens[really_in_image_ind] = $
           map[star_x[really_in_image_ind], star_y[really_in_image_ind]]
        
        rad_to_blank = $
           get_rad_to_blank(mag=star_km[star_ind], intens=intens, band=band)

        for kk = 0, star_ct-1 do begin
           
           this_ra = star_ra[star_ind[kk]]
           this_dec = star_dec[star_ind[kk]]
           this_rad = rad_to_blank[kk]

           dist_from_ctr = sphdist(this_ra, this_dec, this_dat.ra_deg, this_dat.dec_deg, /deg)
           if dist_from_ctr lt 10./3600. then begin
              print, "I consider this star very likely to be the center of the galaxy."
              continue
           endif

           dist = sphdist(ri, di, this_ra, this_dec, /deg)
           blank_ind = where(dist le this_rad, blank_ct)
           if blank_ct gt 0 then begin
              star_mask[blank_ind] = 1B
           endif

        endfor

     endif

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IF REQUESTED, ALSO BLANK AUTOMATICALLY FOUND STARS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  unsharp_mask = mask*0B

  if keyword_set(do_flag_sharp) then begin
     
     thresh = 0.1
     if band eq 'w1' then thresh = 0.05
     if band eq 'w2' then thresh = 0.05
     if band eq 'w3' then thresh = 0.1
     if band eq 'w4' then thresh = 1.0
     if band eq 'nuv' then thresh = 5d-3
     if band eq 'fuv' then thresh = 5d-3

     unsharp = (map - median(map, 11)) 
     unsharp_mask = (unsharp gt thresh)*(rgrid gt fid_rad*inner) + $
                    (unsharp gt 10.*thresh)*(rgrid lt fid_rad*inner)
     unsharp_mask *= (rgrid gt fid_rad)

     ;unsharp = unsharp*(mask ne 10)
     reg = label_region(unsharp_mask)
     n_reg = max(reg)
     print, "Found "+str(n_reg)+" features via unsharp masking."
     if n_reg gt 0 then begin

        unsharp_ra = fltarr(n_reg)*nan
        unsharp_dec = fltarr(n_reg)*nan
        unsharp_intens = fltarr(n_reg)*nan

        for ii = 0, n_reg-1 do begin
           unsharp_intens[ii] = max(map*(reg eq (ii+1)), maxind, /nan)
           unsharp_ra[ii] = ri[maxind]
           unsharp_dec[ii] = di[maxind]        
        endfor

        if n_reg gt 0 then begin        
           rad_to_blank = $
              get_rad_to_blank(intens=unsharp_intens, band=band)
           
           for kk = 0, n_reg-1 do begin
              
              this_ra = unsharp_ra[kk]
              this_dec = unsharp_dec[kk]
              this_rad = rad_to_blank[kk]

              dist = sphdist(ri, di, this_ra, this_dec, /deg)
              blank_ind = where(dist le this_rad, blank_ct)
              if blank_ct gt 0 then begin
                 unsharp_mask[blank_ind] = 1B
              endif

           endfor

        endif

     endif

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COMBINE MASKS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  mask = mask + 1B*star_mask + 2B*unsharp_mask + 3B*gal_mask

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if keyword_set(show) then begin

     !p.multi=0
     
     print, 'PGC'+strcompress(str(pgc_num),/rem)
     print, "Current radius ", fid_rad
     print, "Current position angle ", pa
     print, "Current inclination ", incl

                                ;rms = stddev(map[where(finite(map) and mask ne 10)],/nan)
     rms = mad(map[where(finite(map) and mask ne 10)])
     disp, map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=5, ystyle=5
     contour, mask, lev=[10,100], /overplot, color=cgcolor('blue')
     contour, star_mask, lev=[1], /overplot, color=cgcolor('red')
     contour, gal_mask, lev=[1], /overplot, color=cgcolor('orange')
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
