pro find_point_sources $
   , infile = infile $
   , map = map $
   , hdr = hdr $
   , outfile=outfile $
   , fwhm=beam_deg $
   , exclude_ra=exclude_ra $
   , exclude_dec=exclude_dec $
   , exclude_tol=exclude_tol $
   , n_star=n_star $
   , star_ra=star_ra $
   , star_dec=star_dec $
   , star_intens=star_intens $
   , star_km=star_km $
   , show=show $
   , pause=pause

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

  sz = size(map)
  make_axes, hdr, ri=ri, di=di, ra=ra, da=da
  pix_deg = sphdist(ri[0,0], di[0,0], ri[1,0], di[1,0], /deg)
  if n_elements(beam_deg) eq 0 then $
     beam_deg = sxpar(hdr, 'BMAJ')
  beam_pix = beam_deg/pix_deg
  fwhm = beam_pix

  nan = !values.f_nan

  if n_elements(band) eq 0 then $
     band = 'w1'

  if n_elements(thresh) eq 0 then $
     thresh = 0.1

  if n_elements(sharp) eq 0 then $
     sharplim = [0.2, 1.0]

  if n_elements(round) eq 0 then $
     roundlim = [-1.0, 1.0]
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RUN THE SOURCE FINDER
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  find, map, star_x, star_y, star_flux, star_sharp, star_round $
        , thresh, fwhm, roundlim, sharplim $
        , /silent

  if n_elements(star_x) eq 0 then begin
     print, "No stars found. Returning."
     n_star = 0     
     return
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVERT TO RA, DEC AND EXTRACT INTENSITY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  xyad, hdr, star_x, star_y, star_ra, star_dec

  if n_elements(exclude_ra) gt 0 then begin

     keep = finite(star_ra)
     
     if n_elements(exclude_tol) eq 0 then $
        exclude_tol = exclude_ra*0.0+sxpar(hdr, 'BMAJ')

     for kk = 0, n_elements(exclude_ra)-1 do begin
        dist = sphdist(star_ra, star_dec, exclude_ra[kk], exclude_dec[kk], /deg)
        bad_ind = where(dist le exclude_tol, bad_ct)
        if bad_ct gt 0 then keep[bad_ind] = 0B
     endfor

     ind = where(keep, ct)

     if ct eq 0 then begin
        print, "All stars excluded. Returning."
        n_star = 0
        return
     endif

     star_x = star_x[ind]
     star_y = star_y[ind]
     star_ra = star_ra[ind]
     star_dec = star_dec[ind]

  endif

  n_star = n_elements(star_x)
  print, "Found "+str(n_star)+" stars."

  star_intens = map[round(star_x), round(star_y)]
  
  if band eq 'w1' then begin
     coef = 52e3
  endif else if band eq 'w2' then begin
     coef = 28e3
  endif else if band eq 'w3' then begin
     coef = 4.8e3
  endif else if band eq 'w4' then begin
     coef = 1.5e3
  endif else if band eq 'nuv' then begin
     coef = 38.
  endif else if band eq 'fuv' then begin
     coef = 0.4
  endif else begin
     print, "Band not recognized, using band w1."
     coef = 52e3
  endelse

  star_km = -1.*2.5*alog10(star_intens/coef)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if keyword_set(show) then begin

     !p.multi=0
     loadct, 0

     rms = 0.1

     disp, map, ra, da, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5, color=cgcolor('white',255), /radec

     oplot, star_ra, star_dec, psym=cgsymcat('filledstar') $
            , color=cgcolor('orange'), symsize=2.0
     
     if keyword_set(pause) then ch = get_kbrd(1)

  endif

end
