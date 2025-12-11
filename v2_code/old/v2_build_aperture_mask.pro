pro v2_build_aperture_mask $
   , infile = infile $
   , map = map $
   , hdr = hdr $
   , ra = ra $
   , dec = dec $
   , posang = pa $
   , incl = incl $
   , fid_rad = fid_rad $
   , outfile = outfile $ 
   , show = show $
   , pause = pause

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
  
  nan = !values.f_nan

  min_rad = 30./3600.
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if n_elements(map) eq 0 or n_elements(hdr) eq 0 then begin
     if file_test(infile) eq 0 then begin
        message, 'File not found and map/header not supplied.', /info
        return
     endif     
     map = readfits(infile, hdr, /silent)  
  endif
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD THE COORD MAP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

; Build the coordinates
  
  make_axes, hdr, ri=ri, di=di

; Default to face on
  
  if n_elements(incl) eq 0 or n_elements(pa) eq 0 then begin
     pa = 0.
     incl = 0.
  endif

  if finite(incl) eq 0 or finite(pa) eq 0 then begin
     pa = 0.
     incl = 0.0           
  endif

; Cap the inclination at 60 for aperture purposes
  
  if incl gt 60 then begin
     incl = 60.
  endif

; Default to the center of the image
  
  if n_elements(ra) eq 0 or n_elements(dec) eq 0 then begin
     ra = mean(ri)
     dec = mean(di)
  endif

  if finite(ra) eq 0 or finite(dec) eq 0 then begin
     ra = mean(ri)
     dec = mean(di)
  endif

; Work out the galactocentric radius from the orientation and center
  
  gal_vec = [pa, incl, ra, dec]
  deproject, ri, di, gal_vec, rgrid=rgrid

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; NOTE THE FIDUCIAL RADIUS IN THE HEADER
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Pick the radius to define the area around the galaxy. This will be
; avoided in a background subtraction.

  if n_elements(fid_rad) eq 0 then begin
     fid_rad = min_rad
  endif
 
  if fid_rad lt min_rad or finite(fid_rad) eq 0 then begin
     fid_rad = min_rad
  endif

  sxaddpar, hdr, 'PA', pa, 'Adopted position angle [deg]'
  sxaddpar, hdr, 'INCL', incl, 'Adopted inclination [deg]'
  sxaddpar, hdr, 'FIDRAD', fid_rad, 'Fiducial radius [deg]'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if keyword_set(show) then begin

     !p.multi=0
     
     rms = mad(map[where(finite(map))])
     disp, map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=5, ystyle=5
     contour, rgrid, lev=[1,2,3,4]*fid_rad $
              , /overplot, color=cgcolor('blue')
     if keyword_set(pause) then $
        ch = get_kbrd(1)

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if n_elements(outfile) gt 0 then begin
     sxaddpar, hdr, 'BUNIT', 'DEG'
     writefits, outfile, rgrid, hdr
  endif

end
