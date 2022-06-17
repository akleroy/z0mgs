pro build_aperture_mask $
   , pgc = pgc_num $
   , galdata = this_dat $
   , infile = infile $
   , map = map $
   , hdr = hdr $
   , outfile = outfile $
   , mask = mask $
   , force_rad = override_rad $
   , force_pa = override_pa $
   , force_incl = override_incl $
   , override_pgc_list = override_pgc_list $
   , override_rad_list = override_rad_list $
   , override_pa_list = override_pa_list $
   , override_incl_list = override_incl_list $
   , show = show $
   , pause = pause

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

  if n_elements(override_rad) eq 0 then $
     override_rad = nan
  if n_elements(override_pa) eq 0 then $
     override_pa = nan
  if n_elements(override_incl) eq 0 then $
     override_incl = nan

  if n_elements(override_pgc_list) eq 0 then begin
     if file_test('custom_mask_specs.txt') then begin
        override_file_found = 1B
        readcol, 'custom_mask_specs.txt', format='L,F,F,F', comment='#' $
                 , override_pgc_list, override_rad_list, override_pa_list, override_incl_list
        ind = where(override_pgc_list eq pgc_num, ct)
        if ct eq 1 then begin
           if n_elements(override_rad) eq 0 then $
              override_rad = override_rad_list[ind]
           if n_elements(override_pa) eq 0 then $
              override_pa = override_pa_list[ind]
           if n_elements(override_incl) eq 0 then $
              override_incl = override_incl_list[ind]
        endif
     endif else begin
        override_file_found = 0B
     endelse
  endif
  
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
  mask = finite(map)*0B
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD THE RADIUS MAP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

; Work out coordinates

  make_axes, hdr, ri=ri, di=di

; Set inclination and position angle, default to face on lacking this
; information. Even for edge on systems, assume an inclination of ~60
; to avoid issues with the thickness of the disk.

  pa = this_dat.posang_deg
  incl = this_dat.incl_deg        

  if finite(pa) then $
     sxaddpar, hdr, 'GBPA', pa, 'Database position angle [deg]'

  if finite(incl) then $
     sxaddpar, hdr, 'GBINCL', incl, 'Database inclination [deg]'

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

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; NOTE THE FIDUCIAL RADIUS IN THE HEADER
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

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

  sxaddpar, hdr, 'PA', pa, 'Adopted position angle [deg]'
  sxaddpar, hdr, 'INCL', incl, 'Adopted inclination [deg]'

  sxaddpar, hdr, 'FIDRAD', fid_rad, 'Fiducial radius [deg]'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if keyword_set(show) then begin

     !p.multi=0
     
     print, 'PGC'+strcompress(str(pgc_num),/rem)
     print, "Current radius ", fid_rad
     print, "Current position angle ", pa
     print, "Current inclination ", incl

     rms = mad(map[where(finite(map) and mask ne 10)])
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
