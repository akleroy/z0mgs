pro build_unsharp_mask $
   , pgc=pgc_num $
   , galdata=this_dat $
   , wise1=wise1 $
   , hdr=hdr $
   , mask=output_mask $
   , star_ra=unsharp_ra $
   , star_dec=unsharp_dec $
   , star_intens=unsharp_intens

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
; IDENTIFY THE REGION NEAR THE GALAXY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  sz = size(wise1)
  mask = finite(wise1)*0B

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
; IF REQUESTED, ALSO BLANK AUTOMATICALLY FOUND STARS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  output_mask = mask*0B

  map = wise1
  thresh = 0.1
   
  unsharp = (map - median(map, 5)) 
  unsharp_mask = (unsharp gt thresh)*(rgrid gt fid_rad) + $
                 (unsharp gt 10.*thresh)*(rgrid lt fid_rad)
  unsharp_mask *= (rgrid gt fid_rad)
  
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
              output_mask[blank_ind] = 1B
           endif
           
        endfor
        
     endif
     
  endif

  mask = unsharp_mask

end
