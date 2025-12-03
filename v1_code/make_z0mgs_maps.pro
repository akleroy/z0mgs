pro make_z0mgs_maps $
   , do_m31 = do_m31 $
   , just_pgc=just_pgc $
   , skip_pgc=skip_pgc $   
   , skip_name=skip_name $   
   , start_name=start_name $
   , start_num=start_num $   
   , pause=pause $
   , bands=user_bands $
   , skip_user_mask = skip_user_mask

; Ported back from PHANGS survey work. Process the delivered z0mgs
; GALEX and WISE data and masks to create intensity images to be used
; calculating SFR and Mstar.

; &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
;
; Known problems (from PHANGS survey work)
;
; - a few galaxies still have centers getting flagged. Work on this
;
; - at W1/W2 stars overwhelm some images. Drop to a median profile
;   in the worst cases.
;
; - PSF / Saturation biggest concern at w3/w4. Flagged in catalog, but
;   probably need something for PHANGS, since the density of these
;   sources are highest for us.
;
; - FUV looks quite good overall, a few ring artifacts
;
; - NUV has some star problems, m33 has very weird artifacts in NUV
; 
; &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

; &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
; COMPILE DIRECTORIES AND LIST OF TARGETS
; &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

  tab = mrdfits('../../measurements/delivery_index_gauss15.fits',1,h)
  dat = gal_data(pgc=tab.pgc)
  n_gals = n_elements(tab)
  
  dir_in = '../../delivery/'
  dir_out = '../../dr1_sfrmstar/'
  imdir = dir_out+'qa/'
  user_mask_dir = dir_out+'user_masks/'

  bands = ['w1','w2','w3','w4','nuv','fuv']
  n_bands = 6

; will need to move this to a list
  w1w2_failures = $
     ['ngc1809','eso097-013','ngc2283','ngc2566','ngc3059' $
      ,'ngc5530','ngc6300','ngc6744','ngc6946']

; NGC4945 - need to deal with edge on case

; &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
; LOOP OVER GALAXIES
; &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

  if n_elements(start_num) eq 0 then $
     start_num = 0
  
  first = 1
  for ii = start_num, n_gals-1 do begin

     this_tab = tab[ii]
     this_dat = dat[ii]
     
     if keyword_set(do_m31) eq 0 then $
        if dat[ii].pgc eq 2557 then $
           continue

     phangs_name = $
        strlowcase(strcompress(dat[ii].objname, /rem))

     pgc_name = strcompress(dat[ii].pgcname, /rem)
     
     if n_elements(just) gt 0 then begin
        if total(just eq phangs_name) eq 0 then $
           continue
     endif
     
     if n_elements(start_name) gt 0 and first eq 1 then begin
        if (total(start_name eq phangs_name) eq 0 or $
            total(start_name eq pgc_name) eq 0) then $           
               continue $
        else $
           first = 0
     endif
     
     if n_elements(just_pgc) gt 0 then begin
        if total(just_pgc eq dat.pgc) eq 0 then $
           continue
     endif
     
     if n_elements(skip_pgc) gt 0 then begin
        if total(skip_pgc eq dat.pgc) gt 0 then $
           continue        
     endif
     
     print, '-----------------------------------------'
     print, 'Processing '+pgc_name+' or '+phangs_name
     print, 'which is galaxy '+str(ii)+' of '+str(n_gals)
     print, '-----------------------------------------'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&
; LOOP OVER RESOLUTIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     for kk = 0, 1 do begin

        if kk eq 1 then begin
           res = 'gauss7p5'
           beam = 7.5
        endif

        if kk eq 0 then begin
           res = 'gauss15'
           beam = 15.
        endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&
; LOOP OVER BANDS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        for jj = 0, n_bands-1 do begin

           if n_elements(user_bands) gt 0 then $
              if total(bands[jj] eq user_bands) eq 0 then $
                 continue

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&
; READ FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; The map
           this_file = $
              dir_in + pgc_name + $
              '_'+bands[jj]+'_'+res+'.fits'
           if file_test(this_file) eq 0 then $
              continue          
           map = readfits(this_file, hdr)

; The star mask
           star_file = $
              dir_in + pgc_name + $
              '_'+bands[jj]+'_'+res+'_stars.fits'
           stars = readfits(star_file, star_hdr)

; Any user mask
           user_mask_file = $
              user_mask_dir + pgc_name + $
              '_'+bands[jj]+'_'+res+'_user_mask.fits'
           has_user_mask = 0B
           if file_test(user_mask_file) then begin
              has_user_mask = 1B
              user_mask = readfits(user_mask_file, user_mask_hdr)
          endif
           
; The galaxy mask
           galaxy_file = $
              dir_in + pgc_name + $
              '_'+res+'_galaxies.fits'
           has_gal_mask = 0B
           if file_test(galaxy_file) then begin
              has_gal_mask = 1B
              galaxies = readfits(galaxy_file, gal_hdr)
           endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD GALACTOCENTRIC AND ON THE SKY RADIUS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

           make_axes, hdr, ri=ri, di=di
           this_incl = dat[ii].incl_deg
           if this_incl gt 80. then $
              this_incl = 80.
           this_pa = dat[ii].posang_deg
           if finite(this_pa) eq 0 then begin
              this_pa = 0.0
              this_incl = 0.0
           endif
           this_ra_ctr = dat[ii].ra_deg
           this_dec_ctr = dat[ii].dec_deg
           deproject, ri, di, [this_pa, this_incl, this_ra_ctr, this_dec_ctr] $
                      , rgrid=rgrid, tgrid=tgrid
           flat_rgrid = sphdist(this_ra_ctr, this_dec_ctr, ri, di, /deg)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BLANK OUT STARS, ALWAYS LEAVE THE CENTER UNMASKED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

           blank_map = map
           blank_ind = $
              where(stars and $
                    flat_rgrid gt 60./3600. and $
                    finite(map))
           blank_map[blank_ind] = !values.f_nan
         
           if has_user_mask and not keyword_set(skip_user_mask)then begin
              user_ind = where(user_mask, user_ct)
              if user_ct gt 0 then $
                 blank_map[user_ind] = !values.f_nan
              blank_ind = where(finite(map) and finite(blank_map) eq 0)
           endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD A RADIAL PROFILE STAYING NEAR THE MAJOR AXIS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

           if this_incl lt 80. then begin
              cos_thresh = 0.0
           endif
           if this_incl ge 80. then begin
              cos_thresh = 0.5
           endif

           this_r25 = dat[ii].r25_deg
           exclusion = (2.0*beam/3600 > 0.1*this_r25)
           mask = abs(cos(tgrid)) gt cos_thresh or rgrid le exclusion
           mask_ind = where(mask, mask_ct)
           if mask_ct gt 0 then begin

              xmin = -0.5 * beam/3600.
              xmax = max(rgrid)
              binsize = beam/3600.
              bins = $
                 bin_data(rgrid[mask_ind], blank_map[mask_ind], /nan $
                          , xmin=xmin, xmax=xmax, binsize=binsize)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INTERPOLATE TO FILL IN BLANKED PIXELS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

              interp_map = blank_map
              interp_map[blank_ind] = $
                 interpol(bins.ymed, bins.xmid, rgrid[blank_ind])
              
              if (total(phangs_name eq w1w2_failures) gt 0 or $
                  total(pgc_name eq w1w2_failures) gt 0) and $              
                 (bands[jj] eq 'w1' or bands[jj] eq 'w2') then begin
                 
                 bins = $
                    bin_data(rgrid, map, /nan $
                             , xmin=xmin, xmax=xmax, binsize=binsize)
                 
                 interp_map = interpol(bins.ymed, bins.xmid, rgrid)
              endif
              
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE A RADIAL PROFILE IMAGE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

              rprof_map = $
                 interpol(bins.ymed, bins.xmid, rgrid)
              resid_map = $
                 interp_map - rprof_map

           endif else begin

              rprof_map = map
              resid_map = map*0.0
              
           endelse
           
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DISPLAY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

           aparm=0.001
           title = pgc_name+'_'+bands[jj]+'_'+res

           !p.multi=[0,3,2]
           loadct, 0
           ;disp, alog10(map), title=title+' original'
           disp, asinh(map/aparm)/asinh(1.0/aparm), title=title+' original'
           contour, stars, lev=[1], /overplot, c_color=getcolor('red')
           if has_user_mask then begin
              contour, user_mask, lev=[1], /overplot, thick=3, c_color=getcolor('yellow')
              contour, user_mask, lev=[1], /overplot, c_color=getcolor('blue')
           endif

           loadct, 0
           ;disp, alog10(blank_map), title=title+' blanked'
           disp, asinh(blank_map/aparm)/asinh(1.0/aparm), title=title+' blanked'
           if has_user_mask then begin           
              contour, user_mask, lev=[1], /overplot, thick=3, c_color=getcolor('yellow')
              contour, user_mask, lev=[1], /overplot, c_color=getcolor('blue')
           endif

           loadct, 0
           ;disp, alog10(interp_map), title=title+' interpolated'
           disp, asinh(interp_map/aparm)/asinh(1.0/aparm), title=title+' interpolated'
           loadct, 0
           disp, asinh(rprof_map/aparm)/asinh(1.0/aparm), title=title+' radial only'
           loadct, 0
           disp, asinh(resid_map/aparm)/asinh(1.0/aparm), title=title+' radial resids'

           ;if has_gal_mask then $
           ;   contour, galaxies, lev=[1], /overplot, c_color=getcolor('yellow')
           !p.multi=[0,1,1]
           
           im = tvrd(true=1)
           image_fname = imdir+title+'.png'
           write_png, image_fname, im

           if keyword_set(pause) then begin
              print, "Hit a key to continue."
              ch = get_kbrd(1)
           endif           

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

           this_out_file = $
              dir_out + pgc_name + $
              '_'+bands[jj]+'_'+res+'.fits'
           writefits, this_out_file, map, hdr

           this_out_file = $
              dir_out + pgc_name + $
              '_'+bands[jj]+'_'+res+'_blanked.fits'
           writefits, this_out_file, blank_map, hdr

; The 15" is too poor in some cases. Skip this band-res and instead
; convolve the interpolated 7.5" to coarser resolution.

           if bands[jj] eq 'w1' and res eq 'gauss15' then $
              continue

           this_out_file = $
              dir_out + pgc_name + $
              '_'+bands[jj]+'_'+res+'_interpol.fits'
           writefits, this_out_file, interp_map, hdr

           this_out_file = $
              dir_out + pgc_name + $
              '_'+bands[jj]+'_'+res+'_rprof.fits'
           writefits, this_out_file, rprof_map, hdr

           this_out_file = $
              dir_out + pgc_name + $
              '_'+bands[jj]+'_'+res+'_resid.fits'
           writefits, this_out_file, resid_map, hdr

           if bands[jj] eq 'w1' and res eq 'gauss7p5' then begin

; Convolve the interpolated map
              
              this_in_file = $
                 dir_out + pgc_name + $
                 '_'+bands[jj]+'_gauss7p5_interpol.fits'
              
              new_out_file =  $
                 dir_out + pgc_name + $
                 '_'+bands[jj]+'_gauss15_interpol.fits'

              conv_with_gauss $
                 , data=this_in_file $
                 , target_beam = 15.*[1,1,0] $
                 , out_data = map_out $
                 , out_hdr = hdr_out

              target_hdr = $
                 headfits(dir_out + pgc_name + $
                          '_'+bands[jj]+'_gauss15.fits')

              hastrom, map_out, hdr_out, target_hdr $
                       , cubic=-0.5, missing=!values.f_nan, interp=2
              writefits, new_out_file, map_out, hdr_out

; Convolve the radial profile map

              this_in_file = $
                 dir_out + pgc_name + $
                 '_'+bands[jj]+'_gauss7p5_rprof.fits'
              
              new_out_file =  $
                 dir_out + pgc_name + $
                 '_'+bands[jj]+'_gauss15_rprof.fits'

              conv_with_gauss $
                 , data=this_in_file $
                 , target_beam = 15.*[1,1,0] $
                 , out_data = map_out $
                 , out_hdr = hdr_out

              target_hdr = $
                 headfits(dir_out + pgc_name + $
                          '_'+bands[jj]+'_gauss15.fits')
              
              hastrom, map_out, hdr_out, target_hdr $
                       , cubic=-0.5, missing=!values.f_nan, interp=2
              writefits, new_out_file, map_out, hdr_out


           endif
           
        endfor
        
     endfor
     
  endfor

end
