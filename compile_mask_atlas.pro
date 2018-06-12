pro compile_mask_atlas $
   , rad = do_rad_mask $
   , gal = do_gal_mask $
   , bright = do_bright_star_mask $
   , unsharp = do_unsharp_mask $
   , show = show $
   , pause = pause $
   , only = only $
   , just = just $
   , tag = tag $
   , start = start_num $
   , stop = stop_num $
   , incremental = incremental

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  unwise_dir = '../unwise/atlas/' 
  galex_dir = '../galex/atlas/'
  out_dir = '../masks/'

  build_galaxy_list $
     , in_dir = in_dir $
     , tag=tag $
     , just=only $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = gal_data $
     , start = start_num $
     , stop = stop_num $
     , exclude = ['PGC17223']
  
  n_pgc = n_elements(pgc_list)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE RADIAL MASKS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  if keyword_set(do_rad_mask) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = strcompress(pgc_list[ii], /rem)
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Radial mask construction for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        for jj = 0, 1 do begin

           if jj eq 0 then res = 'gauss15'
           if jj eq 1 then res = 'gauss7p5'

           aperture_mask_file = out_dir+pgc_name+'_'+res+'_rgrid.fits'
           
           if keyword_set(incremental) then begin
              test = file_search(aperture_test_file, count=test_ct)
              if test_ct gt 0 then continue
           endif
           
           infile = unwise_dir+pgc_name+'_w1_'+res+'_small.fits'
           if file_test(infile) eq 0 then $
              infile = unwise_dir+pgc_name+'_w2_'+res+'_small.fits'
           if file_test(infile) eq 0 then begin
              print, "No WISE1 or WISE2. Stopping."
           endif
           
           build_aperture_mask $
              , pgc = this_dat.pgc $
              , galdata = this_dat $
              , infile = infile $
              , outfile = aperture_mask_file $
              , override_pgc_list = override_pgc_list $
              , override_rad_list = override_rad_list $
              , override_pa_list = override_pa_list $
              , override_incl_list = override_incl_list $              
              , show = show $
              , pause = pause
        
        endfor

     endfor     

  endif  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE GALAXY OVERLAP MASKS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  if keyword_set(do_gal_mask) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = strcompress(pgc_list[ii], /rem)
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Galaxy mask construction for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        for jj = 0, 1 do begin

           if jj eq 0 then res = 'gauss15'
           if jj eq 1 then res = 'gauss7p5'

           galaxy_mask_file = out_dir+pgc_name+'_'+res+'_galaxies.fits'
           
           if keyword_set(incremental) then begin
              test = file_search(galaxy_test_file, count=test_ct)
              if test_ct gt 0 then continue
           endif
           
           infile = unwise_dir+pgc_name+'_w1_'+res+'_small.fits'
           if file_test(infile) eq 0 then $
              infile = unwise_dir+pgc_name+'_w2_'+res+'_small.fits'
           if file_test(infile) eq 0 then begin
              print, "No WISE1 or WISE2. Stopping."
           endif

           build_galaxy_mask $
              , infile = infile $
              , skip_pgc = [this_dat.pgc] $
              , outfile = galaxy_mask_file $
              , galdata = all_gals $
              , rad_to_blank = 1.0 $
              , n_found = n_found $
              , pgc_found = overlap_pgc $
              , show=show $
              , pause=pause
           
;          SAVE AN IDL FILE LISTING THE GALAXIES THAT OVERLAP THIS FIELD
           if n_found gt 0 then begin
              outfile = out_dir+pgc_name+'_galaxies.idl'
              save, file=outfile, overlap_pgc
           endif
        
        endfor

     endfor  

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE BRIGHT STAR MASKS BASED ON 2MASS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  if keyword_set(do_bright_star_mask) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = strcompress(pgc_list[ii], /rem)
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Bright star mask construction for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        for jj = 0, 5 do begin

           atlas_dir = unwise_dir
           if jj eq 0 then band = 'w1'
           if jj eq 1 then band = 'w2'
           if jj eq 2 then band = 'w3'
           if jj eq 3 then band = 'w4'

           if jj eq 4 then band = 'nuv'
           if jj eq 5 then band = 'fuv'
           if jj eq 4 or jj eq 5 then atlas_dir = galex_dir
           
           for mm = 0, 1 do begin

              if mm eq 0 then res = 'gauss15'
              if mm eq 1 then res = 'gauss7p5'
                         
              infile = atlas_dir+pgc_name+'_'+band+'_'+res+'_small.fits'
              if file_test(infile) eq 0 then begin
                 print, "File missing, proceeding ... ", infile
                 continue
              endif
              
              bright_star_mask_file = $
                 out_dir+pgc_name+'_'+band+'_'+res+'_bright_stars.fits'
              
              build_bright_star__mask $
                 , infile = infile $
                 , outfile = bright_star_mask_file $     
                 , band = band $
                 , star_tol = star_tol $
                 , star_ra = star_ra $
                 , star_dec = star_dec $
                 , star_km = star_km $
                 , n_found = n_found $
                 , ra_found = overlap_ra $
                 , dec_found = overlap_dec $
                 , km_found = overlap_km $
                 , show=show $
                 , pause=pause           
              
;             SAVE AN IDL FILE LISTING THE STARS THAT OVERLAP THIS FIELD
              if jj eq 0 then begin
                 if n_found gt 0 then begin
                    outfile = out_dir+pgc_name+'_bright_stars.idl'
                    save, file=outfile, overlap_ra, overlap_dec, overlap_km
                 endif
              endif

           endfor

        endfor
        
     endfor     

  endif  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE UNSHARP MASKS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  if keyword_set(do_unsharp_mask) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = strcompress(pgc_list[ii], /rem)
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Unsharp star mask construction for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

;       FIND THE STARS IN THE WISE1 7.5" IMAGE

        
        
;       LOOP OVER BANDS TO REMOVE THE STARS

        for jj = 0, 5 do begin

           atlas_dir = unwise_dir
           if jj eq 0 then band = 'w1'
           if jj eq 1 then band = 'w2'
           if jj eq 2 then band = 'w3'
           if jj eq 3 then band = 'w4'

           if jj eq 4 then band = 'nuv'
           if jj eq 5 then band = 'fuv'
           if jj eq 4 or jj eq 5 then atlas_dir = galex_dir
           
           for mm = 0, 1 do begin

              if mm eq 0 then res = 'gauss15'
              if mm eq 1 then res = 'gauss7p5'
                         
              infile = atlas_dir+pgc_name+'_'+band+'_'+res+'_small.fits'
              if file_test(infile) eq 0 then begin
                 print, "File missing, proceeding ... ", infile
                 continue
              endif

           endfor
           
        endfor

     endfor

  endif


; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE CUSTOM MASKS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%


end
