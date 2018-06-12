pro analyze_profiles $
   , show = show $
   , pause = pause $
   , just = just $
   , tag = tag $
   , start = start_num $
   , stop = stop_num $
   , galex = do_galex $
   , wise = do_wise $
   , isophot = do_isophot $
   , profile = do_profile

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  unwise_dir = '../unwise/atlas/'
  galex_dir = '../galex/atlas/'
  prof_dir = '../measurements/profiles/'

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  out_dir = '../unwise/atlas/'

  build_galaxy_list $
     , in_dir = in_dir $
     , tag=tag $
     , just=just $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = gal_data $
     , start = start_num $
     , stop = stop_num $
     , exclude = ['PGC17223']
     
  n_pgc = n_elements(pgc_list)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; 
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP AND ANALYZE PROFILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(do_profile) then begin

     if keyword_set(do_wise) then begin
        prof_struct = analyze_profile(infile='', /empty)
        wise_phot = replicate(prof_struct, 4, n_pgc)
     endif

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        print, ''
        print, 'Analyzing radial profiles for for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        if keyword_set(do_wise) then begin
           
           for band = 1, 4 do begin

              prof_file = prof_dir+pgc_name+'_w'+str(band)+'_prof.txt'
              test = file_search(prof_file, count=ct)
              if ct eq 0 then begin
                 print, "File "+prof_file+" not found."
                 continue
              endif
              
              prof_struct = $
                 analyze_profile( $
                 infile=prof_file)
              wise_phot[band-1,ii] = prof_struct

           endfor

        endif

     endfor

     if keyword_set(do_wise) then begin
        save, file='../measurements/wise_phot.idl', wise_phot
     endif

  endif
  
  stop

end
