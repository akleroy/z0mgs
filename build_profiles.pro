pro build_profiles $
   , show = show $
   , pause = pause $
   , just = just $
   , tag = tag $
   , start = start_num $
   , stop = stop_num $
   , galex = do_galex $
   , wise = do_wise

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  prof_dir = '../measurements/profiles/'
  atlas_dir = '../delivery/'
  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'

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

  index = mrdfits('../measurements/delivery_index.fits',1,h)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP AND MEASURE RADIAL PROFILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  for ii = 0, n_pgc-1 do begin
     
     pgc_name = pgc_list[ii]
     this_dat = gal_data[ii]
     
     print, ''
     print, 'Radial profiles for for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
     print, ''

     outfile = prof_dir + pgc_name + '_profiles.fits'

     index_ind = where(index.pgc_name eq pgc_name)
     rings = $
        build_one_profile( $
        pgc_name = pgc_name $
        , gal_data = this_dat $
        , index = index[index_ind] $
        , data_dir = atlas_dir $
        , outfile=outfile $
        , bin_size=7.5 $
        , beam_size=15.0 $
        , /vary_orient $
        , show=show $
        , pause=pause)

  endfor

  !p.multi=0

  stop

end
