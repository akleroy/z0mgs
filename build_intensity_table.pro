pro build_intensity_table $
   , show = show $
   , pause = pause $
   , just = just $
   , tag = tag $
   , start = start_num $
   , stop = stop_num

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  out_dir = '../measurements/samples/'
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
     print, 'Sampling intensity for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
     print, ''

     outfile = out_dir + pgc_name + '_samples.fits'

     index_ind = where(strcompress(index.pgc_name, /rem) $
                       eq pgc_name, index_ct)
     if index_ct eq 0 then begin
        print, "Index mismatch. Stopping."
        stop
     endif
     sample_one_galaxy $
        , pgc_name = pgc_name $
        , gal_data = this_dat $
        , index = index[index_ind] $
        , data_dir = atlas_dir $
        , outfile = outfile $
        , spacing = 7.5 $
        , beam_size = 15.0

  endfor

  !p.multi=0

end
