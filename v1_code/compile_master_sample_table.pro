pro compile_master_sample_table

  in_dir = '../measurements/samples/'

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
     print, 'Reading samples for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
     print, ''

     tab = mrdfits(in_dir + pgc_name + '_samples.fits', 1, h)

     stop

  endfor

  !p.multi=0

end


end
