pro print_missing_From_galbase

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  galex_dir = '../galex/atlas/'
  wise_dir = '../unwise/atlas/'
  out_dir = '../delivery/'

  build_galaxy_list $
     , in_dir = in_dir $
     , tag=tag $
     , just=only $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = dat $
     , start = start_num $
     , stop = stop_num $
     , exclude = ['PGC17223']
  n_pgc = n_elements(pgc_list)

  s = gal_data(pgc=pgc_num)
    
  ind = where(s.pgc eq 0, ct) 
  if ct eq 0 then begin
     print, 'No galaxies missing.'
     return
  endif

  for ii = 0, ct-1 do $
     print, pgc_list[ind[ii]]

end
