pro compare_flags $
   , reason=reason

  if n_elements(reason) eq 0 then $
     reason = 'radius'

  readcol, comment='#', 'akl_galaxies_to_flag.txt', format='A,A' $
           , pgc_akl, flag_akl
  pgc_akl = strcompress(pgc_akl, /rem)
  flag_akl = strcompress(flag_akl, /rem)

  readcol, comment='#', 'mg_galaxies_to_flag.txt', format='A,A' $
           , pgc_mg, flag_mg
  pgc_mg = strcompress(pgc_mg, /rem)
  flag_mg = strcompress(flag_mg, /rem)
  
  ind_akl = where(flag_akl eq reason, akl_ct)
  ind_mg = where(flag_mg eq reason, mg_ct)

  for ii = 0, akl_ct-1 do begin
     mg_agree = total(pgc_akl[ind_akl[ii]] eq pgc_mg[ind_mg]) eq 1 ? "yes" : "no"
     print, pgc_akl[ind_akl[ii]], ' Molly agrees? ', mg_agree
  endfor

  print, '------------------'
  print, '------------------'
  print, '------------------'

  for ii = 0, mg_ct-1 do begin
     akl_agree = total(pgc_mg[ind_mg[ii]] eq pgc_akl[ind_akl]) eq 1 ? "yes" : "no"
     print, pgc_mg[ind_mg[ii]], ' Adam agrees? ', akl_agree
  endfor

  print, '------------------'
  print, reason
  print, '------------------'

  pgc = [pgc_mg[ind_mg], pgc_akl[ind_akl]]
  pgc = pgc[uniq(pgc, sort(pgc))]
  pgc = pgc[sort(pgc)]
  for ii = 0, n_elements(pgc)-1 do begin
     print, pgc[ii]
  endfor
  
  stop

end
