pro z0mgs_list

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  
  flist = file_search(in_dir+'PGC*', count=file_ct)
  pgc_list = strarr(file_ct)
  pgc_num = lonarr(file_ct)
  for ii = 0, file_ct-1 do begin
     pgc_list[ii] = strmid(flist[ii],strlen(in_dir),strlen(flist[ii]))
     pgc_num[ii] = long(strmid(pgc_list[ii],3,strlen(pgc_list[ii])-3))
  endfor
  n_pgc = n_elements(pgc_list)
  
  openw, 1, 'survey_z0mgs.txt'

  for ii = 0, n_pgc-1 do $
     printf, 1, pgc_list[ii]

  close, 1

end
