pro build_galaxy_list $
   , in_dir = in_dir $
   , tag = tag $
   , just = just $
   , pgc_list = out_pgc_list $
   , pgc_num = out_pgc_num $
   , dat = out_gal_data $
   , start = start_num $
   , stop = stop_num $
   , exclude = exclude

  if n_elements(in_dir) eq 0 then begin    
     in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  endif

; BUILD A LIST OF PGC GALAXIES FOR WHICH WE HAVE UNWISE
  flist = file_search(in_dir+'PGC*', count=file_ct)
  pgc_list = strarr(file_ct)
  pgc_num = lonarr(file_ct)
  for ii = 0, file_ct-1 do begin
     pgc_list[ii] = strmid(flist[ii],strlen(in_dir),strlen(flist[ii]))
     pgc_num[ii] = long(strmid(pgc_list[ii],3,strlen(pgc_list[ii])-3))
  endfor
  n_pgc = n_elements(pgc_list)

; TAG
  if n_elements(tag) gt 0 then begin                
     all_data = gal_data(tag=tag)
  endif else begin
     all_data = gal_data(/all)
  endelse

; LOOP OVER AND PARE DOWN 
  first = 1B
  ctr = 0B
  for ii = 0L, n_pgc-1 do begin
     
     counter, ii, n_pgc, 'Building sample '
          
     pgc_name = pgc_list[ii]
     if n_elements(exclude) gt 0 then begin
        if total(pgc_name eq exclude) gt 0 then $
           continue
     endif

     if n_elements(just) gt 0 then $
        if total(pgc_name eq just) eq 0 then $
           continue
     
     if total(all_data.pgc eq pgc_num[ii]) eq 0 then $
        continue

     ctr += 1

     if n_elements(start_num) gt 0 then begin
        if ctr lt start_num then continue
     endif

     if n_elements(stop_num) gt 0 then begin
        if ctr gt stop_num then continue
     endif         

     if first then begin
        first = 0B
        out_pgc_list = [pgc_list[ii]]
        out_pgc_num = [pgc_num[ii]]
     endif else begin
        out_pgc_list = [out_pgc_list,pgc_list[ii]]
        out_pgc_num = [out_pgc_num,pgc_num[ii]]
     endelse
     
  endfor

  out_gal_data = gal_data(out_pgc_list)

end
