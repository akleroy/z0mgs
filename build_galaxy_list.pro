pro build_galaxy_list $
   , in_dir = in_dir $
   , tag = tag $
   , just = just $
   , pgc_list = out_pgc_list $
   , pgc_num = out_pgc_num $
   , dat = out_gal_data $
   , start = start_num $
   , stop = stop_num $
   , exclude = exclude $
   , rebuild = rebuild

  if keyword_set(rebuild) then begin

     if n_elements(in_dir) eq 0 then begin    
        in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
     endif
     
;    BUILD A LIST OF PGC GALAXIES FOR WHICH WE HAVE UNWISE
     flist = file_search(in_dir+'PGC*', count=file_ct)
     pgc_list = strarr(file_ct)
     pgc_num = lonarr(file_ct)
     for ii = 0, file_ct-1 do begin
        pgc_list[ii] = strmid(flist[ii],strlen(in_dir),strlen(flist[ii]))
        pgc_num[ii] = long(strmid(pgc_list[ii],3,strlen(pgc_list[ii])-3))
     endfor
     n_pgc = n_elements(pgc_list)
     
;    GET THE GALBASE ENTRY FOR EACH GALAXY
     all_gal_data = gal_data(pgc_list)     

     save $
        , file='../measurements/z0mgs_list.idl' $
        , pgc_list, pgc_num, all_gal_data

  endif else begin

     restore, '../measurements/z0mgs_list.idl', /v

  endelse

; LOOP OVER AND PARE DOWN 
  ctr = 0B
  n_pgc = n_elements(pgc_list)
  keep = bytarr(n_pgc)
  for ii = 0L, n_pgc-1 do begin
     
     counter, ii, n_pgc, 'Building sample '
          
     pgc_name = pgc_list[ii]
     this_dat = all_gal_data[ii]

     if n_elements(exclude) gt 0 then begin
        if total(pgc_name eq exclude) gt 0 then $
           continue
     endif

     if n_elements(just) gt 0 then $
        if total(pgc_name eq just) eq 0 then $
           continue
     
     if n_elements(tag) gt 0 then begin
        this_tags = strsplit(strcompress(this_dat.tags,/rem), ';', /extract)
        match = 0B
        for jj = 0, n_elements(this_tags)-1 do $
           if this_tags[jj] eq strupcase(tag) then $
              match = 1B
        if match eq 0B then continue
     endif

     ctr += 1

     if n_elements(start_num) gt 0 then begin
        if ctr lt start_num then continue
     endif

     if n_elements(stop_num) gt 0 then begin
        if ctr gt stop_num then continue
     endif         

     keep[ii] = 1B
     
  endfor

  ind = where(keep, ct)
  if ct eq 0 then return
  out_pgc_list = pgc_list[ind]
  out_pgc_num = pgc_num[ind]
  out_gal_data = all_gal_data[ind]

end
