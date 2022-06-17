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

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Build our list and database of galaxy data.
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(rebuild) then begin

     if n_elements(in_dir) eq 0 then begin    
        in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
     endif

;         
;    1. Build a list of galaxies for which we already have unwise data.
;

     flist = file_search(in_dir+'PGC*', count=file_ct)
     pgc_list = strarr(file_ct)
     pgc_num = lonarr(file_ct)
     for ii = 0, file_ct-1 do begin
        pgc_list[ii] = strmid(flist[ii],strlen(in_dir),strlen(flist[ii]))
        pgc_num[ii] = long(strmid(pgc_list[ii],3,strlen(pgc_list[ii])-3))
     endfor

;
;    2. Supplement this with a list of additional targets.
;

     readcol $
        , '../measurements/selection_v2.txt', format='L' $
        , new_pgc_num

;    3. Combine the lists into a list of unique.

     pgc_num = [pgc_num, new_pgc_num]
     pgc_num = pgc_num[sort(pgc_num)]
     pgc_num = pgc_num[uniq(pgc_num)]
     
;    4. Look up the information for each target in the galaxy database.

     all_gal_data = gal_data(pgc=pgc_num, /full)     

;    ... a handful of galaxies have dropped out of our database

     zero_ind = where(all_gal_data.pgc eq 0, zero_ct)
     if zero_ct gt 0 then begin
        nonzero_ind = where(all_gal_data.pgc ne 0, nonzero_ct)
        all_gal_data = all_gal_data[nonzero_ind]
        pgc_num = pgc_num[nonzero_ind]
     endif

;    ... build the string list

     pgc_list = 'PGC'+str(all_gal_data.pgc)

;    5. Save the results as an IDL file

     save $
        , file='../measurements/z0mgs_list.idl' $
        , pgc_list, pgc_num, all_gal_data

  endif else begin

     restore, '../measurements/z0mgs_list.idl', /v

  endelse

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Allow some downselect
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  n_pgc = n_elements(pgc_list)

; initialize a counter and a list of flags
  ctr = 0L
  keep = bytarr(n_pgc)

; loop over targets
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

; Return only the selection

  ind = where(keep, ct)
  if ct eq 0 then return
  out_pgc_list = pgc_list[ind]
  out_pgc_num = pgc_num[ind]
  out_gal_data = all_gal_data[ind]

end
