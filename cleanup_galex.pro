pro cleanup_galex

; We have some "missing" files. Try to refetch these if that's the
; request.

; The missing rrhr files don't unzip.
  scratch_dir = '../galex/scratch_tiles/'
  miss_file_list = file_search(scratch_dir+'*.gz', count=ct)
  
  miss_list = strarr(ct)
  for ii = 0, ct-1 do begin
     miss_list[ii] = strmid(miss_file_list[ii],strlen(scratch_dir))
     print, miss_list[ii]   
  endfor

  ;fetch_all_galex_tiles $
  ;   , only=miss_list $
  ;   , start=0B $
  ;   , stop=200000L $
  ;   , /rrhr

  for ii = 0, ct-1 do begin
     line = 'rm -rf ../galex/scratch_tiles/'+miss_list[ii]
     spawn, line
     line = 'cp ../galex/all_tiles/'+miss_list[ii]+' ../galex/scratch_tiles/.'
     spawn, line
  endfor

  stop

end
