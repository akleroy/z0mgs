pro blank_map_from_list $
   , infile=infile $
   , outfile=outfile $
   , pgc=pgc $
   , band=band $
   , blank_pgc = blank_pgc $
   , blank_band = blank_band $
   , blank_shape = blank_shape $
   , blank_x1 = blank_x1 $
   , blank_y1 = blank_y1 $
   , blank_x2 = blank_x2 $
   , blank_y2 = blank_y2
  
  if n_elements(blank_pgc) eq 0 then begin
     readcol, 'blank_regions.txt', comment='#' $
              , format='L,A,A,F,F,F,F' $
              , blank_pgc, blank_band, blank_shape $
              , blank_x1, blank_y1, blank_x2, blank_y2
  endif

; READ THE DATA
  if n_elements(infile) gt 0 then begin
     if file_test(infile) eq 0 then begin
        print, "File not found."
        return
     endif
     map = readfits(infile, hdr)
  endif

; LOGIC TO MATCH NAME AND BAND
  ind = where(blank_pgc eq pgc, ct)
  if ct eq 0 then begin
     print, "no matches. Output will equal input."
  endif else begin
     make_axes, hdr, ri=ri, di=di
  endelse

; LOGIC TO BLANK  
  for ii = 0, ct-1 do begin

; ... IN A BOX
     if strcompress(blank_shape[ind[ii]], /rem) eq  'box' then begin
        to_blank = where((ri ge blank_x1[ind[ii]]) and $
                         (ri le blank_x2[ind[ii]]) and $
                         (di ge blank_y1[ind[ii]]) and $
                         (di le blank_y2[ind[ii]]), blank_ct)
        if blank_ct gt 0 then $
           map[to_blank] = !values.f_nan
     endif

; ... IN A CIRCLE
     if strcompress(blank_shape[ind[ii]], /rem) eq  'circle' then begin
        dist = sphdist(ri, di, blank_x1[ind[ii]], blank_y1[ind[ii]], /deg)
        to_blank = where(dist le blank_x2[ind[ii]], blank_ct)
        if blank_ct gt 0 then $
           map[to_blank] = !values.f_nan
     endif

  endfor

; WRITE TO DISK
  if n_elements(outfile) gt 0 then begin
     writefits, outfile, map, hdr
  endif

end
