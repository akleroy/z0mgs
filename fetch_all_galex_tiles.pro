pro fetch_all_galex_tiles $
   , start=start $
   , stop=stop $
   , out_dir=out_dir $
   , in_file=in_file $
   , rrhr=rrhr $
   , intbgsub=intbgsub $
   , only=only

; Fetch the all GALEX tiles based off of CAS job output.

  if n_elements(out_dir) eq 0 then $
     out_dir = '/data/tycho/0/leroy.42/allsky/galex/all_tiles/'
  
  if n_elements(in_file) eq 0 then $
     in_file = "MyAllSkyTable_akleroy.csv"

  if n_elements(start) eq 0 then $
     start = 0L
  
  if n_elements(stop) eq 0 then $
     stop = 0L  

; Read the input file

  print, "Reading the input file."

  nlines = file_lines(in_file)
  lines = strarr(nlines)
  openr, 1, in_file
  readf, 1, lines
  close, 1

  first_char = strmid(lines,0,1)
  keep = where(first_char ne '#', keep_ct)
  lines = lines[keep]
  nlines = n_elements(lines)

  if keyword_set(rrhr) then begin
     print, "I am replacing all -int.fits strings with -rrhr.fits"
     for ii = 0, nlines-1 do begin
        this_line = lines[ii]
        if strpos(this_line,'-int.fits') eq -1 then $
           continue
        new_line = mg_streplace(this_line,'-int.fits','-rrhr.fits')
        if n_elements(new_lines) eq 0 then $
           new_lines = new_line $
        else $
           new_lines = [new_lines, new_line]
     endfor
     lines = new_lines
     nlines = n_elements(lines)
  endif

  if keyword_set(intbgsub) then begin
     print, "I am replacing all -int.fits strings with -intbgsub.fits"
     for ii = 0, nlines-1 do begin
        this_line = lines[ii]
        if strpos(this_line,'-int.fits') eq -1 then $
           continue
        new_line = mg_streplace(this_line,'-int.fits','-intbgsub.fits')
        if n_elements(new_lines) eq 0 then $
           new_lines = new_line $
        else $
           new_lines = [new_lines, new_line]
     endfor
     lines = new_lines
     nlines = n_elements(lines)
  endif

  print, "Input file parsed, I found "+str(nlines)+" commands."
  print, "I will fetch from "+str(start)+" to "+str(stop)

; Change diretory and set up the command

  cd, out_dir

  for ii = long(start), long(stop) do begin

     if ii gt nlines then continue

     if n_elements(only) gt 0 then begin
        counter, ii, long(stop), " Testing line "

        use_line = 0B
        for jj = 0, n_elements(only)-1 do begin
           if strpos(lines[ii], only[jj]) ne -1 then begin
              use_line = 1B
              spawn, 'rm -rf '+out_dir+only[jj]
           endif
        endfor
     endif else begin
        use_line = 1B
     endelse

     if use_line then $
        spawn, lines[ii]

  endfor

end
