pro print_mag_to_intens

  readcol, 'mag_to_intens.txt' $
           , format='A,A,A,F' $
           , mag, band, res, s
  s = alog10(s)
  for ii = 0, n_elements(s)-1 do begin
     print, strupcase(mag[ii]), ' ' $
            , strupcase(band[ii]), ' ' $
            , strupcase(res[ii]), ' ' $
            , s[ii]
  endfor

end
