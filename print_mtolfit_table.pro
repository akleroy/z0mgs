pro print_mtolfit_table

  readcol $
     , 'mtol_fits.txt', format='A,F,F,F,F', comment='#' $
     , mtolfit_tag, mtolfit_lo, mtolfit_loval, mtolfit_hi, mtolfit_hival
  mtolfit_tag = strcompress(mtolfit_tag, /rem)
  n_tags = n_elements(mtolfit_tag)

  for ii = 0, n_tags-1 do begin
     slope = $
        (mtolfit_hival[ii] - mtolfit_loval[ii]) / $
        (mtolfit_hi[ii]-mtolfit_lo[ii])
     print, mtolfit_tag[ii]+' '+string(mtolfit_lo[ii])+' '+ $
            string(slope)+' '+string(mtolfit_hi[ii])
  endfor

end
