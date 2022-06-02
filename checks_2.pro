pro checks_2

  readcol, '../measurements/z0mgs_fullsfrmstar.txt' $
           , format='L,A,A,A,F,F,F,F,F,A,F,F,A,F,A' $
           , delim=',', comment='#' $
           , pgc, ngc, ugc, ic $
           , dist, edistdex, logmstar, elogmstar, mtol, methodmtol $
           , logsfr, elogsfr, methodsfr, deltams, flags $
           , /nan

  tab = mrdfits('../measurements/delivery_index_gauss15.fits',1,h)

  ind = where(finite(logmstar) eq 0, ct)
  for ii = 0, ct-1 do begin
     this_pgc = pgc[ind[ii]]
     this_ind = where(tab.pgc eq this_pgc, this_ct)
     if this_ct eq 0 then stop
     print, this_pgc $
            , mtol[ind[ii]] $
            , methodmtol[ind[ii]] $
            , tab[this_ind].has_wise1 $
            , tab[this_ind].flux_wise1 $
            , tab[this_ind].has_fuv $
            , tab[this_ind].flux_fuv $
            , tab[this_ind].has_nuv $
            , tab[this_ind].flux_nuv $
            , tab[this_ind].has_wise3 $
            , tab[this_ind].flux_wise3 $
            , tab[this_ind].has_wise4 $
            , tab[this_ind].flux_wise4 $
            , tab[this_ind].mtol $
            , tab[this_ind].logsfr
            
     
  endfor

  print, "Through list."
  stop

end
