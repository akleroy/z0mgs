readcol, '../measurements/z0mgs_fullsfrmstar.txt' $
         , format='L,A,A,A', delim=',', comment='#' $
         , pgc, ngc, ugc, ic

ngcind = where(strcompress(ngc, /rem) ne '0', ngcct)
pgcngc = (gal_data(ngc[ngcind])).pgc
print, total(pgcngc ne pgc[ngcind])
badngc = where(pgcngc ne pgc[ngcind], badct)
if badct gt 0 then $
   print, ngc[ngcind[badngc]] $
else $
   print, "No mismatches"

ugcind = where(strcompress(ugc, /rem) ne '0', ugcct)
pgcugc = (gal_data(ugc[ugcind])).pgc
print, total(pgcugc ne pgc[ugcind])
badugc = where(pgcugc ne pgc[ugcind], badct)
if badct gt 0 then $
   print, ugc[ugcind[badugc]] $
else $
   print, "No mismatches"

icind = where(strcompress(ic, /rem) ne '0', icct)
pgcic = (gal_data(ic[icind])).pgc
print, total(pgcic ne pgc[icind])
badic = where(pgcic ne pgc[icind], badct)
if badct gt 0 then $
   print, ic[icind[badic]] $
else $
   print, "No mismatches"
