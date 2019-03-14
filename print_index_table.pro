pro print_index_table

  atlas_dir = '../delivery/'
  gaia_dir = '../stars/gaia/'
  
  tab = mrdfits('../measurements/delivery_index_gauss7p5.fits', 1, h)
  n_pgc = n_elements(tab)

  openw, 1, '../measurements/z0mgs_radec.txt'
  
  for ii = 0, n_pgc-1 do begin

     pgc_name = strcompress(tab[ii].pgc_name, /rem)
     line = pgc_name+','+string(tab[ii].ra_deg,format='(F9.5)')+ $
            ','+string(tab[ii].dec_deg,format='(F9.5)')
     printf,1,line

  endfor

  close,1
  
end
