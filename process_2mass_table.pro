pro process_2mass_table
  
  tab = read_ipac_table('../measurements/fp_2mass.fp_psc22192.tbl')

  star_ra = tab.ra
  star_dec = tab.dec
  star_name = tab.designation
  star_km = tab.k_m
   
  save $
     , file='../measurements/2mass_stars.idl' $
     , star_ra, star_dec, star_name, star_km

end
