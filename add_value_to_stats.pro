pro add_value_to_stats

  restore, '../measurements/unwise_stats.idl', /v
  dat = gal_data(stats[*,0,0,0].pgcname)
  glactc, dat.ra_deg, dat.dec_deg, 2000., l, b, 1, /degree

  save $
     , file='../measurements/unwise_stats_with_dat.idl' $
     , stats, dat, l, b 

end
