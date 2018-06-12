pro read_z0mgs_profile
  
  readcol, outfile, format='X,A,A', numline = 10, key, val
  prof_target = strcompress(val[0], /rem)
  prof_band = strcompress(val[1], /rem)
  prof_file = strcompress(val[2], /rem)
  prof_pa = float(val[3])
  prof_incl = float(val[3])
  prof_ra = float(val[4])
  prof_dec = float(val[5])
  prof_srpix = float(val[6])
  prof_rbin = float(val[7])
  prof_r25as = float(val[8])

  readcol, outfile, comment = '#' $
           , prof_rmid, prof_flux, prof_mean, prof_med $
           , prof_counts, prof_mad, prof_std, prof_16th, prof_84th

end
