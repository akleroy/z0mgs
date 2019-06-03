pro print_flag_statistics

  tab = mrdfits('../measurements/delivery_index_gauss7p5.fits', 1, h)

  print, "Overlapping galaxy:"
  print, "... ", total(tab.galaxy_overlap_flag*1.0)
  print, "... ", total(tab.galaxy_overlap_flag*1.0)/total(finite(tab.galaxy_overlap_flag))

  print, "Stars - WISE1:"
  print, "... ", total(tab.star_flag_wise1*1.0)
  print, "... ", total(tab.star_flag_wise1*1.0)/total(tab.has_wise1*1.0)

  print, "Stars - WISE2:"
  print, "... ", total(tab.star_flag_wise2*1.0)
  print, "... ", total(tab.star_flag_wise2*1.0)/total(tab.has_wise2*1.0)

  print, "Stars - WISE3:"
  print, "... ", total(tab.star_flag_wise3*1.0)
  print, "... ", total(tab.star_flag_wise3*1.0)/total(tab.has_wise3*1.0)

  print, "Stars - NUV:"
  print, "... ", total(tab.star_flag_nuv*1.0)
  print, "... ", total(tab.star_flag_nuv*1.0)/total(tab.has_nuv*1.0)

  print, "Stars - FUV:"
  print, "... ", total(tab.star_flag_fuv*1.0)
  print, "... ", total(tab.star_flag_fuv*1.0)/total(tab.has_fuv*1.0)

  print, "Saturation - WISE1:"
  print, "... ", total(tab.sat_effects_wise1*1.0)
  print, "... ", total(tab.sat_effects_wise1*1.0)/total(tab.has_wise1*1.0)

  print, "Saturation - WISE2:"
  print, "... ", total(tab.sat_effects_wise2*1.0)
  print, "... ", total(tab.sat_effects_wise2*1.0)/total(tab.has_wise2*1.0)

  print, "Saturation - WISE3:"
  print, "... ", total(tab.sat_effects_wise3*1.0)
  print, "... ", total(tab.sat_effects_wise3*1.0)/total(tab.has_wise3*1.0)

  tab = mrdfits('../measurements/delivery_index_gauss15.fits', 1, h)

  print, "Saturation - WISE4:"
  print, "... ", total(tab.sat_effects_wise4*1.0)
  print, "... ", total(tab.sat_effects_wise4*1.0)/total(tab.has_wise4*1.0)
  
  print, "Stars - WISE4:"
  print, "... ", total(tab.star_flag_wise4*1.0)
  print, "... ", total(tab.star_flag_wise4*1.0)/total(tab.has_wise4*1.0)

end
