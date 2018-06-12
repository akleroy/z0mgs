pro print_oversamp

  ctr_ra = 30.
  ctr_dec = 30.
  rad_limit = 0.5
  hex_grid $
     , ctr_x = ctr_ra $
     , ctr_y = ctr_dec $
     , spacing = 7.5/3600. $
     , /radec $
     , xout = samp_ra $
     , yout = samp_dec $
     , r_limit = rad_limit $
     , /center

  beam_area = (15./2./3600.)^2*!pi/alog(2)
  n_beams = (rad_limit)^2*!pi/beam_area
  print, "Beams: ", n_beams
  print, "Samples: ", n_elements(samp_ra)
  print, "Oversamp: ", n_elements(samp_ra)/n_beams

end
