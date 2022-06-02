pro write_z0mgs_galbase_tag

  index = mrdfits('../measurements/delivery_index_gauss15.fits',1,hi)
  openw, 1, '~/idl/galbase/gal_data/survey_z0mgs.txt'
  for ii = 0L, n_elements(index)-1 do begin
     printf,1,strcompress(index[ii].pgc_name,/rem)
  endfor
  close, 1

  stop
  
; This still breaks because the PGC name is not an alias for about 200
; of these galaxies...

  index = mrdfits('../measurements/delivery_index_gauss15.fits',1,hi)
  dat1 = gal_data(tag='Z0MGS', /full)
  dat2 = gal_data(pgc=index.pgc, /full)

end
