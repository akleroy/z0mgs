if pgc_name eq 'PGC2557' then begin

  rgrid = readfits(radfile, rhdr)
  tempfile = infile+'.temp.fits'  

  blank_map_from_list $
     , infile=infile $
     , outfile=tempfile $
     , pgc=this_dat.pgc $
     , band=band $
     , blank_pgc = blank_pgc $
     , blank_band = blank_band $
     , blank_shape = blank_shape $
     , blank_x1 = blank_x1 $
     , blank_y1 = blank_y1 $
     , blank_x2 = blank_x2 $
     , blank_y2 = blank_y2

  map = readfits(tempfile, hdr)
  ind = where(rgrid gt 1.5*sxpar(rhdr, 'FIDRAD'))
  map[ind] = !values.f_nan
  writefits, tempfile, map, hdr

  bkfit_galex $
     , mapfile=tempfile $
     , wtfile=wtfile $
     , outfile=outfile $
     , rejfile=rejfile $
     , radfile=radfile $
     , masklist=masklist $
     , band=str(band) $
     , rejected=rejected $
     , show=show $
     , pause=pause $
     , aperture=0.8 $
     , plane=0

  bkgrd = readfits(outfile, bkgrd_hdr)
  ind = where(rgrid gt 1.5*sxpar(rhdr, 'FIDRAD'))
  map[ind] = !values.f_nan
  writefits, outfile, bkgrd, bkgrd_hdr
  
  ;mask = abs(map gt 3d-4)
  ;ind = where(mask eq 0)
  ;medval = median(map[ind])
  ;map -= medval
  ;mask = abs(map gt 3d-4)
  ;ind = where(mask eq 0)
  ;medval = median(map[ind])
  ;map -= medval


endif

if pgc_name eq 'PGC5818' then begin
  
  rgrid = readfits(radfile, rhdr)

  tempfile = infile+'.temp.fits'  

  map = readfits(infile, hdr)
  ind = where(rgrid gt 2.0*sxpar(rhdr, 'FIDRAD'))
  map[ind] = !values.f_nan
  writefits, tempfile, map, hdr
  
  bkfit_galex $
     , mapfile=tempfile $
      , wtfile=wtfile $
      , outfile=outfile $
      , rejfile=rejfile $
      , radfile=radfile $
      , masklist=masklist $
      , band=str(band) $
      , rejected=rejected $
      , show=show $
      , pause=pause $
      , aperture=1.0 $
      , /plane                 

  bkgrd = readfits(outfile, bkgrd_hdr)
  ind = where(rgrid gt 2.0*sxpar(rhdr, 'FIDRAD'))
  map[ind] = !values.f_nan
  writefits, outfile, bkgrd, bkgrd_hdr

endif
