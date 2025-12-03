pro v2_index_sdss_tiles

  do_build_list = 0B
  
  sdss_dir='/data/fourier/leroy.42/allsky/sdss_dr12/frames/'
  tab_dir='../../working_data/sdss/index/'

  file_list = mrdfits(tab_dir+'sdss_frame_list.fits',1)  
  n_files = n_elements(file_list)

  nan = !values.f_nan
  empty = { $
          aklid:-1L, $
          fname:'', $
          frame_id:'', $
          filter:'', $
          exptime:nan, $
          ra:nan, $
          dec:nan, $
          naxis1:-1L, $
          naxis2:-1L, $
          crpix1:-1L, $
          crpix2:-1L $
  }
  tiles = replicate(empty, n_files)

  start_time = systime(/seconds)
  
  for ii = 0, n_files-1 do begin

     if ii mod 100 eq 0 then $
        print, 'Header ', ii, ' out of ', n_files

     hdr = 'dummy'
     hdr = headfits(file_list[ii].fname, compress='bunzip2')

     if n_elements(hdr) eq 1 then begin
        print, "Failed to read. Stopping for debug."
        stop
     endif
     
     tiles[ii].aklid = file_list[ii].id
     tiles[ii].fname = file_list[ii].fname
     tiles[ii].filter = sxpar(hdr,'FILTER')
     tiles[ii].exptime = sxpar(hdr,'EXPTIME')
     tiles[ii].naxis1 = sxpar(hdr,'NAXIS1')
     tiles[ii].naxis2 = sxpar(hdr,'NAXIS2')
     tiles[ii].crpix1 = sxpar(hdr,'CRPIX1')
     tiles[ii].crpix2 = sxpar(hdr,'CRPIX2')
     tiles[ii].ra = sxpar(hdr,'CRVAL1')
     tiles[ii].dec = sxpar(hdr,'CRVAL2')

     if ii mod 1000 eq 0 then begin
        delta_time = systime(/seconds)- start_time
        print, "Seconds elapsed: ", delta_time
        print, "Predicted duration (hours): ", (1.0*n_files)/(1.0*ii)*(delta_time)/3600.
     endif
     
  endfor

  mwrfits, tiles, tab_dir+'sdss_frame_index.fits', /create
  
  stop
  
end
