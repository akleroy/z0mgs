pro z0mgs_photometry $
   , image=image $
   , hdr=hdr $
   , xctr=xctr $
   , yctr=yctr $
   , pa=pa $
   , incl=incl $
   , reject=rejected

; Convert to Jy/pix
  pix_sr = (sxpar(hdr,'CD1_1')*!dtor)^2
  if pix_sr eq 0 then $
     pix_sr = (sxpar(hdr,'CDELT1')*!dtor)^2
  map *= pix_sr*1d6
        
  if n_elements(incl) eq 0 then $
     incl = 0.0

  if finite(incl) eq 0 then $
     incl = 0.0

  if incl gt 80. then $
     incl = 80.
        
  if n_elements(pa) eq 0 then begin
     pa = 0.0
     incl = 0.0
  endif

  if finite(pa) eq 0 then begin
     pa = 0.0
     incl = 0.0
  endif

; Work out galactocentric radius
  make_axes, hdr, ri=ri, di=di
  if n_elements(xctr) eq 0 then $
     xctr = mean(ri, /nan)
  if n_elements(yctr) eq 0 then $
     xctr = mean(di, /nan)

 pos_vec = $
    [pa, incl, xctr, yctr] 

  deproject, ri, di, pos_vec, rgrid=rgrid

; Bin the data
  binsize = 7.5/3600./r25
  bins = $
     bin_data(r25_grid, map $
              , xmin=0.0, xmax=max(rgrid, /nan) $
              , binsize=binsize)

  stop

end
