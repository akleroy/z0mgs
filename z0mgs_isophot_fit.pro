pro z0mgs_isophot_fit $
   , infile = infile $
   , outfile=outfile $
   , outimage=outimage $
   , gal_data = dat $
   , show = show

  area_thresh = 50
  tol = 5.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  map = readfits(infile, hdr)

  rms = sxpar(hdr, 'MADALL', count=rmsct)
  if rmsct eq 0 then begin
     print, 'No noise measurement in header.'
     rms = mad(map)
  endif

  xyad, hdr, 0, 0, ra_1, dec_1
  xyad, hdr, 1, 0, ra_2, dec_2
  pix_scale = sphdist(ra_1, dec_1, ra_2, dec_2, /deg)*3600.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET UP THE FIT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sz = size(map)
  xra = findgen(sz[1]) # (fltarr(sz[2])+1.)
  yra = (fltarr(sz[1])+1.) # findgen(sz[2])

  adxy, hdr, dat.ra_deg, dat.dec_deg, xctr, yctr

  parinfo = { $
            limited:[0,0] $
            , limits:[0.,0.] $
            }
  parinfo_ra = replicate(parinfo, 5)
  parinfo_ra[0].limited = [1,1]
  parinfo_ra[1].limited = [1,1]
  parinfo_ra[2].limited = [1,1]
  parinfo_ra[3].limited = [1,1]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER LEVEL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  min_lev = 3.0*rms
  max_lev = map[xctr, yctr]
  if max_lev lt min_lev then begin
     print, "Not enough bright area found."
     return
  endif

  nlevs = (alog10(max_lev) - alog10(min_lev))/(alog10(sqrt(2.)))
  levs = min_lev*2.0^(0.5*findgen(nlevs))
  
  pa_ra = fltarr(nlevs)*!values.f_nan
  major_ra = fltarr(nlevs)*!values.f_nan
  minor_ra = fltarr(nlevs)*!values.f_nan
  xctr_ra = fltarr(nlevs)*!values.f_nan
  yctr_ra = fltarr(nlevs)*!values.f_nan
  chisq_ra = fltarr(nlevs)*!values.f_nan
  flux_ra = fltarr(nlevs)*!values.f_nan
  area_ra = fltarr(nlevs)*!values.f_nan

  if keyword_set(show) then begin
     loadct, 0
     disp, alog10(map), max=max(alog10(levs), /nan), min=min(alog10(levs), /nan)
  endif

  first = 1B

  if n_elements(outfile) gt 0 then begin
     openw, 1, outfile
     printf, 1, '# Column 1: Isophote in MJy/sr'
     printf, 1, '# Column 2: Fit RA Center in decimal degrees'
     printf, 1, '# Column 3: Fit Dec Center in decimal degrees'
     printf, 1, '# Column 4: Major Axis in arcseconds'
     printf, 1, '# Column 5: Minor Axis in arcseconds'
     printf, 1, '# Column 6: Position angle in radians'
     printf, 1, '# Column 7: Enclosed area in arcseconds^2'
     printf, 1, '# Column 8: Enclosed flux in Jy' 
  endif

  for ii = 0, nlevs-1 do begin

     thresh = levs[ii]
     mask = map gt thresh     
     reg = label_region(mask)
     if reg[xctr, yctr] eq 0 then continue

     mask = mask*(reg eq reg[xctr, yctr])
     area = total(mask)
     if area lt area_thresh then $
        continue

     edge = $
        grow_mask(mask eq 0, iters=2) eq mask     
     edge[0,*] = 0B
     edge[*,0] = 0B
     edge[sz[1]-1,*] = 0B
     edge[*,sz[2]-1] = 0B

     ind = where(edge)
     x = xra[ind]
     y = yra[ind]
     dist = sqrt((x-xctr)^2+(y-yctr)^2)

     mindist = min(dist)
     maxdist = max(dist)
     delta = 15.
     parinfo_ra[2].limits = [xctr-delta,xctr+delta]
     parinfo_ra[3].limits = [yctr-delta,yctr+delta]
     parinfo_ra[0].limits = [mindist, maxdist]
     parinfo_ra[1].limits = [mindist, maxdist]
     guess = [maxdist, mindist, xctr, yctr, dat.posang_deg*!dtor]
     p = mpfitellipse(x, y, guess $
                      , /tilt $
                      , parinfo=parinfo_ra $
                     , /quiet)
     phi = dindgen(101)*2D*!dpi/100
     xm = p[2] + p[0]*cos(phi)*cos(p[4]) + p[1]*sin(phi)*sin(p[4])
     ym = p[3] - p[0]*cos(phi)*sin(p[4]) + p[1]*sin(phi)*cos(p[4])
    
     xyad, hdr, p[2], p[3], ra_fit, dec_fit
     xctr_ra[ii] = ra_fit
     yctr_ra[ii] = dec_fit
     
     offset = sphdist(ra_fit, dec_fit, dat.ra_deg, dat.dec_deg, /deg)*3600.
     if offset gt tol then begin
        print, "Ellipse not centered. Skipping."
        continue
     endif

     major_ra[ii] = p[0]*pix_scale
     minor_ra[ii] = p[1]*pix_scale
     pa_ra[ii] = p[4]
     area_ra[ii] = area*pix_scale^2
     flux_ra[ii] = total(mask*map,/nan)* $
                   (pix_scale/3600.*!dtor)^2*1d6
     
     oplot, xm, ym, color=cgcolor('red')

     line = string(levs[ii], format='(F8.5)')+" "+$
            string(xctr_ra[ii], format='(F9.5)')+" "+$
            string(yctr_ra[ii], format='(F9.5)')+" "+$
            string(major_ra[ii], format='(F8.1)')+" "+$
            string(minor_ra[ii], format='(F8.1)')+" "+$
            string(pa_ra[ii], format='(F6.3)')+" "+$
            string(area_ra[ii], format='(F8.1)')+" "+$
            string(flux_ra[ii], format='(F9.5)')
     if n_elements(outfie) gt 0 then $
        printf,1, line
     print, line

  endfor
  
  if n_elements(outfile) gt 0 then $
     close, 1

  if n_elements(outimage) gt 0 then begin
     im = tvrd(true=1)
     write_png, outimage, im
  endif

end
