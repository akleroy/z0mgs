pro z0mgs_image_slice $
   , infile = infile $
   , outfile=outfile $
   , outimage=outimage $
   , gal_data = dat $
   , show = show

  angle_ra = findgen(40)*5.-10.
  n_slice = n_elements(angle_ra)
  width = 15.0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  map = readfits(infile, hdr)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET UP THE EXERCISE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  make_axes, hdr, rimg=rimg, dimg=dimg

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER SLICE ANGLE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  r25_val = fltarr(n_slice)
  
  for ii = 0, n_slice-1 do begin
     angle = angle_ra[ii]
     deproject $
        , rimg, dimg, [angle, 0.0, dat.ra_deg, dat.dec_deg] $
        , xgrid = major, ygrid=minor
     major *= 3600.
     minor *= 3600.
     ind = where(abs(minor) lt width)
     x = major[ind]
     y = map[ind]
     
     extent = (3.0*dat.r25_deg*3600)

     prof_1 = bin_data(x, y, /nan $
                       , xmin=-7.5/2., xmax=extent $
                       , binsize=7.5)
     prof_2 = bin_data(-1*x, y, /nan $
                       , xmin=-7.5/2., xmax=extent $
                       , binsize=7.5)
     nbins = n_elements(prof_1)
     if ii eq 0 then $
        prof_ra = fltarr(n_slice, nbins)
     min_prof = (prof_1.ymed < prof_2.ymed)
     prof_ra[ii,*] = min_prof
  endfor
  
  rank_ra = prof_ra*0.0
  for ii = 0, nbins-1 do begin
     s = sort(prof_ra[*,ii])
     ind = lindgen(n_slice)
     rank = rank_ra[*,ii]*0.0
     rank[s] = ind
     rank_ra[*,ii] = rank
  endfor

  crosscorr = rank_ra*0.0
  for ii = 0, n_slice-1 do begin
     fid_peak = angle_ra[ii]
     model = abs((angle_ra - fid_peak)) < (abs(angle_ra - (fid_peak+180)) < $
                                           abs(angle_ra - (fid_peak-180)))
     model_ra = crosscorr
     for jj = 0, nbins-1 do $
        model_ra[*,jj] = model

     crosscorr[ii,*] = total(model_ra*rank_ra, 1, /nan)
  endfor     

;  loadct, 0
;  circle, /fill
;  plot, angle_ra, total(rank_ra[*,0.15*nbins:0.35*nbins], 2, /nan) $
;        , psym=-8, ystyle=16

  stop

end
