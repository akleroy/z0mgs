pro z0mgs_profile $
   , infile=infile $
   , outfile=outfile $
   , dat=dat $
   , binsize=binsize $
   , band=band $
   , show=show $
   , pause=pause

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  test = file_search(infile, count=file_ct)
  if file_ct eq 0 then begin
     message, "File not found.", /info
     return
  endif
  map = readfits(infile, hdr, /silent)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD ASTROMETRY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sr_per_pix = get_pixel_scale(hdr)*!dtor
  make_axes, hdr, ri=ri, di=di
  pa = dat.posang_deg
  incl = dat.incl_deg
  if finite(pa) eq 0 or finite(incl) eq 0. then begin
     pa = 0.0
     incl = 0.0
  endif

  if incl gt 80. then begin
     incl = 80.
  endif

  pos_vec = $
     [pa, incl, dat.ra_deg, dat.dec_deg]
  deproject, ri, di, pos_vec, rgrid=rgrid, tgrid=tgrid
  rgrid *= 3600.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BIN THE IMAGE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  rbin = 7.5
  rmin = -0.5*rbin
  delta = (max(di, /nan) - min(di, /nan))*3600./2.0
  rmax = delta

  bins = bin_data(rgrid, map, /nan $
                  , xmin=rmin, xmax=rmax, binsize=rbin $
                  , perc=0.16)
  n_bins = n_elements(bins)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  get_lun, lun
  openw, lun, outfile

  printf, lun, '# Target: ' + strcompress(dat.pgcname, /rem)
  printf, lun, '# Band: ' + strcompress(band, /rem)
  printf, lun, '# File: ' + strcompress(infile, /rem)
  printf, lun, '# Position_angle_deg: ' + string(pa, format='(F8.3)')
  printf, lun, '# Inclination_deg: ' + string(incl, format='(F8.3)')
  printf, lun, '# R.A._center_deg: ' + string(dat.ra_deg, format='(F10.5)')
  printf, lun, '# Dec._center_deg: ' + string(dat.dec_deg, format='(F10.5)')
  printf, lun, '# Steradian_per_pix: ' + string(sr_per_pix, format='(F14.10)')
  printf, lun, '# Binsize_arcsec: ' + string(rbin, format='(F5.2)')
  printf, lun, '# Galaxy_r25_arcsec: '+ string(dat.r25_deg*3600., format='(F6.1)')
  printf, lun, '# Column_1: Radius at middle of bin'
  printf, lun, '# Column_2: Sum in bin (MJy/sr*pix)'
  printf, lun, '# Column_3: Mean in bin (MJy/sr)'
  printf, lun, '# Column_4: Median in bin (MJy/sr)'
  printf, lun, '# Column_5: Counts in bin'
  printf, lun, '# Column_6: Robust scatter in bin (MJy/sr)'
  printf, lun, '# Column_7: Scatter in bin (MJy/sr)'
  printf, lun, '# Column_8: 16th percentile value (MJy/sr)'
  printf, lun, '# Column_9: 84th percentile value (MJy/sr)'

  for ii = 0, n_bins-1 do begin
     line = ''
     line += string(bins[ii].xmid, format='(F8.2)')
     line += ' '+string(bins[ii].ysum, format='(F12.6)')
     line += ' '+string(bins[ii].ymean, format='(F12.6)')
     line += ' '+string(bins[ii].ymed, format='(F12.6)')
     line += ' '+string(bins[ii].counts, format='(I8)')
     line += ' '+string(bins[ii].ymad, format='(F12.6)')
     line += ' '+string(bins[ii].ystd, format='(F12.6)')
     line += ' '+string(bins[ii].ylo_perc, format='(F12.6)')
     line += ' '+string(bins[ii].yhi_perc, format='(F12.6)')
     printf, lun, line
  endfor

  close, lun
  free_lun, lun

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF DESIRED (TESTING READING)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(show) then begin

     readcol, outfile, format='X,A,A', numline = 10, key, val, /silent
     prof_target = strcompress(val[0], /rem)
     prof_band = strcompress(val[1], /rem)
     prof_file = strcompress(val[2], /rem)
     prof_pa = float(val[3])
     prof_incl = float(val[4])
     prof_ra = float(val[5])
     prof_dec = float(val[6])
     prof_srpix = float(val[7])
     prof_rbin = float(val[8])
     prof_r25as = float(val[9])

     readcol, outfile, comment = '#' $
              , prof_rmid, prof_flux, prof_mean, prof_med $
              , prof_counts, prof_mad, prof_std, prof_16th, prof_84th $
              , /silent

     plot, prof_rmid, prof_mean $
           , /nodata
     oplot, prof_rmid, prof_mean, color=cgcolor('cyan'), ps=1
     oplot, prof_rmid, prof_med, color=cgcolor('magenta'), ps=6

  endif

end
