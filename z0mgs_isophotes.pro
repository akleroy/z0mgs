pro z0mgs_isophotes $
   , infile = infile $
   , outfile=outfile $   
   , dat=dat $
   , band=band $
   , show = show $
   , pause = pause

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  test = file_search(infile, count=file_ct)
  if file_ct eq 0 then begin
     message, "File not found.", /info
     return
  endif
  map = readfits(infile, hdr, /silent)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD ASTROMETRY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sr_per_pix = get_pixel_scale(hdr)*!dtor
  adxy, hdr, dat.ra_deg, dat.dec_deg, xctr, yctr
  sz = size(map)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DETERMINE LEVELS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  rms = sxpar(hdr, 'MADALL', count=rmsct)
  if rmsct eq 0 or rms lt 0. or finite(rms) eq 0 then begin
     print, 'No noise measurement in header.'
     rms = mad(map)
  endif

  min_lev = 3.0*rms
  max_lev = map[xctr, yctr]
  if max_lev lt min_lev then begin
     print, "Not enough bright area found."
     return
  endif

  nlevs = ceil((alog10(max_lev) - alog10(min_lev))/(alog10(sqrt(2.))))
  levs = min_lev*2.0^(0.5*findgen(nlevs))

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER LEVELS AND CALCULATE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  flux_ra = fltarr(nlevs)*!values.f_nan
  area_ra = fltarr(nlevs)*!values.f_nan
  perim_ra = fltarr(nlevs)*!values.f_nan

  if keyword_set(show) then begin
     !p.multi=[0, 2, 1]
     loadct, 0
     disp, alog10(map), xstyle=1, ystyle=1
  endif

  for ii = 0, nlevs-1 do begin

     thresh = levs[ii]
     mask = map gt thresh     
     reg = label_region(mask, /ulong)
     if reg[xctr, yctr] eq 0 then continue

     mask = mask*(reg eq reg[xctr, yctr])
     area = total(mask)
     flux = total(mask*map, /nan)

     edge = $
        grow_mask(mask eq 0, iters=2) eq mask     
     edge[0,*] = 0B
     edge[*,0] = 0B
     edge[sz[1]-1,*] = 0B
     edge[*,sz[2]-1] = 0B
     perim = total(edge)

     area_ra[ii] = area
     perim_ra[ii] = perim
     flux_ra[ii] = flux

     if keyword_set(show) then begin
        contour, mask, lev=[1], /overplot, color=cgcolor('red')     
     endif

  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  get_lun, lun
  openw, lun, outfile

  printf, lun, '# Target: ' + strcompress(dat.pgcname, /rem)
  printf, lun, '# Band: ' + strcompress(band, /rem)
  printf, lun, '# File: ' + strcompress(infile, /rem)
  printf, lun, '# Lowest contour: ' + string(min(levs), format='(F8.3)')
  printf, lun, '# Steradian_per_pix: ' + string(sr_per_pix, format='(F10.8)')
  printf, lun, '# Galaxy_r25_arcsec: '+ string(dat.r25_deg*3600., format='(F6.1)')
  printf, lun, '# Column_1: Isophote (MJy/sr)'
  printf, lun, '# Column_2: Area (pix^2)'
  printf, lun, '# Column_3: Perimeter (pix)'
  printf, lun, '# Column_4: Flux (MJy/sr*pix)'

  for ii = 0, nlevs-1 do begin
     line = ''
     line += string(levs[ii], format='(F12.6)')
     line += ' '+string(area_ra[ii], format='(I8)')
     line += ' '+string(perim_ra[ii], format='(I8)')
     line += ' '+string(flux_ra[ii], format='(F12.6)')
     printf, lun, line
  endfor

  close, lun
  free_lun, lun

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF DESIRED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(show) then begin
     plot, levs, area_ra, /xlo, /ylo
     if keyword_set(pause) then begin
        print, "Key to continue."
        ch = get_kbrd(1)
     endif
  endif

end
