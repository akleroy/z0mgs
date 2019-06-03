pro print_gswlc_grid $
   , gridname=gridname

  @constants.bat
  lsun_3p4 = 1.83d18
  restore, '../gswlc/gswlc_data.idl'

  valid_grids = ['ssfr_mstar', 'nuvw1_w3w1', 'fuvw1_w4w1']

  if n_elements(gridname) eq 0 then begin
     print, "Defaulting to SFR/Mstar vs Mstar."
     gridname = 'ssfr_mstar'
  endif

  if total(gridname eq valid_grids) eq 0 then begin
     print, "Invalid grid name."
     stop
  endif

  plot, findgen(10), xtitle='!6'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOAD THE GRIDS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  counts = readfits('../measurements/'+gridname+'_count_grid.fits', ct_hdr)
  mtol_grid = readfits('../measurements/'+gridname+'_mtol_grid.fits', mtol_hdr)
  cfuvw4_grid = readfits('../measurements/'+gridname+'_cfuvw4_grid.fits', cfuvw4_hdr)
  cnuvw4_grid = readfits('../measurements/'+gridname+'_cnuvw4_grid.fits', cnuvw4_hdr)
  cjustw4_grid = readfits('../measurements/'+gridname+'_cjustw4_grid.fits', cjustw4_hdr)
  cfuvw3_grid = readfits('../measurements/'+gridname+'_cfuvw3_grid.fits', cfuvw3_hdr)
  cnuvw3_grid = readfits('../measurements/'+gridname+'_cfuvw3_grid.fits', cnuvw3_hdr)
  cjustw3_grid = readfits('../measurements/'+gridname+'_cjustw3_grid.fits', cjustw3_hdr)

  sz = size(mtol_grid)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE COLUMNS AND GRID
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if gridname eq 'ssfr_mstar' then begin
     ext = 'ssfrmstar'
     xtitle = "log10 SFR / M_* [yr^-1]"
     ytitle = "log10 Stellar Mass from CIGALE [M_sun]"
  endif

  if gridname eq 'nuvw1_w3w1' then begin
     ext = 'nuvw1w3'
     xtitle = "log10 NUV/WISE1"
     ytitle = "log10 WISE3/WISE1"
  endif

  if gridname eq 'fuvw1_w4w1' then begin
     ext = 'fuvw1w4'
     xtitle = "log10 FUV/WISE1"
     ytitle = "log10 WISE4/WISE1"
  endif

  x_axis = sxpar(mtol_hdr,'CRVAL1')+(findgen(sz[1]))*sxpar(mtol_hdr,'CDELT1')
  delta_x = x_axis[1] - x_axis[0]

  y_axis = sxpar(mtol_hdr,'CRVAL2')+(findgen(sz[2]))*sxpar(mtol_hdr,'CDELT2')
  delta_y = y_axis[1] - y_axis[0]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PRINT THE TEXT FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  outfile = '../measurements/fullgrid_'+gridname+'.txt'
  openw, 1, outfile
  printf, 1, '# Column 1: x axis = '+xtitle
  printf, 1, '# Column 2: y axis = '+ytitle
  printf, 1, '# Column 3: log10 galaxies in cell'
  printf, 1, '# Column 4: median WISE1 mass to light ratio [M_sun/L_sun]'
  printf, 1, '# Column 5: log10 C for WISE4 in FUV+WISE4 [M_sun/yr / (erg/s)]'
  printf, 1, '# Column 6: log10 C for WISE4 in NUV+WISE4 [M_sun/yr / (erg/s)]'
  printf, 1, '# Column 7: log10 C for WISE4 in WISE4 only [M_sun/yr / (erg/s)]'
  printf, 1, '# Column 8: log10 C for WISE3 in FUV+WISE3 [M_sun/yr / (erg/s)]'
  printf, 1, '# Column 9: log10 C for WISE3 in NUV+WISE3 [M_sun/yr / (erg/s)]'
  printf, 1, '# Column 10: log10 C for WISE3 in WISE3 only [M_sun/yr / (erg/s)]'

  for ii = 0, sz[1]-1 do begin
     for jj = 0, sz[2]-1 do begin

        line = ''
        line += string(x_axis[ii],format='(F6.2)')+' '
        line += string(y_axis[jj],format='(F6.2)')+' '
        line += string(alog10(counts[ii,jj]),format='(F5.1)')+' '
        line += string(mtol_grid[ii,jj],format='(F6.2)')+' '
        line += string(alog10(cfuvw4_grid[ii,jj]),format='(F6.2)')+' '
        line += string(alog10(cnuvw4_grid[ii,jj]),format='(F6.2)')+' '
        line += string(alog10(cjustw4_grid[ii,jj]),format='(F6.2)')+' '
        line += string(alog10(cfuvw3_grid[ii,jj]),format='(F6.2)')+' '
        line += string(alog10(cnuvw3_grid[ii,jj]),format='(F6.2)')+' '
        line += string(alog10(cjustw3_grid[ii,jj]),format='(F6.2)')

        printf, 1, line

     endfor
  endfor

  close, 1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PRINT THE LATEX STUB
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  outfile = '../measurements/stub_'+gridname+'.tex'
  openw, 1, outfile

  nlines = 0
  print_total = 10

  for ii = 0, sz[1]-1 do begin
     for jj = 0, sz[2]-1 do begin

        if nlines ge print_total then continue
        nlines += 1

        line = ''
        line += '$'+string(x_axis[ii],format='(F6.2)')+'$ & $'
        line += string(y_axis[jj],format='(F6.2)')+'$ & $'
        if counts[ii,jj] eq 0 or finite(alog10(counts[ii,jj])) eq 0 then begin
           for kk = 0, 6 do line += '\ldots$ & $'
           line += '\ldots$ \\'
        endif else begin
           line += string(alog10(counts[ii,jj]),format='(F5.1)')+'$ & $'
           line += string(mtol_grid[ii,jj],format='(F6.2)')+'$ & $'
           line += string(alog10(cfuvw4_grid[ii,jj]),format='(F6.2)')+'$ & $'
           line += string(alog10(cnuvw4_grid[ii,jj]),format='(F6.2)')+'$ & $'
           line += string(alog10(cjustw4_grid[ii,jj]),format='(F6.2)')+'$ & $'
           line += string(alog10(cfuvw3_grid[ii,jj]),format='(F6.2)')+'$ & $'
           line += string(alog10(cnuvw3_grid[ii,jj]),format='(F6.2)')+'$ & $'
           line += string(alog10(cjustw3_grid[ii,jj]),format='(F6.2)')+'$ \\'
        endelse

        printf, 1, line

     endfor
  endfor

  close, 1

end
