pro plot_gswlc_grids $
   , gridname=gridname

  @constants.bat
  lsun_3p4 = 1.83d18
  restore, '../gswlc/gswlc_data.idl'

  plot, findgen(10), xtitle='!6'

  valid_grids = ['ssfr_mstar', 'nuvw1_w3w1', 'fuvw1_w4w1']

  if n_elements(gridname) eq 0 then begin
     print, "Defaulting to SFR/Mstar vs Mstar."
     gridname = 'ssfr_mstar'
  endif

  if total(gridname eq valid_grids) eq 0 then begin
     print, "Invalid grid name."
     stop
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOAD THE GRIDS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  counts = readfits('../measurements/'+gridname+'_count_grid.fits', ct_hdr)
  mtol_grid = readfits('../measurements/'+gridname+'_mtol_grid.fits', mtol_hdr)
  cfuv_grid = readfits('../measurements/'+gridname+'_cfuv_grid.fits', cfuv_hdr)
  cfuvw4_grid = readfits('../measurements/'+gridname+'_cfuvw4_grid.fits', cfuvw4_hdr)
  cjustw4_grid = readfits('../measurements/'+gridname+'_cjustw4_grid.fits', cjustw4_hdr)
  cfuvw3_grid = readfits('../measurements/'+gridname+'_cfuvw3_grid.fits', cfuvw3_hdr)
  cjustw3_grid = readfits('../measurements/'+gridname+'_cjustw3_grid.fits', cjustw3_hdr)

  sz_grid = size(mtol_grid)

  if gridname eq 'ssfr_mstar' then begin
     ext = 'ssfrmstar'
     ytitle = "!6log!d10!n SFR / M!d*!n [yr!u-1!n]"
     xtitle = "!6log!d10!n Stellar Mass from CIGALE [M!d!9n!n!6]"
  endif

  if gridname eq 'nuvw1_w3w1' then begin
     ext = 'nuvw1w3'
     xtitle = "!6log!d10!n NUV/WISE1"
     ytitle = "!6log!d10!n WISE3/WISE1"
  endif

  if gridname eq 'fuvw1_w4w1' then begin
     ext = 'fuvw1w4'
     xtitle = "!6log!d10!n FUV/WISE1"
     ytitle = "!6log!d10!n WISE4/WISE1"
  endif

  x_axis = sxpar(mtol_hdr,'CRVAL1')+(findgen(sz_grid[1]))*sxpar(mtol_hdr,'CDELT1')
  delta_x = x_axis[1] - x_axis[0]

  y_axis = sxpar(mtol_hdr,'CRVAL2')+(findgen(sz_grid[2]))*sxpar(mtol_hdr,'CDELT2')
  delta_y = y_axis[1] - y_axis[0]

  if keyword_set(do_contour) then begin
     vec = counts[where(finite(counts))]
     vec = vec[sort(vec)]
     cdf = total(vec, /cumul)/total(vec,/nan)
     clevs = [interpol(vec, cdf, 0.05) $
              , interpol(vec, cdf, 0.16)]
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MASS TO LIGHT GRID
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  psfile = '../plots/mtol_grid_'+ext+'.eps'
  pnfile = '../plots/mtol_grid_'+ext+'.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  !p.multi=0

  viridis
  disp, alog10(mtol_grid), x_axis, y_axis $
        , xtitle=xtitle $
        , ytitle=ytitle $
        , charsize=1.25, charthick=3 $
        , xstyle=1, ystyle=1 $
        , position=[0.2, 0.2, 0.925, 0.75] $
        , color=cgcolor('black',255) $
        , missing = cgcolor('white',254), reserve=5 $
        , min=-1.25, max=-0.25

  for ii = 0, sz_grid[1]-2 do $
     oplot, (x_axis[ii]+delta_x*0.5)*[1,1], [-100,100], lines=1, color=cgcolor('black')
  for ii = 0, sz_grid[2]-2 do $
     oplot, [-100,100], (y_axis[ii]+delta_y*0.5)*[1,1], lines=1, color=cgcolor('black')

  for ii = 0, sz_grid[1]-1 do $
     for jj = 0, sz_grid[2]-1 do $
        if finite(counts[ii,jj]) then $
           xyouts, x_axis[ii], y_axis[jj] $
                   , string(alog10(counts[ii,jj]),format='(F3.1)') $
                   , align=0.5, charsize=0.75, color=cgcolor('white')
  
  if keyword_set(do_contour) then $
     contour, alog10(counts), x_axis, y_axis $
              , /overplot, color=cgcolor('red') $
              , lev=alog10(clevs), /close

  !p.charthick=3
  viridis
  cgColorbar, Range=[-1.25, -0.25] $
              , Position=[0.2, 0.925, 0.925, 0.975] $
              , xtickformat='(F5.2)' $
              , title='!6log!d10!n WISE1 !7T!6!d*!n [M!d!9n!6!n / L!d!9n!6!n]' $
              , annotatecolor=cgcolor("black", 255) $
              , charthick = 3, charsize=1.25 $
              , textthick=3
  !p.charthick=1
  
  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FUV+W4
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  psfile = '../plots/cfuv_grid_'+ext+'.eps'
  pnfile = '../plots/cfuv_grid_'+ext+'.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  !p.multi=0

  viridis
  disp, alog10(cfuv_grid), x_axis, y_axis $
        , xtitle=xtitle $
        , ytitle=ytitle $
        , charsize=1.25, charthick=3 $
        , xstyle=1, ystyle=1 $
        , position=[0.2, 0.2, 0.925, 0.75] $
        , color=cgcolor('black',255) $
        , missing = cgcolor('white',254), reserve=5 $
        , min=-44.0, max=-43.0

  for ii = 0, sz_grid[1]-2 do $
     oplot, (x_axis[ii]+delta_x*0.5)*[1,1], [-100,100], lines=1, color=cgcolor('black')
  for ii = 0, sz_grid[2]-2 do $
     oplot, [-100,100], (y_axis[ii]+delta_y*0.5)*[1,1], lines=1, color=cgcolor('black')

  for ii = 0, sz_grid[1]-1 do $
     for jj = 0, sz_grid[2]-1 do $
        if finite(counts[ii,jj]) then $
           xyouts, x_axis[ii], y_axis[jj] $
                   , string(alog10(counts[ii,jj]),format='(F3.1)') $
                   , align=0.5, charsize=0.75, color=cgcolor('white')
  
  !p.charthick=3
  viridis
  cgColorbar, Range=[-44.0, -43.0] $
              , Position=[0.2, 0.925, 0.925, 0.975] $
              , xtickformat='(F5.1)' $
              , title='!6log!d10!n C FUV term [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
              , annotatecolor=cgcolor("black", 255) $
              , charthick = 3, charsize=1.25 $
              , textthick=3
  !p.charthick=1
  
  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile


; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FUV+W4
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  psfile = '../plots/cfuvw4_grid_'+ext+'.eps'
  pnfile = '../plots/cfuvw4_grid_'+ext+'.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  !p.multi=0

  viridis
  disp, alog10(cfuvw4_grid), x_axis, y_axis $
        , xtitle=xtitle $
        , ytitle=ytitle $
        , charsize=1.25, charthick=3 $
        , xstyle=1, ystyle=1 $
        , position=[0.2, 0.2, 0.925, 0.75] $
        , color=cgcolor('black',255) $
        , missing = cgcolor('white',254), reserve=5 $
        , min=-43.25, max=-42.25

  for ii = 0, sz_grid[1]-2 do $
     oplot, (x_axis[ii]+delta_x*0.5)*[1,1], [-100,100], lines=1, color=cgcolor('black')
  for ii = 0, sz_grid[2]-2 do $
     oplot, [-100,100], (y_axis[ii]+delta_y*0.5)*[1,1], lines=1, color=cgcolor('black')

  for ii = 0, sz_grid[1]-1 do $
     for jj = 0, sz_grid[2]-1 do $
        if finite(counts[ii,jj]) then $
           xyouts, x_axis[ii], y_axis[jj] $
                   , string(alog10(counts[ii,jj]),format='(F3.1)') $
                   , align=0.5, charsize=0.75, color=cgcolor('white')
  
  !p.charthick=3
  viridis
  cgColorbar, Range=[-43.25, -42.25] $
              , Position=[0.2, 0.925, 0.925, 0.975] $
              , xtickformat='(F5.1)' $
              , title='!6log!d10!n C WISE4+FUV [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
              , annotatecolor=cgcolor("black", 255) $
              , charthick = 3, charsize=1.25 $
              , textthick=3
  !p.charthick=1
  
  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FUV+W3
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  psfile = '../plots/cfuvw3_grid_'+ext+'.eps'
  pnfile = '../plots/cfuvw3_grid_'+ext+'.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  !p.multi=0

  viridis
  disp, alog10(cfuvw3_grid), x_axis, y_axis $
        , xtitle=xtitle $
        , ytitle=ytitle $
        , charsize=1.25, charthick=3 $
        , xstyle=1, ystyle=1 $
        , position=[0.2, 0.2, 0.925, 0.75] $
        , color=cgcolor('black',255) $
        , missing = cgcolor('white',254), reserve=5 $
        , min=-43.25, max=-42.25

  for ii = 0, sz_grid[1]-2 do $
     oplot, (x_axis[ii]+delta_x*0.5)*[1,1], [-100,100], lines=1, color=cgcolor('black')
  for ii = 0, sz_grid[2]-2 do $
     oplot, [-100,100], (y_axis[ii]+delta_y*0.5)*[1,1], lines=1, color=cgcolor('black')

  for ii = 0, sz_grid[1]-1 do $
     for jj = 0, sz_grid[2]-1 do $
        if finite(counts[ii,jj]) then $
           xyouts, x_axis[ii], y_axis[jj] $
                   , string(alog10(counts[ii,jj]),format='(F3.1)') $
                   , align=0.5, charsize=0.75, color=cgcolor('white')
  
  !p.charthick=3
  viridis
  cgColorbar, Range=[-43.25, -42.25] $
              , Position=[0.2, 0.925, 0.925, 0.975] $
              , xtickformat='(F5.1)' $
              , title='!6log!d10!n C WISE3+FUV [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
              , annotatecolor=cgcolor("black", 255) $
              , charthick = 3, charsize=1.25 $
              , textthick=3
  !p.charthick=1
  
  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; JUST WISE 4
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  psfile = '../plots/cjustw4_grid_'+ext+'.eps'
  pnfile = '../plots/cjustw4_grid_'+ext+'.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  !p.multi=0

  viridis
  disp, alog10(cjustw4_grid), x_axis, y_axis $
        , xtitle=xtitle $
        , ytitle=ytitle $
        , charsize=1.25, charthick=3 $
        , xstyle=1, ystyle=1 $
        , position=[0.2, 0.2, 0.925, 0.75] $
        , color=cgcolor('black',255) $
        , missing = cgcolor('white',254), reserve=5 $
        , min=-43.0, max=-42.0

  for ii = 0, sz_grid[1]-2 do $
     oplot, (x_axis[ii]+delta_x*0.5)*[1,1], [-100,100], lines=1, color=cgcolor('black')
  for ii = 0, sz_grid[2]-2 do $
     oplot, [-100,100], (y_axis[ii]+delta_y*0.5)*[1,1], lines=1, color=cgcolor('black')

  for ii = 0, sz_grid[1]-1 do $
     for jj = 0, sz_grid[2]-1 do $
        if finite(counts[ii,jj]) then $
           xyouts, x_axis[ii], y_axis[jj] $
                   , string(alog10(counts[ii,jj]),format='(F3.1)') $
                   , align=0.5, charsize=0.75, color=cgcolor('white')
  
  !p.charthick=3
  viridis
  cgColorbar, Range=[-43.0, -42.0] $
              , Position=[0.2, 0.925, 0.925, 0.975] $
              , xtickformat='(F5.1)' $
              , title='!6log!d10!n C WISE4 only [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
              , annotatecolor=cgcolor("black", 255) $
              , charthick = 3, charsize=1.25 $
              , textthick=3
  !p.charthick=1
  
  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; JUST WISE 3
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  psfile = '../plots/cjustw3_grid_'+ext+'.eps'
  pnfile = '../plots/cjustw3_grid_'+ext+'.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  !p.multi=0

  viridis
  disp, alog10(cjustw3_grid), x_axis, y_axis $
        , xtitle=xtitle $
        , ytitle=ytitle $
        , charsize=1.25, charthick=3 $
        , xstyle=1, ystyle=1 $
        , position=[0.2, 0.2, 0.925, 0.75] $
        , color=cgcolor('black',255) $
        , missing = cgcolor('white',254), reserve=5 $
        , min=-43.25, max=-42.25

  for ii = 0, sz_grid[1]-2 do $
     oplot, (x_axis[ii]+delta_x*0.5)*[1,1], [-100,100], lines=1, color=cgcolor('black')
  for ii = 0, sz_grid[2]-2 do $
     oplot, [-100,100], (y_axis[ii]+delta_y*0.5)*[1,1], lines=1, color=cgcolor('black')

  for ii = 0, sz_grid[1]-1 do $
     for jj = 0, sz_grid[2]-1 do $
        if finite(counts[ii,jj]) then $
           xyouts, x_axis[ii], y_axis[jj] $
                   , string(alog10(counts[ii,jj]),format='(F3.1)') $
                   , align=0.5, charsize=0.75, color=cgcolor('white')
  
  !p.charthick=3
  viridis
  cgColorbar, Range=[-43.25, -42.25] $
              , Position=[0.2, 0.925, 0.925, 0.975] $
              , xtickformat='(F5.1)' $
              , title='!6log!d10!n C WISE3 only [M!d!9n!6!n yr!u-1!n (erg s!u-1!n)!u-1!n]' $
              , annotatecolor=cgcolor("black", 255) $
              , charthick = 3, charsize=1.25 $
              , textthick=3
  !p.charthick=1
  
  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile

  stop

end
