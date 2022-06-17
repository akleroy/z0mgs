pro plot_distribution

  dir = '/data/kant/0/leroy.42/allsky/z0mgs/'
  tab = mrdfits(dir+'../measurements/delivery_index_gauss15.fits',1,h)
  s = gal_data(pgc=tab.pgc)
  
  x = s.dist_mpc*cos(!dtor*s.dec_deg)*cos(!dtor*s.ra_deg)
  y = s.dist_mpc*cos(!dtor*s.dec_deg)*sin(!dtor*s.ra_deg)
  z = s.dist_mpc*sin(!dtor*s.dec_deg)

  bin_mpc = 0.1
  grid_xy = $
     grid_data(x, y $
               , xmin=-50. $
               , xmax=50. $
               , ymin=-50. $
               , ymax=50. $
               , binsize_x=bin_mpc $
               , binsize_y=bin_mpc $
               , xaxis_out = xaxis_xy $
               , yaxis_out = yaxis_xy $
               , /nan)
  grid_xy /= (bin_mpc)^2.

  grid_xz = $
     grid_data(x, z $
               , xmin=-50. $
               , xmax=50. $
               , ymin=-50. $
               , ymax=50. $
               , binsize_x=bin_mpc $
               , binsize_y=bin_mpc $
               , xaxis_out = xaxis_xz $
               , yaxis_out = yaxis_xz $
               , /nan)
  grid_xz /= (bin_mpc)^2.

  kern = psf_gaussian(npix=41, fwhm=21., /norm)
  smooth_xy = convol(grid_xy*1.0, kern)
  smooth_xz = convol(grid_xz*1.0, kern)

; PLOT

  psfile = '../plots/distrib_dist_xy.eps'
  pnfile = '../plots/distrib_dist_xy.png'
  ps, /ps, /def, file=psfile, xs=5, ys=6 $
      , /color, /encaps

  loadct, 3
  reversect
  disp, smooth_xy $
        , xrange=[-50, 50] $
        , yrange=[-50, 50] $
        , thick=5, xthick=5, ythick=5 $
        , charsize=1.0, charthick=3 $
        , xtitle='!8x!6 [Mpc]', ytitle='!8y!6 [Mpc]' $
        , xstyle=1, ystyle=1 $
        , color=cgcolor('black'), reserve=5 $
        , min=0, max=10. $
        , /sq, position=[0.15, 0.15, 0.95, 0.85]

  contour $
     , smooth_xy, xaxis_xy, yaxis_xy $
     , lev=[0.25, 1.0, 4.], /overplot $
     , color=cgcolor('charcoal')
  contour $
     , smooth_xy, xaxis_xy, yaxis_xy $
     , lev=[16.], /overplot $
     , color=cgcolor('lightgray')

  hour_angle = findgen(24)
  xfid = findgen(101)
  for ii = 0, n_elements(hour_angle)-1 do $
     oplot, xfid*cos(!dtor*hour_angle[ii]*15.) $
            , xfid*sin(!dtor*hour_angle[ii]*15.) $
            , color=cgcolor('black') $
            , thick=1, lines=1
  
  hour_label = findgen(4)
  for ii = 0, n_elements(hour_label)-1 do $
     xyouts, [40.]*cos(!dtor*hour_label[ii]*6.*15.) $
             , [40.]*sin(!dtor*hour_label[ii]*6.*15.) $
             , strcompress(string(hour_label[ii]*6. $
                                  ,format='(3I)'),/rem)+'h' $
             , color=cgcolor('lightgray') $
             , align=0.5, charthick=7, charsize=1.0

  hour_label = findgen(4)
  for ii = 0, n_elements(hour_label)-1 do $
     xyouts, [40.]*cos(!dtor*hour_label[ii]*6.*15.) $
             , [40.]*sin(!dtor*hour_label[ii]*6.*15.) $
             , strcompress(string(hour_label[ii]*6. $
                                  ,format='(3I)'),/rem)+'h' $
             , color=cgcolor('black') $
             , align=0.5, charthick=1, charsize=1.0

  !p.charthick=3
  loadct, 3
  reversect
  cgColorbar, Range=[0., 10.] $
              , Position=[0.15, 0.95, 0.95, 0.975] $
              , xtickformat='(F5.2)' $
              , title="!6Density [Galaxies Mpc!u-2!n]" $
              , annotatecolor=cgcolor("black", 255) $
              , charthick = 3, charsize=1.0 $
              , textthick=3
  !p.charthick=1
  loadct, 0

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

  psfile = '../plots/distrib_dist_xz.eps'
  pnfile = '../plots/distrib_dist_xz.png'
  ps, /ps, /def, file=psfile, xs=5, ys=6., /color, /encaps

  loadct, 3
  reversect
  disp, smooth_xz $
        , xrange=[-50, 50] $
        , yrange=[-50, 50] $
        , thick=5, xthick=5, ythick=5 $
        , charsize=1.0, charthick=3 $
        , xtitle='!8x!6 [Mpc]', ytitle='!8z!6 [Mpc]' $
        , xstyle=1, ystyle=1 $
        , color=cgcolor('black'), reserve=5 $
        , min=0, max=10. $
        , /sq, position=[0.15, 0.15, 0.95, 0.85]

  contour $
     , smooth_xz, xaxis_xz, yaxis_xz $
     , lev=[0.25, 1.0, 4.], /overplot $
     , color=cgcolor('charcoal')
  contour $
     , smooth_xz, xaxis_xz, yaxis_xz $
     , lev=[16.], /overplot $
     , color=cgcolor('lightgray')

  dec = findgen(24.)*15.
  dec_for_string = (dec ge 0 and dec le 90)*dec + $
                   (dec gt 90 and dec le 180)*(180-dec) + $
                   (dec gt 180 and dec le 270)*(180-dec) + $
                   (dec gt 270 and dec le 360)*(dec-360)
  dec_string = strarr(n_elements(dec))
  for ii = 0, n_elements(dec)-1 do $
     dec_string[ii] = strcompress(string(dec_for_string[ii] $
                                         ,format='(3I)'),/rem)+ $
     '!ud!n'
  xfid = findgen(101)
  for ii = 0, n_elements(dec)-1 do $
     oplot, xfid*cos(!dtor*dec[ii]) $
            , xfid*sin(!dtor*dec[ii]) $
            , color=cgcolor('black') $
            , thick=1, lines=1

  for ii = 0, n_elements(dec)-1, 6 do $
     xyouts, [40.]*cos(!dtor*dec[ii]) $
             , [40.]*sin(!dtor*dec[ii]) $
             , dec_string[ii] $
             , color=cgcolor('lightgray') $
             , align=0.5, charthick=7, charsize=1.0

  for ii = 0, n_elements(dec)-1, 6 do $
     xyouts, [40.]*cos(!dtor*dec[ii]) $
             , [40.]*sin(!dtor*dec[ii]) $
             , dec_string[ii] $
             , color=cgcolor('black') $
             , align=0.5, charthick=1, charsize=1.0

  !p.charthick=3
  loadct, 3
  reversect
  cgColorbar, Range=[0., 10.] $
              , Position=[0.15, 0.95, 0.95, 0.975] $
              , xtickformat='(F5.2)' $
              , title="!6Density [Galaxies Mpc!u-2!n]" $
              , annotatecolor=cgcolor("black", 255) $
              , charthick = 3, charsize=1.0 $
              , textthick=3
  !p.charthick=1
  loadct, 0

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

  stop

end
