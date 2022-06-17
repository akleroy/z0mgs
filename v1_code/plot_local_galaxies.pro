pro plot_local_galaxies

  @constants.bat

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOAD THE POPULATION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  unwise_dir = '../unwise/atlas/'
  galex_dir = '../galex/atlas/'
  prof_dir = '../measurements/profiles/'

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  out_dir = '../unwise/atlas/'

  build_galaxy_list $
     , in_dir = in_dir $
     , tag=tag $
     , just=just $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = gal_data $
     , start = start_num $
     , stop = stop_num $
     , exclude = ['PGC17223']
  
  n_pgc = n_elements(pgc_list)

  mstar = gal_data.lum_w1*12./lsun

  ind1 = where(mstar ge 1d10 and mstar lt 2d10)
  ind2 = where(mstar ge 2d10 and mstar lt 4d10)
  ind3 = where(mstar gt 4d10)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DISTANCE AND SIZE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
 
  psfile = '../plots/local_size_vs_dist.eps'
  ps, /def, /ps, xs=8, ys=8, /color, /encaps $
      , file=psfile

  xlo = 0.
  xhi = 100.
  ylo = 10.
  yhi = 1d3
  loadct, 0
  plot, [0], [0], /nodata $
        , yrange=[ylo, yhi], /ylog $
        , xrange=[xlo, xhi] $
        , ytitle='!6Radius !8r!d25!n!6 ["]' $
        , xtitle='!6Distance [Mpc]' $
        , xthick=5, ythick=5 $
        , charthick=5, charsize=1.5
  
  yhi = 15.*3./2.
  polyfill, [xlo, xhi, xhi, xlo, xlo] $
            , [ylo, ylo, yhi, yhi, ylo] $
            , color=cgcolor('lightsalmon')

  ylo = 15.*3./2.
  yhi = 15.*10./2.
  polyfill, [xlo, xhi, xhi, xlo, xlo] $
            , [ylo, ylo, yhi, yhi, ylo] $
            , color=cgcolor('lightyellow')

  ylo = 15.*10./2.
  yhi = 1d3
  polyfill, [xlo, xhi, xhi, xlo, xlo] $
            , [ylo, ylo, yhi, yhi, ylo] $
            , color=cgcolor('lightseagreen')

  oplot $
     , gal_data[ind1].dist_mpc, gal_data[ind1].r25_deg*3600. $     
     , ps=cgsymcat('filledcircle'), color=cgcolor('lightgray') $
     , symsize=1.0

  oplot $
     , gal_data[ind2].dist_mpc, gal_data[ind2].r25_deg*3600. $     
     , ps=cgsymcat('filledcircle'), color=cgcolor('gray') $
     , symsize=1.0

  oplot $
     , gal_data[ind3].dist_mpc, gal_data[ind3].r25_deg*3600. $     
     , ps=cgsymcat('filledcircle'), color=cgcolor('darkgray') $
     , symsize=1.0

  fid = findgen(1000)/10.
  oplot, fid, 1./(fid*1d3)/!dtor*3600., lines=2, thick=3
  oplot, fid, 2./(fid*1d3)/!dtor*3600., lines=2, thick=3
  oplot, fid, 4./(fid*1d3)/!dtor*3600., lines=2, thick=3
  oplot, fid, 8./(fid*1d3)/!dtor*3600., lines=2, thick=3
  oplot, fid, 16./(fid*1d3)/!dtor*3600., lines=2, thick=3
  oplot, fid, 32./(fid*1d3)/!dtor*3600., lines=2, thick=3
  oplot, fid, 64./(fid*1d3)/!dtor*3600., lines=2, thick=3

  oplot, fid, fid*0.0+15., lines=0, thick=10
  xyouts, [80], [16.], '!6Resolution (FWHM)', align=0.5, charthick=3, charsize=1.25

  al_legend $
     , /top, /right $
     , box=1, clear=1 $
     , background=cgcolor('white') $
     , charsize=1.5, charthick=3 $
     , ['!6M!d*!n > 4E10 M!d!9n!6!n' $
        ,'!6M!d*!n > 2E10 M!d!9n!6!n' $
        ,'!6M!d*!n >  1E10 M!d!9n!6!n' $
       , '1,2,4,...64 kpc'] $
     , psym=cgsymcat('filledcircle')*[1,1,1,0], lines=[0,0,0,2] $
     , color=[cgcolor('darkgray'), cgcolor('gray'), cgcolor('lightgray'), cgcolor('black')] $
     , pspacing=1.0

  xlo = 0.
  xhi = 100.
  ylo = 10.
  yhi = 1d3
  loadct, 0
  plot, [0], [0], /nodata $
        , yrange=[ylo, yhi], /ylog $
        , xrange=[xlo, xhi] $
        , ytitle='!6Radius !8r!d25!n!6 ["]' $
        , xtitle='!6Distance [Mpc]' $
        , xthick=5, ythick=5 $
        , charthick=5, charsize=1.5 $
        , /noerase

  ps, /xw
  spawn, 'evince '+psfile+' &'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; THE RA AND DEC VIEW
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  psfile = '../plots/local_equatorial.eps'
  ps, /def, /ps, xs=8, ys=8, /color, /encaps $
      , file=psfile

  plot, [0], [0], /nodata $
        , xrange=[0., 360.], xstyle=1 $
        , yrange=[0., 40.] $
        , xtitle='!6Right Ascenscion' $
        , ytitle='!6Distance [Mpc]' $
        , xthick=5, ythick=5 $
        , charthick=5, charsize=1.5

  for ii = 0, 360, 15. do $
     oplot, ii*[1,1], [0., 100], lines=1

  for ii = 15., 345, 30. do $
     xyouts, ii, 2.0, strcompress(string(round(ii / 15.), format='(2I)')+'h',/rem), align=0.5


  x = gal_data.ra_deg
  y = gal_data.dist_mpc
  oplot, x[ind1], y[ind1], psym=cgsymcat('filledcircle'), symsize=1.5, color=cgcolor('lightgray')
  oplot, x[ind2], y[ind2], psym=cgsymcat('filledcircle'), symsize=1.5, color=cgcolor('gray')
  oplot, x[ind3], y[ind3], psym=cgsymcat('filledcircle'), symsize=1.5, color=cgcolor('darkgray')

  ps, /xw
  spawn, 'evince '+psfile+' &'

  psfile = '../plots/local_horizontal.eps'
  ps, /def, /ps, xs=8, ys=8, /color, /encaps $
      , file=psfile

  plot, [0], [0], /nodata $
        , xrange=[-90., 90.], xstyle=1 $
        , yrange=[0., 40.] $
        , xtitle='!6Declination' $
        , ytitle='!6Distance [Mpc]' $
        , xthick=5, ythick=5 $
        , charthick=5, charsize=1.5

  x = gal_data.dec_deg
  y = gal_data.dist_mpc
  oplot, x[ind1], y[ind1], psym=cgsymcat('filledcircle'), symsize=1.5, color=cgcolor('lightgray')
  oplot, x[ind2], y[ind2], psym=cgsymcat('filledcircle'), symsize=1.5, color=cgcolor('gray')
  oplot, x[ind3], y[ind3], psym=cgsymcat('filledcircle'), symsize=1.5, color=cgcolor('darkgray')

  ps, /xw
  spawn, 'evince '+psfile+' &'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DISTRIBUTION VS DIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  psfile = '../plots/local_count_vs_dist.eps'
  ps, /def, /ps, xs=8, ys=8 $
      , /color, /encaps $
      , file=psfile

  plot, [0], [0], /nodata $
        , xrange=[0., 40.], xstyle=1 $
        , yrange=[1., 1d4], /ylo $
        , ytitle='!6Cumulative Count' $
        , xtitle='!6Distance [Mpc]' $
        , xthick=5, ythick=5 $
        , charthick=5, charsize=1.5

  fid = findgen(1000)/10.
  vol = 4./3.*!pi*fid^3.
  oplot, fid, vol/(4./3.*!pi*10.^3.)*10., line=2, thick=10

  vec = gal_data[ind1].dist_mpc
  vec = vec[sort(vec)]
  cdf = total(vec*0.+1., /cumul)
  oplot, vec, cdf, color=cgcolor('lightgray'), thick=15

  vec = gal_data[ind2].dist_mpc
  vec = vec[sort(vec)]
  cdf = total(vec*0.+1., /cumul)
  oplot, vec, cdf, color=cgcolor('gray'), thick=15

  vec = gal_data[ind3].dist_mpc
  vec = vec[sort(vec)]
  cdf = total(vec*0.+1., /cumul)
  oplot, vec, cdf, color=cgcolor('darkgray'), thick=15

  vec = gal_data[where(mstar lt 5d9 and mstar gt 1d9)].dist_mpc
  vec = vec[sort(vec)]
  cdf = total(vec*0.+1., /cumul)
  oplot, vec, cdf, color=cgcolor('firebrick'), thick=15, lines=0
  fid = findgen(1000)/10.
  vol = 4./3.*!pi*fid^3.
  oplot, fid, vol/(4./3.*!pi*10.^3.)*50., line=2, thick=10, color=cgcolor('salmon')

  ps, /xw
  spawn, 'evince '+psfile+' &'

  stop

end
