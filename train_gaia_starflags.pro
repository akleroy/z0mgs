pro train_gaia_starflags

  restore, '../measurements/gaia_stars_gauss7p5.idl', /v

  w1_bins = bin_data(gmag, w1, /nan, xmin=10, xmax=19., binsize=0.5)
  w2_bins = bin_data(gmag, w2, /nan, xmin=10, xmax=19., binsize=0.5)
  w3_bins = bin_data(gmag, w3, /nan, xmin=10, xmax=19., binsize=0.5)
  nuv_bins = bin_data(gmag, nuv, /nan, xmin=10, xmax=19., binsize=0.5)
  fuv_bins = bin_data(gmag, fuv, /nan, xmin=10, xmax=19., binsize=0.5)

  ind = where(w1_bins.xmid ge 10.5 and w1_bins.xmid le 16)
  w1_flux = 10.^(-1.0*w1_bins.xmid/2.5)
  w2_flux = 10.^(-1.0*w2_bins.xmid/2.5)
  w3_flux = 10.^(-1.0*w3_bins.xmid/2.5)
  nuv_flux = 10.^(-1.0*nuv_bins.xmid/2.5)
  fuv_flux = 10.^(-1.0*fuv_bins.xmid/2.5)

  w1_rat = median((w1_bins.ymed / w1_flux)[ind])
  w2_rat = median((w2_bins.ymed / w2_flux)[ind])
  w3_rat = median((w3_bins.ymed / w3_flux)[ind])
  nuv_rat = median((nuv_bins.ymed / nuv_flux)[ind])
  fuv_rat = median((fuv_bins.ymed / fuv_flux)[ind])

  print, 'WISE1 ratio : ', w1_rat
  print, 'WISE2 ratio : ', w2_rat
  print, 'WISE3 ratio : ', w3_rat
  print, 'NUV ratio : ', nuv_rat
  print, 'FUV ratio : ', fuv_rat

  psfile = '../plots/gaia_stacks.eps'
  pnfile = '../plots/gaia_stacks.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  plot, [0], [0], yrange=[-6,2.], xrange=[10., 20.] $
        , xthick=5, ythick=5 $
        , xtitle='!6GAIGA G Magnitude' $
        , ytitle='log!d10!n Atlas Intensity [MJy/sr]' $
        , charthick=3, charsize=1.25

  fid = findgen(101)

  loadct, 0
  for ii = -100, 100 do $
     oplot, ii*1.0*[1,1], [-100,100], lines=1
  for ii = -100, 100 do $
     oplot, [-100,100], ii*1.0*[1,1], lines=1

  oplot, w1_bins.xmid, alog10(w1_bins.ymed) $
         , ps=cgsymcat('filledcircle'), symsize=1.5 $
         , color=cgcolor('firebrick')
  oplot, fid, alog10(10.^(-1.*fid/2.5)*w1_rat), lines=0 $
         , color=cgcolor('firebrick')
  oplot, [0, 1d6], alog10(3e-3*5.)*[1,1], lines=0 $
         , color=cgcolor('firebrick'), thick=10

  oplot, w2_bins.xmid, alog10(w2_bins.ymed) $
        , ps=cgsymcat('filledcircle'), symsize=1.5 $
         , color=cgcolor('salmon')
  oplot, fid, alog10(10.^(-1.*fid/2.5)*w2_rat), lines=0 $
         , color=cgcolor('salmon')
  oplot, [0, 1d6], alog10(3e-3*5.)*[1,1], lines=2 $
         , color=cgcolor('salmon'), thick=10

  oplot, w3_bins.xmid, alog10(w3_bins.ymed) $
         , ps=cgsymcat('filledcircle'), symsize=1.5 $
         , color=cgcolor('goldenrod')
  oplot, fid, alog10(10.^(-1.*fid/2.5)*w3_rat), lines=0 $
         , color=cgcolor('goldenrod')
  oplot, [0, 1d6], alog10(0.1*5.)*[1,1], lines=2 $
         , color=cgcolor('goldenrod'), thick=10

  oplot, nuv_bins.xmid, alog10(nuv_bins.ymed) $
         , ps=cgsymcat('filledcircle'), symsize=1.5 $
         , color=cgcolor('dodgerblue')
  oplot, fid, alog10(10.^(-1.*fid/2.5)*nuv_rat), lines=0 $
         , color=cgcolor('dodgerblue')
  oplot, [0, 1d6], alog10(3e-4*5.)*[1,1], lines=0 $
          , color=cgcolor('dodgerblue'), thick=10

  oplot, fuv_bins.xmid, alog10(fuv_bins.ymed) $
         , ps=cgsymcat('filledcircle'), symsize=1.5 $
         , color=cgcolor('purple')
  oplot, fid, alog10(10.^(-1.*fid/2.5)*fuv_rat), lines=0 $
         , color=cgcolor('purple')
  oplot, [0, 1d6], alog10(3e-4*5.)*[1,1], lines=2 $
         , color=cgcolor('purple'), thick=10

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile+' &'
  
     
  stop

end
