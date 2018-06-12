pro plot_one_galex_galaxy

  tab = mrdfits('../measurements/delivery_index.fits', 1, h)

  ind = where(tab.has_fuv and tab.has_nuv and $
              tab.time_fuv lt 200. and tab.time_nuv lt 200. and $
              tab.gb_deg, ct)
  pick = ind[round(randomu(seed)*ct)]
    
  dir = '../delivery/'
  nuv = readfits(dir+strcompress(tab[pick].pgc_name+'_nuv.fits', /rem), hdr)
  nuv_rej = readfits(dir+strcompress(tab[pick].pgc_name+'_nuv_rejected.fits', /rem), hdr)
  fuv = readfits(dir+strcompress(tab[pick].pgc_name+'_fuv.fits', /rem), hdr)
  fuv_rej = readfits(dir+strcompress(tab[pick].pgc_name+'_fuv_rejected.fits', /rem), hdr)

  fuv_ind = where(fuv_rej lt 0.1)
  nuv_ind = where(nuv_rej lt 0.1)

  nuv_bins = $
     bin_data(nuv[nuv_ind]*1d3, nuv[nuv_ind]*0.0+1.0, /nan $
              , xmin=-5, xmax=5, binsize=0.1)
  
  fuv_bins = $
     bin_data(fuv[fuv_ind]*1d3, fuv[fuv_ind]*0.0+1.0, /nan $
              , xmin=-5, xmax=5, binsize=0.1)

  psfile = '../plots/galex_noise_example.eps'
  
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile

  plot $
     , [0], [0], /nodata $
     , xtitle='!6Intensity [10!u-3!n MJy/sr]' $
     , ytitle='!6log!d10!n Pixels' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[-3.,3], yrange=[0., 1.5*max(nuv_bins.counts)] $
     , /xstyle
  
  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-1d6, 1d6], lines=1, color=cgcolor('charcoal')
  
  for ii = -100, 100 do $
     oplot, [-1d6, 1d6], ii*2d3*[1,1], lines=1, color=cgcolor('charcoal')
  
  
  histplot $
     , nuv_bins.xmid, (nuv_bins.counts > 0.) $
     , /overplot $
     , lthick=3, /nobar, /fline, forient=45 $
     , fcolor=cgcolor('salmon') $
     , lcolor=cgcolor('firebrick')

  histplot $
     , fuv_bins.xmid, (fuv_bins.counts > 0.) $
     , /overplot $
     , lthick=3, /nobar, /fline, forient=-45. $
     , fcolor=cgcolor('dodgerblue') $
     , lcolor=cgcolor('navy')

  fuv_time_string = strcompress(string(round(tab[pick].time_fuv), format='(4I)'),/rem)+'s'
  nuv_time_string = strcompress(string(round(tab[pick].time_nuv), format='(4I)'),/rem)+'s'

  al_legend, /top, /right, clear=1, box=1, background=cgcolor('lightgray') $
             , ['!6FUV '+fuv_time_string,'!6NUV '+nuv_time_string] $
             , lines=-99, charsize=1.25, charthick=3 $
             , textcolor=[cgcolor('dodgerblue'), cgcolor('firebrick')]

  al_legend, /top, /left, clear=1, box=1, background=cgcolor('lightgray') $
             , ['!6'+tab[pick].pgc_name], lines=-99, charsize=1.25, charthick=3

  ps, /xw
  spawn, 'evince '+psfile+' &'
  
  stop

end
