pro plot_galex_time

  tab = mrdfits('../measurements/delivery_index.fits',1,h)
  b = tab.gb_deg

  plot, findgen(10), title='!6Test'
  
  xmin = 1.0
  xmax = 6.0
  binsize = 0.1

  vec = tab[where(tab.has_fuv)].time_fuv
  bins_fuv = $
     bin_data(alog10(vec), vec*0.0+1.0 $
              , xmin=xmin, xmax=xmax, binsize=binsize, /nan)

  vec = tab[where(tab.has_nuv)].time_nuv
  bins_nuv = $
     bin_data(alog10(vec), vec*0.0+1.0 $
              , xmin=xmin, xmax=xmax, binsize=binsize, /nan)
  
  psfile = '../plots/galex_time.eps'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Effective Integration [s]' $
     , ytitle='!6log!d10!n Number of Images' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[xmin, xmax], yrange=[0., 5.]
  
  for ii = -100, 100 do $
     oplot, ii*binsize*10.*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')

  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.25*[1,1], lines=1, color=cgcolor('charcoal')
  
  histplot $
     , bins_nuv.xmid, (alog10(bins_nuv.counts) > (0.)) $
     , /overplot $
     , lthick=3, /nobar, /fill $
     , fcolor=cgcolor('salmon') $
     , lcolor=cgcolor('firebrick')

  histplot $
     , bins_fuv.xmid, (alog10(bins_fuv.counts) > (0.)) $
     , /overplot $
     , lthick=3, /nobar, /fline, forient=45 $
     , fcolor=cgcolor('royalblue') $
     , lcolor=cgcolor('royalblue')
  
  histplot $
     , bins_fuv.xmid, (alog10(bins_fuv.counts) > (0.)) $
     , /overplot $
     , lthick=3, /nobar, /fline, forient=-45 $
     , fcolor=cgcolor('royalblue') $
     , lcolor=cgcolor('royalblue')
  
  al_legend $
     , /top, /right $
     , box=1, clear=1 $
     , background=cgcolor('lightgray') $
     , charsize=1.75, charthick=3 $
     , lines=-99 $
     , ['NUV','FUV'] $
     , textcolor=[cgcolor('salmon'), cgcolor('royalblue')]

  ps, /xw
  spawn, 'evince '+psfile+' &'
  
  stop

end
