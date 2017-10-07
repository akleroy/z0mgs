pro compare_wise_s4g

  tab = mrdfits('../measurements/s4g_comp.fits',1,h)

  bins = bin_data(alog10(tab.i1), tab.w1 $
                  , xmin=-1.5, xmax=3., binsize=0.05, /nan)

  bins_nobk = bin_data(alog10(tab.i1), tab.w1+tab.bk_w1 $
                       , xmin=-1.5, xmax=3., binsize=0.05, /nan)

  bins_s4gbk = bin_data(alog10(tab.i1-tab.bk_i1), tab.w1+tab.bk_w1 $
                        , xmin=-1.5, xmax=3., binsize=0.05, /nan)

  plot, findgen(10), xtitle='!6'

  psfile = '../plots/s4g_vs_w1.eps'
  ps, /def, /ps, xs=8, ys=8, /color, /encaps $
      , file=psfile

  ploterror, bins.xmid, alog10(bins.ymed), bins.ymad_log $
             , ps=cgsymcat('opencircle') $
             , color=cgcolor('black') $
             , xtitle='!6log!d10!n IRAC 1 from S4G [MJy/sr]' $
             , ytitle='!6log!d10!n WISE 1 BINS from Z0MGS [MJy/sr]' $
             , xthick=5, ythick=5, charthick=3

  oploterror, bins_nobk.xmid, alog10(bins_nobk.ymed), bins_nobk.ymad_log $
             , ps=cgsymcat('opencircle') $
             , color=cgcolor('red')

  oploterror, bins_s4gbk.xmid, alog10(bins_s4gbk.ymed), bins_s4gbk.ymad_log $
             , ps=cgsymcat('filledstar') $
             , color=cgcolor('blue')

  equality, color=cgcolor('blue')

  al_legend, /top, /left, lines=-99, textcolor=[cgcolor('black'), cgcolor('red'), cgcolor('blue')] $
             , ['!6Both as is.','!6No background on WISE','!6Background from S4G']

  ps, /xw

  bins = bin_data(alog10(tab.i2), tab.w2 $
                  , xmin=-2, xmax=2.5, binsize=0.05, /nan)

  bins_nobk = bin_data(alog10(tab.i2), tab.w2+tab.bk_w2 $
                       , xmin=-2, xmax=2.5, binsize=0.05, /nan)

  bins_s4gbk = bin_data(alog10(tab.i2-tab.bk_i2), tab.w2+tab.bk_w2 $
                        , xmin=-2, xmax=2.5, binsize=0.05, /nan)

  psfile = '../plots/s4g_vs_w2.eps'
  ps, /def, /ps, xs=8, ys=8, /color, /encaps $
      , file=psfile

  ploterror, bins.xmid, alog10(bins.ymed), bins.ymad_log $
             , ps=cgsymcat('opencircle') $
             , color=cgcolor('black') $
             , xtitle='!6log!d10!n IRAC 2 from S4G [MJy/sr]' $
             , ytitle='!6log!d10!n WISE 2 BINS from Z0MGS [MJy/sr]' $
             , xthick=5, ythick=5, charthick=3

  oploterror, bins_nobk.xmid, alog10(bins_nobk.ymed), bins_nobk.ymad_log $
             , ps=cgsymcat('opencircle') $
             , color=cgcolor('red')

  oploterror, bins_s4gbk.xmid, alog10(bins_s4gbk.ymed), bins_s4gbk.ymad_log $
             , ps=cgsymcat('filledstar') $
             , color=cgcolor('blue')

  equality, color=cgcolor('blue')

  al_legend, /top, /left, lines=-99, textcolor=[cgcolor('black'), cgcolor('red'), cgcolor('blue')] $
             , ['!6Both as is.','!6No background on WISE','!6Background from S4G']

  ps, /xw

  stop

end
