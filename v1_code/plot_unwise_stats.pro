pro plot_unwise_stats $
   , wise=wise $
   , galex=galex

  restore, '../measurements/unwise_stats.idl', /v
  dat = gal_data(stats[*,0,0,0].pgcname)
  glactc, dat.ra_deg, dat.dec_deg, 2000., l, b, 1, /degree

  plot, findgen(10), title='!6Test'

  for this_band = 0, 3 do begin

     if this_band eq 0 then begin
        xmin = -3.0
        xmax = -1.0
        binsize = 0.01
     endif
     if this_band eq 1 then begin
        xmin = -3.0
        xmax = -1.0
        binsize = 0.01
     endif
     if this_band eq 2 then begin
        xmin = -2.5
        xmax = -0.5
        binsize = 0.01
     endif
     if this_band eq 3 then begin
        xmin = -1.5
        xmax = 0.0
        binsize = 0.01
     endif

     this_med = (stats[*,2,this_band,0].med)
     this_mean = (stats[*,2,this_band,0].mean)
     this_rms =  (stats[*,2,this_band,0].rms)
     this_std =  (stats[*,2,this_band,0].std)
     
     
  endfor

  stop

end
