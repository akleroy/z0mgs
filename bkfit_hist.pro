function bkfit_hist $
   , map=map $   
   , mask=mask $
   , rejected=rejected $
   , thresh=rej_thresh $
   , show=show

  if n_elements(rej_thresh) then $
     rej_thresh = 0.1

  rejected = mask eq 10 or $
             finite(map) eq 0 or $
             rejected ge rej_thresh

  fit_ind = where(rejected eq 0, fit_ct)
  if fit_ct eq 0 then return, !values.f_nan
  
  rms = mad(map[fit_ind])
  bins = bin_data(map[fit_ind], map[fit_ind]*0.0+1.0 $
                  , xmin=-10.*rms, xmax=10.*rms $
                  , binsize=0.1*rms, /nan)
  dummy = max(median(bins.counts,3), /nan, maxind)
  
  if keyword_set(show) then begin
     fit = gaussfit(bins.xmid, bins.counts, coefs, nterms=3)
     plot, bins.xmid, bins.counts, ps=10
     oplot, bins.xmid, fit, color=cgcolor('red')
     oplot, coefs[1]*[1,1], [-1d6, 1d6], color=cgcolor('red')
     oplot, bins[maxind].xmid*[1,1], [-1d6, 1d6], color=cgcolor('yellow')
  endif

  return, bins[maxind].xmid ; coefs[1]

end
