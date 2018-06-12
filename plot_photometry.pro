pro plot_photometry

  us = mrdfits('../measurements/z0mgs_photometry.fits', 1, h)
  band = strcompress(us.band, /rem)
  bands = ['FUV','NUV','WISE1','WISE2','WISE3','WISE4']
  
  for ii = 0, 5 do begin
     
     ind = where(band eq bands[ii])
     this = us[ind]
     
     fid = us.mean[1]
     
     for jj = 0, 5 do begin
        rat = us.mean[jj] / fid
        bins = bin_data(alog10(fid), alog10(rat) $
                       , xmin=-5.0 $
                       , xmax=-1.0 $
                       , binsize=0.1, /nan)
        plot, alog10(fid), alog10(rat), ps=3 $
              , xrange=[-5, 1], yrange=[-1., 1.]
        oploterror, bins.xmid, bins.ymed, bins.ymed-bins.ylo_perc, /lo $
                    , color=cgcolor('red'), ps=cgsymcat('filledcircle')
        oploterror, bins.xmid, bins.ymed, bins.ymed-bins.yhi_perc, /hi $
                    , color=cgcolor('red'), ps=cgsymcat('filledcircle')
        ch = get_kbrd(1)
     endfor

     stop

  endfor

end
