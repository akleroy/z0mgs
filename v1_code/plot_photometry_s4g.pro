pro plot_photometry_s4g $
   , inspect=inspect

  plot, findgen(10), xtitle='!6'

  ipac = read_ipac_table('~/idl/galbase/gal_data/s4g_ipac_table.txt')
  z0mgs = mrdfits('../measurements/delivery_index_gauss7p5.fits',1,h)
  pgc_s4g= get_pgc(ipac.object)
  z0mgs_ind = lonarr(n_elements(pgc_s4g))
  for ii = 0, n_elements(pgc_s4g)-1 do $
     z0mgs_ind[ii] = where(z0mgs.pgc eq pgc_s4g[ii])

  s4g_irac1 = alog10(3631*10.^(-1.*ipac.mag1/2.5))
  es4g_irac1 = ipac.emag1/2.5
  us_wise1 = alog10(z0mgs[z0mgs_ind].flux_wise1)
  eus_wise1 = alog10((z0mgs[z0mgs_ind].std_flux_wise1+z0mgs[z0mgs_ind].flux_wise1)/z0mgs[z0mgs_ind].flux_wise1)

  s4g_irac2 = alog10(3631*10.^(-1.*ipac.mag2/2.5))
  es4g_irac2 = ipac.emag2/2.5
  us_wise2 = alog10(z0mgs[z0mgs_ind].flux_wise2)
  eus_wise2 = alog10((z0mgs[z0mgs_ind].std_flux_wise2+z0mgs[z0mgs_ind].flux_wise2)/z0mgs[z0mgs_ind].flux_wise2)

  resid2 = (us_wise2 - s4g_irac2)
  resid1 = (us_wise1 - s4g_irac1) 

  print, "Robust and STDDEV for WISE1: '", mad(resid1), stddev(resid1,/nan)
  print, "Robust and STDDEV for WISE2: '", mad(resid2), stddev(resid2,/nan)

  psfile = '../plots/s4g_integrated_wise1.eps'
  pnfile = '../plots/s4g_integrated_wise1.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  loadct, 0
  plot, [0], [0] $
        , ytitle='!6log!d10!n This Atlas WISE1 Flux [Jy]' $
        , xtitle='!6log!d10!n S4G IRAC1 Flux [Jy]' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xstyle=1, ystyle=1 $
        , xrange=[-3, 1], yrange=[-3,1]
  for ii = -20, 10 do $
     oplot, ii*0.5*[1,1], [-10,10], lines=1
  for ii = -20, 10 do $
     oplot, [-10,10], ii*0.5*[1,1], lines=1

  oploterror, s4g_irac1, us_wise1, es4g_irac1, eus_wise1, ps=cgsymcat('filledcircle'), symsize=0.25, /nohat $
              , color=cgcolor('charcoal')

  fid = findgen(101)-50
  oplot, fid, fid, thick=3, color=cgcolor('firebrick')

  al_legend, /bottom, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Offset: '+string(median(resid1), format='(F6.3)'), $
              'Scatter: '+string(mad(resid1), format='(F6.3)')] $
             , lines=-99, charthick=3, charsize=1.25

  al_legend, /top, /left, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Munoz Mateos','  et al. (2015)'] $
             , lines=-99, charthick=3, charsize=1.25

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

  psfile = '../plots/s4g_integrated_wise2.eps'
  pnfile = '../plots/s4g_integrated_wise2.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  loadct, 0
  plot, [0], [0] $
        , ytitle='!6log!d10!n This Atlas WISE2 Flux [Jy]' $
        , xtitle='!6log!d10!n S4G IRAC2 Flux [Jy]' $
        , xthick=5, ythick=5, charthick=3, charsize=1.25 $
        , xstyle=1, ystyle=1 $
        , xrange=[-3, 1], yrange=[-3,1]
  for ii = -20, 10 do $
     oplot, ii*0.5*[1,1], [-10,10], lines=1
  for ii = -20, 10 do $
     oplot, [-10,10], ii*0.5*[1,1], lines=1

  oploterror, s4g_irac2, us_wise2, es4g_irac2, eus_wise2, ps=cgsymcat('filledcircle'), symsize=0.25, /nohat $
              , color=cgcolor('charcoal')

  fid = findgen(101)-50
  oplot, fid, fid, thick=3, color=cgcolor('firebrick')

  al_legend, /bottom, /right, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Offset: '+string(median(resid2), format='(F6.3)'), $
             'Scatter: '+string(mad(resid2), format='(F6.3)')] $
             , lines=-99, charthick=3, charsize=1.25

  al_legend, /top, /left, box=1, clear=1, background=cgcolor('lightgray') $
             , ['Munoz Mateos','  et al. (2015)'] $
             , lines=-99, charthick=3, charsize=1.25

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile  

  if keyword_set(inspect) then begin

;  ind = where(resid1 lt -0.15, ct)
     ind = where(resid1 gt 0.3, ct)
     print, ct, " bad outliers high."
     for ii = 0, ct-1 do begin
        mapfile = '../delivery/PGC'+strcompress(str(pgc_s4g[ind[ii]]),/rem)+'_w1_gauss7p5.fits'
        if file_test(mapfile) eq 0 then continue
        map = readfits(mapfile, hdr)
        stars = readfits('../delivery/PGC'+strcompress(str(pgc_s4g[ind[ii]]),/rem)+'_w1_gauss7p5_found_stars.fits', hdr)
        bright = readfits('../delivery/PGC'+strcompress(str(pgc_s4g[ind[ii]]),/rem)+'_w1_gauss7p5_bright_stars.fits', hdr)
        loadct, 0
        disp, alog10(map), min=-2., max=1.0, /sq, /xs, /ys
        contour, bright eq 1 or bright eq 11, /overp, lev=[1], color=cgcolor('yellow'), thick=10
        contour, stars eq 1 or stars eq 11, /overp, lev=[1], color=cgcolor('red')
        ch = get_kbrd(1)
     endfor

  endif

  stop
  
end
