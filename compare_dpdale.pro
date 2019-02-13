pro compare_dpdale $
   , inspect=inspect

  data_dir = '../measurements/'

  band = ['w4','w3','w2','w1','nuv','fuv']
  name = ['WISE4','WISE3','WISE2','WISE1','NUV','FUV']
  color = ['red','salmon','goldenrod','lightseagreen','dodgerblue','orchid']

  tab = mrdfits('../measurements/delivery_index_gauss7p5.fits',1,h)

  readcol, data_dir+'dale_pgc_foradam.dat', skip=1, format='L,F,F,F,F,F,F' $
           , dale_pgc, dale_fuv, dale_nuv, dale_w1, dale_w2, dale_w3, dale_w4
  n_dale = n_elements(dale_pgc)
  dale_ind = lonarr(n_dale)
  use_dale = bytarr(n_dale)*0B
  for ii = 0, n_dale-1 do begin
     dale_ind[ii] = where(dale_pgc[ii] eq tab.pgc, ct)
     use_dale[ii] = ct gt 0
  endfor

  readcol, data_dir+'dustpedia_pgc_foradam.dat', skip=1, format='L,F,F,F,F,F,F' $
           , dp_pgc, dp_fuv, dp_nuv, dp_w1, dp_w2, dp_w3, dp_w4
  n_dp = n_elements(dp_pgc)
  dp_ind = lonarr(n_dp)
  use_dp = bytarr(n_dp)*0B
  for ii = 0, n_dp-1 do begin
     dp_ind[ii] = where(dp_pgc[ii] eq tab.pgc, ct)
     use_dp[ii] = ct gt 0
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE A FLUX-FLUX PLOT FOR EACH BAND
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  for ii = 0, 5 do begin

     if ii eq 0 then begin
        range = [-3,1.5]
        band_str = 'wise1'

        keep = where(dale_ind ne -1)
        this_dale = dale_w1[keep]
        this_dale_us = tab[dale_ind[keep]].flux_wise1
        this_dale_eus = tab[dale_ind[keep]].std_flux_wise1        

        keep = where(dp_ind ne -1)
        this_dp = dp_w1[keep]
        this_dp_us = tab[dp_ind[keep]].flux_wise1
        this_dp_eus = tab[dp_ind[keep]].std_flux_wise1        

     endif

     if ii eq 1 then begin
        range = [-3,1.5]
        band_str = 'wise2'

        keep = where(dale_ind ne -1)
        this_dale = dale_w2[keep]
        this_dale_us = tab[dale_ind[keep]].flux_wise2
        this_dale_eus = tab[dale_ind[keep]].std_flux_wise2

        keep = where(dp_ind ne -1)
        this_dp = dp_w2[keep]
        this_dp_us = tab[dp_ind[keep]].flux_wise2
        this_dp_eus = tab[dp_ind[keep]].std_flux_wise2

     endif

     if ii eq 2 then begin
        range = [-3,1.5]
        band_str = 'wise3'

        keep = where(dale_ind ne -1)
        this_dale = dale_w3[keep]
        this_dale_us = tab[dale_ind[keep]].flux_wise3
        this_dale_eus = tab[dale_ind[keep]].std_flux_wise3

        keep = where(dp_ind ne -1)
        this_dp = dp_w3[keep]
        this_dp_us = tab[dp_ind[keep]].flux_wise3
        this_dp_eus = tab[dp_ind[keep]].std_flux_wise3

     endif
     
     if ii eq 3 then begin
        range = [-3,1.5]
        band_str = 'wise4'

        keep = where(dale_ind ne -1)
        this_dale = dale_w4[keep]
        this_dale_us = tab[dale_ind[keep]].flux_wise4
        this_dale_eus = tab[dale_ind[keep]].std_flux_wise4

        keep = where(dp_ind ne -1)
        this_dp = dp_w4[keep]
        this_dp_us = tab[dp_ind[keep]].flux_wise4
        this_dp_eus = tab[dp_ind[keep]].std_flux_wise4

     endif

     if ii eq 4 then begin
        range = [-4,-0.5]
        band_str = 'fuv'

        keep = where(dale_ind ne -1)
        this_dale = dale_fuv[keep]
        this_dale_us = tab[dale_ind[keep]].flux_fuv
        this_dale_eus = tab[dale_ind[keep]].std_flux_fuv

        keep = where(dp_ind ne -1)
        this_dp = dp_fuv[keep]
        this_dp_us = tab[dp_ind[keep]].flux_fuv
        this_dp_eus = tab[dp_ind[keep]].std_flux_fuv

     endif

     if ii eq 5 then begin
        range = [-4,-0.5]
        band_str = 'nuv'

        keep = where(dale_ind ne -1)
        this_dale = dale_nuv[keep]
        this_dale_us = tab[dale_ind[keep]].flux_nuv
        this_dale_eus = tab[dale_ind[keep]].std_flux_nuv

        keep = where(dp_ind ne -1)
        this_dp = dp_nuv[keep]
        this_dp_us = tab[dp_ind[keep]].flux_nuv
        this_dp_eus = tab[dp_ind[keep]].std_flux_nuv

     endif

     resid_dp = alog10(this_dp_us) - alog10(this_dp)
     resid_dale = alog10(this_dale_us) - alog10(this_dale)

     loadct, 0 
     multiplot, /reset
     
     psfile = '../plots/z0mgs_phot_check'+band_str+'.eps'
     pnfile = '../plots/z0mgs_phot_check'+band_str+'.png'
     ps, /def, /ps, xs=5, ys=5, /color, /encaps $
         , file=psfile
     
     loadct, 0
     plot, [0], [0] $
           , ytitle='!6log!d10!n This Atlas '+strupcase(band_str)+' Flux [Jy]' $
           , xtitle='!6log!d10!n Comparison Flux [Jy]' $
           , xthick=5, ythick=5, charthick=3, charsize=1.25 $
           , xstyle=1, ystyle=1 $
           , xrange=range, yrange=range
     for kk = -20, 10 do $
        oplot, kk*0.5*[1,1], [-10,10], lines=1
     for kk = -20, 10 do $
        oplot, [-10,10], kk*0.5*[1,1], lines=1
     
     oploterror, alog10(this_dp), alog10(this_dp_us), alog10((this_dp_eus+this_dp_us)/this_dp_us) $
                 , ps=cgsymcat('filledcircle'), symsize=0.35, /nohat $
                 , color=cgcolor('cornflowerblue')

     oploterror, alog10(this_dale), alog10(this_dale_us), alog10((this_dale_eus+this_dale_us)/this_dale_us) $
                 , ps=cgsymcat('filledcircle'), symsize=0.75, /nohat $
                 , color=cgcolor('black')

     fid = findgen(101)-50
     oplot, fid, fid, thick=3, color=cgcolor('firebrick')

     al_legend, /top, /left, box=1, clear=1, background=cgcolor('lightgray') $
                , ['D17 Offset: '+string(median(resid_dale), format='(F6.3)'), $
                   'D17 Scatter: '+string(mad(resid_dale), format='(F6.3)'), $
                   'C18 Offset: '+string(median(resid_dp), format='(F6.3)'), $
                   'C18 Scatter: '+string(mad(resid_dp), format='(F6.3)')] $
                , lines=-99, charthick=3, charsize=1.25

     al_legend, /bottom, /right, box=1, clear=1, background=cgcolor('lightgray') $
                , ['Dale et al. 2017', 'Clark et al. 2018'] $               
                , psym=cgsymcat('filledcircle'), charthick=3, charsize=1.25 $
                , color=[cgcolor('black'), cgcolor('cornflowerblue')]

     ps, /xw
     spawn, 'evince '+psfile+' &'
     spawn, 'convert -density 300x300 '+psfile+' '+pnfile
     
     if keyword_set(inspect) and band_str eq 'fuv' then begin

;  ind = where(resid1 lt -0.15, ct)
        ind = where(resid_dp lt -0.3, ct)
        print, ct, " bad outliers low."
        for jj = 0, ct-1 do begin
           mapfile = '../delivery/PGC'+strcompress(str(dp_pgc[ind[jj]]),/rem)+'_fuv_gauss7p5.fits'
           if file_test(mapfile) eq 0 then continue
           map = readfits(mapfile, hdr)
           rad =  readfits('../delivery/PGC'+strcompress(str(dp_pgc[ind[jj]]),/rem)+'_gauss7p5_rgrid.fits', rhdr)
           ;bright = readfits('../delivery/PGC'+strcompress(str(dp_pgc[ind[jj]]),/rem)+'_fuv_gauss7p5_bright_stars.fits', hdr)
           loadct, 0
           ;disp, alog10(map), min=-4., max=0.0, /sq, /xs, /ys
           disp, map, min=-5d-3, max=5d-3, /sq, /xs, /ys
           contour, rad, /overp, lev=sxpar(rhdr,'FIDRAD')*[1.25,2.], color=cgcolor('yellow'), thick=3
           ch = get_kbrd(1)
        endfor

     endif

  endfor

  stop

end
