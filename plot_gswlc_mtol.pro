pro plot_gswlc_mtol

  thresh = 10
  ston = 3.0

  @constants.bat
  lsun_3p4 = 1.83d18
  restore, '../gswlc/gswlc_data.idl'
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE SUBSET OF DATA TO USE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  mtol_ind = where(gws_logmstar gt 0 $
                   and w1_lum gt 0 $
                   and gws_w1 gt ston*gws_ew1 $
                   and gws_w3 gt ston*gws_ew3 $
                   and gws_w4 gt ston*gws_ew4 $
                   and gws_nuv gt ston*gws_enuv $
                   and gws_fuv gt ston*gws_efuv $
                   and gws_flagsed eq 0 $
                  )

  mass_bins_min = [9., 9.5, 10., 10.5, 11.]
  mass_bins_max = [9.5, 10., 10.5, 11., 11.5]
  n_mass_bins = n_elements(mass_bins_min)
  mass_bins_color = reverse(['red','salmon','goldenrod','lightseagreen','dodgerblue'])

  readcol $
     , 'mtol_fits.txt', format='A,F,F,F,F', comment='#' $
     , mtolfit_tag, mtolfit_lo, mtolfit_loval, mtolfit_hi, mtolfit_hival
  mtolfit_tag = strcompress(mtolfit_tag, /rem)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT HISTOGRAMS FOR BLUE AND RED GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  xmin = 0.
  xmax = 1.0
  binsize = 0.025
  
  ssfr = gws_ssfr[mtol_ind]
  vec = mtol_w1[mtol_ind]

  vec = vec[where(ssfr lt -11.)]
  vec = vec[sort(vec)]
  n = n_elements(vec)
  red_med = vec[long(0.5*n)]
  red_lo = vec[long(0.16*n)]
  red_hi = vec[long(0.84*n)]

  red_bins = bin_data(vec, vec*0.0+1.0 $
                      , xmin=xmin, xmax=xmax, binsize=binsize, /nan)
  
  ssfr = gws_ssfr[mtol_ind]
  vec = mtol_w1[mtol_ind]
  vec = vec[where(ssfr gt -11.)]
  vec = vec[sort(vec)]
  n = n_elements(vec)
  blue_med = vec[long(0.5*n)]
  blue_lo = vec[long(0.16*n)]
  blue_hi = vec[long(0.84*n)]

  blue_bins = bin_data(vec, vec*0.0+1.0 $
                       , xmin=xmin, xmax=xmax, binsize=binsize, /nan)

  psfile = '../plots/gswlc_mtol_hist.eps'
  pnfile = '../plots/gswlc_mtol_hist.png'
  ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
      , file=psfile
    
  loadct, 0
  plot $
     , [0], [0], /nodata $
     , xtitle='!6WISE1 !7T!6!d*!n [M!d!9n!6!n / L!d!9n!6!n]' $
     , ytitle='!6Fraction of Galaxies' $
     , xthick=5, ythick=5, charthick=3, charsize=1.25 $
     , xrange=[xmin, xmax], yrange=[0., 0.25] $
     , /xstyle

  for ii = -100, 100 do $
     oplot, ii*binsize*10.*[1,1], [-10, 10], lines=1, color=cgcolor('charcoal')

  for ii = -100, 100 do $
     oplot, [-10, 10], ii*0.05*[1,1], lines=1, color=cgcolor('charcoal')

  histplot $
     , red_bins.xmid, ((red_bins.counts/total(red_bins.counts,/nan) > (0.))) $
     , /overplot $
     , lthick=3, /fill $
     , fcolor=cgcolor('salmon') $
     , lcolor=cgcolor('firebrick')

  histplot $
     , blue_bins.xmid, ((blue_bins.counts/total(blue_bins.counts,/nan) > (0.))) $
     , /overplot $
     , lthick=3, /fill, /fline, forient=45, /nobars $
     , fcolor=cgcolor('cornflowerblue') $
     , lcolor=cgcolor('navy')

  al_legend, /top, /right, lines=-99, charsize=1.25, charthick=3 $
             , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
             , textcolor=[cgcolor('firebrick')] $
             , ['!6low sSFR' $
                , '!6med: '+string(red_med,format='(F4.2)') $
                , '!6('+string(red_lo,format='(F4.2)')+'-'+'!6'+string(red_hi,format='(F4.2)')+')']

  al_legend, /top, /left, lines=-99, charsize=1.25, charthick=3 $
             , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
             , textcolor=[cgcolor('navy')] $
             , ['!6high sSFR' $
                , '!6med: '+string(blue_med,format='(F4.2)') $
                , '!6('+string(blue_lo,format='(F4.2)')+'-'+'!6'+string(blue_hi,format='(F4.2)')+')']

  ps, /xw
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile   
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MASS TO LIGHT RATIO VS QUANTITIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  plot, findgen(101), xtitle='!6'

  n_quant = 10

  for ii = 0, n_quant-1 do begin

     if ii ne 2 and ii ne 3 then continue

     if ii eq 0 then begin
        x = gws_logmstar
        xtitle = '!6log!d10!n Stellar Mass from CIGALE [M!d!9n!6!n]'
        tag = 'mstar'
        xmin = 8.
        xmax = 12.
        binmin = 8.5
        binmax = 11.75
        binsize = 0.1
     endif

     if ii eq 1 then begin
        x = gws_ssfr
        xtitle = '!6log!d10!n SFR/M!d*!n from CIGALE [yr!u-1!n]'
        tag = 'ssfr'
        xmin = -12.0
        xmax = -9.0
        binmin = -11.75
        binmax = -8.5
        binsize = 0.1
     endif

     if ii eq 2 then begin
        x = alog10(ssfr_like_fuvw4)
        xtitle = '!6log!d10!n SFR(FUV+W4)/W1 Lum. [M!d!9n!6!n yr!u-1!n L!d!9n!6!n!u-1!n]'
        tag = 'ssfrlike'
        xmin = -12.5
        xmax = -9.5
        binmin = -12.5
        binmax = -10.0
        binsize = 0.1
     endif

     if ii eq 3 then continue
     ;if ii eq 3 then begin
     ;   x = alog10(ssfr_like_nuvw4)
     ;   xtitle = '!6log!d10!n SFR(NUV+W4)/W1 Lum. [M!d!9n!6!n yr!u-1!n L!d!9n!6!n!u-1!n]'
     ;   tag = 'ssfrlike'
     ;   xmin = -12.5
     ;   xmax = -9.5
     ;   binmin = -12.5
     ;   binmax = -10.0
     ;   binsize = 0.1
     ;endif

     if ii eq 4 then begin
        x = w2w1
        xtitle = '!6log!d10!n W2-to-W1'
        tag = 'w2w1'
        xmin = -0.40
        xmax = 0.05
        binmin = -0.30
        binmax = -0.10
        binsize = 0.025
     endif

     if ii eq 5 then begin
        x = w4w1
        xtitle = '!6log!d10!n W4-to-W1'
        tag = 'w4w1'
        xmin = -1.5
        xmax = 1.25
        binmin = -1.25
        binmax = 1.0
        binsize = 0.1
     endif

     if ii eq 6 then begin
        x = w3w1
        xtitle = '!6log!d10!n W3-to-W1'
        tag = 'w3w1'
        xmin = -1.5
        xmax = 1.25
        binmin = -1.25
        binmax = 1.0
        binsize = 0.1
     endif

     if ii eq 7 then begin
        x = fuvw1
        xtitle = '!6log!d10!n FUV-to-W1'
        tag = 'fuvw1'
        xmin = -3.5
        xmax = 0.0
        binmin = -3.0
        binmax = -0.25
        binsize = 0.1
     endif

     if ii eq 8 then begin
        x = nuvw1
        xtitle = '!6log!d10!n NUV-to-W1'
        tag = 'nuvw1'
        xmin = -3.5
        xmax = 0.0
        binmin = -3.0
        binmax = -0.25
        binsize = 0.1
     endif

     if ii eq 9 then begin
        x = alog10(w1_lum/lsun_3p4)
        xtitle = '!6log!d10!n W1 Luminosity [L!d!9n!6!n]'
        tag = 'wise1'
        xmin = 8.5
        xmax = 12.5
        binmin = 9.0
        binmax = 12.0
        binsize = 0.1
     endif

     y = mtol_w1
     mtol_grid = grid_data(x[mtol_ind], y[mtol_ind], /nan $
                           , xaxis_out = xaxis_mtol, yaxis_out = yaxis_mtol $
                           , ymin=0.0, ymax=0.75, binsize_y=0.01 $
                           , xmin=xmin, xmax=xmax, binsize_x=(xmax-xmin)/100.0)

     psfile = '../plots/gswlc_mtol_'+tag+'.eps'
     pnfile = '../plots/gswlc_mtol_'+tag+'.png'
     ps, /def, /ps, xs=5, ys=3.5, /color, /encaps $
         , file=psfile
     
     loadct, 0
     reversect

     disp, alog10(mtol_grid) $
           , xaxis_mtol, yaxis_mtol, /xstyle, /ystyle $
           , xthick=5, ythick=5 $
           , xtitle=xtitle $
           , ytitle='!6WISE1 !7T!6!d*!n [M!d!9n!6!n / L!d!9n!6!n]' $
           , charsize=1.25, charthick=3 $
           , reserve=50, color=cgcolor('black',255) $
           , min=0.0, max=3.0

     if tag eq 'w2w1' then $
        for jj = -100, 100 do $
           oplot, jj*0.1*[1,1], [0, 10], lines=1, color=cgcolor('black') $
     else $
        for jj = -100, 100 do $
           oplot, jj*0.5*[1,1], [0, 10], lines=1, color=cgcolor('black')
     for jj = -100, 100 do $
        oplot, [-100, 100], jj*0.2*[1,1], lines=1, color=cgcolor('black')
     
     contour, alog10(mtol_grid), xaxis_mtol, yaxis_mtol $
              , lev=(findgen(8)+1.)*0.5, /overplot, color=cgcolor('black')
     
;    SPECIAL FOR THIS ONE BECAUSE OF THE HISTORIC FITS

     if tag eq 'w2w1' then begin
        xfid = findgen(101)/100. - 1.
        
        m14 = alog10(3.98*(-1.0*2.5*xfid)+0.13)
        yfid = m14

        oplot, xfid, yfid, thick=15, lines=0, color=cgcolor('lightseagreen')
        oplot, xfid, yfid, thick=7, lines=0, color=cgcolor('darkgreen')

        q15 = 0.316*10.^(-1.*xfid) < 0.6
        yfid = q15

        oplot, xfid, yfid, thick=15, lines=0, color=cgcolor('cornflowerblue')
        oplot, xfid, yfid, thick=7, lines=0, color=cgcolor('royalblue')

        al_legend, /bottom, /left, lines=-99, charsize=1.0, charthick=3 $
                   , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
                   , textcolor=[cgcolor('royalblue'), cgcolor('darkgreen')] $
                   , ['!6Meidt et al. 2014', '!6Querejeta et al. 2015']

     endif

     tag_ind = where(tag eq mtolfit_tag, tag_ct)
     if tag_ct eq 1 then begin
        tag_ind = tag_ind[0]
        xfid = findgen(10001)/10000.*40 + (-20.)
        pivot1 = mtolfit_lo[tag_ind]
        pivot2 = mtolfit_hi[tag_ind]
        zero_pt = mtolfit_loval[tag_ind]
        delta = mtolfit_hival[tag_ind] - mtolfit_loval[tag_ind]
        yfid = zero_pt + $
               (xfid gt pivot1 and xfid le pivot2)*(xfid-pivot1)/(pivot2-pivot1)*(delta) + $
               (xfid gt pivot2)*delta
        oplot, xfid, yfid, thick=15, lines=0, color=cgcolor('black')
        oplot, xfid, yfid, thick=7, lines=0, color=cgcolor('white')
     endif
     
     for jj = 0, n_mass_bins-1 do begin

;       DON'T PLOT FOR WISE1
        if tag eq 'wise1' then continue

        this_data = where(gws_logmstar[mtol_ind] ge mass_bins_min[jj] and $
                          gws_logmstar[mtol_ind] le mass_bins_max[jj], this_ct)
        if this_ct lt 100 then continue

        this_bins = bin_data(x[mtol_ind[this_data]], y[mtol_ind[this_data]] $
                             , xmin=binmin, xmax=binmax, binsize=binsize, /nan)
        plot_ind = where(this_bins.counts ge 100)
        shift = 0.05*binsize
        oploterror, this_bins[plot_ind].xmid+shift*(jj-3) $
                    , this_bins[plot_ind].ymed, this_bins[plot_ind].ymad_log $
                    , color=cgcolor(mass_bins_color[jj]) $
                    , psym=cgsymcat('filledsquare'), symsize=0.75 $
                    , errthick=3, /nohat        
     endfor

     ps, /xw
     spawn, 'evince '+psfile+' &'
     spawn, 'convert -density 300x300 '+psfile+' '+pnfile

;     stop

  endfor

  stop

end
