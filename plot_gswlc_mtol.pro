pro plot_gswlc_mtol

  @constants.bat
  lsun_3p4 = 1.83d18
  thresh = 10

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PREPARE THE GWSLC MEASUREMENTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  gws = mrdfits('~/idl/galbase/gal_data/'+$
                'hlsp_gswlc_galex-sdss-wise_multi_x1_multi_v1_cat.fits' $
                , 1, h_gws)
  readcol, '/data/kant/0/leroy.42/allsky/gswlc/galex_unwise_fluxes_GSWLC-X.dat' $
           , id, gws_fuv, gws_efuv, gws_nuv, gws_enuv, gws_w1_nm, gws_ew1 $
           , gws_w2_nm, gws_ew2, gws_w3_nm, gws_ew3, gws_w4_nm, gws_ew4 $
           , format='L,F,F,F,F,F,F,F,F,F,F,F,F'

  gws_w1 = 3631.*10^(-0.4*(22.5+2.683))*gws_w1_nm
  gws_w2 = 3631*10^(-0.4*(22.5+3.319))*gws_w2_nm
  gws_w3 = 3631*10^(-0.4*(22.5+5.242))*gws_w3_nm
  gws_w4 = 3631*10^(-0.4*(22.5+6.604))*gws_w4_nm

  nu_fuv = c/(154.d-9*1d2)
  nu_nuv = c/(231.d-9*1d2)
  nu_w1 = c/(3.4d-6*1d2)
  nu_w2 = c/(4.5d-6*1d2)
  nu_w3 = c/(12.d-6*1d2)
  nu_w4 = c/(22.d-6*1d2)

  fuv_lum = gws_fuv/1d3*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2
  nuv_lum = gws_nuv/1d3*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2
  w1_lum = gws_w1*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2
  w2_lum = gws_w2*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2
  w3_lum = gws_w3*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2
  w4_lum = gws_w4*1d-23*4.*!pi*(gws.z*c/1d5/70.*1d6*pc)^2

  sfr_fuv_ke12 = $
     lum_to_sfr(band='FUV', cal='KE12', lum=fuv_lum*nu_fuv)
  sfr_nuv_ke12 = $
     lum_to_sfr(band='NUV', cal='KE12', lum=nuv_lum*nu_nuv)
  sfr_w3_j13 = $
     lum_to_sfr(band='WISE3', cal='J13', lum=w3_lum*nu_w3)
  sfr_w4_j13 = $
     lum_to_sfr(band='WISE4', cal='J13', lum=w4_lum*nu_w4)

  sfr_fuvw4_ke12 = $
     lum_to_sfr(band='FUV', cal='KE12' $
                , lum=(fuv_lum*nu_fuv + 3.89*w4_lum*nu_w4))  
  sfr_nuvw3 = $
     sfr_nuv_ke12+sfr_w3_j13

  mtol_w1 = 10.^(gws.logmstar - alog10(w1_lum/lsun_3p4))
  gws_ssfr = gws.logsfrsed-gws.logmstar
  ssfr_like = sfr_nuvw3 / (w1_lum / lsun_3p4)
  w2w1 = alog10(w2_lum/w1_lum)
  w3w1 = alog10(w3_lum/w1_lum)
  nuvw1 = alog10(nuv_lum/w1_lum)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE SUBSET OF DATA TO USE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sane_ind = where(gws.logmstar gt 0 $
                   and w1_lum gt 0 $
                   and gws_w1_nm gt 3.*gws_ew1 $
                   and gws_w3_nm gt 3.*gws_ew3 $
                   and gws_nuv gt 3.*gws_enuv $
                   and mtol_w1 gt 0.02 and mtol_w1 lt 1.0 $
                  )

  mass_bins_min = [9., 9.5, 10., 10.5, 11.]
  mass_bins_max = [9.5, 10., 10.5, 11., 11.5]
  n_mass_bins = n_elements(mass_bins_min)
  mass_bins_color = reverse(['red','salmon','goldenrod','lightseagreen','dodgerblue'])

  stop

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MASS TO LIGHT RATIO VS QUANTITIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  plot, findgen(101), xtitle='!6'

  n_quant = 7

  for ii = 0, n_quant-1 do begin

     if ii eq 0 then begin
        x = gws.logmstar
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
        x = alog10(ssfr_like)
        xtitle = '!6log!d10!n SFR(NUV+W3)/W1 Lum.  [M!d!9n!6!n yr!u-1!n L!d!9n!6!n!u-1!n]'
        tag = 'ssfrlike'
        xmin = -12.5
        xmax = -9.5
        binmin = -12.5
        binmax = -10.0
        binsize = 0.1
     endif

     if ii eq 3 then begin
        x = w2w1
        xtitle = '!6log!d10!n W2-to-W1'
        tag = 'w2w1'
        xmin = -0.40
        xmax = 0.05
        binmin = -0.30
        binmax = -0.10
        binsize = 0.025
     endif

     if ii eq 4 then begin
        x = w3w1
        xtitle = '!6log!d10!n W3-to-W1'
        tag = 'w3w1'
        xmin = -1.5
        xmax = 1.25
        binmin = -1.25
        binmax = 1.0
        binsize = 0.1
     endif

     if ii eq 5 then begin
        x = nuvw1
        xtitle = '!6log!d10!n NUV-to-W1'
        tag = 'nuvw1'
        xmin = -3.5
        xmax = 0.0
        binmin = -3.0
        binmax = -0.25
        binsize = 0.1
     endif

     if ii eq 6 then begin
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
     mtol_grid = grid_data(x[sane_ind], y[sane_ind], /nan $
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

     for jj = -100, 100 do $
        oplot, jj*0.5*[1,1], [0, 10], lines=1, color=cgcolor('black')
     for jj = -100, 100 do $
        oplot, [-100, 100], jj*0.2*[1,1], lines=1, color=cgcolor('black')
          
     contour, alog10(mtol_grid), xaxis_mtol, yaxis_mtol $
              , lev=(findgen(8)+1.)*0.5, /overplot, color=cgcolor('black')
     
     ;if ii eq 2 then begin
     ;   xfid = findgen(101)/100.*5 + (-13.)
     ;   pivot1 = -11.5
     ;   pivot2 = -10.
     ;   zero_pt = 0.475
     ;   delta = -0.30
     ;   yfid = zero_pt + $
     ;          (xfid gt pivot1 and xfid le pivot2)*(xfid-pivot1)/(pivot2-pivot1)*(delta) + $
     ;          (xfid gt pivot2)*delta
     ;   oplot, xfid, yfid, thick=10, lines=0, color=cgcolor('black')
     ;   oplot, xfid, yfid, thick=3, lines=0, color=cgcolor('white')
     ;   al_legend, /top, /right $
     ;              , box=1, clear=1, background=cgcolor('lightgray') $
     ;              , lines=2, thick=5, ['Prescription'] $
     ;              , pspacing=1. $
     ;              , charthick=3, charsize=1.25
     ;endif
  
     if ii eq 3 then begin
        xfid = findgen(101)/100. - 1.
        
        m14 = alog10(3.98*(-1.0*2.5*xfid)+0.13)
        yfid = m14

        oplot, xfid, yfid, thick=15, lines=0, color=cgcolor('lightseagreen')
        oplot, xfid, yfid, thick=7, lines=0, color=cgcolor('darkgreen')

        ;yfid = 10.^(-0.339*((-1*2.5*xfid+alog10(280.9/179.9)))-0.336)
        q15 = 0.316*10.^(-1.*xfid) < 0.6
        yfid = q15

        oplot, xfid, yfid, thick=15, lines=0, color=cgcolor('cornflowerblue')
        oplot, xfid, yfid, thick=7, lines=0, color=cgcolor('royalblue')

        al_legend, /bottom, /left, lines=-99, charsize=1.0, charthick=3 $
                   , box=1, outline=cgcolor('black'), background=cgcolor('lightgray') $
                   , textcolor=[cgcolor('royalblue'), cgcolor('darkgreen')] $
                   , ['!6Meidt et al. 2014', '!6Querejeta et al. 2015']

;        al_legend, /top, /right $
;                   , box=1, clear=1, background=cgcolor('lightgray') $
;                   , color=[cgcolor('black'), cgcolor('royalblue')] $
;                   , lines=[0,0], thick=10, ['Meidt et al. (2014)', 'Querejeta et al. (2015)'] $
;                   , pspacing=1. $
;                   , charthick=3, charsize=1.25
     endif

    if ii eq 4 then begin
        xfid = findgen(101)/100.*5 + (-2.)
        pivot1 = 0.1
        pivot2 = 0.85
        zero_pt = 0.525
        delta = -0.30
        yfid = zero_pt + $
               (xfid gt pivot1 and xfid le pivot2)*(xfid-pivot1)/(pivot2-pivot1)*(delta) + $
               (xfid gt pivot2)*delta
        oplot, xfid, yfid, thick=15, lines=0, color=cgcolor('black')
        oplot, xfid, yfid, thick=7, lines=0, color=cgcolor('white')
     endif

    if ii eq 5 then begin
        xfid = findgen(1001)/1000.*10 + (-5.)
        pivot1 = -1.9
        pivot2 = -0.4
        zero_pt = 0.525
        delta = -0.30
        yfid = zero_pt + $
               (xfid gt pivot1 and xfid le pivot2)*(xfid-pivot1)/(pivot2-pivot1)*(delta) + $
               (xfid gt pivot2)*delta
        oplot, xfid, yfid, thick=15, lines=0, color=cgcolor('black')
        oplot, xfid, yfid, thick=7, lines=0, color=cgcolor('white')
     endif
     
    if ii eq 6 then begin

       ;xfid = findgen(101)/10.*10.+3.
       ;yfid = 10.^(-1.84)*(10.^(xfid))^(0.10)
       ;oplot, xfid, yfid, thick=15, lines=0, color=cgcolor('black')
       ;oplot, xfid, yfid, thick=7, lines=0, color=cgcolor('white')

       ;xfid = findgen(101)/10.*10.+3.
       ;yfid = 10.^(-1.09)*(10.^(xfid))^(0.05)
       ;oplot, xfid, yfid, thick=15, lines=0, color=cgcolor('firebrick')
       ;oplot, xfid, yfid, thick=7, lines=0, color=cgcolor('salmon')

    endif

     for jj = 0, n_mass_bins-1 do begin
        if ii eq 6 then continue

        this_data = where(gws[sane_ind].logmstar ge mass_bins_min[jj] and $
                          gws[sane_ind].logmstar le mass_bins_max[jj], this_ct)
        if this_ct lt 100 then continue

        this_bins = bin_data(x[sane_ind[this_data]], y[sane_ind[this_data]] $
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
