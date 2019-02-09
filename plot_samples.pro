pro plot_samples

  restore, '../measurements/z0mgs_samples.idl', /v

  ind = where(sample.pgc ne 0L and $
              sample.wise1 gt 1d-2 and $
              sample.wise3 gt 0.1 and $
              sample.nuv gt 5d-4)

  w3w1 = alog10(sample[ind].wise3/sample[ind].wise1)
  nuvw1 = alog10(sample[ind].nuv/sample[ind].wise1)

  grid = $
     grid_data(w3w1, nuvw1, /nan $
               , xaxis_out = x_axis, yaxis_out = y_axis $
               , xmin=-1., xmax=1., binsize_x=0.01 $
               , ymin=-3, ymax=0, binsize_y=0.01)
  grid /= (0.01*0.01)
  psf = psf_gaussian(fwhm=5, npix=21, /norm)
  grid = convol(grid, psf, /edge_zero)

  psfile = '../plots/nuvw1_vs_w3w1.eps'
  pnfile = '../plots/nuvw1_vs_w3w1.png'
  ps, /def, /ps, xs=8, ys=8, /color, /encaps $
      , file=psfile

  ;viridis
  loadct, 0
  reversect
  minval = 3.0
  maxval = 6.5
  disp, alog10(grid), x_axis, y_axis $
        , reserve=5, color=cgcolor('black',255) $
        , xthick=5, ythick=5 $
        , charthick=3, charsize=1.5 $
        , xtitle='!6log!d10!n WISE3-to-WISE1' $
        , ytitle='!6log!d10!n NUV-to-WISE1' $
        , max=maxval, min=minval $
        , position=[0.2, 0.2, 0.95, 0.75]

  levs = 3.0+alog10(2.^(findgen(16)))
  contour, alog10(grid), x_axis, y_axis $
           , lev=levs $
           , /overplot, c_color=cgcolor('white')

  for ii = -100, 100 do $
     oplot, ii*0.5*[1,1], [-100, 100], lines=1, color=cgcolor('gray')

  for ii = -100, 100 do $
     oplot, [-100, 100], ii*0.5*[1,1], lines=1, color=cgcolor('gray')

  xyouts, [-0.75], [-2.75], align=0.5, color=cgcolor('firebrick') $
          , ['!8Old stars!6'], charsize=1.5, charthick=7
  xyouts, [-0.75], [-2.75], align=0.5, color=cgcolor('salmon') $
          , ['!8Old stars!6'], charsize=1.5, charthick=3

  xyouts, [0.5], [-2.75], align=0.5, color=cgcolor('firebrick') $
          , ['!8Dusty!6'], charsize=1.5, charthick=7
  xyouts, [0.5], [-2.75], align=0.5, color=cgcolor('salmon') $
          , ['!8Dusty!6'], charsize=1.5, charthick=3

  xyouts, [-0.75], [-0.25], align=0.5, color=cgcolor('firebrick') $
          , ['!8Dust Free!6'], charsize=1.5, charthick=7
  xyouts, [-0.75], [-0.25], align=0.5, color=cgcolor('salmon') $
          , ['!8Dust Free!6'], charsize=1.5, charthick=3

  xyouts, [0.5], [-0.25], align=0.5, color=cgcolor('firebrick') $
          , ['!8Star Formation!6'], charsize=1.5, charthick=7
  xyouts, [0.5], [-0.25], align=0.5, color=cgcolor('salmon') $
          , ['!8Star Formation!6'], charsize=1.5, charthick=3

  ;viridis
  loadct, 0
  reversect
  cgColorBar, range=[minval, maxval] $
              , position = [0.2, 0.925, 0.95, 0.95] $
              , title='!6log!d10!n Density of Data [lines of sight dex!u-2!n]' $
              , textthick=3, charsize=1.5, xthick=5 $
              , color=cgcolor('black',255)

  ps, /xw
  spawn, 'evince '+psfile+' &'

; NOW LOOK AT BAND-TO-BAND CORRELATION

  ind = where(sample.pgc ne 0L and $
              finite(sample.wise1) and $
              finite(sample.wise2) and $
              finite(sample.wise3) and $
              finite(sample.wise4) and $
              finite(sample.nuv) and $
              finite(sample.fuv), ct)

  sig_ind = $
     where(sample.pgc ne 0L and $
           sample.wise1 gt 1d-2 and $
           sample.wise2 gt 2d-2 and $
           sample.wise3 gt 1d-1 and $
           sample.wise4 gt 5d-1 and $
           sample.nuv gt 5d-4 and $
           sample.fuv gt 5d-4 $
           , ct)
  
  label_num = 0.5+lindgen(6)
  label_vec = ['FUV','NUV', 'WISE1', 'WISE2', 'WISE3', 'WISE4'] 

  rcorr_grid = fltarr(6, 6)*!values.f_nan
  sig_rcorr_grid = fltarr(6, 6)*!values.f_nan
  for ii = 0, 5 do begin

     if ii eq 0 then x = sample[ind].fuv
     if ii eq 1 then x = sample[ind].nuv
     if ii eq 2 then x = sample[ind].wise1
     if ii eq 3 then x = sample[ind].wise2
     if ii eq 4 then x = sample[ind].wise3
     if ii eq 5 then x = sample[ind].wise4

     if ii eq 0 then sigx = sample[sig_ind].fuv
     if ii eq 1 then sigx = sample[sig_ind].nuv
     if ii eq 2 then sigx = sample[sig_ind].wise1
     if ii eq 3 then sigx = sample[sig_ind].wise2
     if ii eq 4 then sigx = sample[sig_ind].wise3
     if ii eq 5 then sigx = sample[sig_ind].wise4

     for jj = 0, 5 do begin

        if jj eq 0 then y = sample[ind].fuv
        if jj eq 1 then y = sample[ind].nuv
        if jj eq 2 then y = sample[ind].wise1
        if jj eq 3 then y = sample[ind].wise2
        if jj eq 4 then y = sample[ind].wise3
        if jj eq 5 then y = sample[ind].wise4

        if jj eq 0 then sigy = sample[sig_ind].fuv
        if jj eq 1 then sigy = sample[sig_ind].nuv
        if jj eq 2 then sigy = sample[sig_ind].wise1
        if jj eq 3 then sigy = sample[sig_ind].wise2
        if jj eq 4 then sigy = sample[sig_ind].wise3
        if jj eq 5 then sigy = sample[sig_ind].wise4
        
        rcorr_grid[ii,jj] = (r_correlate(x, y))[0]
        sig_rcorr_grid[ii,jj] = (r_correlate(sigx, sigy))[0]

     endfor

  endfor

  psfile = '../plots/band_vs_band.eps'
  pnfile = '../plots/band_vs_band.png'
  ps, /def, /ps, xs=8, ys=8, /color, /encaps $
      , file=psfile

  viridis

  disp, sig_rcorr_grid $
        , xstyle=5, ystyle=5 $
        , min=0.5, max=1.1 $
        , title='Rank Correlation Between Bands' $
        , reserve=5, color=cgcolor('black', 255) $
        , charthick=3, charsize=1.5 $
        , position=[0.15, 0.15, 0.95, 0.95]

  for ii = 0, 5 do $
     xyouts $
     , -0.5, 0.5+ii, align=0.5, orient=90. $
     , label_vec[ii], color=cgcolor('black') $
     , charsize=1.5, charthick=3

  for ii = 0, 5 do $
     xyouts $
     , 0.5+ii, -0.5, align=0.5, orient=0. $
     , label_vec[ii], color=cgcolor('black') $
     , charsize=1.5, charthick=3
  
  for ii = 0, 5 do $
     for jj = 0, 5 do $
        xyouts, 0.5+ii, 0.5+jj, align=0.5 $
                , '('+string(sig_rcorr_grid[ii,jj], format='(F5.2)')+')'  $
                , charsize=1.5, charthick=3, color=cgcolor('white')

  ps, /xw
  spawn, 'evince '+psfile+' &'


  stop

end
