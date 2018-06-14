pro plot_star_stacks

  restore, '../measurements/star_stacks/star_stack_w1_gauss7p5.idl', /v
  w1_val = val
  w1_mag = mag
  w1_bins = bin_data(alog10(w1_val), w1_mag, /nan $
                     , xmin=-1.0, xmax=2.5, binsize=0.1)
  w1_coef = median(w1_val/(10.^(-1.*w1_mag/2.5)))

  sz = size(stack)
  norm = stack
  for ii = 0, sz[3]-1 do begin
     norm[*,*,ii] = stack[*,*,ii]/val[ii]     
  endfor
  w1_beam = median(norm,dim=3)

  restore, '../measurements/star_stacks/star_stack_w2_1.idl', /v
  w2_val = val
  w2_mag = mag
  w2_bins = bin_data(alog10(w2_val), w2_mag, /nan $
                     , xmin=-1.0, xmax=2.5, binsize=0.1)
  w2_coef = median(w2_val/(10.^(-1.*w2_mag/2.5)))

  sz = size(stack)
  norm = stack
  for ii = 0, sz[3]-1 do begin
     norm[*,*,ii] = stack[*,*,ii]/val[ii]     
  endfor
  w2_beam = median(norm,dim=3)

  restore, '../measurements/star_stacks/star_stack_w3_1.idl', /v
  w3_val = val
  w3_mag = mag
  w3_bins = bin_data(alog10(w3_val), w3_mag, /nan $
                     , xmin=-1.0, xmax=2.5, binsize=0.1)
  w3_coef = median(w3_val/(10.^(-1.*w3_mag/2.5)))

  sz = size(stack)
  norm = stack
  for ii = 0, sz[3]-1 do begin
     norm[*,*,ii] = stack[*,*,ii]/val[ii]     
  endfor
  w3_beam = median(norm,dim=3)

  restore, '../measurements/star_stacks/star_stack_w4_1.idl', /v
  w4_val = val
  w4_mag = mag
  w4_bins = bin_data(alog10(w4_val), w4_mag, /nan $
                     , xmin=-1.0, xmax=2.5, binsize=0.1)
  w4_coef = median(w4_val/(10.^(-1.*w4_mag/2.5)))

  sz = size(stack)
  norm = stack
  for ii = 0, sz[3]-1 do begin
     norm[*,*,ii] = stack[*,*,ii]/val[ii]     
  endfor
  w4_beam = median(norm,dim=3)

  restore, '../measurements/star_stacks/star_stack_nuv_1.idl', /v
  nuv_val = val
  nuv_mag = mag
  nuv_bins = bin_data(alog10(nuv_val), nuv_mag, /nan $
                     , xmin=-1.0, xmax=2.5, binsize=0.1)
  nuv_coef = median(nuv_val/(10.^(-1.*nuv_mag/2.5)))

  sz = size(stack)
  norm = stack
  for ii = 0, sz[3]-1 do begin
     norm[*,*,ii] = stack[*,*,ii]/val[ii]     
  endfor
 nuv_beam = median(norm,dim=3)

  restore, '../measurements/star_stacks/star_stack_fuv_1.idl', /v
  fuv_val = val
  fuv_mag = mag
  fuv_bins = bin_data(alog10(fuv_val), fuv_mag, /nan $
                     , xmin=-1.0, xmax=2.5, binsize=0.1)
  fuv_coef = median(fuv_val/(10.^(-1.*fuv_mag/2.5)))

  sz = size(stack)
  norm = stack
  for ii = 0, sz[3]-1 do begin
     norm[*,*,ii] = stack[*,*,ii]/val[ii]     
  endfor
  fuv_beam = median(norm,dim=3)

  x = findgen(sz[1])
  x -= mean(x)
  xgrid = x # (fltarr(sz[2])+1.)
  y = findgen(sz[2])
  y -= mean(y)
  ygrid =(fltarr(sz[1])+1.) # (y)
  r = sqrt(xgrid^2+ygrid^2)
  pix = 2.75

  plot, findgen(10)
  psfile = '../plots/mag_vs_intens_wise.eps'
  pnfile = '../plots/mag_vs_intens_wise.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Intensity [MJy sr!u-1!n]' $
     , ytitle='!6log!d10!n Ks (2MASS) Magnitude' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[-1.0, 2.5], yrange=[0., 11.] $
     , /xstyle

  oplot, alog10(w4_val), w4_mag, ps=1, symsize=0.5, color=cgcolor('lightseagreen')
  oplot, alog10(w3_val), w3_mag, ps=1, symsize=0.5, color=cgcolor('goldenrod')
  oplot, alog10(w2_val), w2_mag, ps=1, symsize=0.5, color=cgcolor('salmon')
  oplot, alog10(w1_val), w1_mag, ps=1, symsize=0.5, color=cgcolor('firebrick')

  fid = 10^(findgen(101)/100.*5.-2.)

  oplot, alog10(fid) $
         , -1.*alog10(1./w1_coef*fid)*2.5, thick=3, lines=2, color=cgcolor('black')
  oplot, alog10(fid) $
         , -1.*alog10(1./w2_coef*fid)*2.5, thick=3, lines=2, color=cgcolor('black')
  oplot, alog10(fid) $
         , -1.*alog10(1./w3_coef*fid)*2.5, thick=3, lines=2, color=cgcolor('black')
  oplot, alog10(fid) $
         , -1.*alog10(1./w4_coef*fid)*2.5, thick=3, lines=2, color=cgcolor('black')

  ps, /xw
  spawn, 'evince '+psfile+' &'     

  plot, findgen(10)
  psfile = '../plots/mag_vs_intens_galex.eps'
  pnfile = '../plots/mag_vs_intens_galex.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  plot $
     , [0], [0], /nodata $
     , xtitle='!6log!d10!n Intensity [MJy sr!u-1!n]' $
     , ytitle='!6log!d10!n Ks (2MASS) Magnitude' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[-5.0, 1.0], yrange=[0., 11.] $
     , /xstyle

  oplot, alog10(nuv_val), nuv_mag, ps=1, symsize=0.5, color=cgcolor('dodgerblue')
  oplot, alog10(fuv_val), fuv_mag, ps=1, symsize=0.5, color=cgcolor('orchid')

  fid = 10^(findgen(101)/100.*10.-5.)

  oplot, alog10(fid) $
         , -1.*alog10(1./nuv_coef*fid)*2.5, thick=3, lines=2, color=cgcolor('black')
  oplot, alog10(fid) $
         , -1.*alog10(1./fuv_coef*fid)*2.5, thick=3, lines=2, color=cgcolor('black')

  ps, /xw
  spawn, 'evince '+psfile+' &'

  plot, findgen(10)
  loadct, 0
  psfile = '../plots/wise_beam.eps'
  pnfile = '../plots/wise_beam.png'
  ps, /def, /ps, xs=5, ys=5, /color, /encaps $
      , file=psfile
  
  plot $
     , [0], [0], /nodata $
     , ytitle='!6log!d10!n Max |Intensity| [Fraction of Peak]' $
     , xtitle='!6Radius ["]' $
     , xthick=5, ythick=5, charthick=3, charsize=1.5 $
     , xrange=[0., 200.], yrange=[-5., 0.] $
     , /xstyle

;  oplot, r*pix, alog10(beam), ps=1
  
  fid = findgen(201)
  oplot, fid, alog10(exp(-1.*(fid)^2/2./(15./2.354)^2)), color=cgcolor('black') $
         , thick=10

  oplot, fid, alog10(1d-3*exp(-1.*(fid)^2/2./(125./2.354)^2)), color=cgcolor('black') $
         , thick=10

  bins = bin_data(r*pix, abs(w1_beam), xmin=-1, xmax=200., binsize=2.0, /nan)
  oplot, bins.xmid, alog10(bins.ymax), color=cgcolor('firebrick'), thick=5, lines=0
;  oplot, bins.xmid, alog10(bins.ymed), color=cgcolor('firebrick'), thick=5, lines=0

  bins = bin_data(r*pix, abs(w2_beam), xmin=-1, xmax=200., binsize=2.0, /nan)
  oplot, bins.xmid, alog10(bins.ymax), color=cgcolor('salmon'), thick=5, lines=0
;  oplot, bins.xmid, alog10(bins.ymed), color=cgcolor('salmon'), thick=5, lines=0

  bins = bin_data(r*pix, abs(w3_beam), xmin=-1, xmax=200., binsize=2.0, /nan)
  oplot, bins.xmid, alog10(bins.ymax), color=cgcolor('goldenrod'), thick=5, lines=0
;  oplot, bins.xmid, alog10(bins.ymed), color=cgcolor('goldenrod'), thick=5, lines=0

  bins = bin_data(r*pix, abs(w4_beam), xmin=-1, xmax=200., binsize=2.0, /nan)
  oplot, bins.xmid, alog10(bins.ymax), color=cgcolor('lightseagreen'), thick=5, lines=0
;  oplot, bins.xmid, alog10(bins.ymed), color=cgcolor('lightseagreen'), thick=5, lines=0

  bins = bin_data(r*pix, abs(nuv_beam), xmin=-1, xmax=200., binsize=2.0, /nan)
  oplot, bins.xmid, alog10(bins.ymax), color=cgcolor('dodgerblue'), thick=5, lines=0
;  oplot, bins.xmid, alog10(bins.ymed), color=cgcolor('dodgerblue'), thick=5, lines=0

  bins = bin_data(r*pix, abs(fuv_beam), xmin=-1, xmax=200., binsize=2.0, /nan)
  oplot, bins.xmid, alog10(bins.ymax), color=cgcolor('orchid'), thick=5, lines=0
;  oplot, bins.xmid, alog10(bins.ymed), color=cgcolor('orchid'), thick=5, lines=0

  ps, /xw
  spawn, 'evince '+psfile+' &'       

  print, "Best Fit Coefficients (Note Scatter for UV):"
  print, 'WISE 1: ', w1_coef
  print, 'WISE 2: ', w2_coef
  print, 'WISE 3: ', w3_coef
  print, 'WISE 4: ', w4_coef
  print, 'NUV: ', nuv_coef
  print, 'FUV: ', fuv_coef

  stop

end
