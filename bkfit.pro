function bkfit $
   , map=map $   
   , mask=mask $
   , method=method $
   , niter=niter $
   , thresh=thresh $
   , rejected=rejected $
   , coefs=coefs

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; TUNING PARAMETERS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(niter) eq 0 then $
     niter = 5
  if n_elements(thresh) eq 0 then $
     thresh = 3.0

  sz = size(map)
  x = findgen(sz[1]) # (fltarr(sz[2])+1.0)
  y = (fltarr(sz[1])+1.0) #  findgen(sz[2])
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; OUTLIER REJECTING MEDIAN FIT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  rejected = mask eq 10 or finite(map) eq 0
  for ii = 0, niter-1 do begin
     bkind = where(rejected eq 0, ct)
     if ct eq 0 then continue
     vec = map[bkind]
     rms = mad(vec)
     med = median(vec)
     bad_ind = where(abs(vec-med) gt thresh*rms, bad_ct)
     if bad_ct gt 0 then $
        rejected[bkind[bad_ind]] = 1B
  endfor

  bksub = map - med
  coefs = med
  
  if method eq 'MEDIAN' then begin
     return, bksub
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLANE FIT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = 0, niter-1 do begin
     bkind = where(rejected eq 0, ct)
     if ct eq 0 then continue
     this_x = x[bkind]
     this_y = y[bkind]
     vec = map[bkind]
     
     coefs = planefit(this_x,this_y,vec,0., yfit)
     resid = vec - yfit
     rms = mad(resid)
     bad_ind = where(abs(resid) gt thresh*rms, bad_ct)
     if bad_ct gt 0 then $
        rejected[bkind[bad_ind]] = 1B     
  endfor

  fit = coefs[0] + coefs[1]*x + coefs[2]*y
  bksub = map-fit

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  loadct, 0
  !p.multi=0
  disp, bksub, max=rms*5., min=-5.*rms
  contour, rejected, /overplot, lev=[1], color=cgcolor('red')

  return, bksub

end
