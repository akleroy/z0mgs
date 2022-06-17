function bkfit $
   , map=map $   
   , mask=mask $
   , method=method $
   , niter=niter $
   , thresh=thresh $
   , rejected=rejected $
   , coefs=coefs $
   , show=show $
   , pause=pause

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

  aperture = mask eq 100 and finite(map)
  rejected = mask*0B

  for ii = 0, niter-1 do begin
     bkind = where(aperture, ct)
     if ct eq 0 then begin
        message, 'No background pixels', /info        
        return, map
        continue
     endif
     vec = map[bkind]
     rms = mad(vec)
     med = median(vec)
     bad_ind = where(abs(vec-med) gt thresh*rms, bad_ct)
     if bad_ct gt 0 then $
        aperture[bkind[bad_ind]] = 0B    
  endfor  

  vec_med = mean(vec,/nan)
  vec_rms = mad(vec)

  bksub = map - vec_med ;vec_mean ;med
  coefs = med
  rejected = (mask ne 10 and (abs(bksub) gt thresh*rms)) or ((mask mod 10) eq 1)

  if method eq 'MEDIAN' then begin
     if keyword_set(show) then begin
        loadct, 0
        !p.multi=0
        disp, bksub, max=rms*3., min=-3.*rms
        contour, rejected, /overplot, lev=[1], color=cgcolor('red')
        contour, mask, /overplot, lev=[10,100], color=cgcolor('blue')
        if keyword_set(pause) then ch = get_kbrd(1)
     endif

     return, bksub
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLANE FIT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = 0, niter-1 do begin
     bkind = where(aperture, ct)
     if ct eq 0 then continue
     this_x = x[bkind]
     this_y = y[bkind]
     vec = map[bkind]
     
     coefs = planefit(this_x,this_y,vec,0., yfit)
     resid = vec - yfit
     rms = mad(resid)
     bad_ind = where(abs(resid) gt thresh*rms, bad_ct)
     if bad_ct gt 0 then $
        aperture[bkind[bad_ind]] = 0B
  endfor

  fit = coefs[0] + coefs[1]*x + coefs[2]*y
  bksub = map-fit

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(show) then begin
     loadct, 0
     !p.multi=0
     disp, bksub, max=rms*3., min=-3.*rms
     contour, rejected, /overplot, lev=[1], color=cgcolor('red')
     contour, mask, /overplot, lev=[10,100], color=cgcolor('blue')
     if keyword_set(pause) then ch = get_kbrd(1)
  endif

  return, bksub

end
