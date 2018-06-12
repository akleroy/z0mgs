function bkfit_galex $
   , map=map $   
   , mask=mask $
   , wt=wt $
   , niter=niter $
   , thresh=thresh $
   , rejected=rejected $
   , coefs=coefs $
   , pause=pause $
   , show=show

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
; MAKE A COPY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  noiselike = 1./sqrt(wt/median(wt))

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; OUTLIER REJECTING MEDIAN FIT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  aperture = mask eq 100 and finite(map)
  rejected = mask*0B

  for ii = 0, niter-1 do begin
     bkind = where(aperture, ct)
     if ct eq 0 then continue
     vec = map[bkind]
     wtvec = noiselike[bkind]
     rms = mad(vec/wtvec)
     med = median(vec)
     bad_ind = where(abs(vec-med)/wtvec gt thresh*rms, bad_ct)
     if bad_ct gt 0 then $
        aperture[bkind[bad_ind]] = 0B    
  endfor  
  vec_mean = mean(vec,/nan)
  vec_med = mean(vec,/nan)
  vec_rms = mad(vec)

  bksub = map - vec_med ;vec_mean ;med
  coefs = med
  rejected = mask ne 10 and (abs(bksub) gt thresh*rms)
  
  if keyword_set(show) then begin
     loadct, 0
     !p.multi=0
     rms = stddev(bksub,/nan)
     disp, bksub, max=vec_rms*3., min=-3.*vec_rms, /sq
     contour, rejected, /overplot, lev=[1], color=cgcolor('red')
     contour, mask, /overplot, lev=[10,100], color=cgcolor('blue')
     if keyword_set(pause) then begin
        print, "Hit a key."
        ch = get_kbrd(1)
     endif
  endif

  return, bksub

end
