function bkfit_galex $
   , map=map $   
   , mask=mask $
   , wt=wt $
   , kernel=kernel $
   , method=method $
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
  if n_elements(kernel) eq 0 then $
     kernel=3

  sz = size(map)
  x = findgen(sz[1]) # (fltarr(sz[2])+1.0)
  y = (fltarr(sz[1])+1.0) #  findgen(sz[2])
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE A COPY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  orig = map
  map = smooth(map, kernel, /nan)
  noiselike = 1./sqrt(wt/median(wt))

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; OUTLIER REJECTING MEDIAN FIT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  rejected = mask eq 10 or finite(map) eq 0
  for ii = 0, niter-1 do begin
     bkind = where(rejected eq 0, ct)
     if ct eq 0 then continue
     vec = map[bkind]
     wtvec = noiselike[bkind]
     rms = mad(vec/noiselike)
     med = median(vec)
     bad_ind = where(abs(vec-med)/wtvec gt thresh*rms, bad_ct)
     if bad_ct gt 0 then $
        rejected[bkind[bad_ind]] = 1B
  endfor

  bksub = map - med
  coefs = med
  
  if method eq 'MEDIAN' then begin
     if keyword_set(show) then begin
        loadct, 0
        !p.multi=0
        disp, bksub, max=rms*5., min=-5.*rms
        contour, rejected, /overplot, lev=[1], color=cgcolor('red')
        if keyword_set(pause) then begin
           print, "Hit a key."
           ch = get_kbrd(1)
        endif
     endif

     return, bksub
  endif

end
