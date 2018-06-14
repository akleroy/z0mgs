pro build_unsharp_mask $
   , infile = infile $
   , map = map $
   , hdr = hdr $
   , outfile=outfile $
   , n_star=n_reg $
   , star_ra=unsharp_ra $
   , star_dec=unsharp_dec $
   , star_intens=unsharp_intens $
   , star_km=unsharp_km $
   , show=show $
   , pause=pause

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if n_elements(map) eq 0 or n_elements(hdr) eq 0 then begin
     if file_test(infile) eq 0 then begin
        message, 'Target file not found and map and header not supplied.', /info
        return
     endif     
     map = readfits(infile, hdr, /silent)  
  endif

  sz = size(map)
  make_axes, hdr, ri=ri, di=di, ra=ra, da=da
  
  nan = !values.f_nan

  if n_elements(band) eq 0 then $
     band = 'w1'

  if n_elements(thresh) eq 0 then $
     thresh = 0.1

  if n_elements(filter) eq 0 then $
     filter = 5

; I think that we expect a peak to filter ratio of about 4.44 for a
; Gaussian with a filter width 5 pixels and 7.5" FWHM and 2.75" pixels.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RUN THE FILTERING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  med_map = median(map, filter)
  
  unsharp_diff = map - med_map

  unsharp_rat = map / med_map

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  unsharp_mask = unsharp_diff gt thresh

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LABEL REGIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  reg = label_region(unsharp_mask)
  n_reg = max(reg)
  print, "Found "+str(n_reg)+" features via unsharp masking."

  if n_reg gt 0 then begin
     
     unsharp_ra = fltarr(n_reg)*nan
     unsharp_dec = fltarr(n_reg)*nan
     unsharp_intens = fltarr(n_reg)*nan
     unsharp_km = fltarr(n_reg)*nan

     for ii = 0, n_reg-1 do begin
        mask = reg eq (ii+1)
        ind = where(mask, ct)
        if ct eq 0 then continue
        
        unsharp_intens[ii] = max(map[ind], /nan)
        unsharp_ra[ii] = total(ri[ind]*map[ind],/nan)/total(map[ind],/nan)
        unsharp_dec[ii] = total(di[ind]*map[ind],/nan)/total(map[ind],/nan)

     endfor
     
     if band eq 'w1' then begin
        coef = 52e3
     endif else if band eq 'w2' then begin
        coef = 28e3
     endif else if band eq 'w3' then begin
        coef = 4.8e3
     endif else if band eq 'w4' then begin
        coef = 1.5e3
     endif else if band eq 'nuv' then begin
        coef = 38.
     endif else if band eq 'fuv' then begin
        coef = 0.4
     endif else begin
        print, "Band not recognized, using band w1."
        coef = 52e3
     endelse

     unsharp_km = -1.*2.5*alog10(unsharp_intens/coef)
    
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if keyword_set(show) then begin

     !p.multi=0
     loadct, 0

;     if total(finite(map)) lt 10 then $ 
;        rms = 0.0 $
;     else $
;        rms = mad(map[where(finite(map))])
     rms = 0.1

     disp, map, ra, da, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5, color=cgcolor('white',255), /radec

     contour, unsharp_mask, ra, da, lev=[1], /overplot, color=cgcolor('red', 255)

     oplot, unsharp_ra, unsharp_dec, psym=cgsymcat('filledstar') $
            , color=cgcolor('orange'), symsize=1.5

     if keyword_set(pause) then ch = get_kbrd(1)

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if n_elements(outfile) gt 0 then begin
     sxaddpar, hdr, 'BUNIT', 'MASK'
     writefits, outfile, unsharp_mask, hdr
  endif

end
