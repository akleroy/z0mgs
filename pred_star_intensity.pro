function pred_star_intensity $
   , mag, gaia=gaia, band=band, res=res
  
  if n_elements(res) eq 0 then $
     res = 'gauss7p5'

  coef = !values.f_nan

  if keyword_set(gaia) then begin

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GAIA G MAG
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if res eq 'gauss15' then begin
        print, "Resolution 15as not implemented for GAIA."
        return, !values.f_nan
     endif
     
     if band eq 'w1' then begin
        coef = 829964.d
     endif

     if band eq 'w2' then begin
        coef = 450745.d
     endif

     if band eq 'w3' then begin
        coef = 77349.7d
     endif

     if band eq 'nuv' then begin
        coef = 13067.6d
     endif

     if band eq 'fuv' then begin
        coef = 27.065d
     endif
     
  endif else begin

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; 2MASS KS MAG
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if band eq 'w1' then begin
        if res eq 'gauss7p5' then begin           
           coef = 164.9d3
        endif else if res eq 'gauss15' then begin
           coef = 52d3
        endif else begin
           return, !values._fnan
        endelse
     endif

     if band eq 'w2' then begin
        if res eq 'gauss7p5' then begin           
           coef = 91.8d3
        endif else if res eq 'gauss15' then begin
           coef = 27.4d3
        endif else begin
           return, !values._fnan
        endelse
     endif

     if band eq 'w3' then begin
        if res eq 'gauss7p5' then begin           
           coef = 16.8d3
        endif else if res eq 'gauss15' then begin
           coef = 4.76d3
        endif else begin
           return, !values._fnan
        endelse
     endif
     
     if band eq 'w4' then begin
        if res eq 'gauss7p5' then begin           
           coef = !values.f_nan
        endif else if res eq 'gauss15' then begin
           coef = 1.55d3
        endif else begin
           return, !values._fnan
        endelse
     endif

     if band eq 'nuv' then begin
        if res eq 'gauss7p5' then begin           
           coef = 202.d
        endif else if res eq 'gauss15' then begin
           coef = 69.7d
        endif else begin
           return, !values._fnan
        endelse
     endif
     
     if band eq 'fuv' then begin
        if res eq 'gauss7p5' then begin           
           coef = 1.35d
        endif else if res eq 'gauss15' then begin
           coef = 0.6d
        endif else begin
           return, !values._fnan
        endelse
     endif

  endelse

  intens = 10.^(-1.0d*mag/2.5)*coef

  return, intens

end
