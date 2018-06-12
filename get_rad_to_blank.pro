function get_rad_to_blank $
   , mag=mag $
   , intens=intens $
   , band=band $
   , fwhm=fwhm $
   , use_error=use_error

  if n_elements(band) eq 0 then begin
     message, 'Defaulting to WISE1.', /info
     band = 'w1'
  endif

  if band eq 'w1' then begin
     coef = 52e3
     thresh = 10.^(-2.5)
  endif else if band eq 'w2' then begin
     coef = 28e3
     thresh = 10.^(-2.4)
  endif else if band eq 'w3' then begin
     coef = 4.8e3
     thresh = 10.^(-1.75)
  endif else if band eq 'w4' then begin
     coef = 1.5e3
     thresh = 10.^(-0.75)
  endif else if band eq 'nuv' then begin
     coef = 38.
     thresh = 10.^(-3.5)
  endif else if band eq 'fuv' then begin
     coef = 0.4
     thresh = 10.^(-3.5)
  endif else begin
     print, "Band not recognized, using band w1."
     coef = 52e3
  endelse

  if n_elements(mag) gt 0 then begin
     pred_intens = 10.^(-1.*mag/2.5)*coef
  endif

  if n_elements(intens) eq 0 then $
     intens = pred_intens

  missing = where(finite(intens) eq 0, missing_ct)
  if missing_ct gt 0 and n_elements(mag) gt 0 then $
     intens[missing] = pred_intens[missing]

  if n_elements(fwhm) eq 0 then begin
     fwhm = 15./3600.
  endif

  error_beam = 125./3600.

  rad_to_blank = fltarr(n_elements(intens))
  
  ind = where(intens gt thresh and intens/thresh le 1000., ct)
  if ct gt 0 then $
     rad_to_blank[ind] = ((fwhm/2.354)*sqrt(alog(intens/thresh)*2.0d))[ind]

  if keyword_set(use_error) then begin
     ind = where(intens/thresh ge 1000., ct)
     if ct gt 0 then $
        rad_to_blank[ind] = ((error_beam/2.354)*sqrt(alog(intens/1d3/thresh)*2.0d))[ind]
  endif

  return, rad_to_blank

end
