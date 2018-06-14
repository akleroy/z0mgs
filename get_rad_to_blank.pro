function get_rad_to_blank $
   , mag=mag $
   , intens=intens $
   , band=band $
   , fwhm=fwhm $
   , error_rad=error_rad $
   , minpeak=minpeak $
   , thresh=thresh

  if n_elements(fwhm) eq 0 then begin
     fwhm = 15./3600.
  endif

  if n_elements(band) eq 0 then begin
     message, 'Defaulting to WISE1.', /info
     band = 'w1'
  endif

  if band eq 'w1' then begin
     coef = 52e3
     if fwhm eq 7.5/3600. then begin
        coef = 168e3
     endif
     if n_elements(thresh) eq 0 then begin
        thresh = 0.01
     endif
     if n_elements(minpeak) eq 0 then begin
        minpeak = 0.1
     endif
  endif else if band eq 'w2' then begin
     coef = 28e3
     if fwhm eq 7.5/3600. then begin
        coef = 27.7e3
     endif
     if n_elements(thresh) eq 0 then begin
        thresh = 0.01
     endif
     if n_elements(minpeak) eq 0 then begin
        minpeak = 0.1
     endif
  endif else if band eq 'w3' then begin
     coef = 4.8e3
     if fwhm eq 7.5/3600. then begin
        coef = 16.2e3
     endif
     if n_elements(thresh) eq 0 then begin
        thresh = 0.05
     endif
     if n_elements(minpeak) eq 0 then begin
        minpeak = 0.5
     endif
  endif else if band eq 'w4' then begin
     coef = 1.5e3
     if n_elements(thresh) eq 0 then begin
        thresh = 0.5
     endif
     if n_elements(minpeak) eq 0 then begin
        minpeak = 2.0
     endif
  endif else if band eq 'nuv' then begin
     coef = 38.
     if fwhm eq 7.5/3600. then begin
        coef = 139.
     endif
     if n_elements(thresh) eq 0 then begin
        thresh = 0.001
     endif
     if n_elements(minpeak) eq 0 then begin
        minpeak = 0.01
     endif
  endif else if band eq 'fuv' then begin
     coef = 0.4
     if fwhm eq 7.5/3600. then begin
        coef = 1.23
     endif
     if n_elements(thresh) eq 0 then begin
        thresh = 0.001
     endif
     if n_elements(minpeak) eq 0 then begin
        minpeak = 0.01
     endif
  endif else begin
     print, "Band not recognized, using band w1."
     coef = 52e3
     if n_elements(thresh) eq 0 then begin
        thresh = 0.01
     endif
     if n_elements(minpeak) eq 0 then begin
        minpeak = 0.1
     endif
  endelse

  if n_elements(mag) gt 0 then begin
     pred_intens = 10.^(-1.*mag/2.5)*coef
  endif

  if n_elements(intens) eq 0 then $
     intens = pred_intens

  missing = where(finite(intens) eq 0, missing_ct)
  if missing_ct gt 0 and n_elements(mag) gt 0 then $
     intens[missing] = pred_intens[missing]

  rad_to_blank = ((fwhm/2.354)*sqrt(alog(intens/thresh)*2.0d))
  too_low = where(intens le minpeak, low_ct)
  if low_ct gt 0 then $
     rad_to_blank[too_low] = !values.f_nan

  error_beam = 125./3600.
  error_rad = rad_to_blank*!values.f_nan
  ind = where(intens/thresh ge 1000., ct)
  error_rad[ind] = ((error_beam/2.354)*sqrt(alog(intens/1d3/thresh)*2.0d))[ind]

  return, rad_to_blank

end
