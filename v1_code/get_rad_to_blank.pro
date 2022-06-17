function get_rad_to_blank $
   , intens=intens $
   , fwhm=fwhm $
   , thresh=thresh

  if n_elements(fwhm) eq 0 then begin
     fwhm = 7.5/3600.
  endif

  rad_to_blank = intens*0.0
  ind = where(intens gt thresh, ct)
  if ct gt 0 then $
     rad_to_blank[ind] = (((fwhm/2.354)*sqrt(alog(intens/thresh)*2.0d)))[ind]

  return, rad_to_blank

end
