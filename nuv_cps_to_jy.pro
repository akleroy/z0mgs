function nuv_cps_to_jy $
   , cps

;+
;
; Counts per second to jy transformation for GALEX FUV band. Based on 
;
; http://galexgi.gsfc.nasa.gov/docs/galex/FAQ/counts_background.html
;
; mAB = -2.5 log10 CPS + 20.08
;
; and mAB = -2.5 log10 (f_nu / 3630.8)
;
;-

  return, cps*3.37289e-05

end
