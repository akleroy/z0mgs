function fuv_cps_to_jy $
   , cps

;+
;
; Counts per second to jy transformation for GALEX FUV band. Based on 
;
; http://galexgi.gsfc.nasa.gov/docs/galex/FAQ/counts_background.html
;
; mAB = -2.5 log10 CPS + 18.82
;
; and mAB = -2.5 log10 (f_nu / 3630.8)
;
;-

  return, cps*0.000107647

end
