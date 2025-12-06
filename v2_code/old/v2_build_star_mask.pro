pro v2_build_star_mask $
   , map = image_file $
   , pred = pred_file $
   , outfile = outfile $
   , thresh = thresh $
   , use_rms = use_rms $   
   , show=show $
   , pause=pause

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Initialize
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if n_elements(pred_file) eq 0 then begin
     print, "prediction file required."
     return
  endif

  if n_elements(image_file) eq 0 then begin
     print, "image file required."
     return
  endif  
  
  if file_test(pred_file) eq 0 then begin
     print, pred_file, " not found."
     return
  endif

  if file_test(image_file) eq 0 then begin
     print, image_file, " not found."
     return
  endif
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Read the map and prediction and make the mask
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  image = readfits(infile, hdr_image)
  pred = readfits(pred_file, hdr_pred)

; set the scaling multiplied by the threshold used to mask the
; prediction mask. If use_rms is called then the rms of the image is
; calculated and multiplied by the threshold
  
  scaling = 1.0
  if keyword_set(use_rms) then begin

     fin_ind = where(finite(image), fin_ct)
     if fin_ct gt 10 then begin        
        scaling = mad(map[fin_ind])     
     endelse
     
  endif

; The mask is just where the prediction is greater than the thresh
; times scaling
  
  mask = pred ge thresh*scaling

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Write to disk
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if n_elements(outfile) gt 0 then begin

     
     
  endif
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Show the results  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(show) then begin
     
     !p.multi=[0,2,2]

     intens_map = readfits(native_res_pred_file, intens_hdr)
     
     if total(finite(map)) lt 10 then $ 
        rms = 0.0 $
     else $
        rms = mad(map[where(finite(map))])
     
     loadct, 0
     disp, map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, intens_map, lev=[1.*rms], color=cgcolor('red'), /overplot
     
     loadct, 0
     disp, intens_map, max=5*rms, min=-5.*rms $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, intens_map, lev=[1.*rms], color=cgcolor('red'), /overplot
     
     loadct, 0
     disp, map-intens_map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, intens_map, lev=[1.*rms], color=cgcolor('red'), /overplot
     

     if keyword_set(pause) then begin
        ch = get_kbrd(1)
     endif
     
     !p.multi=0
     
  endif

  
end
