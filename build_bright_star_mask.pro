pro build_bright_star_mask $
   , pgcname=pgcname $
   , galdata=this_data $
   , band=band $
   , res=res $
   , infile=infile $
   , predfile=predfile $
   , outfile=outfile $
   , show=show $
   , pause=pause $
   , thresh=thresh $
   , sat=sat_thresh $
   , pad=pad

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS AND READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  nan = !values.f_nan

  gaia_dir = '../stars/gaia/'

  if n_elements(res) eq 0 then  begin
     res = 'gauss7p5'
  endif
  
  valid_res = ['gauss7p5', 'gauss15']
  if total(res eq valid_res) eq 0 then begin
     print, "Defaulting to 7.5as Gaussian"
     res = 'gauss7p5'
  endif

  valid_bands = ['w1','w2','w3','w4','nuv','fuv']
  if total(band eq valid_bands) eq 0 then begin
     print, "Defaulting to WISE1"
     band = 'w1'
  endif

  if file_test(infile) eq 0 then begin
     message, 'Target file not found and map and header not supplied.', /info
     return
  endif     
  map = readfits(infile, hdr, /silent)  

  if file_test(predfile) eq 0 then begin
     message, 'Prediction file not found and is required.', /info
     return
  endif     
  pred = readfits(predfile, predhdr, /silent)  

; check that they have the same astrometry?
    
  error_beam = (sxpar(hdr, 'BMAJ')*3600. > 7.5/3600.)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; THRESHOLDS / TUNING PARAMETERS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(thresh) eq 0 then begin

     if res eq 'gauss7p5' then begin
        if band eq 'w1' then thresh = 0.5d-2
        if band eq 'w2' then thresh = 0.5d-2
        if band eq 'w3' then thresh = 1d-2
        if band eq 'w4' then thresh = !values.f_nan
        if band eq 'nuv' then thresh = 5d-4
        if band eq 'fuv' then thresh = 5d-4
     endif else begin
        if band eq 'w1' then thresh = 1d-2
        if band eq 'w2' then thresh = 1d-2
        if band eq 'w3' then thresh = 1d-2
        if band eq 'w4' then thresh = 5d-1
        if band eq 'nuv' then thresh = 5d-4
        if band eq 'fuv' then thresh = 5d-4
     endelse
     
  endif

  if n_elements(sat_thresh) eq 0 then begin

     if band eq 'w1' then sat_thresh = 1d2
     if band eq 'w2' then sat_thresh = 1d2
     if band eq 'w3' then sat_thresh = 3d2
     if band eq 'w4' then sat_thresh = 3d2
     if band eq 'nuv' then sat_thresh = 1d6
     if band eq 'fuv' then sat_thresh = 1d6

     if band eq 'w1' then sat_edge = 1d1
     if band eq 'w2' then sat_edge = 1d1
     if band eq 'w3' then sat_edge = 3d1
     if band eq 'w4' then sat_edge = 3d1
     if band eq 'nuv' then sat_edge = 1d6
     if band eq 'fuv' then sat_edge = 1d6

  endif

  if n_elements(pad) eq 0 then $
     pad = 2

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE MASK AND WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  star_mask = pred ge thresh

  sat_mask = pred gt sat_thresh
  if total(sat_mask) gt 0 then begin
     print, "... likely saturation detected."

;    This is the most dangerous step because it interacts with the image ...
;     sat_mask = grow_mask(sat_mask, constraint=map gt sat_edge)
     conv_with_gauss $
        , data=pred $
        , hdr=predhdr $
        , target_beam=error_beam*[1,1,0] $
        , out_data = smoothpred $
        , /quiet

;     sat_mask = grow_mask(sat_mask, iters = 10)
     sat_mask = grow_mask(sat_mask, constraint=smoothpred gt thresh)
     star_mask = star_mask or sat_mask     
     ;pause = 0B
     ;show = 0B
  endif ;else begin
     ;pause = 0B
     ;show = 0B
     
;  endelse

  supersat_mask = pred gt sat_thresh*10.
  if total(supersat_mask) gt 0 then begin
     print, "... likely SUPER saturated stars in image."

     conv_with_gauss $
        , data=pred*sat_mask $
        , hdr=predhdr $
        , target_beam=5.*error_beam*[1,1,0] $
        , out_data = smoothpred $
        , /quiet

     supersat_mask = grow_mask(supersat_mask, constraint=smoothpred gt thresh)
     star_mask = star_mask or supersat_mask     
     ;pause = 1B
     ;show = 1B
  endif

  if pad gt 0 then $
     star_mask = grow_mask(star_mask, iters=pad)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if keyword_set(show) then begin

     !p.multi=[0,2,2]

     if total(finite(map)) lt 10 then $ 
        rms =0.0 $
     else $
        rms = mad(map[where(finite(map))])

     loadct, 0
     disp, map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, star_mask, lev=[1,11], /overplot, color=cgcolor('blue', 255)

     loadct, 0
     disp, map, max=0.25, min=-0.01 $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, star_mask, lev=[1,11], /overplot, color=cgcolor('blue', 255)

     loadct, 0
     disp, map*(star_mask ne 1) $
           , max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5

     loadct, 0
     disp, map*(star_mask ne 1) $
           , max=0.25, min=-0.01 $
           , /sq, xstyle=1, ystyle=1, reserve=5

     if keyword_set(pause) then begin
        ch = get_kbrd(1)
     endif

     !p.multi=0

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if n_elements(outfile) gt 0 then begin
     sxaddpar, hdr, 'BUNIT', 'MASK'
     writefits, outfile, star_mask, hdr
  endif

end
