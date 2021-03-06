pro bkfit_unwise $
   , mapfile=infile $  
   , map=map $
   , hdr=hdr $
   , outfile=outfile $  
   , radfile=radfile $
   , fidrad=fidrad_in $
   , masklist=masklist $
   , band=band $
   , rejfile=rejfile $
   , rejected=rejected $
   , bksub=bksub $
   , aperture_scale=aperture_scale $
   , show=show $
   , pause=pause $
   , plane=plane
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

   if n_elements(band) eq 0 then $
      band = 'w1'
 
  if n_elements(map) eq 0 or n_elements(hdr) eq 0 then begin
     if file_test(infile) eq 0 then begin
        message, 'Target file not found and map and header not supplied.', /info
        return
     endif     
     map = readfits(infile, hdr, /silent)  
  endif

  mask = finite(map) eq 0
  n_masks = n_elements(masklist)
  for ii = 0, n_masks-1 do begin
     if file_test(masklist[ii]) eq 0 then continue
     this_mask = readfits(masklist[ii], this_mask_hdr, /silent)
     if total(size(this_mask,/dim) ne size(map,/dim)) gt 0 then begin
        message, "Mask wrong size. This shouldn't happen. Stopping.", /info
        stop
     endif
     this_mask = (this_mask mod 10) gt 0
     mask = mask or this_mask
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; TUNING PARAMETERS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(niter) eq 0 then $
     niter = 5
  if n_elements(thresh) eq 0 then $
     thresh = 2.0
  if n_elements(aperture_scale) eq 0 then $
     aperture_scale = 2.0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MASK OUT REGION NEAR THE GALAXY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  aperture_mask = mask*0B
  if n_elements(radfile) gt 0 then begin
     if file_test(radfile) then begin
        rgrid = readfits(radfile, rhdr, /silent)
        if total(size(rgrid,/dim) ne size(map,/dim)) gt 0 then begin
           message, "Radius grid wrong size. This shouldn't happen. Stopping.", /info
           stop
        endif
        if n_elements(fidrad_in) then begin
           fidrad = fidrad_in
        endif else begin
           fidrad = sxpar(rhdr, 'FIDRAD',missing=0.0)
        endelse
        aperture_mask = rgrid lt fidrad*aperture_scale
     endif
  endif
  mask = mask or aperture_mask
  mask_frac = (total(mask)*1. - total(aperture_mask)*1.) $
              /(n_elements(mask)*1. - total(aperture_mask)*1.0)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEAL WITH THE CASE WHERE ALL DATA ARE MASKED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if mask_frac gt 0.99 then begin
     message, 'More than 99% of the image appears masked. I will unsmask all but the aperture.', /info
     message, 'This galaxy should be flagged during delivery.', /info
     mask = aperture_mask
  endif

  if total(mask eq 0 and finite(map)) eq 0 then begin
     message, "There isn't enough room in the image to fit a background. Returning.", /info
     return
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; REMOVE A BACKGROUND
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  aperture = mask eq 0

  for ii = 0, niter-1 do begin

     bkind = where(aperture, ct)
     if ct eq 0 then begin
        message, 'No background pixels', /info        
        continue
     endif

     vec = map[bkind]
     rms = mad(vec)
     med = median(vec)

     bad_ind = where(abs(vec-med) gt thresh*rms, bad_ct)
     if bad_ct gt 0 then $
        aperture[bkind[bad_ind]] = 0B

  endfor 

  mask = aperture eq 0
  bksub = map - med

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLANE FIT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if band eq 'w3' or band eq 'w4' or keyword_set(plane) then begin

     sz = size(map)
     x = findgen(sz[1]) # (fltarr(sz[2])+1.0)
     y = (fltarr(sz[1])+1.0) #  findgen(sz[2])  
     
     for ii = 0, niter-1 do begin

        bkind = where(aperture, ct)
        if ct eq 0 then continue
        this_x = x[bkind]
        this_y = y[bkind]
        vec = map[bkind]
        
        coefs = planefit(this_x,this_y,vec,0., yfit)
        resid = vec - yfit
        rms = mad(resid)
        bad_ind = where(abs(resid) gt thresh*rms, bad_ct)
        if bad_ct gt 0 then $
           aperture[bkind[bad_ind]] = 0B

     endfor
     
     fit = coefs[0] + coefs[1]*x + coefs[2]*y
     bksub = map-fit

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RECENTER HISTOGRAM
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if rms lt 1d-6 or finite(rms) eq 0 then begin
     print, "Noise too low. Stopping."
  endif

  fit_ind = where(mask eq 0, fit_ct)     
  if abs(median(bksub[fit_ind])) gt 5.*rms then begin
     print, "Background subtraction has failed. Stopping."
     stop
  endif

  fit_ind = where(mask eq 0, fit_ct)     
  bins = bin_data(bksub[fit_ind], bksub[fit_ind]*0.0+1.0 $
                  , xmin=-5.*rms, xmax=5.*rms $
                  , binsize=0.05*rms, /nan)     
  hist = convol(bins.counts*1.0, psf_gaussian(npix=21,fwhm=7,ndim=1),/nan)  
  maxval = max(hist, maxind, /nan)
  bksub = bksub - bins[maxind].xmid

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; STATS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
 
  std = stddev(bksub[where(aperture_mask eq 0)],/nan)
  rms = mad(bksub[fit_ind])
  rej_frac = (total(mask)*1. - total(aperture_mask)*1.)/(n_elements(mask)*1. - total(aperture_mask)*1.)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(show) then begin

     !p.multi=[0,2,2]
     
     loadct, 0     
     disp, map, max=rms*3., min=-3.*rms, title=band
     contour, mask, /overplot, lev=[1], color=cgcolor('blue')
     
     loadct, 0
     disp, bksub, max=rms*3., min=-3.*rms, title=band
     contour, mask, /overplot, lev=[1], color=cgcolor('blue')

     fit_ind = where(mask eq 0, fit_ct)     
     bins = bin_data(map[fit_ind], map[fit_ind]*0.0+1.0 $
                     , xmin=-5.*rms, xmax=5.*rms $
                     , binsize=0.05*rms, /nan)     
     hist = convol(bins.counts*1.0, psf_gaussian(npix=21,fwhm=7,ndim=1),/nan)

     plot, bins.xmid, hist, ps=10, title=band
     oplot, 0.0*[1,1], [-1d6, 1d6], color=cgcolor('red')

     fit_ind = where(mask eq 0, fit_ct)     
     bins = bin_data(bksub[fit_ind], bksub[fit_ind]*0.0+1.0 $
                     , xmin=-5.*rms, xmax=5.*rms $
                     , binsize=0.05*rms, /nan)     
     hist = convol(bins.counts*1.0, psf_gaussian(npix=21,fwhm=7,ndim=1),/nan)
  
     plot, bins.xmid, hist, ps=10, title=band
     oplot, 0.0*[1,1], [-1d6, 1d6], color=cgcolor('red')

     if keyword_set(pause) then ch = get_kbrd(1)
     !p.multi=0
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if n_elements(outfile) gt 0 then begin
     hdr_copy = hdr
     sxaddpar, hdr_copy, 'RMS', rms
     sxaddpar, hdr_copy, 'STDDEV', std
     sxaddpar, hdr_copy, 'MASKFRAC', mask_frac
     sxaddpar, hdr_copy, 'REJFRAC', rej_frac 
     writefits, outfile, bksub, hdr_copy
  endif

  if n_elements(rejfile) gt 0 then begin
     hdr_copy = hdr
     sxaddpar, hdr_copy, 'RMS', rms
     sxaddpar, hdr_copy, 'STDDEV', std     
     sxaddpar, hdr_copy, 'MASKFRAC', mask_frac
     sxaddpar, hdr_copy, 'REJFRAC', rej_frac 
     writefits, rejfile, rejected, hdr_copy
  endif

end
