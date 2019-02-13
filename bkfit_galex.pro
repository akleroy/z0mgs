pro bkfit_galex $
   , mapfile=infile $  
   , map=map $
   , hdr=hdr $
   , wtfile=wtfile $
   , outfile=outfile $  
   , radfile=radfile $
   , fidrad=fidrad_in $
   , masklist=masklist $
   , band=band $
   , rejfile=rejfile $
   , rejected=rejected $
   , bkgrd=bkgrd_out $   
   , aperture_scale=aperture_scale $
   , plane=plane $
   , show=show $
   , pause=pause
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  noise_floor = 5e-5

   if n_elements(band) eq 0 then $
      band = 'nuv'
 
  if n_elements(map) eq 0 or n_elements(hdr) eq 0 then begin
     if file_test(infile) eq 0 then begin
        message, 'Target file not found and map and header not supplied.', /info
        return
     endif     
     map = readfits(infile, hdr, /silent)  
  endif

  wt = finite(map)*1.0
  if n_elements(wtfile) ne 0 then begin
     if file_test(wtfile) eq 0 then begin
        message, 'Target weight file not found and map and header not supplied.', /info
        return
     endif     
     wt = readfits(wtfile, wt_hdr, /silent)  
  endif

  mask = finite(map) eq 0
  n_masks = n_elements(masklist)
  for ii = 0, n_masks-1 do begin
     if file_test(masklist[ii]) eq 0 then continue
     this_mask = readfits(masklist[ii], this_mask_hdr, /silent)
     this_mask = (this_mask mod 10) gt 0
     mask = mask or this_mask
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE A BACKGROUND
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  sz = size(map)
  bkgrd = 0.0*fltarr(sz[1], sz[2])

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
; CREATE A VECTOR TO NORMALIZE FOR EXPOSURE TIME
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  noiselike = 1./sqrt(wt/median(wt)) 

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MASK OUT REGION NEAR THE GALAXY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  aperture_mask = mask*0B
  if n_elements(radfile) gt 0 then begin
     if file_test(radfile) then begin
        rgrid = readfits(radfile, rhdr, /silent)
        if total(size(rgrid, /dim) ne size(map, /dim)) gt 0 then begin
           print, "Mismatched size. Returning."
           return
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
     wtvec = noiselike[bkind]
     rms = mad(vec) > noise_floor
     med = median(vec)
     bad_ind = where(abs(vec-med) gt thresh*rms/wtvec, bad_ct)
     if bad_ct gt 0 then $
        aperture[bkind[bad_ind]] = 0B    

  endfor 

  mask = aperture eq 0
  bkgrd = med
  bksub = map - bkgrd

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLANE FIT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  fit_a_plane = 0B
  if keyword_set(plane) and mask_frac lt 0.75 then begin

     x = findgen(sz[1]) # (fltarr(sz[2])+1.0)
     y = (fltarr(sz[1])+1.0) #  findgen(sz[2])  
     
     for ii = 0, niter-1 do begin

        bkind = where(aperture, ct)
        if ct eq 0 then continue
        this_x = x[bkind]
        this_y = y[bkind]
        vec = map[bkind]        
        wtvec = noiselike[bkind]

        coefs = planefit(this_x,this_y,vec,0., yfit)
        resid = vec - yfit
        rms = mad(resid) > noise_floor
        bad_ind = where(abs(resid) gt thresh*rms*wtvec, bad_ct)
        if bad_ct gt 0 then $
           aperture[bkind[bad_ind]] = 0B

     endfor
     
     bkgrd = coefs[0] + coefs[1]*x + coefs[2]*y
     fit_a_plane = 1B
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RECENTER HISTOGRAM
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if rms lt 1d-6 or finite(rms) eq 0 then begin
     print, "Noise too low. Stopping."
     stop
  endif

  fit_ind = where(mask eq 0, fit_ct)
  if abs(median((map - bkgrd)[fit_ind])) gt 5.*rms then begin
     print, "Background subtraction has failed. Stopping."
     stop
  endif

  bins = bin_data((map-bkgrd)[fit_ind], (map-bkgrd)[fit_ind]*0.0+1.0 $
                  , xmin=-5.*rms, xmax=5.*rms $
                  , binsize=0.05*rms, /nan)
  hist = convol(bins.counts*1.0, psf_gaussian(npix=21,fwhm=7,ndim=1),/nan)  
  maxval = max(hist, maxind, /nan)
  bkgrd = bkgrd + bins[maxind].xmid

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; STATS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
 
  std = stddev((map - bkgrd)[where(aperture_mask eq 0)],/nan)
  rms = mad((map - bkgrd)[fit_ind])
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
     disp, (map-bkgrd), max=rms*3., min=-3.*rms, title=band
     contour, mask, /overplot, lev=[1], color=cgcolor('blue')
     
     fit_ind = where(mask eq 0, fit_ct)          
     bins = bin_data(map[fit_ind], map[fit_ind]*0.0+1.0 $
                     , xmin=-5.*rms, xmax=5.*rms $
                     , binsize=0.05*rms, /nan)     
     hist = convol(bins.counts*1.0, psf_gaussian(npix=21,fwhm=7,ndim=1),/nan)

     plot, bins.xmid, hist, ps=10, title=band
     oplot, 0.0*[1,1], [-1d6, 1d6], color=cgcolor('red')

     fit_ind = where(mask eq 0, fit_ct)     
     bins = bin_data((map-bkgrd)[fit_ind], (map-bkgrd)[fit_ind]*0.0+1.0 $
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
     if finite(rms) eq 0 then rms = -1.
     sxaddpar, hdr_copy, 'RMS', rms
     if finite(std) eq 0 then std = -1.
     sxaddpar, hdr_copy, 'STDDEV', std
     sxaddpar, hdr_copy, 'MASKFRAC', mask_frac
     sxaddpar, hdr_copy, 'REJFRAC', rej_frac 
     sxaddpar, hdr_copy, 'FITPLANE', fit_a_plane
     writefits, outfile, bkgrd, hdr_copy
  endif

  if n_elements(rejfile) gt 0 then begin
     hdr_copy = hdr
     if finite(rms) eq 0 then rms = -1.
     sxaddpar, hdr_copy, 'RMS', rms
     if finite(std) eq 0 then std = -1.
     sxaddpar, hdr_copy, 'STDDEV', std     
     sxaddpar, hdr_copy, 'MASKFRAC', mask_frac
     sxaddpar, hdr_copy, 'REJFRAC', rej_frac 
     writefits, rejfile, rejected, hdr_copy
  endif

end
