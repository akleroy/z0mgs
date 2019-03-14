function lookup_mtol $
   , nuvw1 = nuvw1 $
   , w3w1 = w3w1 $   
   , lumw1 = lumw1 $
   , ssfr = ssfr $
   , unc = unc $
   , dir = dir

  if n_elements(dir) eq 0 then $
     dir = '/data/kant/0/leroy.42/allsky/measurements/'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FIGURE OUT WHAT WE'RE DOING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  has_nuvw1 = n_elements(nuvw1) gt 0
  has_w3w1 = n_elements(w3w1) gt 0
  has_lumw1 = n_elements(lumw1) gt 0
  has_ssfr = n_elements(ssfr) gt 0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WE HAVE SPECIFIC STAR FORMATION RATE - JUST USE THIS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if has_ssfr then begin
     pivot1 = -10.9
     pivot2 = -9.5
     zero_pt = 0.485
     delta = -0.285
     mtol = zero_pt + $
            (ssfr gt pivot1 and ssfr le pivot2)*(ssfr-pivot1)/(pivot2-pivot1)*(delta) + $
            (ssfr gt pivot2)*ssfr
     unc = !values.f_nan*mtol
     return, mtol
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SINGLE COLOR LOOKUP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; WE ONLY HAVE NUV/W1

  if has_nuvw1 and (has_w3w1 eq 0) and (has_lumw1 eq 0) then begin     
     pivot1 = -1.9
     pivot2 = -0.4
     zero_pt = 0.485
     delta = -0.285
     mtol = zero_pt + $
            (nuvw1 gt pivot1 and nuvw1 le pivot2)*(nuvw1-pivot1)/(pivot2-pivot1)*(delta) + $
            (nuvw1 gt pivot2)*delta
     unc = !values.f_nan*mtol
     return, mtol
  endif

; WE ONLY HAVE W3/W1

  if has_w3w1 and (has_nuvw1 eq 0) and (has_lumw1 eq 0) then begin     
     pivot1 = 0.1
     pivot2 = 0.85
     zero_pt = 0.485
     delta = -0.285
     mtol = zero_pt + $
            (w3w1 gt pivot1 and w3w1 le pivot2)*(w3w1-pivot1)/(pivot2-pivot1)*(delta) + $
            (w3w1 gt pivot2)*delta
     unc = !values.f_nan*mtol
     return, mtol
  endif

; WE ONLY HAVE LUMW1

  if has_lumw1 and (has_w3w1 eq 0) and (has_nuvw1 eq 0)then begin     
     mtol = 0.5 + 0.0*lumw1
     unc = !values.f_nan*mtol
     return, mtol
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; TWO DIMENSIONAL GRIDS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if has_w3w1 and has_lumw1 and (has_nuvw1 eq 0) then begin

     gridfile = dir+'w3w1_lumw1_mtol_grid.fits'
     rmsfile = dir+'w3w1_lumw1_madmtol_grid.fits'
     if file_test(gridfile) eq 0 or file_test(rmsfile) eq 0 then begin
        print, gridfile, ' or ', rmsfile, ' not found. Returning'
        return, !values.f_nan
     endif

     grid = $
        readfits(gridfile, hdr)
     rms = $
        readfits(rmsfile, rms_hdr)

     sz = size(grid)
     xaxis = findgen(sz[1])
     yaxis = findgen(sz[2])

     w3w1_axis = ((findgen(sz[1])+1.)-sxpar(hdr,'CRPIX1'))*sxpar(hdr, 'CDELT1')+sxpar(hdr,'CRVAL1')
     lumw1_axis = ((findgen(sz[2])+1.)-sxpar(hdr,'CRPIX2'))*sxpar(hdr, 'CDELT2')+sxpar(hdr,'CRVAL2')

     copy_lumw1 = (lumw1 < max(lumw1_axis)) > min(lumw1_axis)

     x = interpol(xaxis, w3w1_axis, w3w1)
     y = interpol(yaxis, lumw1_axis, copy_lumw1)

     mtol = !values.f_nan * w3w1
     unc = !values.f_nan * w3w1

     ind = where(x ge 0 and y ge 0  $
                 and x lt sz[1] and y lt sz[2], ct)    
     if ct eq 0 then return, mtol
     mtol[ind] = interpolate(grid, x[ind], y[ind])
     unc[ind] = interpolate(rms, x[ind], y[ind])     

     return, mtol
     
  endif

  if has_w3w1 and has_nuvw1 and (has_lumw1 eq 0) then begin

     gridfile = dir+'nuvw1_w3w1_mtol_grid.fits'
     rmsfile = dir+'nuvw1_w3w1_madmtol_grid.fits'
     if file_test(gridfile) eq 0 or file_test(rmsfile) eq 0 then begin
        print, gridfile, ' or ', rmsfile, ' not found. Returning'
        return, !values.f_nan
     endif

     grid = $
        readfits(gridfile, hdr)
     rms = $
        readfits(rmsfile, rms_hdr)

     sz = size(grid)
     xaxis = findgen(sz[1])
     yaxis = findgen(sz[2])

     nuvw1_axis = ((findgen(sz[1])+1.)-sxpar(hdr,'CRPIX1'))*sxpar(hdr, 'CDELT1')+sxpar(hdr,'CRVAL1')
     w3w1_axis = ((findgen(sz[2])+1.)-sxpar(hdr,'CRPIX2'))*sxpar(hdr, 'CDELT2')+sxpar(hdr,'CRVAL2')

     x = interpol(xaxis, nuvw1_axis, nuvw1)
     y = interpol(yaxis, w3w1_axis, w3w1)

     mtol = !values.f_nan * w3w1
     unc = !values.f_nan * w3w1

     ind = where(x ge 0 and y ge 0  $
                 and x lt sz[1] and y lt sz[2], ct)    
     if ct eq 0 then return, mtol
     mtol[ind] = interpolate(grid, x[ind], y[ind])
     unc[ind] = interpolate(rms, x[ind], y[ind])     

     return, mtol

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOKUP IN THE CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if has_w3w1 and has_nuvw1 and has_lumw1 then begin

     gridfile = dir+'mtol_cube.fits'
     rmsfile = dir+'madmtol_cube.fits'
     if file_test(gridfile) eq 0 or file_test(rmsfile) eq 0 then begin
        print, gridfile, ' or ', rmsfile, ' not found. Returning'
        return, !values.f_nan
     endif

     grid = $
        readfits(gridfile, hdr)
     rms = $
        readfits(rmsfile, rms_hdr)

     sz = size(grid)
     xaxis = findgen(sz[1])
     yaxis = findgen(sz[2])
     zaxis = findgen(sz[3])

     nuvw1_axis = ((findgen(sz[1])+1.)-sxpar(hdr,'CRPIX1'))*sxpar(hdr, 'CDELT1')+sxpar(hdr,'CRVAL1')
     w3w1_axis = ((findgen(sz[2])+1.)-sxpar(hdr,'CRPIX2'))*sxpar(hdr, 'CDELT2')+sxpar(hdr,'CRVAL2')
     lumw1_axis = ((findgen(sz[3])+1.)-sxpar(hdr,'CRPIX3'))*sxpar(hdr, 'CDELT3')+sxpar(hdr,'CRVAL3')

     copy_lumw1 = (lumw1 < max(lumw1_axis)) > min(lumw1_axis)

     x = interpol(xaxis, nuvw1_axis, nuvw1)
     y = interpol(yaxis, w3w1_axis, w3w1)
     z = interpol(zaxis, lumw1_axis, copy_lumw1)

     mtol = !values.f_nan*x
     unc = !values.f_nan*x

     ind = where(x ge 0 and y ge 0 and z ge 0 $
                 and x lt sz[1] and y lt sz[2] and z lt sz[3], ct)    
     if ct eq 0 then return, mtol
     mtol[ind] = interpolate(grid, x[ind], y[ind], z[ind])
     unc[ind] = interpolate(rms, x[ind], y[ind], z[ind])

  endif

  return, mtol

end
