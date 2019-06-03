function lookup_mtol $
   , nuvw1 = nuvw1 $
   , fuvw1 = fuvw1 $  
   , w3w1 = w3w1 $  
   , w4w1 = w4w1 $
   , ssfrlike = ssfrlike $
   , truessfr =truessfr $
   , tag = tag $
   , dir = dir

  low_bound = 0.2
  high_bound = 0.5

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DIRECTORY FOR GRIDS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(dir) eq 0 then $
     dir = '/data/kant/0/leroy.42/allsky/'
  grid_dir = dir+'measurements/'
  fit_dir = dir+'z0mgs/'

  readcol $
     , fit_dir+'mtol_fits.txt', format='A,F,F,F,F', comment='#' $
     , fit_tag, mtolfit_lo, mtolfit_loval, mtolfit_hi, mtolfit_hival
  fit_tag = strcompress(fit_tag, /rem)
  n_tags = n_elements(fit_tag)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FIGURE OUT WHAT WE'RE DOING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  has_nuvw1 = n_elements(nuvw1) gt 0
  has_fuvw1 = n_elements(fuvw1) gt 0
  has_w3w1 = n_elements(w3w1) gt 0
  has_w4w1 = n_elements(w4w1) gt 0
  has_truessfr = n_elements(truessfr) gt 0
  has_ssfrlike = n_elements(ssfrlike) gt 0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WE HAVE SPECIFIC STAR FORMATION RATE - JUST USE THIS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if has_truessfr then begin
     
     tag_ind = (where(fit_tag eq 'ssfr'))[0]
     xin = truessfr

     pivot1 = mtolfit_lo[tag_ind]
     pivot2 = mtolfit_hi[tag_ind]
     zero_pt = mtolfit_loval[tag_ind]
     delta = mtolfit_hival[tag_ind] - mtolfit_loval[tag_ind]
     
     mtol = zero_pt + $
            (xin gt pivot1 and xin le pivot2)*(xin-pivot1)/(pivot2-pivot1)*(delta) + $
            (xin gt pivot2)*delta
     unc = !values.f_nan*mtol

     return, mtol
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; USE A ONE DIMENSIONAL COLOR
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  found_color = 0B

  if has_ssfrlike then begin     
     tag_ind = (where(fit_tag eq 'ssfrlike'))[0]
     xin = ssfrlike
     found_color = 1B
  endif else if has_w4w1 then begin
     tag_ind = (where(fit_tag eq 'w4w1'))[0]
     xin = w4w1
     found_color = 1B
  endif else if has_fuvw1 then begin
     tag_ind = (where(fit_tag eq 'fuvw1'))[0]
     xin = fuvw1
     found_color = 1B
  endif else if has_nuvw1 then begin
     tag_ind = (where(fit_tag eq 'nuvw1'))[0]
     xin = nuvw1
     found_color = 1B
  endif else if has_w3w1 then begin
     tag_ind = (where(fit_tag eq 'w3w1'))[0]
     xin = w3w1
     found_color = 1B
  endif

  if found_color then begin
     pivot1 = mtolfit_lo[tag_ind]
     pivot2 = mtolfit_hi[tag_ind]
     zero_pt = mtolfit_loval[tag_ind]
     delta = mtolfit_hival[tag_ind] - mtolfit_loval[tag_ind]
     
     mtol = zero_pt + $
            (xin gt pivot1 and xin le pivot2)*(xin-pivot1)/(pivot2-pivot1)*(delta) + $
            (xin gt pivot2)*delta
     unc = !values.f_nan*mtol
     
     return, mtol
  endif  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WE DID NOT FIND ANYTHING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  return, 0.35d

end
