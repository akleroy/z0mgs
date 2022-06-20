function v2_calc_tile_overlap $
   , ra_ctr=ra_ctr $
   , dec_ctr=dec_ctr $
   , pad=pad $
   , max_ra=max_ra $
   , min_ra=min_ra $
   , max_dec=max_dec $
   , min_dec=min_dec

  mean_dec = (min_dec+max_dec)*0.5
  
; Find all tiles in the declination range
  overlap = ((min_dec-pad) lt dec_ctr) and $
            ((max_dec+pad) gt dec_ctr)
  
; Trap high/low declination case. NB as far as I can tell only one
; tile lies at |dec| > 85 so don't stress about this case.

  if dec_ctr+pad gt 85. then begin
     overlap = mean_dec gt 85.
     return, overlap
  endif

  if dec_ctr-pad lt 85. then begin
     overlap = mean_dec lt 85.
     return, overlap
  endif
  
  ra_pad = pad / cos(!dtor*mean_dec)

; MERIDIAN CASES
  merid = where(max_ra lt min_ra)
  overlap[merid] = overlap[merid] and $
                   ((((min_ra-ra_pad) lt ra_ctr) or $
                     ((max_ra+ra_pad) gt ra_ctr)))[merid]

; BORING CASE
  normal = where(max_ra gt min_ra)  
  overlap[normal] = overlap[normal] and $
                    ((((min_ra-ra_pad) lt ra_ctr) and $
                      ((max_ra+ra_pad) gt ra_ctr)))[normal]
  return, overlap

end
