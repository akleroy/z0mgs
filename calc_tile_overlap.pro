function calc_tile_overlap $
   , ra_ctr=ra_ctr $
   , dec_ctr=dec_ctr $
   , pad=pad $
   , min_ra=min_ra $
   , max_ra=max_ra $
   , min_dec=min_dec $
   , max_dec=max_dec

  overlap = ((min_dec-pad) lt dec_ctr) and $
            ((max_dec+pad) gt dec_ctr)
  
; TRAP HIGH LATITUDE CASE AND (I GUESS) TOSS BACK ALL TILES. DO BETTER LATER.
  mean_dec = (min_dec+max_dec)*0.5    
  if abs(dec_ctr)+pad gt 88. then begin
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
