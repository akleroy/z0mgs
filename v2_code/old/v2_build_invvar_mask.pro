pro v2_build_invvar_mask $
   , infile=infile_invvar $
   , outfile=outfile_mask $
   , show=show $
   , pause=pause

  radius = 3
  
  invvar = readfits(infile_invvar, hdr, /silent)

  rms = stddev(invvar,/nan)
  med = median(invvar)
  resid = (invvar-med)/rms
  
  mask = grow_mask(abs(resid) gt 5., radius=radius)

  if keyword_set(show) then begin
     loadct, 3
     disp, resid, xs=1, ys=1
     contour, mask, lev=[1], /overplot, c_color=cgcolor('green')
  endif

  if keyword_set(pause) then begin
     print, "Key to continue."
     ch = get_kbrd(1)
  endif

  writefits, outfile_mask, mask, hdr
  
end
