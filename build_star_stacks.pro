pro build_star_stacks $
   , band=band $
   , start=start $
   , stop=stop

  if n_elements(band) eq 0 then band='w1'

  tag = 1
  outfile = '../measurements/star_stacks/star_stack_'+band+'_'+ $
            strcompress(str(tag),/rem)+'.idl'
  print, band, outfile

  restore, '../measurements/2mass_stars.idl', /v
      
  unwise_dir = '../unwise/atlas/' 
  galex_dir = '../galex/atlas/' 

  if band eq 'w1' or band eq 'w2' or band eq 'w3' or band eq 'w4' then $
     dir = unwise_dir $
  else $
     dir = galex_dir

  index = mrdfits('../measurements/delivery_index.fits',1,h)
  pgc_list = index.pgc
  n_pgc = n_elements(index)
    
  if n_elements(start) eq 0 then start=0
  if n_elements(stop) eq 0 then stop=n_pgc-1

  for ii = start, stop-1 do begin

     print, ii, ' / ', stop

     pgc_string = strcompress(str(index[ii].pgc),/rem)

     mask_file = dir+'PGC'+pgc_string+'_mask.fits'
     if file_test(mask_file) eq 0 then begin
        print, 'PGC'+strcompress(str(index[ii].pgc),/rem)+' lacks a mask. Write this down!'
        continue
     endif
     mask = readfits(mask_file, mask_hdr)

     if band eq 'w1' and index[ii].has_wise1 eq 0 then continue
     if band eq 'w2' and index[ii].has_wise2 eq 0 then continue
     if band eq 'w3' and index[ii].has_wise3 eq 0 then continue
     if band eq 'w4' and index[ii].has_wise4 eq 0 then continue
     if band eq 'nuv' and index[ii].has_nuv eq 0 then continue
     if band eq 'fuv' and index[ii].has_fuv eq 0 then continue

     if band eq 'w1' or band eq 'w2' or band eq 'w3' or band eq 'w4' then $        
        map = readfits(dir+'PGC'+pgc_string+'_'+band+'_bksub_gauss15.fits', map_hdr)
     if band eq 'nuv' or band eq 'fuv' then $
        map = readfits(dir+'PGC'+pgc_string+'_'+band+'_gauss15_bksub.fits', map_hdr)
     
     adxy, map_hdr, star_ra, star_dec, star_x, star_y
     sz = size(map)
     in_im = where(star_x ge 50 and star_y ge 50 and $
                   star_x lt sz[1]-51 and star_y lt sz[2]-51 $
                   , star_ct)
     if star_ct eq 0 then continue
     
     for jj = 0, star_ct-1 do begin
        xlo = star_x[in_im[jj]]-50
        xhi = star_x[in_im[jj]]+50
        ylo = star_y[in_im[jj]]-50
        yhi = star_y[in_im[jj]]+50
        cutout = map[xlo:xhi, ylo:yhi]

        if ii mod 1000 eq 0 and ii ne 0 then begin
           save, file=outfile $
                 , stack, mag, val
           stack = [cutout]
           mag = [star_km[in_im[jj]]]
           val = [cutout[50,50]]        
           tag += 1
           outfile = '../measurements/star_stacks/star_stack_'+band+'_'+ $
                     strcompress(str(tag),/rem)+'.idl'           
        endif else if  n_elements(stack) eq 0 then begin
           stack = [cutout]
           mag = [star_km[in_im[jj]]]
           val = [cutout[50,50]]
        endif else begin
           stack = [[[stack]],[[cutout]]]
           mag = [mag, star_km[in_im[jj]]]
           val = [val, cutout[50,50]]
        endelse
     endfor 

     if ((ii mod 100) eq 0) and $
        ((ii mod 1000) ne 0) then begin
        disp, median(stack, dim=3), max=0.01
     endif
     
  endfor

  save, file=outfile $
        , stack, mag, val

  
end
