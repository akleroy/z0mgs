pro mg_eyeball_sample $
   , just=just $
   , start=start $
   , stop=stop $
   , user=user $
   , high=high $
   , norej=norej $
   , big=big
      
  atlas_dir = '../delivery/'  
  index = mrdfits('../measurements/delivery_index.fits',1,h)
  if keyword_set(big) then begin
     s = gal_data(pgc=index.pgc)
     ind = where(s.r25_deg gt 100./3600., ct)
     index = index[ind]
  endif

  ind = sort(index.pgc)
  ind = reverse(ind)
  index = index[ind]

  pgc_list = index.pgc
  n_pgc = n_elements(index)
    
  bands = ['FUV','NUV','WISE1','WISE2','WISE3','WISE4']
  name = ['fuv','nuv','w1','w2','w3','w4']
  n_bands = n_elements(bands)

  if n_elements(start) eq 0 then start=0
  if n_elements(stop) eq 0 then stop=n_pgc-1

  if n_elements(user) eq 0 then user='akl'
  flag_file = '../measurements/flags_'+user+'.idl'
  if file_test(flag_file) eq 0 then begin
     print, 'Re-initializing.'
     flags = bytarr(n_elements(index))*0B     
     save, file=flag_file, pgc_list, flags
  endif else begin
     restore, flag_file, /v
  endelse

  for ii = start, stop-1 do begin

     if n_elements(just) gt 0 then $
        if total(just eq index[ii].pgc) eq 0 then $
           continue
       
     loadct, 0
     !p.multi=[0,3,2]

     mask_file = atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_mask.fits'
     if file_test(mask_file) eq 0 then begin
        print, 'PGC'+strcompress(str(index[ii].pgc),/rem)+' lacks a mask. Write this down!'
        continue
     endif
     mask = readfits(mask_file, mask_hdr)

     for jj = 0, n_bands-1 do begin

        if jj eq 0 then begin
           if index[ii].has_wise1 eq 0 then continue
           w1 = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w1.fits', w1_hdr)
           w1_rej = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w1_rejected.fits')
           rms = index[ii].rms_wise1
           tag = 'WISE1'
           map = w1
           rej = w1_rej           
        endif

        if jj eq 1 then begin
           if index[ii].has_wise2 eq 0 then continue
           w2 = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w2.fits', w2_hdr)
           w2_rej = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w2_rejected.fits')
           rms = index[ii].rms_wise1
           tag = 'WISE2'
           map = w2
           rej = w2_rej           
        endif

        if jj eq 2 then begin
           if index[ii].has_wise3 eq 0 then continue
           w3 = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w3.fits', w3_hdr)
           w3_rej = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w3_rejected.fits')
           rms = index[ii].rms_wise1
           tag = 'WISE3'
           map = w3
           rej = w3_rej           
        endif

        if jj eq 3 then begin
           if index[ii].has_wise4 eq 0 then continue
           w4 = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w4.fits', w4_hdr)
           w4_rej = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w4_rejected.fits')
           rms = index[ii].rms_wise4
           tag = 'WISE4'
           map = w4
           rej = w4_rej           
        endif

        if jj eq 4 then begin
           if index[ii].has_nuv eq 0 then continue
           nuv = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_nuv.fits', nuv_hdr)
           nuv_rej = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_nuv_rejected.fits')
           rms = index[ii].rms_nuv
           rms  = 0.5d-3
           tag = 'NUV'
           map = nuv
           rej = nuv_rej           
        endif

        if jj eq 5 then begin
           if index[ii].has_fuv eq 0 then continue
           fuv = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_fuv.fits', fuv_hdr)
           fuv_rej = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_fuv_rejected.fits')
           rms = index[ii].rms_fuv
           rms = 0.5d-3
           tag = 'FUV'
           map = fuv
           rej = fuv_rej           
        endif

        if keyword_set(norej) then begin
           dispmap = map*(rej eq 0)
        endif else begin
           dispmap = map
        endelse
        if keyword_set(high) then begin
           vec = map[sort(dispmap)] 
           vec = vec[where(finite(vec))]
           n  = n_elements(vec)
           disp, dispmap, /sq, min=vec[round(0.025*n)], max=vec[round(0.975*n)], reserve=5 $
                 , title=tag, color=cgcolor('white',254), /xstyle, /ystyle
        endif else begin
           disp, dispmap, /sq, min=-10.*rms, max=10.*rms, reserve=5 $
                 , title=tag, color=cgcolor('white',254), /xstyle, /ystyle
        endelse        
        contour, rej, lev=[0.25], color=cgcolor('red',255), /overplot
        contour, mask, lev=[10,100], color=cgcolor('blue',255), /overplot
     
     endfor

     !p.multi = 0
     
     print, ''
     print, 'PGC'+str(index[ii].pgc)+' galaxy '+str(ii)
     print, 'Currently flagged? ', flags[ii]
     print, '(F/f) to flag, (-) to go back, (s) to stop, other key to go forward.'
     print, ''
     ch = get_kbrd(1)
     if strupcase(ch) eq 'F' then flags[ii] = 1
     if ch eq '-' then ii = (ii-2) > (-1)
     if strupcase(ch) eq 'S' then begin
        stop
     endif
     save, file=flag_file, pgc_list, flags

  endfor

  save, file=flag_file, pgc_list, flags

end
