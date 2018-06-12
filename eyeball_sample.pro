pro eyeball_sample $
   , just=just $
   , start=start $
   , stop=stop $
   , user=user $
   , big=big

;  restore, '../measurements/2mass_stars.idl', /v
      
  atlas_dir = '../delivery/'  
  index = mrdfits('../measurements/delivery_index.fits',1,h)
  s = gal_data(pgc=index.pgc)
  if keyword_set(big) then begin
     ind = where(s.r25_deg gt 100./3600., ct)
     index = index[ind]
     s = s[ind]
  endif

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
       
     !p.multi=[0,3,2]

     for jj = 0, n_bands-1 do begin

        if jj eq 0 then begin
           if index[ii].has_wise1 eq 0 then continue
           w1 = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w1.fits', w1_hdr, /silent)
           w2 = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w2.fits', w2_hdr, /silent)
           rms = index[ii].rms_wise1
           tag = 'WISE1'
           band = 'w1'
           map = w1
           hdr = w1_hdr
           minval = -1d-1
           maxval = 1d-1       
        endif

        if jj eq 1 then begin
           if index[ii].has_wise2 eq 0 then continue
           w2 = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w2.fits', w2_hdr, /silent)
           rms = index[ii].rms_wise1
           tag = 'WISE2'
           band = 'w2'
           map = w2
           hdr = w2_hdr
           minval = -1d-1
           maxval = 1d-1      
        endif

        if jj eq 2 then begin
           if index[ii].has_wise3 eq 0 then continue
           w3 = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w3.fits', w3_hdr, /silent)
           rms = index[ii].rms_wise1
           tag = 'WISE3'
           band = 'w3'
           map = w3
           hdr = w3_hdr
           minval = -1d-1
           maxval = 1d-1    
        endif

        if jj eq 3 then begin
           if index[ii].has_wise4 eq 0 then continue
           w4 = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_w4.fits', w4_hdr, /silent)
           rms = index[ii].rms_wise4
           tag = 'WISE4'
           band = 'w4'
           map = w4
           hdr = w4_hdr
           minval = -1d0
           maxval = 1d0     
        endif

        if jj eq 4 then begin
           if index[ii].has_nuv eq 0 then continue
           nuv = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_nuv.fits', nuv_hdr, /silent)
           rms = index[ii].rms_nuv
           tag = 'NUV'
           band = 'nuv'
           map = nuv
           hdr = nuv_hdr
           minval = -5d-3
           maxval = 5d-3     
        endif

        if jj eq 5 then begin
           if index[ii].has_fuv eq 0 then continue
           fuv = readfits(atlas_dir+'PGC'+strcompress(str(index[ii].pgc),/rem)+'_fuv.fits', fuv_hdr, /silent)
           rms = index[ii].rms_fuv
           tag = 'FUV'
           band = 'fuv'
           map = fuv
           hdr = fuv_hdr
           minval = -5d-3
           maxval = 5d-3
        endif
        
        loadct, 0
        dispmap = map
        
        if jj eq 0 then begin
           build_unsharp_mask $
              , pgc=index[ii].pgc $
              , galdata=s[ii] $
              , wise1=w1 $
              , hdr=hdr $
              , mask=unsharp_mask $
              , star_ra=unsharp_ra $
              , star_dec=unsharp_dec $
              , star_intens=unsharp_intens           
        endif

        if 1 eq 0 then begin
           build_mask, map=dispmap, hdr=hdr, band=band $
                       , pgc = index[ii].pgc $
                       , galdata=s[ii] $
                       , mask=mask $
                       , /do_flag_gals, /do_flag_stars, /do_flag_sharp $
                       , allgaldata=allgals $
                       , star_ra=star_ra, star_dec=star_dec, star_km=star_km
        endif

        loadct, 0
        disp, dispmap, /sq, min=minval, max=maxval, reserve=5 $
              , title=tag, color=cgcolor('white',254), /xstyle, /ystyle        

        contour, unsharp_mask, lev=[1], /overplot, color=cgcolor('red')
        ;contour, mask eq 10, lev=[1], /overplot, color=cgcolor('blue')
        ;contour, (mask mod 10) gt 0, lev=[1], /overplot, color=cgcolor('red')
 
        ;contour, mask, lev=[10,100], color=cgcolor('blue',255), /overplot
        ;if star_ct gt 0 then begin
        ;   x = star_x[in_im] 
        ;   y = star_y[in_im]
        ;   m = star_km[in_im]
           
        ;   oplot, x, y, color=cgcolor('red') $
        ;          , symsize=2, psym=cgsymcat('filledstar')
        ;   print, m
        ;   ind = where(m lt 8., ct)
        ;   if ct gt 0 then $
        ;      oplot, x[ind], y[ind], color=cgcolor('blue') $
        ;             , symsize=2, psym=cgsymcat('filledstar')
        ;endif

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
     
     ;loadct, 0
     ;disp, w1, /sq, min=-0.1, max=0.5, /xstyle, /ystyle
     ;test = w1*w2*nuv
     ;mask = test gt 10.
     ;mask = grow_mask(mask, iters=10)
     ;contour, mask, /overplot, lev=[1], color=cgcolor('red')
     ;ch = get_kbrd(1)

  endfor

  save, file=flag_file, pgc_list, flags

end
