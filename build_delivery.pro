pro build_delivery $
   , wise=do_wise $
   , galex=do_galex $
   , masks=do_masks $
   , reset=do_reset $
   , index=do_index $
   , just=just $
   , only=only $
   , tag = tag $
   , start = start_num $
   , stop = stop_num
    
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  galex_dir = '../galex/atlas/'
  wise_dir = '../unwise/atlas/'
  mask_dir = '../masks/'
  out_dir = '../delivery/'

  build_galaxy_list $
     , in_dir = in_dir $
     , tag=tag $
     , just=only $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = gal_data $
     , start = start_num $
     , stop = stop_num $
     , exclude = ['PGC17223']
  n_pgc = n_elements(pgc_list)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WIPE EVERYTHING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_reset) then begin

     print, 'Hit y to really wiped the delivery and rebuild.'
     chr = get_kbrd(1)
     if chr eq 'y' then begin
        spawn, 'rm -rf '+out_dir
        spawn, 'mkdir '+out_dir
     endif

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER GALAXIES, REBIN TO ~5.5" PIXELS, WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = 0, n_pgc-1 do begin
     
     pgc_name = strcompress(pgc_list[ii], /rem)
     this_dat = gal_data[ii]
         
     if n_elements(just) gt 0 then $
        if total(pgc_name eq just) eq 0 then $
           continue
     
     print, ''
     print, 'Packaging '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
     print, ''

;    WISE BANDS 1 TO 4 IMAGE
     for jj = 0, 3 do begin
        
        if jj eq 0 then band = 'w1'
        if jj eq 1 then band = 'w2'
        if jj eq 2 then band = 'w3'
        if jj eq 3 then band = 'w4'

        if keyword_set(do_wise) eq 0 then continue
        
        for mm = 0, 1 do begin
           
           if mm eq 0  then begin
              res_str = 'gauss15'
           endif
           
           if mm eq 1 then begin
              if band eq 'w4' then continue
              res_str = 'gauss7p5'
           endif

           infile = wise_dir + pgc_name+'_'+band+'_'+res_str+'_bksub.fits'           
           outfile = out_dir + pgc_name+'_'+band+'_'+res_str+'.fits'           
           if file_test(infile) then begin
              map = readfits(infile, hdr, /silent)
              writefits, outfile, map, hdr
           endif else begin
              print, "Missing ... ", infile
           endelse

        endfor

     endfor

;    GALEX FUV AND NUV IMAGE AND WEIGHT
     for jj = 0, 1 do begin

        if keyword_set(do_galex) eq 0 then continue

        if jj eq 0 then band = 'fuv'
        if jj eq 1 then band = 'nuv'

        for mm = 0, 1 do begin
           
           if mm eq 0  then begin
              res_str = 'gauss15'
           endif
           
           if mm eq 1 then begin
              res_str = 'gauss7p5'
           endif

           infile = galex_dir + pgc_name+'_'+band+'_'+res_str+'_extcorr.fits'           
           outfile = out_dir + pgc_name+'_'+band+'_'+res_str+'.fits'           
           if file_test(infile) then begin
              hdr = headfits(infile)
              if sxpar(hdr, 'SKIP') eq 1 then $
                 continue
              map = readfits(infile, hdr, /silent)
              writefits, outfile, map, hdr
           endif           


           infile = galex_dir + pgc_name+'_'+band+'_weight_'+res_str+'_small.fits'
           outfile = out_dir + pgc_name+'_'+band+'_'+res_str+'_weight.fits'           
           if file_test(infile) then begin
              map = readfits(infile, hdr, /silent)
              writefits, outfile, map, hdr
           endif
           
        endfor
        
     endfor
     
;    MASKS
     for mm = 0, 1 do begin
        
        if keyword_set(do_masks) eq 0 then continue
        
        if mm eq 0  then begin
           res_str = 'gauss15'
        endif
        
        if mm eq 1 then begin
           res_str = 'gauss7p5'
        endif
        
        infile = mask_dir + pgc_name+'_'+res_str+'_rgrid.fits'
        outfile = out_dir + pgc_name+'_'+res_str+'_rgrid.fits'
        if file_test(infile) then begin
           map = readfits(infile, hdr, /silent)
           writefits, outfile, map, hdr
        endif

        infile = mask_dir + pgc_name+'_'+res_str+'_galaxies.fits'
        outfile = out_dir + pgc_name+'_'+res_str+'_galaxies.fits'
        if file_test(infile) then begin
           map = readfits(infile, hdr, /silent)
           writefits, outfile, map, hdr
        endif

        for jj = 0, 5 do begin

           if jj eq 0 then band = 'fuv'
           if jj eq 1 then band = 'nuv'
           if jj eq 2 then band = 'w1'
           if jj eq 3 then band = 'w2'
           if jj eq 4 then band = 'w3'
           if jj eq 5 then begin
              band = 'w4'
              if res_str eq 'gauss7p5' then continue
           endif

           infile = mask_dir + pgc_name+'_'+band+'_'+res_str+'_bright_stars.fits'           
           outfile = out_dir + pgc_name+'_'+band+'_'+res_str+'_bright_stars.fits'            
          if file_test(infile) then begin
              map = readfits(infile, hdr, /silent)
              writefits, outfile, map, hdr
           endif
                 
           infile = mask_dir + pgc_name+'_'+band+'_'+res_str+'_found_stars.fits'           
           outfile = out_dir + pgc_name+'_'+band+'_'+res_str+'_found_stars.fits'           
           if file_test(infile) then begin
              map = readfits(infile, hdr, /silent)
              writefits, outfile, map, hdr
           endif
                      
           infile = mask_dir + pgc_name+'_'+band+'_'+res_str+'_custom.fits'           
           outfile = out_dir + pgc_name+'_'+band+'_'+res_str+'_custom_mask.fits'           
           if file_test(infile) then begin
              map = readfits(infile, hdr, /silent)
              writefits, outfile, map, hdr
           endif

        endfor
        
     endfor

  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INDEX THE ATLAS CONTENTS INTO A TABLE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_index) then begin

     empty_index = index_one_galaxy(/empty)
     index = replicate(empty_index, n_pgc)

     for jj = 0, 1 do begin
        
        if jj eq 0 then res_str = 'gauss7p5'
        if jj eq 1 then res_str = 'gauss15'

        for ii = 0, n_pgc-1 do begin
           
           pgc_name = pgc_list[ii]
           this_dat = gal_data[ii]
           
           print, ''
           print, 'Indexing '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
           print, ''
           
           this_index = $
              index_one_galaxy( $
              pgc_name=pgc_name $
              , gal_data=this_dat $
              , atlas_dir= out_dir $
              , res_str=res_str)
           
           index[ii] = this_index
           
        endfor

        mwrfits, index, '../measurements/delivery_index_'+res_str+'.fits', /create

     endfor

  endif

end
