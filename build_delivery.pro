pro build_delivery $
   , wise=do_wise $
   , galex=do_galex $
   , reset=do_reset $
   , index=do_index $
   , just=just $
   , only=only
    
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  galex_dir = '../galex/atlas/'
  wise_dir = '../unwise/atlas/'
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
     
     pgc_name = strcompresS(pgc_list[ii], /rem)
     this_dat = gal_data[ii]
         
     if n_elements(just) gt 0 then $
        if total(pgc_name eq just) eq 0 then $
           continue
     
     print, ''
     print, 'Packaging '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
     print, ''


;    WISE BANDS 1 TO 4 IMAGE AND REJECTED PIXEL
     for band = 1, 4 do begin
        
        if keyword_set(do_wise) eq 0 then continue

        map_fname = wise_dir + pgc_name + '_w'+str(band)+'_bksub_gauss15.fits'
        if file_test(map_fname) eq 0 then begin
           message, "No map for "+map_fname, /info
           continue
        endif
        mask_fname = wise_dir + pgc_name + '_mask.fits'
        rej_fname = wise_dir + pgc_name + '_w'+str(band)+'_rejected_gauss15.fits'

        map = readfits(map_fname, hdr, /silent)
        mask = readfits(mask_fname, mask_hdr, /silent)
        rejected = readfits(rej_fname, rejected_hdr, /silent)        

        sz = size(mask)        
        ind = where(mask ge 10, ct)
        ind_to_xyv, ind, sz=sz, x=x, y=y

        hextract, map, hdr, min(x), max(x), min(y), max(y)        
        sz = size(map)
        new_sz = [sz[1]/2, sz[2]/2]
        hrebin, map, hdr, out=new_sz

        ;hastrom, mask, mask_hdr, hdr, interp=0
        ;hastrom, rejected, rejected_hdr, hdr $
        ;         , interp=2, cubic=-0.5, missing=!values.f_nan

        ;rejected = grow_mask(rejected, iters=1)

        writefits, out_dir+pgc_name + '_w'+str(band)+'.fits', map, hdr
        ;writefits, out_dir+pgc_name + '_mask.fits', mask, hdr
        ;writefits, out_dir+pgc_name + '_w'+str(band)+'_rejected.fits', rejected, hdr

     endfor

;    GALEX FUV AND NUV IMAGE, WEIGHT, AND REJECTED PIXEL
     for jj = 0, 1 do begin
        if jj eq 0 then band = 'fuv'
        if jj eq 1 then band = 'nuv'

        if keyword_set(do_galex) eq 0 then continue
                
        map_fname = galex_dir + pgc_name + '_'+band+'_extcorr.fits'
        if file_test(map_fname) eq 0 then continue
        ;mask_fname = galex_dir + pgc_name + '_mask.fits'
        ;rej_fname = galex_dir + pgc_name + '_'+band+'_gauss15_rejected.fits'
        wt_fname = galex_dir + pgc_name + '_'+band+'_weight_gauss15_align.fits'

        map = readfits(map_fname, hdr, /silent)
        if sxpar(hdr, 'SKIP') eq 1 then $
           continue
        mask = readfits(mask_fname, mask_hdr, /silent)
        ;rejected = readfits(rej_fname, rejected_hdr, /silent)
        weight = readfits(wt_fname, weight_hdr, /silent)

        ind = where(mask ge 10, ct)
        sz = size(mask)
        ind_to_xyv, ind, sz=sz, x=x, y=y
        hextract, map, hdr, min(x), max(x), min(y), max(y)        
        sz = size(map)
        new_sz = [sz[1]/2, sz[2]/2]
        hrebin, map, hdr, out=new_sz

        ;hastrom, mask, mask_hdr, hdr, interp=0
        hastrom, weight, weight_hdr, hdr $
                 , interp=2, cubic=-0.5, missing=!values.f_nan
        ;hastrom, rejected, rejected_hdr, hdr $
        ;         , interp=2, cubic=-0.5, missing=!values.f_nan

        ;rejected = grow_mask(rejected, iters=1)

        writefits, out_dir+pgc_name + '_'+band+'.fits', map, hdr
        ;writefits, out_dir+pgc_name + '_'+band+'_rejected.fits', rejected, hdr
        writefits, out_dir+pgc_name + '_'+band+'_weight.fits', weight, weight_hdr

     endfor
     
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INDEX THE ATLAS CONTENTS INTO A TABLE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_index) then begin

     empty_index = index_one_galaxy(/empty)
     index = replicate(empty_index, n_pgc)

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
           , atlas_dir= out_dir)
     
        index[ii] = this_index
               
     endfor

     mwrfits, index, '../measurements/delivery_index.fits', /create

  endif

end
