pro build_delivery $
   , wise=do_wise $
   , galex=do_galex $
   , reset=do_reset $
   , index=do_index
    
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
     , just=just $
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
     
     pgc_name = pgc_list[ii]
     this_dat = gal_data[ii]
         
     print, ''
     print, 'Packaging '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
     print, ''


;    WISE BANDS 1 TO 4 IMAGE AND REJECTED PIXEL
     for band = 1, 4 do begin
        
        if keyword_set(do_wise) eq 0 then continue

        map = readfits(wise_dir + pgc_name + '_w'+str(band)+'_gauss15.fits', hdr, /silent)
        sz = size(map)
        new_sz = [sz[1]/2, sz[2]/2]

        hrebin, map, hdr, out=new_sz
        writefits, out_dir+pgc_name + '_w'+str(band)+'.fits', map, hdr

        rejected = readfits(wise_dir + pgc_name + '_w'+str(band)+'_rejected_gauss15.fits', hdr, /silent)
        hrebin, rejected, hdr, out=new_sz
        writefits, out_dir+pgc_name + '_w'+str(band)+'_rejected.fits', map, hdr

     endfor

;    GALEX FUV AND NUV IMAGE, WEIGHT, AND REJECTED PIXEL
     for jj = 0, 1 do begin
        if jj eq 0 then band = 'fuv'
        if jj eq 1 then band = 'nuv'

        if keyword_set(do_galex) eq 0 then continue
                
        infile = galex_dir + pgc_name + '_'+band+'_extcorr.fits'

        test = file_search(infile, count=ct)
        if ct eq 0 then $
           continue

        map = readfits(infile, hdr, /silent)
        if sxpar(hdr, 'SKIP') eq 1 then $
           continue
        sz = size(map)
        new_sz = [sz[1]/2, sz[2]/2]

        hrebin, map, hdr, out=new_sz
        writefits, out_dir+pgc_name + '_'+band+'.fits', map, hdr

        rejected = readfits(galex_dir + pgc_name + '_'+band+'_rejected.fits', hdr, /silent)
        hrebin, rejected, hdr, out=new_sz
        writefits, out_dir+pgc_name + '_'+str(band)+'_rejected.fits', map, hdr        

        rejected = readfits(galex_dir + pgc_name + '_'+band+'_weight_gauss15_align.fits', hdr, /silent)
        hrebin, rejected, hdr, out=new_sz
        writefits, out_dir+pgc_name + '_'+str(band)+'_weight.fits', map, hdr        

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
