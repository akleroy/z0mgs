pro v2_galex_build_atlas $
   , subsample = subsample $
   , just_galaxy = just_galaxy $
   , skip_galaxy = skip_galaxy $   
   , start_galaxy = start_galaxy $
   , stop_galaxy = stop_galaxy $
   , stage = do_stage $
   , bkgrd = do_bkgrd $ 
   , convol = do_convol $
   , show = show $
   , pause = pause $
   , incremental = incremental

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RECURSIVE WRAPPER  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if subsample eq 'all' then begin
     subsample_list = ['localgroup','localvolume','largeleda' $
                       , 'smallleda', 'manga']
     for ii = 0, n_elements(subsample_list)-1 do begin
        this_subsample = subsample_list[ii]
        v2_galex_build_atlas $
           , subsample = this_subsample $
           , just_galaxy = just_galaxy $
           , skip_galaxy = skip_galaxy $   
           , start_galaxy = start_galaxy $
           , stop_galaxy = stop_galaxy $
           , stage = do_stage $
           , bkgrd = do_bkgrd $ 
           , convol = do_convol $
           , show = show $
           , pause = pause $
           , incremental = incremental
     endfor
     return
  endif 
    
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOAD UNWISE META DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  tdir = '../../measurements/'
  
; Default  
  if n_elements(subsample) eq 0 then begin
     subsample = 'largeleda'
  endif

  if subsample eq 'localgroup' then begin
     tab = mrdfits(tdir+'unwise_v2_index_localgroup.fits',1,tab_hdr)
  endif
  
  if subsample eq 'localvolume' then begin
     tab = mrdfits(tdir+'unwise_v2_index_localvolume.fits',1,tab_hdr)
  endif
  
  if subsample eq 'largeleda' then begin
     tab = mrdfits(tdir+'unwise_v2_index_largeleda.fits',1,tab_hdr)
  endif

  if subsample eq 'smallleda' then begin
     tab = mrdfits(tdir+'unwise_v2_index_smallleda.fits',1,tab_hdr)
  endif

  if subsample eq 'manga' then begin
     tab = mrdfits(tdir+'unwise_v2_index_manga.fits',1,tab_hdr)
  endif

  n_gal = n_elements(tab)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  staged_dir = '../../working_data/galex/staged/'+subsample+'/'
  bkgrd_dir = '../../working_data/galex/bkgrd/'+subsample+'/'
  mask_dir = '../../working_data/galex/masks/'+subsample+'/'
  convolved_dir = '../../working_data/galex/convolved/'+subsample+'/'
  final_dir = '../../working_data/galex/final/'+subsample+'/'

  bands = ['nuv','fuv']
  n_bands = n_elements(bands)
  
  galex_tile_index = mrdfits('galex_index_file.fits',1,tih)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER TARGETS  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  first = 1B
  last = 0B
  not_found_ct = 0L
  
  for ii = 0, n_gal-1 do begin

     this_galaxy = strcompress(tab[ii].z0mgs_name,/rem)
     
     if n_elements(just_galaxy) gt 0 then begin
        if total(just_galaxy eq this_galaxy) eq 0 then begin
           continue
        endif
     endif

     if n_elements(skip_galaxy) gt 0 then begin
        if total(skip_galaxy eq this_galaxy) gt 0 then begin
           continue
        endif
     endif

     if n_elements(start_galaxy) gt 0 then begin
        if first eq 1B then begin
           if total(start_galaxy eq this_galaxy) eq 0 then begin
              continue
           endif else begin
              first = 0B
           endelse
        endif
     endif

     if n_elements(stop_galaxy) gt 0 then begin
        if last eq 1B then begin
           continue
        endif
        if last eq 0B then begin
           if total(stop_galaxy eq this_galaxy) eq 0 then begin
              continue
           endif else begin
              last = 1B
           endelse
        endif
     endif

     print, ''
     print, 'Processing galaxy '+str(ii)+' / '+str(n_gal)+' ... '+this_galaxy
     print, ''
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CREATE CUTOUTS FROM TILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if keyword_set(do_stage) then begin

        w1_file = strcompress(tab[ii].w1_fname, /rem)
        
        for vv = 0, n_bands-1 do begin

           ;if vv eq 0 and tab[ii].use_fuv eq 0 then continue
           ;if vv eq 1 and tab[ii].use_nuv eq 0 then continue
           
           this_band = bands[vv]           

           outfile_image = staged_dir+this_galaxy+ $
                         '_'+this_band+'_mjysr.fits'
           outfile_weight = staged_dir+this_galaxy+ $
                            '_'+this_band+'_weight.fits'
           
           if keyword_set(incremental) then begin
              image_present = file_test(outfile_image)
              weight_present = file_test(outfile_weight)
              if image_present gt 0 and $
                 weight_present gt 0 then $
                 continue
           endif

           ra_ctr = tab[ii].ra_ctr
           dec_ctr = tab[ii].dec_ctr
           size_deg = (tab[ii].trc_dec - tab[ii].blc_dec)
           
           v2_extract_galex_stamp $
              , fuv=(this_band eq 'fuv') $
              , ra_ctr=ra_ctr $
              , dec_ctr=dec_ctr $
              , size_deg=size_deg $
              , index=galex_tile_index $
              , /useint $
              , show=show $
              , pause=pause $
              , image=image $
              , weight=weight $
              , hdr=hdr

           writefits, outfile_image, image, hdr
           writefits, outfile_weight, weight, hdr
           
        endfor

     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Mask and fit the background
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Additional star and user masking
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Make convolved masked and unmasked images
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  endfor
     
end
