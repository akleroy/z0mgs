pro v2_unwise_build_atlas $
   , subsample = subsample $
   , just_galaxy = just_galaxy $
   , skip_galaxy = skip_galaxy $   
   , start_galaxy = start_galaxy $
   , stop_galaxy = stop_galaxy $
   , stage = do_stage $
   , gaia = do_gaia $
   , stack = do_stacks $
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
        v2_unwise_build_atlas $
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
; LOAD META DATA
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

  staged_dir = '../../working_data/unwise/staged/'+subsample+'/'
  bkgrd_dir = '../../working_data/unwise/bkgrd/'+subsample+'/'
  mask_dir = '../../working_data/unwise/masks/'+subsample+'/'
  convolved_dir = '../../working_data/unwise/convolved/'+subsample+'/'
  final_dir = '../../working_data/unwise/final/'+subsample+'/'

  bands = ['w1','w2','w3','w4']
  n_bands = n_elements(bands)
  
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
; COPY FILES, CHANGE UNITS, FLAG ON INVVAR TO STAGE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if keyword_set(do_stage) then begin

        w1_file = strcompress(tab[ii].w1_fname, /rem)
        
        for vv = 0, n_bands-1 do begin

           if vv eq 0 and tab[ii].use_w1 eq 0 then continue
           if vv eq 1 and tab[ii].use_w2 eq 0 then continue
           if vv eq 2 and tab[ii].use_w3 eq 0 then continue
           if vv eq 3 and tab[ii].use_w4 eq 0 then continue
           
           this_band = bands[vv]           
           
           infile_img = str_replace(w1_file,'w1', this_band)
           infile_invar = str_replace(infile_img,'-img-m.fits' $
                                      , '-invvar-m.fits')

           test = file_search(infile_img, count=ct)
           if ct eq 0 then begin
              ;printf, 1, pgc_name, band, infile
              print, 'NOT FOUND: ', this_galaxy, this_band, ' ', infile_img
              not_found_ct += 1
              stop
           endif

           outfile_img = staged_dir+this_galaxy+'_'+this_band+'_mjysr.fits'
           outfile_mask = mask_dir+this_galaxy+'_'+this_band+'_artifacts.fits'
           
           if keyword_set(incremental) then begin
              img_present = file_search(outfile_img, count=img_ct)
              mask_present = file_search(outfile_mask, count=mask_ct)              
              if img_present gt 0 and mask_present gt 0 then continue
           endif

           v2_build_invvar_mask $
              , infile=infile_invar $
              , outfile=outfile_mask $
              , show=show $
              , pause=pause
           
           v2_convert_unwise_to_mjysr $
              , infile = infile_img $
              , outfile = outfile_img $
              , band = this_band
           
        endfor

     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Query Gaia for stars
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if keyword_set(do_gaia) then begin

     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Stack stars
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if keyword_set(do_stack_stars) then begin

     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Make star masks
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if keyword_set(do_make_star_masks) then begin

     endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Make galaxy masks
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
     if keyword_set(do_galaxy_masks) then begin

     endif
     
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Mask and fit the background
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if keyword_set(do_bkgrd) then begin

     endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Additional star and user masking
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Make convolved masked and unmasked images
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  endfor
     
end
