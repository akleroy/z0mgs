pro v2_galex_build_atlas $
   , subsample = subsample $
   , just_galaxy = just_galaxy $
   , skip_galaxy = skip_galaxy $   
   , start_galaxy = start_galaxy $
   , stop_galaxy = stop_galaxy $
   , stage = do_stage $
   , galaxy_mask = do_galaxy_mask $
   , coord_mask = do_coord_images $
   , star_pred = do_star_pred $
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
  gaia_dir = '../../working_data/gaia/'+subsample+'/'
  
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
     this_pgc = tab[ii].pgc
     
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

           if strcompress(this_galaxy,/rem) eq 'PGC2557' then begin
              print, "... special case of M31"
              dat = gal_data(this_galaxy)
              print, "... size was: ", size_deg
              size_deg = 3.0*dat.r25_deg
              print, "... adjusted to: ", size_deg
           endif
           
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
; Make galaxy masks
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if keyword_set(do_galaxy_mask) then begin
        
        for vv = 0, n_bands-1 do begin

           this_band = bands[vv]

           print, "... making galaxy mask for ", this_band
           
           template_image_file = $
              staged_dir+this_galaxy+ $
              '_'+this_band+'_mjysr.fits'
           
           galaxy_mask_file = $
              mask_dir+this_galaxy+ $
              '_'+this_band+'_galmask.fits'
           
           v2_build_galaxy_mask $           
              , this_pgc=this_pgc $
              , all_gal_data=all_gal_data $
              , infile=template_image_file $
              , outfile=galaxy_mask_file $
              , show=show $
              , pause=pause

        endfor
        
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Make coord images
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Make radius images
     
     if keyword_set(do_coord_images) then begin

        
        
     endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Make star predictions
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
      
; Makes a predicted star flux file with a Jy value in the nearest
; pixel. Compared to the v1 approach this is a new intermediate data
; product, and uses GAIA DR3. Convolving and masking comes after.
     
     if keyword_set(do_star_pred) then begin

        for vv = 0, n_bands-1 do begin

           this_band = bands[vv]

           psf_dir = $
              '/export/bell-tycho/leroy.42/ellohess/kernels/PSF_FITS_Files/'           
           if this_band eq 'nuv' then $
              psf_file = psf_dir+'PSF_Corrected_GALEX_NUV_added_wing.fits'     
           if this_band eq 'fuv' then $
              psf_file = psf_dir+'PSF_Corrected_GALEX_FUV_added_wing.fits'     

           template_image_file = $
              staged_dir+this_galaxy+ $
              '_'+this_band+'_mjysr.fits'
           
           this_gaia_file = $
              gaia_dir+this_galaxy+ $
              '_gaia_dr3.fits'

           star_flux_file = $
              mask_dir+this_galaxy+ $
              '_'+this_band+'_starflux.fits'

           native_res_pred_file = $
              mask_dir+this_galaxy+ $
              '_'+this_band+'_starintens.fits'
           
           v2_build_star_pred $
              , ra_ctr=tab[ii].ra_ctr $
              , dec_ctr=tab[ii].dec_ctr $
              , band=this_band $
              , infile=template_image_file $
              , gaia_file=this_gaia_file $
              , pred_flux_file=star_flux_file $
              , psf_file=psf_file $
              , pred_image_file=native_res_pred_file $
              , ks_struct=ks_struct $
              , show=show $
              , pause=pause
           
        endfor
        
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Fit and subtract a background
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
     if keyword_set(do_bkgrd) then begin
        
     endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Convolve
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
     if keyword_set(do_convolve) then begin        
        
     endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Extract a subimage
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%       

     if keyword_set(do_extract) then begin
        
     endif
     
  endfor
  
end
