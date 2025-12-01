pro v2_unwise_build_atlas $
   , subsample = subsample $
   , just_galaxy = just_galaxy $
   , skip_galaxy = skip_galaxy $   
   , start_galaxy = start_galaxy $
   , stop_galaxy = stop_galaxy $
   , stage = do_stage $
   , galaxy_mask = do_galaxy_mask $
   , aperture_mask = do_aperture_mask $
   , star_pred = do_star_pred $
   , star_mask = do_star_mask $
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
  dat = gal_data(pgc=tab.pgc, /full)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  staged_dir = '../../working_data/unwise/staged/'+subsample+'/'
  bkgrd_dir = '../../working_data/unwise/bkgrd/'+subsample+'/'
  mask_dir = '../../working_data/unwise/masks/'+subsample+'/'
  convolved_dir = '../../working_data/unwise/convolved/'+subsample+'/'
  final_dir = '../../working_data/unwise/final/'+subsample+'/'
  gaia_dir = '../../working_data/gaia/'+subsample+'/'
  
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
; Make galaxy masks
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Use HyperLEDA's catalog to mask known galaxies in the field
     
     if keyword_set(do_galaxy_mask) then begin
        
        print, "... making galaxy mask"

        template_image_file = 'None'
        
        for vv = 0, n_bands-1 do begin
           
           this_band = bands[vv]

           test_image_file = $
              staged_dir+this_galaxy+ $
              '_'+this_band+'_mjysr.fits'

           if file_test(test_image_file) then $
              template_image_file = test_image_file
        endfor
        
        galaxy_mask_file = $
           mask_dir+this_galaxy+ $
           '_galmask.fits'
        
        v2_build_galaxy_mask $           
           , this_pgc=this_pgc $
           , all_gal_data=all_gal_data $
           , infile=template_image_file $
           , outfile=galaxy_mask_file $
           , show=show $
           , pause=pause
        
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Make aperture masks / radius images
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Make an image of the location of the galaxy
     
     if keyword_set(do_aperture_mask) then begin        

        print, "... making aperture mask"

        template_image_file = 'None'
        
        for vv = 0, n_bands-1 do begin
           
           this_band = bands[vv]

           test_image_file = $
              staged_dir+this_galaxy+ $
              '_'+this_band+'_mjysr.fits'

           if file_test(test_image_file) then $
              template_image_file = test_image_file
        endfor

        galaxy_mask_file = $
           mask_dir+this_galaxy+ $
           '_radius.fits'

        ra = dat[ii].ra_deg
        dec = dat[ii].dec_deg
        pa = dat[ii].posang_deg
        incl = dat[ii].incl_deg
        fid_rad = dat[ii].r25_deg
        
        v2_build_aperture_mask $
           , infile = template_image_file $
           , ra = ra $
           , dec = dec $
           , posang = pa $
           , incl = incl $
           , fid_rad = fid_rad $
           , outfile = galaxy_mask_file $ 
           , show = show $
           , pause = pause
        
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
           if this_band eq 'w1' then $
              psf_file = psf_dir+'PSF_Corrected_WISE_ATLAS_3.4_added_wing.fits'
           if this_band eq 'w2' then $
              psf_file = psf_dir+'PSF_Corrected_WISE_ATLAS_4.6_added_wing.fits'
           if this_band eq 'w3' then $
              psf_file = psf_dir+'PSF_Corrected_WISE_ATLAS_11.6_added_wing.fits'     
           if this_band eq 'w4' then $
              psf_file = psf_dir+'PSF_Corrected_WISE_ATLAS_22.1_added_wing.fits'      

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
; Turn the star predictions into masks
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if keyword_set(do_make_star_masks) then begin
 
        for vv = 0, n_bands-1 do begin

           this_band = bands[vv]

           template_image_file = $
              staged_dir+this_galaxy+ $
              '_'+this_band+'_mjysr.fits'
           
           native_res_pred_file = $
              mask_dir+this_galaxy+ $
              '_'+this_band+'_starintens.fits'

           star_mask_file = $
              mask_dir+this_galaxy+ $
              '_'+this_band+'_starmask.fits'

           
           v2_build_star_mask $
              , map = template_image_file $
              , pred = native_res_pred_file $
              , outfile = star_mask_file $
              , show = show $
              , pause = pause
           
        endfor           
        
     endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Mask and fit the background
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if keyword_set(do_bkgrd) then begin

        for vv = 0, n_bands-1 do begin

           this_band = bands[vv]

           infile = $
              staged_dir+this_galaxy+ $
              '_'+this_band+'_mjysr.fits'

           outfile = $
              staged_dir+this_galaxy+ $
              '_'+this_band+'_bksub.fits'
           
           star_mask_file = $
              mask_dir+this_galaxy+ $
              '_'+this_band+'_starmask.fits'

           gal_mask_file = $
              mask_dir+this_galaxy+ $
              '_'+this_band+'_galmask.fits'

           rad_file = $
              mask_dir+this_galaxy+ $
              '_'+this_band+'_radius.fits'
           
           v2_bkfit_unwise $
              , infile=infile $
              , band=this_band $              
              , outfile=outfile $  
              , radfile=rad_file $
              , masklist=[star_mask_file, gal_mask_file] $
              , show=show $
              , pause=pause $
              , plane=plane
           
        endfor           

        
     endif
          
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Make convolved masked and unmasked images
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if keyword_set(do_convolve) then begin

        for vv = 0, n_bands-1 do begin

           this_band = bands[vv]

           infile = $
              staged_dir+this_galaxy+ $
              '_'+this_band+'_bksub.fits'
           
           star_mask_file = $
              mask_dir+this_galaxy+ $
              '_'+this_band+'_starmask.fits'

           gal_mask_file = $
              mask_dir+this_galaxy+ $
              '_'+this_band+'_galmask.fits'
           
           
           
        endfor           
        
     endif
     
  endfor
     
end
