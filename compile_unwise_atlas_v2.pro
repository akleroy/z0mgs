pro compile_unwise_atlas_v2 $
   , units = do_units $
   , convol = do_convol $
   , extract = do_extract $
   , bksub = do_bksub $
   , show = show $
   , only = only $
   , pause = pause $
   , just = just $
   , tag = tag $
   , start = start_num $
   , stop = stop_num $
   , incremental = incremental

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  mask_dir = '../masks/'
  out_dir = '../unwise/atlas/'

  build_galaxy_list $
     , in_dir = in_dir $
     , tag=tag $
     , just=only $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = gal_data $
     , start = start_num $
     , stop = stop_num $
     , exclude = ['PGC17223','PGC89980']
  
  n_pgc = n_elements(pgc_list)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY FILES AND CHANGE UNITS TO STAGE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_units) then begin
     
     not_found_ct = 0L
     openw, 1, 'filesnotfound.txt'

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Unit conversion for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        for band = 1, 4 do begin
           
           infile = in_dir+pgc_name+'/'+'unwise-'+pgc_name+'-w'+str(band)+'-img-m.fits'

           test = file_search(infile, count=ct)
           if ct eq 0 then begin
              printf, 1, pgc_name, band, infile
              print, 'NOT FOUND: ', pgc_name, band, ' ', infile
              not_found_ct += 1
              continue
           endif

           outfile = out_dir+pgc_name+'_w'+str(band)+'_mjysr.fits'

           if keyword_set(incremental) then begin
              test = file_search(outfile, count=test_ct)
              if test_ct gt 0 then continue
           endif

           convert_unwise_to_mjysr $
              , infile = infile $
              , outfile = outfile $
              , band = band
           
        endfor

     endfor

     close, 1

     print, "Total files not found "+str(not_found_ct)

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RUN THE CONVOLUTIONS TO MATCH BEAMS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  !p.multi=[0,2,1]

  if keyword_set(do_convol) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Convolutions for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        for band = 1, 4 do begin

           if keyword_set(incremental) then begin
              outfile = out_dir+pgc_name+'_w'+str(band)+'_gauss15.fits'
              test_1 = file_search(outfile, count=test_ct_1)

              if band le 3 then begin                 
                 outfile = out_dir+pgc_name+'_w'+str(band)+'_gauss7p5.fits'
                 test_2 = file_search(outfile, count=test_ct_2)
                 found = test_ct_1 gt 0 and test_ct_2 gt 0
              endif else begin
                 found = test_ct_1 gt 0
              endelse

              if found then continue
           endif

           infile = out_dir+pgc_name+'_w'+str(band)+'_mjysr.fits'

           test = file_search(infile, count=ct)
           if ct eq 0 then begin
              message, 'File not found '+infile, /info
              continue
           endif

           test = readfits(infile, hdr)
           if test[0] eq -1 then begin
              message, "Problematic FITS file. Skipping.", /info
              continue
           endif

           outfile = out_dir+pgc_name+'_w'+str(band)+'_gauss15.fits'

           conv_z0mg_galaxy, $
              infile=infile, $
              start_psf='w'+str(band), $
              end_psf='g15', $
              outfile=outfile

           map = readfits(outfile, hdr)
           sxaddpar, hdr, 'BMAJ', 15./3600.
           sxaddpar, hdr, 'BMIN', 15./3600.
           sxaddpar, hdr, 'BMPA', 0.0
           writefits, outfile, map, hdr

           if band le 3 then begin
              outfile = out_dir+pgc_name+'_w'+str(band)+'_gauss7p5.fits'

              conv_z0mg_galaxy, $
                 infile=infile, $
                 start_psf='w'+str(band), $
                 end_psf='g7p5', $
                 outfile=outfile

              map = readfits(outfile, hdr)
              sxaddpar, hdr, 'BMAJ', 7.5/3600.
              sxaddpar, hdr, 'BMIN', 7.5/3600.
              sxaddpar, hdr, 'BMPA', 0.0
              writefits, outfile, map, hdr

           endif
           
        endfor
        
     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; EXTRACT A SMALLER FOOTPRINT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_extract) then begin

     readcol, 'custom_mask_specs.txt', format='L,F,F,F', comment='#' $
                 , override_pgc_list, override_rad_list, override_pa_list, override_incl_list

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Extract smaller footprint for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''
        
        for band = 1, 4 do begin

           for mm = 0, 1 do begin
          
              if mm eq 0 and band eq 4 then continue

              if mm eq 0 and band lt 4 then begin
                 res_str = 'gauss7p5'
                 infile = out_dir+pgc_name+'_w'+str(band)+'_gauss7p5.fits'
                 outfile = out_dir+pgc_name+'_w'+str(band)+'_gauss7p5_small.fits'
                 do_rebin = 0B
              endif else begin
                 infile = out_dir+pgc_name+'_w'+str(band)+'_gauss15.fits'
                 outfile = out_dir+pgc_name+'_w'+str(band)+'_gauss15_small.fits'
                 do_rebin = 0B
              endelse
              
              if file_test(infile) eq 0 then begin
                 message, 'File not found '+infile, /info
                 continue
              endif
              
              if keyword_set(incremental) and file_test(outfile) then begin
                 message, 'Image already in place '+outfile, /info
                 continue
              endif
              
              this_rad = this_dat.r25_deg
              if finite(this_rad) eq 0 then this_rad = 30./3600.
              delta = 4.0*(this_rad > 30./3600.)

              override_ind = where(override_pgc_list eq this_dat.pgc, override_ct)
              if override_ct eq 1 then begin
                 if finite(override_rad_list[override_ind]) then $
                    delta = (override_rad_list[override_ind])[0]*2.0
              endif

              map = readfits(infile, hdr, /silent)
              pix = abs(sxpar(hdr, 'CD1_1'))
              delta_pix = delta / pix

              sz = size(map)
              mean_x = sz[1]/2.0
              mean_y = sz[2]/2.0

              low_x = floor(mean_x - delta_pix) > 0
              low_y = floor(mean_y - delta_pix) > 0
              high_x = ceil(mean_x + delta_pix) < (sz[1]-1)
              high_y = ceil(mean_y + delta_pix) < (sz[2]-1)

              hextract, map, hdr, low_x, high_x, low_y, high_y, /silent
              sz = size(map)
              if do_rebin then begin
                 new_sz = [sz[1]/2, sz[2]/2]
                 hrebin, map, hdr, out=new_sz
              endif
              
              writefits, outfile, map, hdr

           endfor

        endfor
        
     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FIT BACKGROUND
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  !p.multi=0

  if keyword_set(do_bksub) then begin
     
     for ii = 0, n_pgc-1 do begin
        
        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]
        
        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Background fit for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        for band = 1, 4 do begin

           for mm = 0, 1 do begin

              if mm eq 0  then begin
                 res_str = 'gauss15'
              endif
              
              if mm eq 1 then begin
                 if band eq 4 then continue
                 res_str = 'gauss7p5'
              endif

              infile = out_dir+pgc_name+'_w'+str(band)+'_'+res_str+'_small.fits'
              outfile = out_dir+pgc_name+'_w'+str(band)+'_'+res_str+'_bksub.fits'
              rejfile = out_dir+pgc_name+'_w'+str(band)+'_'+res_str+'_rejected.fits'

              radfile = mask_dir+pgc_name+'_'+res_str+'_rgrid.fits'
              galfile = mask_dir+pgc_name+'_'+res_str+'_galaxies.fits'
              brightfile = mask_dir+pgc_name+'_w'+str(band)+'_'+res_str+'_bright_stars.fits'
              foundfile = mask_dir+pgc_name+'_w'+str(band)+'_'+res_str+'_found_stars.fits'
              handfile = mask_dir+pgc_name+'_w'+str(band)+'_'+res_str+'_custom.fits'

              masklist = [galfile, brightfile, foundfile, handfile]

              if file_test(infile) eq 0 then begin
                 message, "File missing. Skipping.", /info
                 continue
              endif
              
              if keyword_set(incremental) and file_test(outfile) then begin
                 message, 'Image already in place '+outfile, /info
                 continue
              endif
              
              bkfit_unwise $
                 , mapfile=infile $
                 , outfile=outfile $
                 , rejfile=rejfile $
                 , radfile=radfile $
                 , masklist=masklist $
                 , band='w'+str(band) $
                 , rejected=rejected $
                 , show=show $
                 , pause=pause

           endfor

        endfor     
        
     endfor

  endif
     
end
