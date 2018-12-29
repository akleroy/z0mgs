pro compile_galex_atlas $
   , cutouts = do_cutouts $
   , inventory=do_inv $
   , convol = do_convol $
   , extract = do_extract $
   , mask = do_mask $
   , bksub = do_bksub $
   , special = do_special $
   , stats = do_stats $
   , extcorr = do_extcorr $
   , show = show $
   , pause = pause $
   , tag = tag $
   , start = start_num $
   , stop = stop_num $
   , incremental = incremental $
   , only = only $
   , just = just

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  out_dir = '../galex/atlas/'
  mask_dir = '../masks/'

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
; CREATE THE CUTOUTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_cutouts) then begin

     !p.multi = [0,4,4]

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Extracting GALEX cutouts '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        found = 0B
        if keyword_set(incremental) then begin
           outfile = out_dir + pgc_name+'_fuv_cutout.fits'
           test = file_search(outfile, count=im_ct)

           outfile = out_dir + pgc_name+'_fuv_weight.fits'
           test = file_search(outfile, count=wt_ct)

           found = (im_ct eq 1 and wt_ct eq 1)
        endif

        if keyword_set(incremental) eq 0 or $
           found eq 0 then begin

           size_deg = (this_dat.r25_deg*6.0) > (655.*2.75/3600.)
           if pgc_name eq 'PGC2557' then $
              size_deg = (this_dat.r25_deg*3.0)

           extract_galex_stamp $
              , /fuv $
              , ra_ctr = this_dat.ra_deg $
              , dec_ctr = this_dat.dec_deg $
              , size_deg = size_deg $
              , index = index_file $
              , image = image $
              , weight = weight $
              , hdr = hdr $
              , /useint $
              , show=show
           
           outfile = out_dir + pgc_name+'_fuv_cutout.fits'
           writefits, outfile, image, hdr

           outfile = out_dir + pgc_name+'_fuv_weight.fits'
           writefits, outfile, weight, hdr

        endif

        found = 0B
        if keyword_set(incremental) then begin
           outfile = out_dir + pgc_name+'_nuv_cutout.fits'
           test = file_search(outfile, count=im_ct)

           outfile = out_dir + pgc_name+'_nuv_weight.fits'
           test = file_search(outfile, count=wt_ct)

           found = (im_ct eq 1 and wt_ct eq 1)
        endif

        if keyword_set(incremental) eq 0 or $
           found eq 0 then begin

           extract_galex_stamp $
              , fuv=0 $
              , ra_ctr = this_dat.ra_deg $
              , dec_ctr = this_dat.dec_deg $
              , size_deg = size_deg $
              , index = index_file $
              , image = image $
              , weight = weight $
              , hdr = hdr $
              , /useint $
              , show=show
           
           outfile = out_dir + pgc_name+'_nuv_cutout.fits'
           writefits, outfile, image, hdr

           outfile = out_dir + pgc_name+'_nuv_weight.fits'
           writefits, outfile, weight, hdr

        endif

     endfor

     !p.multi = 0

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; TAKE AN INVENTORY OF COVERAGE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_inv) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Running GALEX inventory '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        for jj = 0, 1 do begin
           
           if jj eq 0 then band = 'fuv'
           if jj eq 1 then band = 'nuv'
           
           im_file = out_dir+pgc_name+'_'+band+'_cutout.fits'

           if keyword_set(incremental) then begin
              im_hdr = headfits(im_file)
              dummy = sxpar(hdr, 'SKIP', count=count_skip)
              dummy = sxpar(hdr, 'MEANINT', count=count_int)
              dummy = sxpar(hdr, 'FRACCOV', count=count_cov)
              if count_skip ge 1 and count_int ge 1 and count_cov ge 1 then $
                 continue
           endif

           im = readfits(im_file, hdr, /silent)
           
           make_axes, hdr, ri=ri, di=di
           
           pa = this_dat.posang_deg
           incl = this_dat.incl_deg        
           
           if finite(incl) eq 0 then begin
              incl = 0.0           
           endif
           
           if incl gt 60 then begin
              incl = 60.
           endif
           
           if finite(pa) eq 0 then begin
              pa = 0.0
              incl = 0.0
           endif
           
           xctr = this_dat.ra_deg
           yctr = this_dat.dec_deg
           
           gal_vec = [pa, incl, xctr, yctr]
           deproject, ri, di, gal_vec, rgrid=rgrid

           fid_rad = this_dat.r25_deg > 30./3600.
           gal_ind = where(rgrid lt fid_rad, gal_ct)
           if gal_ct eq 0 then begin
              print, "No pixels inside the galaxy for "+pgc_name
           endif

           rrhr_file = out_dir+pgc_name+'_'+band+'_weight.fits'
           test = file_search(rrhr_file, count=ct)
           if ct eq 0 then begin
              print, "Didn't find rrhr file "+rrhr_file
              continue
           endif
           rrhr = readfits(rrhr_file, /silent, rrhr_hdr)

           mean_int = mean(rrhr[gal_ind])
           frac_covered = total(finite(im[gal_ind]) and $
                                (rrhr[gal_ind] gt 0.0))/(1.0d*gal_ct)
           skip = frac_covered lt 0.75

           print, '... INT TIME: '+str(mean_int)+' FRAC: '+ $
                  str(frac_covered)+' SKIP? '+str(1L*skip)

           sxaddpar, hdr, 'SKIP', skip, 'PROCESS / NOT'
           sxaddpar, hdr, 'MEANINT', mean_int, 'SECONDS'
           sxaddpar, hdr, 'FRACCOV', frac_covered, 'SECONDS'
           sxaddpar, rrhr_hdr, 'FRACCOV', frac_covered, 'FRACTION OF GALAXY COVERED'
           
           writefits, rrhr_file, rrhr, rrhr_hdr
           writefits, im_file, im, hdr

        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLUTIONS TO MATCH BEAMS AND ALIGN TO THE WISE ASTROMETRY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  if keyword_set(do_convol) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Convolution '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''
        
        wise_hdr_file = '../unwise/atlas/'+pgc_name+'_w1_gauss15.fits'
        if file_test(wise_hdr_file) eq 0 then $
           wise_hdr_file = '../unwise/atlas/'+pgc_name+'_w2_gauss15.fits'
        if  file_test(wise_hdr_file) eq 0 then begin
           message, 'Missing WISE header information.', /info
           continue
        endif
        target_hdr = headfits(wise_hdr_file)
        
        for jj = 0, 1 do begin

           if jj eq 0 then band = 'fuv'
           if jj eq 1 then band = 'nuv'

           infile = out_dir+pgc_name+'_'+band+'_cutout.fits'
           test = file_search(infile, count=ct)
           if ct eq 0 then begin
              message, 'File not found '+infile, /info
              continue
           endif

           for mm = 0, 3 do begin

              if mm eq 0 then begin
                 infile = out_dir+pgc_name+'_'+band+'_cutout.fits'
                 outfile = out_dir+pgc_name+'_'+band+'_gauss15.fits'
                 alignfile = out_dir+pgc_name+'_'+band+'_gauss15_align.fits'
                 end_psf = 'g15'
                 bmaj = 15./3600.
              endif
              if mm eq 1 then begin
                 infile = out_dir+pgc_name+'_'+band+'_weight.fits'
                 outfile = out_dir+pgc_name+'_'+band+'_weight_gauss15.fits'
                 alignfile = out_dir+pgc_name+'_'+band+'_weight_gauss15_align.fits'
                 end_psf = 'g15'
                 bmaj = 15./3600.
              endif
              if mm eq 2 then begin
                 infile = out_dir+pgc_name+'_'+band+'_cutout.fits'
                 outfile = out_dir+pgc_name+'_'+band+'_gauss7p5.fits'
                 alignfile = out_dir+pgc_name+'_'+band+'_gauss7p5_align.fits'
                 end_psf = 'g7p5'
                 bmaj = 7.5/3600.
              endif
              if mm eq 3 then begin
                 infile = out_dir+pgc_name+'_'+band+'_weight.fits'
                 outfile = out_dir+pgc_name+'_'+band+'_weight_gauss7p5.fits'
                 alignfile = out_dir+pgc_name+'_'+band+'_weight_gauss7p5_align.fits'
                 end_psf = 'g7p5'
                 bmaj = 7.5/3600.
              endif

              test = readfits(infile, hdr)
              if test[0] eq -1 then begin
                 message, "Problematic FITS file. Skipping.", /info
                 continue
              endif

              conv_z0mg_galaxy, $
                 infile=infile, $
                 start_psf=band, $
                 end_psf=end_psf, $
                 outfile=outfile
              
              conv_image = $
                 readfits(outfile, conv_hdr)
              
              hastrom, conv_image, conv_hdr, target_hdr $
                       , cubic=-0.5, interp=2, missing=!values.f_nan
              
              sxaddpar, conv_hdr, 'BMAJ', bmaj
              sxaddpar, conv_hdr, 'BMIN', bmaj
              sxaddpar, conv_hdr, 'BPA', 0.0

              writefits, alignfile, conv_image, conv_hdr

           endfor
           
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
        
        for jj = 0, 1 do begin
           
           if jj eq 0 then band = 'fuv'
           if jj eq 1 then band = 'nuv'

           for mm = 0, 1 do begin
          
              if mm eq 0 then begin
                 infile = out_dir+pgc_name+'_'+band+'_gauss7p5_align.fits'
                 wtfile = out_dir+pgc_name+'_'+band+'_weight_gauss7p5_align.fits'
                 outfile = out_dir+pgc_name+'_'+band+'_gauss7p5_small.fits'
                 wtoutfile = out_dir+pgc_name+'_'+band+'_weight_gauss7p5_small.fits'
                 do_rebin = 0B
              endif else begin
                 infile = out_dir+pgc_name+'_'+band+'_gauss15_align.fits'
                 wtfile = out_dir+pgc_name+'_'+band+'_weight_gauss15_align.fits'
                 outfile = out_dir+pgc_name+'_'+band+'_gauss15_small.fits'
                 wtoutfile = out_dir+pgc_name+'_'+band+'_weight_gauss15_small.fits'
                 do_rebin = 1B
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
              pix = abs(sxpar(hdr, 'CDELT1'))
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
              
              weight = readfits(wtfile, weight_hdr, /silent)
              hastrom, weight, weight_hdr, hdr $
                       , interp=2, cubic=-0.5, missing=!values.f_nan
              writefits, wtoutfile, weight, weight_hdr

           endfor
           
        endfor
        
     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BACKGROUND SUBTRACTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_bksub) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Background subtraction '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        maskfile = out_dir+pgc_name+'_mask.fits'
        test = file_search(maskfile, count=ct)
        if ct eq 0 then $
           continue
        mask = readfits(maskfile, /silent, mask_hdr)
        
        for jj = 0, 1 do begin
           if jj eq 0 then band = 'fuv'
           if jj eq 1 then band = 'nuv'
           
           for mm = 0, 1 do begin
              if mm eq 0 then res_str = 'gauss15'
              if mm eq 1 then res_str = 'gauss7p5'

              infile = out_dir+pgc_name+'_'+band+'_'+res_str+'_small.fits'
              wtfile = out_dir+pgc_name+'_'+band+'_weight_'+res_str+'_small.fits'
              outfile = out_dir+pgc_name+'_'+band+'_'+res_str+'_bksub.fits'
              rejfile = out_dir+pgc_name+'_'+band+'_'+res_str+'_rejected.fits'

              radfile = mask_dir+pgc_name+'_'+res_str+'_rgrid.fits'
              galfile = mask_dir+pgc_name+'_'+res_str+'_galaxies.fits'
              brightfile = mask_dir+pgc_name+'_'+str(band)+'_'+res_str+'_bright_stars.fits'
              foundfile = mask_dir+pgc_name+'_'+str(band)+'_'+res_str+'_found_stars.fits'
              handfile = mask_dir+pgc_name+'_'+str(band)+'_'+res_str+'_custom.fits'
              
              masklist = [galfile, brightfile, foundfile, handfile]

              if file_test(infile) eq 0 then begin
                 message, "File missing. Skipping.", /info
                 continue
              endif

              if keyword_set(incremental) and file_test(outfile) then begin
                 message, 'Image already in place '+outfile, /info
                 continue
              endif

              hdr = headfits(infile, /silent)
              if sxpar(hdr,'SKIP') then begin
                 print, "Skipping "+infile+" based on header."
                 continue
              endif
          
              bkfit_galex $
                 , mapfile=infile $
                 , wtfile=wtfile $
                 , outfile=outfile $
                 , rejfile=rejfile $
                 , radfile=radfile $
                 , masklist=masklist $
                 , band=str(band) $
                 , rejected=rejected $
                 , show=show $
                 , pause=pause $
                 , /plane
              
           endfor
        
        endfor

     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SPECIAL POST- OR RE-PROCESSING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_special) then begin

     !p.multi=0

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii] 
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Special processing for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        for jj = 0, 1 do begin
           if jj eq 0 then begin
              band = 'fuv'
           endif
           if jj eq 1 then begin
              band = 'nuv'
           endif

           for mm = 0, 1 do begin
              if mm eq 0 then res_str = 'gauss15'
              if mm eq 1 then res_str = 'gauss7p5'

              infile = out_dir+pgc_name+'_'+band+'_'+res_str+'_small.fits'
              wtfile = out_dir+pgc_name+'_'+band+'_weight_'+res_str+'_small.fits'
              outfile = out_dir+pgc_name+'_'+band+'_'+res_str+'_bksub.fits'
              rejfile = out_dir+pgc_name+'_'+band+'_'+res_str+'_rejected.fits'

              radfile = mask_dir+pgc_name+'_'+res_str+'_rgrid.fits'
              galfile = mask_dir+pgc_name+'_'+res_str+'_galaxies.fits'
              brightfile = mask_dir+pgc_name+'_'+str(band)+'_'+res_str+'_bright_stars.fits'
              foundfile = mask_dir+pgc_name+'_'+str(band)+'_'+res_str+'_found_stars.fits'
              handfile = mask_dir+pgc_name+'_'+str(band)+'_'+res_str+'_custom.fits'

              @special_galex_processing.pro

           endfor

        endfor

     endfor

  endif
        
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; APPLY EXTINCTION CORRECTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  if keyword_set(do_extcorr) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ''
        print, 'Applying extinction correction for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        ebv_sfd98 = this_dat.av_sfd98/3.1
        afuv = calc_afuv(ebv_sfd98)
        anuv = calc_anuv(ebv_sfd98)

        for jj = 0, 1 do begin
           if jj eq 0 then begin
              band = 'fuv'
              ext = afuv
           endif
           if jj eq 1 then begin
              band = 'nuv'
              ext = anuv
           endif

           for mm = 0, 1 do begin
              if mm eq 0 then res_str = 'gauss15'
              if mm eq 1 then res_str = 'gauss7p5'

              infile =  out_dir+pgc_name+'_'+band+'_'+res_str+'_bksub.fits'
              outfile = out_dir+pgc_name+'_'+band+'_'+res_str+'_extcorr.fits'
              if file_test(infile) eq 0 then $
                 continue
              
              map = readfits(infile, hdr, /silent)
              map *= 10.^(ext/2.5)
              sxaddpar, hdr, 'MWEXT', anuv, 'APPLIED MAG - MILKY WAY EXTINCTION'
              writefits, outfile, map, hdr

           endfor

        endfor

     endfor
     
  endif
     
end
