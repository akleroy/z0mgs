pro compile_galex_atlas $
   , cutouts = do_cutouts $
   , mask = do_mask $
   , inventory=do_inv $
   , bksub = do_bksub $
   , convol = do_convol $
   , stats = do_stats $
   , show = show $
   , just = just $
   , tag = tag $
   , start = start_num $
   , stop = stop_num $
   , incremental = incremental

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  out_dir = '../galex/atlas/'

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
; CREATE THE CUTOUTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_cutouts) then begin

     !p.multi = [0,4,4]

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

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

           extract_galex_stamp $
              , /fuv $
              , ra_ctr = this_dat.ra_deg $
              , dec_ctr = this_dat.dec_deg $
              , size_deg = (this_dat.r25_deg*6.0) > (655.*2.75/3600.) $
              , index = index_file $
              , image = image $
              , weight = weight $
              , hdr = hdr $
              , /useint $
              , /show
           
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
              , size_deg = (this_dat.r25_deg*6.0) > (655.*2.75/3600.) $
              , index = index_file $
              , image = image $
              , weight = weight $
              , hdr = hdr $
              , /useint $
              , /show
           
           outfile = out_dir + pgc_name+'_nuv_cutout.fits'
           writefits, outfile, image, hdr

           outfile = out_dir + pgc_name+'_nuv_weight.fits'
           writefits, outfile, weight, hdr

        endif

     endfor

     !p.multi = 0

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD A PRELIMINARY MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_mask) then begin

     !p.multi=[0,4,4]

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        print, ''
        print, 'Building a basic mask '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        maskfile = out_dir+pgc_name+'_mask.fits'

        if keyword_set(incremental) then begin
           test = file_search(maskfile, count=test_ct)
           if test_ct gt 0 then continue
        endif

        build_galex_mask $
           , pgc = pgc_list[ii] $
           , galdata = this_dat $
           , outfile = maskfile $
           , /show

     endfor     
     
     !p.multi=0

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; TAKE AN INVENTORY OF COVERAGE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_inv) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

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

           maskfile = out_dir+pgc_name+'_mask.fits'
           test = file_search(maskfile, count=ct)
           if ct eq 0 then $
              continue
           mask = readfits(maskfile, /silent, mask_hdr)

           im = readfits(im_file, hdr, /silent)

           gal_ind = where(mask eq 10, gal_ct)
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
           sxaddpar, rrhr_hdr, 'FRACCOV', mean_int, 'FRACTION OF GALAXY COVERED'
           
           writefits, rrhr_file, rrhr, rrhr_hdr
           writefits, im_file, im, hdr

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

        print, ''
        print, 'Background subtraction '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        maskfile = out_dir+pgc_name+'_mask.fits'
        test = file_search(maskfile, count=ct)
        if ct eq 0 then $
           continue
        mask = readfits(maskfile, /silent, mask_hdr)
        
        for band = 1, 2 do begin           

           if band eq 1 then begin
              infile = out_dir+pgc_name+'_fuv_cutout.fits'
              wtfile = out_dir+pgc_name+'_fuv_weight.fits'
              outfile = out_dir+pgc_name+'_fuv_bksub.fits'
              rejfile = out_dir+pgc_name+'_fuv_rejected.fits'
           endif
           if band eq 2 then begin
              infile = out_dir+pgc_name+'_nuv_cutout.fits'
              wtfile = out_dir+pgc_name+'_nuv_weight.fits'
              outfile = out_dir+pgc_name+'_nuv_bksub.fits'
              rejfile = out_dir+pgc_name+'_nuv_rejected.fits'
           endif
           test = file_search(infile, count=ct)
           if ct eq 0 then $
              continue
           map = readfits(infile, hdr, /silent)
           wt = readfits(wtfile, /silent)

           if sxpar(hdr,'SKIP') then begin
              print, "Skipping "+infile
              continue
           endif

           bksub = $
              bkfit_galex( $
              map=map $
              , mask=mask $
              , wt=wt $
              , kernel=5 $
              , rejected=rejected $
              , niter=5 $
              , thresh=3.0 $
              , method='MEDIAN' $
              , coefs=coefs $
              , show=show)

           sxaddpar, hdr, 'BKPLANE0', coefs[0]
           
           writefits, outfile, bksub, hdr
           writefits, rejfile, rejected*1.0, hdr
           
        endfor
        
     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLUTIONS TO MATCH BEAMS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  !p.multi=[0,4,4]

  if keyword_set(do_convol) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        print, ''
        print, 'Convolution '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''
        
        wise_hdr_file = '../unwise/atlas/'+pgc_name+'_w3_gauss15.fits'
        test = file_search(wise_hdr_file, count=hdr_ct)
        if hdr_ct eq 0 then begin
           message, 'Missing WISE header information.', /info
           continue
        endif
        target_hdr = headfits(wise_hdr_file)
        
        for jj = 0, 1 do begin

           if jj eq 0 then band = 'fuv'
           if jj eq 1 then band = 'nuv'

           infile = out_dir+pgc_name+'_'+band+'_bksub.fits'
           test = file_search(infile, count=ct)
           if ct eq 0 then begin
              message, 'File not found '+infile, /info
              continue
           endif

           for mm = 0, 2 do begin

              if mm eq 0 then ext = 'bksub'
              if mm eq 1 then ext = 'weight'
              if mm eq 2 then ext = 'rejected'

              infile = out_dir+pgc_name+'_'+band+'_'+ext+'.fits'

              test = readfits(infile, hdr)
              if test[0] eq -1 then begin
                 message, "Problematic FITS file. Skipping.", /info
                 continue
              endif

              outfile = $
                 out_dir+pgc_name+'_'+$
                 band+'_'+ext+'_gauss15.fits'
              
              conv_z0mg_galaxy, $
                 infile=infile, $
                 start_psf=band, $
                 end_psf='g15', $
                 outfile=outfile
              
              conv_image = $
                 readfits(outfile, conv_hdr)
              
              hastrom, conv_image, conv_hdr, target_hdr $
                       , cubic=-0.5, interp=2, missing=!values.f_nan

              alignfile = $
                 out_dir+pgc_name+'_'+$
                 band+'_'+ext+'_gauss15_align.fits'
              
              writefits, alignfile, conv_image, conv_hdr

              outfile = $
                 out_dir+pgc_name+'_'+$
                 band+'_'+ext+'_gauss7p5.fits'
              
              conv_z0mg_galaxy, $
                 infile=infile, $
                 start_psf=band, $
                 end_psf='g7p5', $
                 outfile=outfile
              
              conv_image = $
                 readfits(outfile, conv_hdr)
              
              hastrom, conv_image, conv_hdr, target_hdr $
                       , cubic=-0.5, interp=2, missing=!values.f_nan

              alignfile = $
                 out_dir+pgc_name+'_'+ $
                 band+'_'+ext+'_gauss7p5_align.fits'
              
              writefits, alignfile, conv_image, conv_hdr

           endfor
           
        endfor

     endfor
     
  endif
  
end
