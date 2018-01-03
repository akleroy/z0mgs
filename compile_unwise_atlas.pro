pro compile_unwise_atlas $
   , units = do_units $
   , mask = do_mask $
   , bksub = do_bksub $
   , convol = do_convol $
   , stats = do_stats $
   , show = show $
   , just = just $
   , tag = tag $
   , start = start_num $
   , stop = stop_num

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  out_dir = '../unwise/atlas/'

  build_galaxy_list $
     , in_dir = in_dir $
     , tag=tag $
     , just=just $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = gal_data $
     , start = start_num $
     , stop = stop_num
  n_pgc = n_elements(pgc_list)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY FILES AND CHANGE UNITS TO STAGE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_units) then begin
     
     not_found_ct = 0L
     openw, 1, 'filesnotfound.txt'

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]

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
; MAKE A BASIC MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=[0,5,5]

  if keyword_set(do_mask) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        print, ''
        print, 'Mask construction for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        maskfile = out_dir+pgc_name+'_mask.fits'
        build_unwise_mask $
           , pgc = pgc_list[ii] $
           , galdata = this_dat $
           , outfile = maskfile $           
           , /show

     endfor     

  endif

  !p.multi=0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RUN A BACKGROUND SUBTRACTION OUTSIDE THE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_bksub) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        print, ''
        print, 'Background subtraction for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''
        
        maskfile = out_dir+pgc_name+'_mask.fits'
        test = file_search(maskfile, count=ct)
        if ct eq 0 then $
           continue
        mask = readfits(maskfile, /silent, mask_hdr)
        
        for band = 1, 4 do begin
           
           infile = out_dir+pgc_name+'_w'+str(band)+'_mjysr.fits'
           test = file_search(infile, count=ct)
           if ct eq 0 then $
              continue
           map = readfits(infile, hdr, /silent)

           bksub = $
              bkfit( $
              map=map $
              , mask=mask $
              , rejected=rejected $
              , niter=5 $
              , thresh=3.0 $
              , method='PLANE' $
              , coefs=coefs)

           sxaddpar, hdr, 'BKPLANE0', coefs[0]
           sxaddpar, hdr, 'BKPLANE1', coefs[1]
           sxaddpar, hdr, 'BKPLANE2', coefs[2]
           
           outfile = out_dir+pgc_name+'_w'+str(band)+'_bksub.fits'
           writefits, outfile, bksub, hdr

           outfile = out_dir+pgc_name+'_w'+str(band)+'_rejected.fits'
           writefits, outfile, rejected*1.0, hdr
           
        endfor
        
     endfor
     
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

        print, ''
        print, 'Convolutions for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        for band = 1, 4 do begin

           infile = out_dir+pgc_name+'_w'+str(band)+'_bksub.fits'

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

           if band le 3 then begin
              outfile = out_dir+pgc_name+'_w'+str(band)+'_gauss7p5.fits'
              conv_z0mg_galaxy, $
                 infile=infile, $
                 start_psf='w'+str(band), $
                 end_psf='g7p5', $
                 outfile=outfile
           endif

;          Convolve the rejected pixel file, too.

           infile = out_dir+pgc_name+'_w'+str(band)+'_rejected.fits'

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

           outfile = out_dir+pgc_name+'_w'+str(band)+'_rejected_gauss15.fits'

           conv_z0mg_galaxy, $
              infile=infile, $
              start_psf='w'+str(band), $
              end_psf='g15', $
              outfile=outfile

           if band le 3 then begin
              outfile = out_dir+pgc_name+'_w'+str(band)+'_rejected_gauss7p5.fits'              
              conv_z0mg_galaxy, $
                 infile=infile, $
                 start_psf='w'+str(band), $
                 end_psf='g7p5', $
                 outfile=outfile
           endif
           
        endfor
        
     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RUN STATISTICS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  if keyword_set(do_stats) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        print, ''
        print, 'Runnings stats for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        for band = 1, 4 do begin

           maskfile = out_dir+pgc_name+'_mask.fits'
           test = file_search(maskfile, count=ct)
           if ct eq 0 then begin
              message, 'Mask not found '+infile, /info
              continue
           endif
           
           for res = 0, 2, do begin

              if res eq 1 and band eq 4 then continue
              
              if res eq 0 then res_str = ''
              if res eq 1 then res_str = '_gauss7p5'
              if res eq 2 then res_str = '_gauss15'
              
              infile = out_dir+pgc_name+'_w'+str(band)+res_str+'.fits'
              rejectfile = out_dir+pgc_name+'_w'+str(band)+'_rejected'+res_str+'.fits'
              outfile = out_dir+pgc_name+'_w'+str(band)+res_str+'.fits'

              test = file_search(infile, count=ct)
              if ct eq 0 then begin
                 message, 'File not found '+infile, /info
                 continue
              endif
           
              z0mgs_stat_image $
                 , infile=infile $
                 , outfile=outfile $
                 , mask=maskfile $
                 , reject=rejectfile $
                 , /print

           endfor

        endfor
        
     endfor
     
  endif

end
