pro compile_unwise_atlas $
   , units = do_units $
   , mask = do_mask $
   , catquery = do_catquery $
   , bksub = do_bksub $
   , convol = do_convol $
   , stats = do_stats $
   , isophot = do_isophot $
   , slice = do_slice $
   , show = show $
   , just = just $
   , tag = tag

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  out_dir = '../unwise/atlas/'

; BUILD A LIST OF PGC GALAXIES THAT WE HAVE RIGHT NOW
  flist = file_search(in_dir+'PGC*', count=file_ct)
  pgc_list = strarr(file_ct)
  pgc_num = lonarr(file_ct)
  for ii = 0, file_ct-1 do begin
     pgc_list[ii] = strmid(flist[ii],strlen(in_dir),strlen(flist[ii]))
     pgc_num[ii] = long(strmid(pgc_list[ii],3,strlen(pgc_list[ii])-3))
  endfor
  n_pgc = n_elements(pgc_list)

; WRITE 

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY FILES AND CHANGE UNITS TO STAGE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_units) then begin

     if n_elements(tag) gt 0 then begin                
        all_data = gal_data(tag=tag)
     endif else begin
        all_data = gal_data(/all)
     endelse
     
     openw, 1, 'filesnotfound.txt'

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Unit conversion '

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        for band = 1, 4 do begin
           
           infile = in_dir+pgc_name+'/'+'unwise-'+pgc_name+'-w'+str(band)+'-img-m.fits'

           test = file_search(infile, count=ct)
           if ct eq 0 then begin
              printf, 1, pgc_name, band, infile
              print, 'NOT FOUND: ', pgc_name, band, ' ', infile
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

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; QUERY IRSA CATALOGS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_catquery) then begin

     if n_elements(tag) gt 0 then begin                
        all_data = gal_data(tag=tag)
     endif else begin
        all_data = gal_data(/all)
     endelse

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Catalog query '

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue
        
        this_dat = all_data[where(all_data.pgc eq pgc_num[ii])]        

        coords = [this_dat.ra_deg, this_dat.dec_deg]

        radius = this_dat.r25_deg*3600.*5.

        outfile = out_dir+'../catalogs/'+pgc_name+'_allwise_cat.txt'

        z0mgs_query_irsa_cat $
           , coords $
           , catalog='allwise_p3as_psd' $
           , radius=radius $
           , radunits='arcsec' $
           , outfile=outfile $
           , query=url $
           , /noread
        
        print, url

     endfor
     
  endif  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE A BASIC MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=[0,5,5]

  if keyword_set(do_mask) then begin

     if n_elements(tag) gt 0 then begin                
        all_data = gal_data(tag=tag)
     endif else begin
        all_data = gal_data(/all)
     endelse

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Preliminary mask construction '

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        this_ind = where(all_data.pgc eq pgc_num[ii], ct)
        if ct eq 0 then continue
        this_dat = all_data[this_ind]

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

     if n_elements(tag) gt 0 then begin                
        all_data = gal_data(tag=tag)
     endif else begin
        all_data = gal_data(/all)
     endelse

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Background subtraction '

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        this_ind = where(all_data.pgc eq pgc_num[ii], ct)
        if ct eq 0 then continue
        this_dat = all_data[this_ind]

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

     if n_elements(tag) gt 0 then begin                
        all_data = gal_data(tag=tag)
     endif else begin
        all_data = gal_data(/all)
     endelse

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Convolution '

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        this_ind = where(all_data.pgc eq pgc_num[ii], ct)
        if ct eq 0 then continue
        this_dat = all_data[this_ind]

        pgc_name = pgc_list[ii]

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

     if n_elements(tag) gt 0 then begin                
        all_data = gal_data(tag=tag)
     endif else begin
        all_data = gal_data(/all)
     endelse

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Statistics '

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        this_ind = where(all_data.pgc eq pgc_num[ii], ct)
        if ct eq 0 then continue
        this_dat = all_data[this_ind]

        print, pgc_name
        for band = 1, 4 do begin

           infile = out_dir+pgc_name+'_w'+str(band)+'_gauss15.fits'
           maskfile = out_dir+pgc_name+'_mask.fits'
           rejectfile = out_dir+pgc_name+'_w'+str(band)+'_rejected_gauss15.fits'
           outfile = out_dir+pgc_name+'_w'+str(band)+'_gauss15.fits'

           test = file_search(infile, count=ct)
           if ct eq 0 then begin
              message, 'File not found '+infile, /info
              continue
           endif

           test = file_search(maskfile, count=ct)
           if ct eq 0 then begin
              message, 'Mask not found '+infile, /info
              continue
           endif
           
           z0mgs_stat_image $
              , infile=infile $
              , outfile=outfile $
              , mask=maskfile $
              , reject=rejectfile $
              , /print

           if band le 3 then begin

              infile = out_dir+pgc_name+'_w'+str(band)+'_gauss7p5.fits'
              maskfile = out_dir+pgc_name+'_mask.fits'
              rejectfile = out_dir+pgc_name+'_w'+str(band)+'_rejected_gauss7p5.fits'
              outfile = out_dir+pgc_name+'_w'+str(band)+'_gauss7p5.fits'

              test = file_search(infile, count=ct)
              if ct eq 0 then begin
                 message, 'File not found '+infile, /info
                 continue
              endif

              test = file_search(maskfile, count=ct)
              if ct eq 0 then begin
                 message, 'Mask not found '+infile, /info
                 continue
              endif
              
              z0mgs_stat_image $
                 , infile=infile $
                 , outfile=outfile $
                 , mask=maskfile $
                 , reject=rejectfile $
                 , /print

           endif

        endfor
        
     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ISOPHOTAL FITTING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  if keyword_set(do_isophot) then begin

     if n_elements(tag) gt 0 then begin                
        all_data = gal_data(tag=tag)
     endif else begin
        all_data = gal_data(/all)
     endelse

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Isophotal fitting '

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        this_ind = where(all_data.pgc eq pgc_num[ii], ct)
        if ct eq 0 then continue
        this_dat = all_data[this_ind]

        band = 1
        infile = out_dir+pgc_name+'_w'+str(band)+'_gauss7p5.fits'
        z0mgs_isophot_fit $
           , outfile = '../measurements/'+pgc_name+'_isofit.txt' $
           , outimage = '../measurements/'+pgc_name+'_isofit.png' $
           , infile = infile $
           , gal_data = this_dat $
           , /show

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SLICES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  if keyword_set(do_slice) then begin

     if n_elements(tag) gt 0 then begin                
        all_data = gal_data(tag=tag)
     endif else begin
        all_data = gal_data(/all)
     endelse

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Slice construction '

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        this_ind = where(all_data.pgc eq pgc_num[ii], ct)
        if ct eq 0 then continue
        this_dat = all_data[this_ind]

        band = 1
        infile = out_dir+pgc_name+'_w'+str(band)+'_gauss15.fits'
        z0mgs_image_slice $
           , outfile = '../measurements/'+pgc_name+'_slice.txt' $
           , outimage = '../measurements/'+pgc_name+'_slice.png' $
           , infile = infile $
           , gal_data = this_dat $
           , /show        

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PHOTOMETRY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

end
