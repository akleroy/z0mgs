pro compile_atlas $
   , units = do_units $
   , mask = do_mask $
   , catquery = do_catquery $
   , bksub = do_bksub $
   , convol = do_convol $
   , show = show $
   , just = just

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
     
     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Unit conversion '

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        for band = 1, 4 do begin
           
           infile = in_dir+pgc_name+'/'+'unwise-'+pgc_name+'-w'+str(band)+'-img-m.fits'

           outfile = out_dir+pgc_name+'_w'+str(band)+'_mjysr.fits'

           convert_unwise_to_mjysr $
              , infile = infile $
              , outfile = outfile $
              , band = band
           
        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; QUERY IRSA CATALOGS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_catquery) then begin

     all_data = gal_data(/all)

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

     all_data = gal_data(/all)

     for ii = 0, n_pgc-1 do begin

        ;counter, ii, n_pgc, 'Unit conversion '
        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        print, ii*1.0/n_pgc, ':', pgc_name                
        this_dat = all_data[where(all_data.pgc eq pgc_num[ii])]

        build_unwise_mask $
           , pgc = pgc_list[ii] $
           , galdata = this_dat $
           , outfile = outfile $
           , /show

     endfor     

  endif

  !p.multi=0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RUN A BACKGROUND SUBTRACTION OUTSIDE THE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=[0,4,5]

  if keyword_set(do_bksub) then begin

     all_data = gal_data(/all)

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Background subtraction '

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

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

           bkind = where(mask eq 100)
           bklev = median(map[bkind])
           rms = mad(map[bkind])
           map -= bklev
           sxaddpar, hdr, 'NOISE', rms
           sxaddpar, hdr, 'BKGRD', bklev

           outfile =  out_dir+pgc_name+'_w'+str(band)+'_bksub.fits'

           writefits, outfile, map, hdr

           if keyword_set(show) then begin
              disp, map, /sq, xstyle=5, ystyle=5, min=-3.*rms, max=3.*rms
           endif
           
        endfor
        
     endfor
          
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RUN THE CONVOLUTIONS TO MATCH BEAMS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

end
