pro query_unwise $
   , ra=ra_str $
   , dec=dec_str $
   , outdir=outdir $
   , outname=outname $
   , outroot=outroot $
   , nodownload=nodownload $
   , process=process $
   , target_hdr=target_hdr $
   , make_hdr=make_hdr $
   , pause=pause $
   , size_str=size_str $
   , mac=mac

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; PRE-PROCESSING
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=  
  
; PROCESS THE INPUT STRINGS TO FIT IN A URL  
  ra_str = strcompress(ra_str,/rem)
  dec_str = strcompress(dec_str,/rem)

  if n_elements(size_str) eq 0 then $
     size_str='100'
  
; DEFAULT OUTPUT NAME  
  if n_elements(outname) eq 0 then $
     outname='unwise.fits.tar.gz'

; CALIBRATION TO GO FROM VEGAS TO ABMAG
  w1_vtoab = 2.683
  w2_vtoab = 3.319
  w3_vtoab = 5.242
  w4_vtoab = 6.604

; NORMALIZATION OF UNITY IN VEGAS MAG
  norm_mag = 22.5
  pix_as = 2.75
  
; COUNTS -> JY CONVERSION
  w1_to_mjysr = 10.^((norm_mag+w1_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2
  w2_to_mjysr = 10.^((norm_mag+w2_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2
  w3_to_mjysr = 10.^((norm_mag+w3_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2
  w4_to_mjysr = 10.^((norm_mag+w4_vtoab)/(-2.5))*3631.0d/1d6/(pix_as/3600.*!dtor)^2
  
; OUTPUT DIRECTORY
  if n_elements(outdir) eq 0 then $
     outdir='./'

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; QUERY FROM THE WEB
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; MAC-CENTRIC COMMAND
  if keyword_set(mac) then begin
     command = $
        'curl -o '+outdir+outname+' -O ' + $
        '"http://unwise.me/cutout_fits?version=allwise&ra='+ $
        ra_str+"&dec="+dec_str+ $
        '&size='+size_str+'&bands=1234"'
  endif else begin
     command = $
        'wget '+ $
        '-O "'+outdir+outname+'" '+ $
        '"http://unwise.me/cutout_fits?version=allwise&ra='+ $
        ra_str+"&dec="+dec_str+ $
        '&size='+size_str+'&bands=1234"'
  endelse

; CALL IT IF NOT TOLD NOT TO  
  print, command
  if keyword_set(nodownload) eq 0 then $
     spawn, command

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; PROCESS THE DATA
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; ... IF NOT, THEN RETURN  
  if keyword_set(process) eq 0 then begin
     return
  endif

; ... UNTAR  
  spawn, 'tar -zxvf "'+outdir+outname+'" -C '+outdir, tar_file_list
  n_tar = n_elements(tar_file_list)

; ... MAKE A HEADER IF REQUESTED
  if keyword_set(make_hdr) then begin
     pix_len = long(size_str)
     pix_scale = 2.75/3600.
     mkhdr, target_hdr, 2, [pix_len, pix_len]
     sxaddpar, target_hdr, 'NAXIS1', pix_len
     sxaddpar, target_hdr, 'NAXIS2', pix_len
     sxaddpar, target_hdr, 'CTYPE1', 'RA---TAN'
     sxaddpar, target_hdr, 'CRVAL1', double(ra_str)
     sxaddpar, target_hdr, 'CRPIX1', (pix_len/2.)*1.
     sxaddpar, target_hdr, 'CDELT1', -1.0*pix_scale
     sxaddpar, target_hdr, 'CTYPE2', 'DEC--TAN'
     sxaddpar, target_hdr, 'CRVAL2', double(dec_str)
     sxaddpar, target_hdr, 'CRPIX2', (pix_len/2.)*1.
     sxaddpar, target_hdr, 'CDELT2', pix_scale
     sxaddpar, target_hdr, 'EQUINOX', 2000.
  endif

; ... SET UP THE OUTPUT
  make_axes, target_hdr, ri=ri_targ, di=di_targ
  sz_out = size(ri_targ)

; ... LOOP OVER BANDS  
  for band = 1, 4 do begin

     band_file = bytarr(n_tar)
     for kk = 0, n_tar-1 do begin
        if strpos(tar_file_list[kk],'w'+str(band)+'-img-m.fits') ne -1 then $
           band_file[kk] = 1B
     endfor
     flist = outdir+tar_file_list[where(band_file, count)]

     if count eq 0 then begin
        message, 'No result for band '+str(band)+' and '+outname, /info
     endif

     outim = ri_targ*!values.f_nan
     for jj = 0, count-1 do begin
        
        im = readfits(flist[jj], hdr)
        make_axes, hdr, ri=ri, di=di
        
        if band eq 1 then begin
           im *= w1_to_mjysr
        endif
        if band eq 2 then begin
           im *= w2_to_mjysr
        endif
        if band eq 3 then begin
           im *= w3_to_mjysr
        endif
        if band eq 4 then begin
           im *= w4_to_mjysr
        endif
        
        adxy, target_hdr, ri, di, x, y
        
        in_image = (x gt 0 and x lt (sz_out[1]-1)) and $
                   (y gt 0 and y lt (sz_out[2]-1))
        if total(in_image) eq 0 then begin
           print, "No overlap. Proceeding."
           continue
        endif

        hastrom, im, hdr, target_hdr $
                 , interp=1, missing=!values.f_nan
        
        useful = where(finite(im), usect)
        outim[useful] = im[useful]

     endfor

     disp, outim
     if keyword_set(pause) then begin
        print, "Band "+str(band)+" hit a key to continue..."
        ch = get_kbrd(1)
     endif

     sxaddpar, hdr, 'BUNIT', 'MJY/SR'
     writefits, outdir+outroot+'_wise_band'+str(band)+'.fits' $
                , outim, hdr
  endfor

;  stop
  
end
