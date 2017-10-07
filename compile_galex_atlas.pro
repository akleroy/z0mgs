pro compile_galex_atlas $
   , convol = do_convol $
   , show = show $
   , just = just

  in_dir = '/data/tycho/0/lewis.1590/galbase_allsky/cutouts/'
  out_dir = '../galex/atlas/'

; BUILD A LIST OF PGC GALAXIES THAT WE HAVE RIGHT NOW
  pgc_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  flist = file_search(pgc_dir+'PGC*', count=file_ct)
  pgc_list = strarr(file_ct)
  pgc_num = lonarr(file_ct)
  for ii = 0, file_ct-1 do begin
     pgc_list[ii] = strmid(flist[ii],strlen(pgc_dir),strlen(flist[ii]))
     pgc_num[ii] = long(strmid(pgc_list[ii],3,strlen(pgc_list[ii])-3))
  endfor
  n_pgc = n_elements(pgc_list)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RUN THE CONVOLUTIONS TO MATCH BEAMS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  !p.multi=[0,3,1]

  if keyword_set(do_convol) then begin

     all_data = gal_data(/all)

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Convolution '

        pgc_name = pgc_list[ii]        

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        this_data = all_data[where(all_data.pgc eq pgc_num[ii])]

        wise_hdr_file = '../unwise/atlas/'+pgc_name+'_w3_gauss15.fits'
        test = file_search(wise_hdr_file, count=hdr_ct)
        if hdr_ct eq 0 then $
           continue
        target_hdr = headfits(wise_hdr_file)
              
        for jj = 0, 1 do begin

           if jj eq 0 then begin
              infile = in_dir+strcompress(this_data.name,/rem)+'_FUV.FITS'
              outfile = out_dir+pgc_name+'_fuv_gauss15.fits'
              out_align = out_dir+pgc_name+'_fuv_align.fits'
           endif

           if jj eq 1 then begin
              infile = in_dir+strcompress(this_data.name,/rem)+'_NUV.FITS'
              outfile = out_dir+pgc_name+'_nuv_gauss15.fits'
              out_align = out_dir+pgc_name+'_nuv_align.fits'
           endif

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

           conv_z0mg_galaxy, $
              infile=infile, $
              start_psf='fuv', $
              end_psf='g15', $
              outfile=outfile
           
           before = readfits(infile)
           after = readfits(outfile, after_hdr)
           
           loadct, 33
           disp, before, /sq, max=1d-2, min=0
           disp, after, /sq, max=1d-2, min=0        
           
           hastrom, after, after_hdr, target_hdr $
                    , cubic=-0.5, interp=2, missing=!values.f_nan

           disp, after, /sq, max=1d-2, min=0        

           writefits, out_align, after, after_hdr

        endfor

     endfor
     
  endif
  
  
end
