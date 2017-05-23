pro compile_galex_index

  hr_ra = lindgen(24)

  data_dir = '../galex/sorted_tiles/'
  
  flist = file_search(data_dir+'*/*-int.fits', count=count)

  mean_ra = fltarr(count)*!values.f_nan
  mean_dec = fltarr(count)*!values.f_nan

  empty_field = {fname: '' $
                 , flagfile: '' $
                 , rrhrfile: '' $
                 , fuv: 0B $
                 , nuv: 0B $
                 , mean_ra: !values.f_nan $
                 , mean_dec: !values.f_nan $
                 , max_ra: !values.f_nan $
                 , min_ra: !values.f_nan $
                 , max_dec: !values.f_nan $
                 , min_dec: !values.f_nan $
                }

  fields = replicate(empty_field, count)

  for ii = 0, count-1 do begin

     counter, ii, count, ' out of '

     intname = strmid(flist[ii], strlen(data_dir))
     h = headfits(data_dir+intname)

     flagname = mg_streplace(intname,'-int.fits','-flags.fits')     
     rrhrname = mg_streplace(intname,'-int.fits','-rrhr.fits')     
   
     nuv = (strpos(intname,'-nd-') ne -1)
     fuv = (strpos(intname,'-fd-') ne -1)

     if nuv and fuv then begin
        message, 'Should not be both NUV and FUV.', /info
        stop
     endif

     nx = sxpar(h,'NAXIS1')
     ny = sxpar(h,'NAXIS2')

     xyad, h, 0, 0, r1, d1
     xyad, h, 0, ny-1, r2, d2
     xyad, h, nx-1, 0, r3, d3
     xyad, h, nx-1, ny-1, r4, d4

     ri = [r1, r2, r3, r4]
     di = [d1, d2, d3, d4]

     xyad, h, nx/2, ny/2, mean_ra, mean_dec

     max_ra = max(ri)
     min_ra = min(ri)
     max_dec = max(di)
     min_dec = min(di)

     fields[ii].fname = intname
     fields[ii].flagfile = flagname
     fields[ii].rrhrfile = rrhrname

     fields[ii].fuv = fuv
     fields[ii].nuv = nuv

     fields[ii].mean_ra = mean_ra
     fields[ii].mean_dec = mean_dec

;    Deal with the wrapping of RA
     if (max_ra - min_ra) gt 180. then begin
        temp = max_ra
        max_ra = min_ra
        min_ra = temp
     endif

     fields[ii].min_ra = min_ra
     fields[ii].max_ra = max_ra

     fields[ii].min_dec = min_dec
     fields[ii].max_dec = max_dec

  endfor

  mwrfits, fields, 'galex_index_file.fits', /create

  stop

end
