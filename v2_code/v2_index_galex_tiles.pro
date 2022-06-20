pro v2_index_galex_tiles

  root_dir = '/data/fourier/leroy.42/allsky/all_galex_tiles/'

  empty_field = {tile_dir: '' $
                 , int_file: '' $
                 , flag_file: '' $
                 , rrhr_file: '' $
                 , bgsub_file: '' $
                 , int_found: 0B $
                 , flag_found: 0B $
                 , rrhr_found: 0B $
                 , bgsub_found: 0B $
                 , fuv: 0B $
                 , nuv: 0B $
                 , mean_ra: !values.f_nan $
                 , mean_dec: !values.f_nan $
                 , max_ra: !values.f_nan $
                 , min_ra: !values.f_nan $
                 , max_dec: !values.f_nan $
                 , min_dec: !values.f_nan $
                }

  for kk = 0, 1 do begin
     
     if kk eq 0 then begin
        tile_dir = root_dir+'fd/'
     endif

     if kk eq 1 then begin
        tile_dir = root_dir+'nd/'
     endif

     flist = file_search(tile_dir+'*-int.fits.gz', count=count)
     print, 'I found ', count, ' tiles.'
     
     mean_ra = fltarr(count)*!values.f_nan
     mean_dec = fltarr(count)*!values.f_nan
     
     fields = replicate(empty_field, count)
     fields.tile_dir = tile_dir
     
     for ii = 0, count-1 do begin

        counter, ii, count, ' out of '

        intname = strmid(flist[ii], strlen(tile_dir))
        h = headfits(tile_dir+intname)

        flagname = mg_streplace(intname,'-int.fits.gz' $
                                ,'-flags.fits.gz')     
        rrhrname = mg_streplace(intname,'-int.fits.gz' $
                                ,'-rrhr.fits.gz')     
        bgsubname = mg_streplace(intname,'-int.fits.gz' $
                                 ,'-intbgsub.fits.gz')     
        
        nuv = (strpos(intname,'-nd-') ne -1)
        fuv = (strpos(intname,'-fd-') ne -1)

        if nuv and fuv then begin
           message, 'Should not be both NUV and FUV.', /info
           stop
        endif

        fields[ii].int_file = intname
        fields[ii].flag_file = flagname
        fields[ii].rrhr_file = rrhrname
        fields[ii].bgsub_file = bgsubname

        fields[ii].int_found = file_test(tile_dir+intname)
        fields[ii].flag_found = file_test(tile_dir+flagname)
        fields[ii].rrhr_found = file_test(tile_dir+rrhrname)
        fields[ii].bgsub_found = file_test(tile_dir+bgsubname)
        
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
     
     endfor

     if kk eq 0 then $
        dbase = fields $
     else $
        dbase = [dbase, fields]
  endfor
  
  mwrfits, dbase, 'galex_index_file.fits', /create

  stop

end
