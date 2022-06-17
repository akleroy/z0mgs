pro print_unwise_patch

;+
;
;
;-

  readcol $
     , 'filesnotfound.txt', format='A,I' $
     , pgc_name, band
  n_calls = n_elements(pgc_name)
  pgc_num = lonarr(n_calls)
  for ii = 0, n_calls-1 do $
     pgc_num[ii] = (strmid(pgc_name[ii],3,strlen(pgc_name[ii])-1))

  data = gal_data(/all)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE THE CALL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Pixel scale of unwise in arcseconds.

  pix_scale = 2.75

; Loop over galaxies and write calls to extract the data.

  base_call = 'python -u unwise_coadd.py'
  band1_ext = ' --band 1'
  band2_ext = ' --band 2'
  band3_ext = ' --band 3 --bgmatch --medfilt 0'
  band4_ext = ' --band 4 --bgmatch --medfilt 0'
  
  out_dir = 'z0mgs/'
  
  openw, unit, 'unwise_makeup_call.txt', /get_lun

  for ii = 0, n_calls-1 do begin
     
     ind = where(data.pgc eq pgc_num[ii])
     this_data = data[ind]
     
     pgc_name = 'PGC'+str(this_data.pgc)
     ra_str = str(this_data.ra_deg)
     dec_str = str(this_data.dec_deg)

     this_r25 = this_data.r25_deg*3600.
     if finite(this_r25) eq 0 then $     
        size_str = '655' $
     else begin
        target_size = 6.*this_r25/pix_scale
        if target_size gt 655 then begin
           print, "Big galaxy: ", this_data.name
           size_str = str(long(round(target_size)))
        endif else begin
           size_str = '655'
        endelse
     endelse

     if band[ii] eq 1 then this_band_ext = band1_ext
     if band[ii] eq 2 then this_band_ext = band2_ext
     if band[ii] eq 3 then this_band_ext = band3_ext
     if band[ii] eq 4 then this_band_ext = band4_ext

     call = base_call + ' -o '+out_dir+ $
            ' --size '+size_str+ $
            ' --ra '+ra_str+ $
            ' --dec '+dec_str+ $
            ' --name '+pgc_name+ $
            this_band_ext
     
     printf, unit, call
     
  endfor

  close, unit

end
