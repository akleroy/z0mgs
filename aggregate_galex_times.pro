pro aggregate_galex_times $
   , just = just $
   , tag = tag $
   , start = start_num $
   , stop = stop_num
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  pgc_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  atlas_dir = '../galex/atlas/'
  out_dir = '../measurements/'

  build_galaxy_list $
     , in_dir = pgc_dir $
     , tag=tag $
     , just=just $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = gal_data $
     , start = start_num $
     , stop = stop_num
  n_pgc = n_elements(pgc_list)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  nan = !values.f_nan
  empty_cov = $
     { $
     pgcname:'' $
     , band:'' $
     , skip_fuv:nan $
     , fuv_cov:nan $  
     , fuv_time:nan $
     , skip_nuv:nan $
     , nuv_cov:nan $  
     , nuv_time:nan $
     }
  cov = replicate(empty_cov,n_pgc)
  
  for ii = 0, n_pgc-1 do begin
     
     counter, ii, n_pgc, 'Parsing header '     

     pgc_name = pgc_list[ii]        

     for jj = 0, 1 do begin        
        
        if jj eq 0 then band = 'fuv'
        if jj eq 1 then band = 'nuv'
        
        infile = atlas_dir+pgc_name+'_'+band+'_cutout.fits'
        test = file_search(infile, count=ct)
        if ct eq 0 then $
           continue
        hdr = headfits(infile, /silent)

        cov[ii].pgcname=pgc_name
        cov[ii].band = strupcase(band)
        if band eq 'fuv' then begin
           cov[ii].skip_fuv = sxpar(hdr,'SKIP')
           cov[ii].fuv_time = sxpar(hdr,'MEANINT')
           cov[ii].fuv_cov = sxpar(hdr,'FRACCOV')
        endif else begin
           cov[ii].skip_nuv = sxpar(hdr,'SKIP')
           cov[ii].nuv_time = sxpar(hdr,'MEANINT')
           cov[ii].nuv_cov = sxpar(hdr,'FRACCOV')
        endelse

     endfor
     
  endfor

  mwrfits, cov, '../measurements/galex_cov.fits', /create
  save, cov, file='../measurements/galex_cov.idl'

  stop

end
