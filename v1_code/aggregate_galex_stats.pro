pro aggregate_galex_stats $
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
  empty_stat = $
     { $
     pgcname:'' $
     , band:'' $
     , max:nan $  
     , med:nan $
     , mean:nan $
     , rms:nan $
     , std:nan $
     , rejfrac:nan $
     , flatrms:nan $
     , flatstd:nan $
     , meanwt:nan $
     , medwt:nan $
     }
  stats = replicate(empty_stat,n_pgc,2)

  for ii = 0, n_pgc-1 do begin
     
     counter, ii, n_pgc, 'Parsing header '     

     pgc_name = pgc_list[ii]        

     for jj = 0, 1 do begin

        if jj eq 0 then band = 'fuv'
        if jj eq 1 then band = 'nuv'

        infile = atlas_dir+pgc_name+'_'+str(band)+'_bksub.fits'
        test = file_search(infile, count=ct)
        if ct eq 0 then $
           continue
        hdr = headfits(infile, /silent)
        
        ext = 'ALL'
        stats[ii,jj].pgcname = pgc_name
        stats[ii,jj].band = band
        stats[ii,jj].max = sxpar(hdr,'MAX'+ext)
        stats[ii,jj].med = sxpar(hdr,'MED'+ext)
        stats[ii,jj].mean = sxpar(hdr,'MEAN'+ext)
        stats[ii,jj].rms = sxpar(hdr,'MAD'+ext)
        stats[ii,jj].std = sxpar(hdr,'STD'+ext)
        stats[ii,jj].rejfrac = sxpar(hdr,'REJ'+ext)
        stats[ii,jj].flatrms = sxpar(hdr,'FLATMADA')
        stats[ii,jj].flatstd = sxpar(hdr,'FLATSTDA')
        stats[ii,jj].meanwt = sxpar(hdr,'MEANWTAL')
        stats[ii,jj].medwt = sxpar(hdr,'MEDWTALL')

     endfor

  endfor

  mwrfits, stats, '../measurements/galex_stats.fits', /create
  save, stats, file='../measurements/galex_stats.idl'

  stop

end
