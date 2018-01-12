pro aggregate_stats $
   , just = just $
   , tag = tag $
   , start = start_num $
   , stop = stop_num
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  pgc_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  atlas_dir = '../unwise/atlas/'
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
     , band:-1 $
     , region:'' $
     , res:'' $
     , max:nan $  
     , med:nan $
     , mean:nan $
     , rms:nan $
     , std:nan $
     , rejfrac:nan $
     }
  stats = replicate(empty_stat,n_pgc,3,4,5)

  stats[*,0,*,*].res = 'NATIVE'
  stats[*,1,*,*].res = '7P5'
  stats[*,2,*,*].res = '15'

  stats[*,*,0,*].band = 1
  stats[*,*,1,*].band = 2
  stats[*,*,2,*].band = 3
  stats[*,*,3,*].band = 4

  stats[*,*,*,0].region = 'ALL'
  stats[*,*,*,1].region = 'Q1'
  stats[*,*,*,2].region = 'Q2'
  stats[*,*,*,3].region = 'Q3'
  stats[*,*,*,4].region = 'Q4'

  for ii = 0, n_pgc-1 do begin
     
     counter, ii, n_pgc, 'Parsing header '     

     pgc_name = pgc_list[ii]        

     for res = 0, 2 do begin

        if res eq 0 then res_str = '_bksub'
        if res eq 1 then res_str = '_gauss7p5'
        if res eq 2 then res_str = '_gauss15'

        for band = 1, 4 do begin        

           infile = atlas_dir+pgc_name+'_w'+str(band)+res_str+'.fits'
           test = file_search(infile, count=ct)
           if ct eq 0 then $
              continue
           hdr = headfits(infile, /silent)

           for reg = 0, 4 do begin

              if reg eq 0 then begin
                 ext = 'ALL'
              endif
              if reg eq 1 then begin
                 ext = 'Q1'
              endif
              if reg eq 2 then begin
                 ext = 'Q2'
              endif
              if reg eq 3 then begin
                 ext = 'Q3'
              endif
              if reg eq 4 then begin
                 ext = 'Q4'
              endif

              stats[ii,res,band-1,reg].pgcname = pgc_name
              stats[ii,res,band-1,reg].band = band
              stats[ii,res,band-1,reg].region = ext
              stats[ii,res,band-1,reg].res = res_str
              stats[ii,res,band-1,reg].max = sxpar(hdr,'MAX'+ext)
              stats[ii,res,band-1,reg].med = sxpar(hdr,'MED'+ext)
              stats[ii,res,band-1,reg].mean = sxpar(hdr,'MEAN'+ext)
              stats[ii,res,band-1,reg].rms = sxpar(hdr,'MAD'+ext)
              stats[ii,res,band-1,reg].std = sxpar(hdr,'STD'+ext)
              stats[ii,res,band-1,reg].rejfrac = sxpar(hdr,'REJ'+ext)

           endfor

        endfor

     endfor

  endfor

  mwrfits, stats, '../measurements/unwise_stats.fits', /create
  save, stats, file='../measurements/unwise_stats.idl'

  stop

end
