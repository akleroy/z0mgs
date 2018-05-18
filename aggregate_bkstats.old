pro aggregate_bkstats $
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

  coefs = fltarr(3, 4, n_pgc)*!values.f_nan
  imsize = fltarr(2, 4, n_pgc)

  openw, 1, out_dir+'unwise_bkgrd.txt'
  printf, 1, '# column 1: PGC NAME'
  printf, 1, '# column 2: WISE BAND'
  printf, 1, '# column 3: x size'
  printf, 1, '# column 4: y size'
  printf, 1, '# column 5: coefs 0'
  printf, 1, '# column 6: coefs 1'
  printf, 1, '# column 7: coefs 2'

  for ii = 0, n_pgc-1 do begin
     
     counter, ii, n_pgc, 'Parsing header '

     pgc_name = pgc_list[ii]        
     for band = 0, 3 do begin        
        infile = atlas_dir+pgc_name+'_w'+str(band+1)+'_bksub.fits'
        test = file_search(infile, count=ct)
        if ct eq 0 then $
           continue
        hdr = headfits(infile, /silent)
        coefs[0,band,ii] = sxpar(hdr, 'BKPLANE0')
        coefs[1,band,ii] = sxpar(hdr, 'BKPLANE1')
        coefs[2,band,ii] = sxpar(hdr, 'BKPLANE2')
        imsize[0,band,ii] = sxpar(hdr, 'NAXIS1')
        imsize[1,band,ii] = sxpar(hdr, 'NAXIS2')

        line = $
           pgc_name+' '+$
           str(band+1)+' '+$
           str(imsize[0,band,ii])+' '+$
           str(imsize[1,band,ii])+' '+$
           str(coefs[0,band,ii])+' '+$
           str(coefs[1,band,ii])+' '+$
           str(coefs[2,band,ii])

        printf, 1, line
        
     endfor     

  endfor

  close, 1

  readcol, '../measurements/unwise_bkgrd.txt', format='A,I,I,I,F,F,F' $
           , gal, band, x, y, c0, c1, c2
  stop

end
