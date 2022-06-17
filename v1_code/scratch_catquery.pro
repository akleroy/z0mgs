   , catquery = do_catquery $

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; QUERY IRSA CATALOGS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_catquery) then begin

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Catalog query '

        pgc_name = pgc_list[ii]
        
        this_dat = gal_data[ii]        
        
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
