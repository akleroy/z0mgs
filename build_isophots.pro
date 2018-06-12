; ... not functional right now. Stashing code.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP AND MEASURE ISOPHOTAL PROFILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  !p.multi=[0,2,1]

  if keyword_set(do_isophot) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        print, ''
        print, 'Isophotes for for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        if keyword_set(do_wise) then begin
           
           for band = 1, 4 do begin

              infile = unwise_dir+pgc_name+'_w'+str(band)+'_gauss15.fits'
              outfile = prof_dir+pgc_name+'_w'+str(band)+'_isophot.txt'

              z0mgs_isophotes $
                 , infile=infile $
                 , outfile=outfile $
                 , dat=this_dat $
                 , band='WISE'+str(band) $
                 , show=show $
                 , pause=pause
              
           endfor

        endif

        if keyword_set(do_galex) then begin

           infile = galex_dir+pgc_name+'_fuv_gauss15_align.fits'
           outfile = prof_dir+pgc_name+'_fuv_isohpot.txt'

           z0mgs_isophotes $
              , infile=infile $
              , outfile=outfile $
              , dat=this_dat $
              , band='FUV' $
              , show=show $
              , pause=pause

           infile = galex_dir+pgc_name+'_nuv_gauss15_align.fits'
           outfile = prof_dir+pgc_name+'_nuv_isophot.txt'

           z0mgs_isophotes $
              , infile=infile $
              , outfile=outfile $
              , dat=this_dat $
              , band='NUV' $
              , show=show $
              , pause=pause
              
        endif

     endfor

  endif

  !p.multi=0
