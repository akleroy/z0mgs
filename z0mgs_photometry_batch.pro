pro z0mgs_photometry_batch $
   , show = show $
   , pause = pause $
   , just = just $
   , tag = tag $
   , start = start_num $
   , stop = stop_num $
   , galex = do_galex $
   , wise = do_wise $
   , isophot = do_isophot $
   , profile = do_profile

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  unwise_dir = '../unwise/atlas/'
  galex_dir = '../galex/atlas/'
  prof_dir = '../measurements/profiles/'

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  out_dir = '../unwise/atlas/'

  build_galaxy_list $
     , in_dir = in_dir $
     , tag=tag $
     , just=just $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = gal_data $
     , start = start_num $
     , stop = stop_num $
     , exclude = ['PGC17223']
     
  n_pgc = n_elements(pgc_list)

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
           outfile = prof_dir+pgc_name+'_fuv_prof.txt'

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

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP AND MEASURE RADIAL PROFILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  !p.multi=[0,1,4]

  if keyword_set(do_profile) then begin

     for ii = 0, n_pgc-1 do begin

        pgc_name = pgc_list[ii]
        this_dat = gal_data[ii]

        print, ''
        print, 'Radial profiles for for '+str(ii)+' / '+str(n_pgc)+' ... '+pgc_name
        print, ''

        if keyword_set(do_wise) then begin
           
           for band = 1, 4 do begin

              infile = unwise_dir+pgc_name+'_w'+str(band)+'_gauss15.fits'
              outfile = prof_dir+pgc_name+'_w'+str(band)+'_prof.txt'

              z0mgs_profile $
                 , infile=infile $
                 , outfile=outfile $
                 , dat=this_dat $
                 , binsize=7.5 $
                 , band='WISE'+str(band) $
                 , show=show
              
           endfor

        endif

        if keyword_set(do_galex) then begin

           infile = galex_dir+pgc_name+'_fuv_gauss15_align.fits'
           outfile = prof_dir+pgc_name+'_fuv_prof.txt'

           z0mgs_profile $
              , infile=infile $
              , outfile=outfile $
              , dat=this_dat $
              , binsize=7.5 $
              , band='FUV' $
              , show=show

           infile = galex_dir+pgc_name+'_nuv_gauss15_align.fits'
           outfile = prof_dir+pgc_name+'_nuv_prof.txt'

           z0mgs_profile $
              , infile=infile $
              , outfile=outfile $
              , dat=this_dat $
              , binsize=7.5 $
              , band='NUV' $
              , show=show
              
        endif

     endfor

  endif

  !p.multi=0

end
