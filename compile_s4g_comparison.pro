pro compile_s4g_comparison

; Convolve S4G maps to Z0MGs resolution appropriate for comparison to
; understand how best to interpret the unWISE measurements.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SETUP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  index_dir = '/data/tycho/0/leroy.42/ellohess/code/index/'
  out_dir = '../cutouts/s4g/'
  atlas_dir = '../unwise/atlas/'

  readcol $
     , index_dir + 'processed_irac.txt' $
     , comment='#', format='A,A,A,A' $
     , gal, survey, band, fname

  gdata = gal_data(gal)

  ind = where(survey eq 's4g_release', n_gals)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER ALL GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = 0, n_gals-1 do begin
     
;    CHECK IF THE GALAXY IS IN THE ATLAS
     pgcname = 'PGC' + str(gdata[ind[ii]].pgc)
     test_file = file_search(atlas_dir + $
                             pgcname+'_w1_gauss15.fits' $
                            , count=ct)
     if ct eq 0 then continue
     target_hdr = headfits(test_file)

;    CONVOLVE THE GALAXY TO 15"
     outfile_name = out_dir + pgcname + $
                    '_' + band[ind[ii]] + '_gauss15.fits'

     conv_with_gauss $
        , data=index_dir+fname[ind[ii]] $
        , target_beam=15.*[1,1,0] $
        , out_data=out_data $
        , out_hdr=out_hdr

     !p.multi=[0,2,1]
     loadct, 0
     disp, out_data, /sq, title=fname[ind[ii]]

; ALIGN TO THE WISE ASTROMETRY
     hastrom, out_data, out_hdr, target_hdr $
              , interp=2, cubic=-0.5, missing=!values.f_nan

     disp, out_data, /sq, title=outfile_name

; WRITE TO DISK
     writefits, outfile_name, out_data, out_hdr

  endfor

end
