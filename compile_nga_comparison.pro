pro compile_nga_comparison $
   , convolve=do_convolve $
   , bksub=do_bksub $
   , sample=do_sample $
   , overwrite=overwrite

; Convolve NGA maps to Z0MGs resolution appropriate for comparison to
; understand how best to interpret the unWISE measurements.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SETUP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  index_dir = '/data/tycho/0/leroy.42/ellohess/code/index/'
  out_dir = '../cutouts/nga/'
  atlas_dir = '../delivery/'
  galex_dir = '../galex/atlas/'
  ;mask_dir = '../masks/'

  readcol $
     , index_dir + 'processed_galexatlas.txt' $
     , comment='#', format='A,A,A,A' $
     , gal, survey, band, fname

  gdata = gal_data(gal, /full)
  ind = where(survey eq 'nga_release', n_gals)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLVE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_convolve) then begin

     loadct, 33
     !p.multi=[0,2,2]

     for ii = 0, n_gals-1 do begin
        
        for jj = 0, 1 do begin

           if jj eq 0 then begin
              target_res = 7.5
              res_str = 'gauss7p5'
           endif
           if jj eq 1 then begin
              target_res = 15.
              res_str = 'gauss15'
           endif

;    CHECK IF THE GALAXY IS IN THE ATLAS
           pgcname = 'PGC' + str(gdata[ind[ii]].pgc)
           test_file = file_search(atlas_dir + $
                                   pgcname+'_nuv_'+res_str+'.fits' $
                                   , count=ct)
           if ct eq 0 then continue
           target_hdr = headfits(test_file)

;    CONVOLVE THE GALAXY TO 15"
           outfile_name = out_dir + pgcname + $
                          '_' + band[ind[ii]] + '_'+res_str+'.fits'

           if keyword_set(overwrite) eq 0 then $
              if file_test(outfile_name) eq 1 then $
                 continue
           
           conv_with_gauss $
              , data=index_dir+fname[ind[ii]] $
              , target_beam=target_res*[1,1,0] $
              , out_data=out_data $
              , out_hdr=out_hdr
           
; ALIGN TO THE WISE ASTROMETRY
           hastrom, out_data, out_hdr, target_hdr $
                    , interp=2, cubic=-0.5, missing=!values.f_nan
           
           disp, alog10(out_data), /sq, title=outfile_name
           
; WRITE TO DISK
           writefits, outfile_name, out_data, out_hdr
           
        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BACKGROUND SUBTRACTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_bksub) then begin

     loadct, 0
     !p.multi=0

     for ii = 0, n_gals-1 do begin
        
        for jj = 0, 1 do begin

           if jj eq 0 then begin
              target_res = 7.5
              res_str = 'gauss7p5'
           endif
           if jj eq 1 then begin
              target_res = 15.
              res_str = 'gauss15'
           endif

;    CHECK IF THE GALAXY IS IN THE ATLAS
           pgcname = 'PGC' + str(gdata[ind[ii]].pgc)
           print, pgcname

           infile = out_dir + pgcname + $
                    '_' + band[ind[ii]] + '_'+res_str+'.fits'
           
           if file_test(infile) eq 0 then continue           

           radfile = atlas_dir + pgcname + $
                     '_'+res_str+'_rgrid.fits'

           outfile = out_dir + pgcname + $
                     '_' + band[ind[ii]] + '_'+res_str+'_bksub.fits'
           
           show = 0B
           if (ii mod 25) eq 0 then show = 1B
           
           bkfit_galex $
              , mapfile=infile $
              , radfile=radfile $
              , bkgrd=bkgrd $
              , aperture=1.0 $
              , show=show
           
           map = readfits(infile, hdr)
           map = map - bkgrd
           writefits, outfile, map, hdr
           
        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD S4G SAMPLING FOR COMPARISON
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_sample) then begin
     
     loadct, 0
     !p.multi=0

     for jj = 0, 1 do begin
        
        if jj eq 0 then begin
           target_res = 7.5
           res_str = 'gauss7p5'
        endif
        if jj eq 1 then begin
           target_res = 15.
           res_str = 'gauss15'
        endif
        
        fid_n = 1e7
        nan = !values.f_nan
        nga_comp = fltarr(fid_n)*nan
        z0mgs_comp = nga_comp
        noise_comp = nga_comp
        band_comp = replicate('', fid_n)
        gal_comp = replicate('', fid_n)
        tracker = 0L

        for ii = 0, n_gals-1 do begin
           
           counter, ii, n_gals, ' out of '

           pgcname = 'PGC' + str(gdata[ind[ii]].pgc)
           
           if pgcname eq 'PGC2557' or pgcname eq 'PGC5818' then continue

           infile_name = out_dir + pgcname + $
                         '_' + band[ind[ii]] + '_'+res_str+'_bksub.fits'
           
           if file_test(infile_name) eq 0 then begin
              print, "Missing comparison NGA image.", infile_name
              continue
           endif
           
           map_nga = readfits(infile_name, nga_hdr, /silent)
           
           if band[ind[ii]] eq 'fuv' then begin
              map_z0mgs = readfits(atlas_dir+pgcname+'_fuv_'+res_str+'.fits', atlas_hdr, /silent)
              this_band = 'FUV'
              ext_fac = 10.^(sxpar(atlas_hdr, 'MWEXT')/2.5)
           endif else begin
              map_z0mgs = readfits(atlas_dir+pgcname+'_nuv_'+res_str+'.fits', atlas_hdr, /silent)
              this_band = 'NUV'
              ext_fac = 10.^(sxpar(atlas_hdr, 'MWEXT')/2.5)
           endelse

           map_nga *= ext_fac

           radfile = atlas_dir + pgcname + $
                     '_'+res_str+'_rgrid.fits'
           rgrid = readfits(radfile, rhdr, /silent)
           
           noise = sxpar(atlas_hdr,'RMS')
           samp_ind = where(rgrid le sxpar(rhdr,'FIDRAD'), samp_ct)
           if samp_ct eq 0 then continue
                      
           nga_comp[tracker:(tracker+samp_ct-1)] = map_nga[samp_ind]
           z0mgs_comp[tracker:(tracker+samp_ct-1)] = map_z0mgs[samp_ind]
           noise_comp[tracker:(tracker+samp_ct-1)] = noise
           band_comp[tracker:(tracker+samp_ct-1)] = this_band
           gal_comp[tracker:(tracker+samp_ct-1)] = pgcname

           tracker += samp_ct

           ;if (ii mod 100) eq 0 then begin
           ;   plot, z0mgs_comp, nga_comp, ps=3, /xlo, /ylo $
           ;         , xrange=[1d-4,1d2], yrange=[1d-4,1d2]
           ;endif

        endfor     

        save, file='../measurements/nga_z0mgs_comp_'+res_str+'.idl' $
              , nga_comp, z0mgs_comp, noise_comp, band_comp, gal_comp

     endfor

     stop
     
  endif

end
