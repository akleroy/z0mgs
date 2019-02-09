pro compile_s4g_comparison $
   , overwrite=overwrite $
   , convolve=do_convolve $
   , bksub=do_bksub $
   , sample=do_sample $
   , start=start $
   , stop=stop   
 
; Convolve S4G maps to Z0MGs resolution appropriate for comparison to
; understand how best to interpret the unWISE measurements.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SETUP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  index_dir = '/data/tycho/0/leroy.42/ellohess/code/index/'
  out_dir = '../cutouts/s4g/'
  atlas_dir = '../delivery/'
  unwise_dir = '../unwise/atlas/'
  mask_dir = '../masks/'

  readcol $
     , index_dir + 'processed_irac.txt' $
     , comment='#', format='A,A,A,A' $
     , gal, survey, band, fname

  gdata = gal_data(gal)

  ind = where(survey eq 's4g_release', n_gals)
  
  if n_elements(start) eq 0 then start = 0
  if n_elements(stop) eq 0 then stop = n_gals-1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLVE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_convolve) then begin

     loadct, 33
     !p.multi=[0,2,2]

     for ii = start, stop do begin
        
        for jj = 0, 1 do begin

           if jj eq 0 then begin
              res_str = 'gauss15'
              target_res = 15.
           endif
           if jj eq 1 then begin
              res_str = 'gauss7p5'
              target_res = 7.5
           endif

;          CHECK IF THE GALAXY IS IN THE ATLAS
           pgcname = 'PGC' + str(gdata[ind[ii]].pgc)
           test_file = file_search(atlas_dir + $
                                   pgcname+'_w1_'+res_str+'.fits' $
                                   , count=ct)
           if ct eq 0 then continue
           target_hdr = headfits(test_file)

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
           
;          ALIGN TO THE WISE ASTROMETRY
           hastrom, out_data, out_hdr, target_hdr $
                    , interp=2, cubic=-0.5, missing=!values.f_nan
           
           disp, alog10(out_data), /sq, title=outfile_name
           
;          WRITE TO DISK
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

     for ii = start, stop do begin
        
        counter, ii, n_gals, ' out of '
        
        for jj = 0, 1 do begin

           if jj eq 0 then begin
              target_res = 7.5
              res_str = 'gauss7p5'
           endif
           if jj eq 1 then begin
              target_res = 15.
              res_str = 'gauss15'
           endif

           pgcname = 'PGC' + str(gdata[ind[ii]].pgc)
           print, pgcname

           infile = out_dir + pgcname + $
                    '_' + band[ind[ii]] + '_'+res_str+'.fits'
           
           if file_test(infile) eq 0 then continue           

           radfile = mask_dir + pgcname + $
                     '_'+res_str+'_rgrid.fits'

           outfile = out_dir + pgcname + $
                     '_' + band[ind[ii]] + '_'+res_str+'_bksub.fits'
           
           show = 0B
           if (ii mod 25) eq 0 then show = 1B

           bkfit_unwise $
              , mapfile=infile $
              , radfile=radfile $
              , outfile=outfile $
              , aperture=1.0 $
              , show=show

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

        fid_n = 3e7
        nan = !values.f_nan
        s4g_comp = fltarr(fid_n)*nan
        z0mgs_comp = s4g_comp
        noise_comp = s4g_comp
        band_comp = replicate('', fid_n)
        gal_comp = replicate('', fid_n)
        tracker = 0L

        for ii = start, stop do begin        

           counter, ii, n_gals, ' out of '
           
           pgcname = 'PGC' + str(gdata[ind[ii]].pgc)
        
           infile_name = out_dir + pgcname + $
                         '_' + band[ind[ii]] + '_'+res_str+'_bksub.fits'

           if file_test(infile_name) eq 0 then continue
        
           map_s4g = readfits(infile_name, s4g_hdr, /silent)

           if band[ind[ii]] eq 'irac1' then begin
              map_z0mgs = readfits(atlas_dir+pgcname+'_w1_'+res_str+'.fits', atlas_hdr, /silent)
              this_band = 'WISE1'
           endif else if band[ind[ii]] eq 'irac2' then begin
              map_z0mgs = readfits(atlas_dir+pgcname+'_w2_'+res_str+'.fits', atlas_hdr, /silent)
              this_band = 'WISE2'
           endif else begin
              continue
           endelse
        
           radfile = mask_dir + pgcname + $
                     '_'+res_str+'_rgrid.fits'
           rgrid = readfits(radfile, rhdr, /silent)
           
           noise = sxpar(atlas_hdr,'RMS')
           samp_ind = where(rgrid le sxpar(rhdr,'FIDRAD'), samp_ct)
           if samp_ct eq 0 then continue
        
           s4g_comp[tracker:(tracker+samp_ct-1)] = map_s4g[samp_ind]
           z0mgs_comp[tracker:(tracker+samp_ct-1)] = map_z0mgs[samp_ind]
           noise_comp[tracker:(tracker+samp_ct-1)] = noise
           band_comp[tracker:(tracker+samp_ct-1)] = this_band
           gal_comp[tracker:(tracker+samp_ct-1)] = pgcname
           
           tracker += samp_ct
           
        endfor     
        
        save, file='../measurements/s4g_z0mgs_comp_'+res_str+'.idl' $
              , s4g_comp, z0mgs_comp, noise_comp, band_comp, gal_comp

     endfor

  endif

end
