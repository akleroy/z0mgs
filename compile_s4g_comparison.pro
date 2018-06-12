pro compile_s4g_comparison $
   , overwrite=overwrite $
   , convolve=do_convolve $
   , bksub=do_bksub $
   , sample=do_sample

; Convolve S4G maps to Z0MGs resolution appropriate for comparison to
; understand how best to interpret the unWISE measurements.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SETUP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  index_dir = '/data/tycho/0/leroy.42/ellohess/code/index/'
  out_dir = '../cutouts/s4g/'
  atlas_dir = '../delivery/'
  unwise_dir = '../unwise/atlas/'

  readcol $
     , index_dir + 'processed_irac.txt' $
     , comment='#', format='A,A,A,A' $
     , gal, survey, band, fname

  gdata = gal_data(gal)

  ind = where(survey eq 's4g_release', n_gals)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLVE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_convolve) then begin

     loadct, 33
     !p.multi=[0,2,2]

     for ii = 0, n_gals-1 do begin
        
;    CHECK IF THE GALAXY IS IN THE ATLAS
        pgcname = 'PGC' + str(gdata[ind[ii]].pgc)
        test_file = file_search(atlas_dir + $
                                pgcname+'_w1.fits' $
                                , count=ct)
        if ct eq 0 then continue
        target_hdr = headfits(test_file)

;    CONVOLVE THE GALAXY TO 15"
        outfile_name = out_dir + pgcname + $
                       '_' + band[ind[ii]] + '_gauss15.fits'

        if keyword_set(overwrite) eq 0 then $
           if file_test(outfile_name) eq 1 then $
              continue
        
        conv_with_gauss $
           , data=index_dir+fname[ind[ii]] $
           , target_beam=15.*[1,1,0] $
           , out_data=out_data $
           , out_hdr=out_hdr
        
; ALIGN TO THE WISE ASTROMETRY
        hastrom, out_data, out_hdr, target_hdr $
                 , interp=2, cubic=-0.5, missing=!values.f_nan
        
        disp, alog10(out_data), /sq, title=outfile_name
        
; WRITE TO DISK
        writefits, outfile_name, out_data, out_hdr
        
     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BACKGROUND SUBTRACTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_bksub) then begin

     loadct, 0
     !p.multi=0

     for ii = 0, n_gals-1 do begin
        
        counter, ii, n_gals, ' out of '
        
;    CHECK IF THE GALAXY IS IN THE ATLAS
        pgcname = 'PGC' + str(gdata[ind[ii]].pgc)

        infile_name = out_dir + pgcname + $
                      '_' + band[ind[ii]] + '_gauss15.fits'

        test_file = file_search(infile_name $
                                , count=ct)
        
        if ct eq 0 then continue
        map = readfits(infile_name, hdr, /silent)

        maskfile = unwise_dir+pgcname+'_mask.fits'
        test = file_search(maskfile, count=ct)
        if ct eq 0 then $
           continue
        mask = readfits(maskfile, /silent, mask_hdr)
        hastrom, mask, mask_hdr, hdr, interp=0

        bksub = $
           bkfit( $
           map=map $
           , mask=mask $
           , rejected=rejected $
           , niter=5 $
           , thresh=3.0 $
           , method='MEDIAN' $
           , coefs=coefs)
        
        sxaddpar, hdr, 'BKPLANE0', coefs[0]

        outfile = out_dir + pgcname + $
                  '_' + band[ind[ii]] + '_bksub.fits'
        
        writefits, outfile, bksub, hdr

        if ii mod 25 eq 0 then begin
           disp, bksub, /sq, max=mad(bksub)*5., min=-1.*mad(bksub)*5.
           contour, rejected, lev=[1], /overplot, color=cgcolor('red')
        endif

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD S4G SAMPLING FOR COMPARISON
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_sample) then begin
     
     loadct, 0
     !p.multi=0

     for ii = 0, n_gals-1 do begin
        
        counter, ii, n_gals, ' out of '

        pgcname = 'PGC' + str(gdata[ind[ii]].pgc)
        
        infile_name = out_dir + pgcname + $
                      '_' + band[ind[ii]] + '_bksub.fits'

        if file_test(infile_name) eq 0 then continue
        
        map_s4g = readfits(infile_name, s4g_hdr, /silent)

        if band[ind[ii]] eq 'irac1' then begin
           map_z0mgs = readfits(atlas_dir+pgcname+'_w1.fits', atlas_hdr, /silent)
           this_band = 'WISE1'
        endif else if band[ind[ii]] eq 'irac2' then begin
           map_z0mgs = readfits(atlas_dir+pgcname+'_w2.fits', atlas_hdr, /silent)
           this_band = 'WISE2'
        endif else begin
           continue
        endelse

        maskfile = unwise_dir+pgcname+'_mask.fits'
        test = file_search(maskfile, count=ct)
        if ct eq 0 then $
           continue
        mask = readfits(maskfile, /silent, mask_hdr)
        hastrom, mask, mask_hdr, atlas_hdr, interp=0
        
        off_ind = where(mask eq 10, off_ct)
        offset_1 = median(map_s4g[off_ind] - map_z0mgs[off_ind])
        sixlin, map_z0mgs[off_ind], map_s4g[off_ind], a, sa, b, sb

        off_ind = where(mask eq 100, off_ct)
        offset_2 = median(map_s4g[off_ind] - map_z0mgs[off_ind])

        noise = sxpar(atlas_hdr,'MADALL')
        samp_ind = where(mask eq 10 , samp_ct)
        if samp_ct eq 0 then continue
        
        offset = median(map_s4g[samp_ind] - map_z0mgs[samp_ind])

        if n_elements(s4g_comp) eq 0 then begin
           s4g_comp = map_s4g[samp_ind]
           z0mgs_comp = map_z0mgs[samp_ind]
           slope_comp = replicate(b[2], samp_ct)
           inter_comp = replicate(a[2], samp_ct)
           offset1_comp = replicate(offset_1, samp_ct)
           offset2_comp = replicate(offset_2, samp_ct)
           noise_comp = replicate(noise, samp_ct)
           band_comp = replicate(this_band, samp_ct)
           gal_comp = replicate(pgcname, samp_ct)
        endif else begin
           s4g_comp = [s4g_comp, map_s4g[samp_ind]]
           z0mgs_comp = [z0mgs_comp, map_z0mgs[samp_ind]]
           slope_comp = [slope_comp, replicate(b[2], samp_ct)]
           inter_comp = [inter_comp, replicate(a[2], samp_ct)]
           offset1_comp = [offset1_comp, replicate(offset_1, samp_ct)]
           offset2_comp = [offset2_comp, replicate(offset_2, samp_ct)]
           noise_comp = [noise_comp, replicate(noise, samp_ct)]
           band_comp = [band_comp, replicate(this_band, samp_ct)]
           gal_comp = [gal_comp, replicate(pgcname, samp_ct)]
        endelse

        if ii mod 100 eq 0 then begin
           plot, s4g_comp, z0mgs_comp, /xlo, /ylo, ps=3
           equality, color=cgcolor('red')
        endif

     endfor     

     save, file='../measurements/s4g_z0mgs_comp.idl' $
           , s4g_comp, z0mgs_comp, offset1_comp, offset2_comp $
           , slope_comp, inter_comp $
           , noise_comp, band_comp, gal_comp

  endif

end
