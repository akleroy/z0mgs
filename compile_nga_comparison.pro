pro compile_nga_comparison $
   , convolve=do_convolve $
   , bksub=do_bksub $
   , sample=do_sample

; Convolve NGA maps to Z0MGs resolution appropriate for comparison to
; understand how best to interpret the unWISE measurements.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SETUP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  index_dir = '/data/tycho/0/leroy.42/ellohess/code/index/'
  out_dir = '../cutouts/nga/'
  atlas_dir = '../delivery/'
  galex_dir = '../galex/atlas/'

  readcol $
     , index_dir + 'processed_galexatlas.txt' $
     , comment='#', format='A,A,A,A' $
     , gal, survey, band, fname

  gdata = gal_data(gal)

  ind = where(survey eq 'nga_release', n_gals)
  
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
                                pgcname+'_nuv.fits' $
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
        
;    CHECK IF THE GALAXY IS IN THE ATLAS
        pgcname = 'PGC' + str(gdata[ind[ii]].pgc)

        infile_name = out_dir + pgcname + $
                      '_' + band[ind[ii]] + '_gauss15.fits'

        test_file = file_search(infile_name $
                                , count=ct)
        
        if ct eq 0 then continue
        map = readfits(infile_name, hdr, /silent)

        maskfile = galex_dir+pgcname+'_mask.fits'
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
        
        map_nga = readfits(infile_name, nga_hdr, /silent)

        if band[ind[ii]] eq 'fuv' then begin
           map_z0mgs = readfits(atlas_dir+pgcname+'_fuv.fits', atlas_hdr, /silent)
           this_band = 'FUV'
           ext_fac = 10.^(sxpar(atlas_hdr, 'AFUV')/2.5)
        endif else begin
           map_z0mgs = readfits(atlas_dir+pgcname+'_nuv.fits', atlas_hdr, /silent)
           this_band = 'NUV'
           ext_fac = 10.^(sxpar(atlas_hdr, 'ANUV')/2.5)
        endelse

        map_nga *= ext_fac

        maskfile = galex_dir+pgcname+'_mask.fits'
        test = file_search(maskfile, count=ct)
        if ct eq 0 then $
           continue
        mask = readfits(maskfile, /silent, mask_hdr)
        hastrom, mask, mask_hdr, atlas_hdr, interp=0

        off_ind = where(mask eq 10, off_ct)
        offset_1 = median(map_nga[off_ind] - map_z0mgs[off_ind])
        sixlin, map_z0mgs[off_ind], map_nga[off_ind], a, sa, b, sb
;        plot,  map_z0mgs[off_ind], map_nga[off_ind], ps=1        
;        xfid = findgen(10000)-1000.
;        oplot, xfid, b[2]*xfid+a[2], lines=0, color=cgcolor('red')

        off_ind = where(mask eq 100, off_ct)
        offset_2 = median(map_nga[off_ind] - map_z0mgs[off_ind])

        noise = sxpar(atlas_hdr,'MADALL')
        samp_ind = where(mask eq 10, samp_ct)
        if samp_ct eq 0 then continue

        if n_elements(nga_comp) eq 0 then begin
           nga_comp = map_nga[samp_ind]
           z0mgs_comp = map_z0mgs[samp_ind]
           slope_comp = replicate(b[2], samp_ct)
           inter_comp = replicate(a[2], samp_ct)
           offset1_comp = replicate(offset_1, samp_ct)
           offset2_comp = replicate(offset_2, samp_ct)
           noise_comp = replicate(noise, samp_ct)
           band_comp = replicate(this_band, samp_ct)
           gal_comp = replicate(pgcname, samp_ct)
        endif else begin
           nga_comp = [nga_comp, map_nga[samp_ind]]
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
;           plot, nga_comp, z0mgs_comp, /xlo, /ylo, ps=3
           !p.multi=[0,2,1]
           plot, offset1_comp, inter_comp, ps=1
           equality, color=cgcolor('red')
           plot, offset2_comp, inter_comp, ps=1
           equality, color=cgcolor('red')
           !p.multi=0
        endif

     endfor     

     save, file='../measurements/nga_z0mgs_comp.idl' $
           , nga_comp, z0mgs_comp, offset1_comp, offset2_comp, slope_comp, inter_comp $
           , noise_comp, band_comp, gal_comp

     stop

  endif

end
