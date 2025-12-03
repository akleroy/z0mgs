pro make_z0mgs_sfr_mstar $
   , do_m31 = do_m31 $
   , just=just $
   , start=start $
   , pause=pause

; Convert processed z0mgs GALEX and WISE data into maps of SFR and
; Mstar. 

; Known problems
;
; Make sure that delete / deprecate the gauss7p5 cases with WISE4 in
; them. They were blank.
;
; Consider to build mass-to-light ratios at lower resolution, and
; possibly smooth them.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COMPILE DIRECTORIES AND LIST OF TARGETS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  @key_dir.pro
  @key_target.pro

  z0mgs_dir_in = '../delivery/'
  z0mgs_dir_out = '../working_data/z0mgs/'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  first = 1
  for ii = 0, n_gals-1 do begin

     if keyword_set(do_m31) eq 0 then $
        if dat[ii].pgc eq 2557 then $
           continue

     phangs_name = $
        strlowcase(strcompress(dat[ii].objname, /rem))

     if n_elements(just) gt 0 then begin
        if total(just eq phangs_name) eq 0 then $
           continue
     endif
     
     if n_elements(start) gt 0 and first eq 1 then begin
        if total(start eq phangs_name) eq 0 then $
           continue $
        else $
           first = 0
     endif

     print, '-----------------------------------------'
     print, 'Processing '+phangs_name
     print, '-----------------------------------------'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&
; LOOP OVER RESOLUTIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     for kk = 0, 1 do begin

        if kk eq 0 then begin
           ;print, "No 7.5arcsec for now."
           ;continue
           res = 'gauss7p5'
           beam = 7.5
        endif

        if kk eq 1 then begin
           res = 'gauss15'
           beam = 15.
        endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        
        w1_fname = z0mgs_dir_in+phangs_name+'_w1_'+res+'_interpol.fits'
        has_w1 = file_test(w1_fname)
        if has_w1 then $
           w1 = readfits(w1_fname, w1_hdr) $
        else begin
           print, "Can't proceed without w1."
           continue
        endelse

        w2_fname = z0mgs_dir_in+phangs_name+'_w2_'+res+'_interpol.fits'
        has_w2 = file_test(w2_fname)
        if has_w2 then $
           w2 = readfits(w2_fname, w2_hdr) $
        else $
           w2 = !values.f_nan*w2

        w3_fname = z0mgs_dir_in+phangs_name+'_w3_'+res+'_interpol.fits'
        has_w3 = file_test(w3_fname)
        if has_w3 then $
           w3 = readfits(w3_fname, w3_hdr) $
        else $
           w3 = !values.f_nan*w3

        w4_fname = z0mgs_dir_in+phangs_name+'_w4_'+res+'_interpol.fits'
        has_w4 = file_test(w4_fname)
        if has_w4 then $
           w4 = readfits(w4_fname, w4_hdr) $
        else $
           w4 = !values.f_nan*w1

        nuv_fname = z0mgs_dir_in+phangs_name+'_nuv_'+res+'_interpol.fits'
        has_nuv = file_test(nuv_fname)
        if has_nuv then $
           nuv = readfits(nuv_fname, nuv_hdr) $
        else $
           nuv = !values.f_nan*w1

        fuv_fname = z0mgs_dir_in+phangs_name+'_fuv_'+res+'_interpol.fits'
        has_fuv = file_test(fuv_fname)
        if has_fuv then $
           fuv = readfits(fuv_fname, fuv_hdr) $
        else $
           fuv = !values.f_nan*w1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CALCULATE STAR FORMATION RATE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        sfr_struct = $
           int_to_sfr( $
           fuv=fuv, nuv=nuv $
           , w1=w1, w2=w2, w3=w3, w4=w4)
        
        if has_w4 and has_fuv then begin
           best_sfr = sfr_struct.fuvw4 
           wise_used = 'w4'
        endif else if has_w4 and has_nuv then begin
           best_sfr = sfr_struct.nuvw4
           wise_used = 'w4'
        endif else if has_w4 then begin
           best_sfr = sfr_struct.w4only
           wise_used = 'w4'
        endif else if has_w3 and has_fuv then begin
           best_sfr = sfr_struct.fuvw3
           wise_used = 'w3'
        endif else if has_w3 and has_nuv then begin
           best_sfr = sfr_struct.nuvw3
           wise_used = 'w3'
        endif else if has_w3 then begin
           best_sfr = sfr_struct.w3only
           wise_used = 'w3'
        endif else begin
           print, "No data for "+phangs_name+". Skipping."
           continue
        endelse
        
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE A RADIAL PROFILE OF MASS-TO-LIGHT RATIO
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Calculate RA and DEC on the WISE1 grid.

        make_axes, w1_hdr, ri=ri, di=di

; Get the orientation of the galaxy

        this_incl = dat[ii].incl_deg
        if this_incl gt 80. then $
           this_incl = 80.

        this_pa = dat[ii].posang_deg
        if finite(this_pa) eq 0 then begin
           this_pa = 0.0
           this_incl = 0.0
        endif

        this_ra_ctr = dat[ii].ra_deg
        this_dec_ctr = dat[ii].dec_deg

; Calculate galactocentric radius and distance on the sky

        deproject, ri, di, [this_pa, this_incl, this_ra_ctr, this_dec_ctr] $
                   , rgrid=rgrid, tgrid=tgrid
        flat_rgrid = sphdist(this_ra_ctr, this_dec_ctr, ri, di, /deg)

; Exclude the center of the galaxy and the minor axis from the calculation

        this_r25 = dat[ii].r25_deg
        exclusion = (2.0*beam/3600 > 0.1*this_r25)
        mask = abs(cos(tgrid)) gt 0.25 or rgrid le exclusion
        ind = where(mask)

; Bin the radial profiles of WISE1 and SFR

        xmin = -0.5 * beam/3600.
        xmax = max(rgrid)
        binsize = beam/3600.

        w1_bins = $
           bin_data(rgrid[ind], w1[ind], /nan $
                    , xmin=xmin, xmax=xmax, binsize=binsize)

        sfr_bins = $
           bin_data(rgrid[ind], best_sfr[ind], /nan $
                    , xmin=xmin, xmax=xmax, binsize=binsize)
        
; Estimate the mass-to-light ratio based on SFR/WISE1

; 30Jul20 - modify this to use explicit W4-to-W1 recipe?

        mtol_prof = int_to_mtol(sfr=sfr_bins.ymean, w1=w1_bins.ymean)

; Convert the mass-to-light ratio into a map in locations with
; significant detections.

        outer =  (w1_bins.ymean gt 5.*sxpar(w1_hdr, 'RMS')) and $
                 (w1_bins.xmid le this_r25) and $
                 (w1_bins.xmid ge 0.5*this_r25)
        outer_ct = total(outer)

        thresh = (w1_bins.ymean gt 5.*sxpar(w1_hdr, 'RMS')) and $
                 (w1_bins.xmid le this_r25)
        thresh_ct = total(thresh)

        if outer_ct gt 0 then $
           med_ind = where(outer eq 1, ct) $
        else $
           med_ind = where(thresh eq 1, ct)
        median_mtol = median(mtol_prof[med_ind])

        empty_ind = where(thresh eq 0, empty_ct)
        if ct gt 0 then mtol_prof[empty_ind] = median_mtol
        
; Smooth a bit?

        ;kern = [0.25, 0.5, 0.25]
        target_pc = 2.0e3/2.354
        sig_deg = target_pc/(dat[ii].dist_mpc*1d6)/!dtor
        sig_beams = sig_deg/(beam/3600.)

        xfid = findgen(ceil(sig_beams*2.354*3)+1)
        xfid = xfid-mean(xfid)
        xfid *= binsize
        kern = exp(-1.*xfid^2/2./sig_deg^2)
        kern = kern/total(kern)
        mtol_prof = convol(mtol_prof, kern, /edge_truncate)
        mtol_prof = convol(mtol_prof, kern, /edge_truncate)

; Grid the mass-to-light ratio back onto a radial grid.

        mtol_map = interpol(mtol_prof, w1_bins.xmid, rgrid)
           
; Make a mass map

        mass_map = 3.3d2*(mtol_map/0.5) * w1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW M/L MAP AND PROFILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        !p.multi=[0,2,1]
        loadct, 0
        disp, mtol_map, /sq, min=0., max=0.7
        loadct, 0
        plot, w1_bins.xmid, mtol_prof, ps=1 $
              , xrange=[-0.1*this_r25, 2.0*this_r25]
        ind = where(thresh eq 1, ct)        
        oplot, w1_bins[ind].xmid, mtol_prof[ind], lines=0, color=cgcolor('red')
        !p.multi=0

        if keyword_set(pause) then begin
           print, "Hit a key to continue."
           ch = get_kbrd(1)
        endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE RESULTS TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        if wise_used eq 'w4' then $
           sfr_hdr = w4_hdr
        if wise_used eq 'w3' then $
           sfr_hdr = w3_hdr

        this_out_file = $
           z0mgs_dir_out + phangs_name + $
           '_mtol_'+res+'.fits'
        hdr = w1_hdr
        sxaddpar, hdr, 'BUNIT', 'WISEMTOL'
        writefits, this_out_file, mtol_map, hdr

        this_out_file = $
           z0mgs_dir_out + phangs_name + $
           '_mstar_'+res+'.fits'
        hdr = w1_hdr
        sxaddpar, hdr, 'BUNIT', 'MSUN/PC^2'
        
        writefits, this_out_file, mass_map, hdr

        if has_fuv then begin
           this_out_file = $
              z0mgs_dir_out + phangs_name + $
              '_sfrfuvonly_'+res+'.fits'
           hdr = sfr_hdr
           sxaddpar, hdr, 'BUNIT', 'MSUN/YR/KPC2'
           writefits, this_out_file, sfr_struct.fuvonly, hdr
        endif

        if has_nuv then begin
           this_out_file = $
              z0mgs_dir_out + phangs_name + $
              '_sfrnuvonly_'+res+'.fits'
           hdr = sfr_hdr
           sxaddpar, hdr, 'BUNIT', 'MSUN/YR/KPC2'
           writefits, this_out_file, sfr_struct.nuvonly, hdr
        endif

        if has_w3 then begin
           this_out_file = $
              z0mgs_dir_out + phangs_name + $
              '_sfrw3only_'+res+'.fits'
           hdr = sfr_hdr
           sxaddpar, hdr, 'BUNIT', 'MSUN/YR/KPC2'
           writefits, this_out_file, sfr_struct.w3only, hdr
        endif

        if has_w4 then begin
           this_out_file = $
              z0mgs_dir_out + phangs_name + $
              '_sfrw4only_'+res+'.fits'
           hdr = sfr_hdr
           sxaddpar, hdr, 'BUNIT', 'MSUN/YR/KPC2'
           writefits, this_out_file, sfr_struct.w4only, hdr
        endif

        if has_fuv and has_w4 then begin
           this_out_file = $
              z0mgs_dir_out + phangs_name + $
              '_sfrfuvw4_'+res+'.fits'
           hdr = sfr_hdr
           sxaddpar, hdr, 'BUNIT', 'MSUN/YR/KPC2'
           writefits, this_out_file, sfr_struct.fuvw4, hdr
        endif

        if has_nuv and has_w4 then begin
           this_out_file = $
              z0mgs_dir_out + phangs_name + $
              '_sfrnuvw4_'+res+'.fits'
           hdr = sfr_hdr
           sxaddpar, hdr, 'BUNIT', 'MSUN/YR/KPC2'
           writefits, this_out_file, sfr_struct.nuvw4, hdr
        endif

        if has_fuv and has_w3 then begin
           this_out_file = $
              z0mgs_dir_out + phangs_name + $
              '_sfrfuvw3_'+res+'.fits'
           hdr = sfr_hdr
           sxaddpar, hdr, 'BUNIT', 'MSUN/YR/KPC2'
           writefits, this_out_file, sfr_struct.fuvw3, hdr
        endif

        if has_nuv and has_w3 then begin
           this_out_file = $
              z0mgs_dir_out + phangs_name + $
              '_sfrnuvw3_'+res+'.fits'
           hdr = sfr_hdr
           sxaddpar, hdr, 'BUNIT', 'MSUN/YR/KPC2'
           writefits, this_out_file, sfr_struct.nuvw3, hdr
        endif

     endfor

  endfor

end
