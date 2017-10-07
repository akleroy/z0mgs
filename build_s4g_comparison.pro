pro build_s4g_comparison

  atlas_dir = '../unwise/atlas/'
  s4g_dir = '../cutouts/s4g/'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GET A LIST OF S4G GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  s4g_dat = gal_data(tag='S4G')
  n_gal = n_elements(s4g_dat)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE THE HEADER
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  nan = !values.f_nan
  empty_entry = $
     { $
     gal: '' $
     , ra_deg: nan $
     , dec_deg: nan $
     , dist_mpc: nan $
     , posang_deg: nan $
     , incl_deg: nan $
     , beam_as: nan $
     , rgal_as: nan $
     , rgal_kpc: nan $
     , rgal_r25: nan $
     , w1: nan $
     , e_w1: nan $
     , bk_w1: nan $
     , w2: nan $
     , e_w2: nan $
     , bk_w2: nan $
     , w3: nan $
     , e_w3: nan $
     , bk_w3: nan $
     , w4: nan $
     , e_w4: nan $
     , bk_w4: nan $
     , i1: nan $
     , e_i1: nan $
     , bk_i1: nan $
     , i2: nan $
     , e_i2: nan $
     , bk_i2: nan $
     , st: nan $
     , e_st: nan $
     , ns: nan $
     , e_ns: nan $
     }

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = 0, n_gal-1 do begin

     ;if ii gt 10 then continue

     pgcname = 'PGC' + str(s4g_dat[ii].pgc)
     
;      CHECK IF WE HAVE BOTH UNWISE AND S4G
     s4g_file = file_search(s4g_dir+pgcname+'_irac1_gauss15.fits' $
                            , count=s4g_ct)
     atlas_file = file_search(atlas_dir+pgcname+ $
                              '_w1_gauss15.fits' $
                              , count=atlas_ct)
     if atlas_ct eq 0 or s4g_ct eq 0 then $
        continue
     
;    LOAD IN THE DATA
     w1 = readfits(atlas_dir+pgcname+'_w1_gauss15.fits', w1_hdr)
     w2 = readfits(atlas_dir+pgcname+'_w2_gauss15.fits', w2_hdr)
     w3 = readfits(atlas_dir+pgcname+'_w3_gauss15.fits', w3_hdr)
     w4 = readfits(atlas_dir+pgcname+'_w4_gauss15.fits', w4_hdr)

     mask = readfits(atlas_dir+pgcname+'_mask.fits', w4_hdr)

     i1_fname = s4g_dir+pgcname+'_irac1_gauss15.fits'
     i1 = readfits(s4g_dir+pgcname+'_irac1_gauss15.fits', i1_hdr)
     bkind = where(mask eq 100 and finite(i1), bkct)
     if bkct gt 10 then $
        i1_bkgrd = median(i1[bkind]) $
     else $
        i1_bkgrd = !values.f_nan
        
     i2 = readfits(s4g_dir+pgcname+'_irac2_gauss15.fits', i2_hdr)
     if bkct gt 10 then $
        i2_bkgrd = median(i2[bkind]) $
     else $
        i2_bkgrd = !values.f_nan

     st = readfits(s4g_dir+pgcname+'_stellar_gauss15.fits', st_hdr)
     found_st = st[0] ne -1
     ns = readfits(s4g_dir+pgcname+'_nonstellar_gauss15.fits', ns_hdr)
     found_ns = ns[0] ne -1

;    BUILD A SET OF SAMPLING POINTS
     pa = s4g_dat[ii].posang_deg
     if finite(pa) eq 0 then begin
        pa = 0.0 
        incl = 0.0
     endif else begin
        incl = (s4g_dat[ii].incl_deg < 80.)
     endelse 

     make_axes, w1_hdr, ri=ri, di=di
     deproject $
        , s4g_dat[ii].ra_deg, s4g_dat[ii].dec_deg $
        , [pa, incl $
           ,s4g_dat[ii].ra_deg,s4g_dat[ii].dec_deg] $
        , rgrid=rgal_deg, tgrid=theta_rad, /vector  

     mask_hdr = w3_hdr
     mask = finite(i1) and finite(i2) and $
            rgal_deg lt 1.5*s4g_dat[ii].r25_deg

     target_res_as = 15.
     spacing = target_res_as / 3600. / 2.0
     make_sampling_points $
        , ra_ctr = s4g_dat[ii].ra_deg $
        , dec_ctr = s4g_dat[ii].dec_deg $
        , max_rad = 1.5*s4g_dat[ii].r25_deg $
        , spacing = spacing $
        , mask = mask $
        , hdr_mask = mask_hdr $
        , overlay = i1_fname $
        , /show $
        , samp_ra = samp_ra $
        , samp_dec = samp_dec
     if n_elements(samp_ra) le 1 then begin
        message, 'Skipping this galaxy '+pgcname, /info
        continue
     endif

     adxy, w3_hdr, samp_ra, samp_dec, samp_x, samp_y
     n_pts = n_elements(samp_x)

;    CARRY OUT A COUPLE OF NOISE ESTIMATES
     rms_w1 = mad(w1[where(mask eq 0)])
     rms_w2 = mad(w2[where(mask eq 0)])
     rms_w3 = mad(w3[where(mask eq 0)])
     rms_w4 = mad(w4[where(mask eq 0)])

     rms_i1 = mad(i1[where(mask eq 0)])
     rms_i2 = mad(i2[where(mask eq 0)])
     if found_ns then $
        rms_ns = mad(ns[where(mask eq 0)])
     if found_st then $
        rms_st = mad(st[where(mask eq 0)])

;    RECORD THESE INTO THE MINI-DATABASE     
     this_dat = replicate(empty_entry, n_pts)
     
     this_dat.gal = pgcname
     this_dat.ra_deg = samp_ra
     this_dat.dec_deg = samp_dec
     this_dat.posang_deg = s4g_dat[ii].posang_deg
     this_dat.incl_deg = s4g_dat[ii].incl_deg
     this_dat.beam_as = 15.
     this_dat.rgal_as = rgal_deg[samp_x, samp_y]*3600.
     this_dat.rgal_kpc = rgal_deg[samp_x, samp_y]*!dtor*s4g_dat[ii].dist_mpc
     this_dat.rgal_r25 = rgal_deg[samp_x, samp_y]/s4g_dat[ii].r25_deg

     this_dat.w1 = w1[samp_x, samp_y]
     this_dat.w2 = w2[samp_x, samp_y]
     this_dat.w3 = w3[samp_x, samp_y]
     this_dat.w4 = w4[samp_x, samp_y]

     this_dat.e_w1 = rms_w1
     this_dat.e_w2 = rms_w2
     this_dat.e_w3 = rms_w3
     this_dat.e_w4 = rms_w4

     this_dat.bk_w1 = sxpar(w1_hdr, 'BKGRD')
     this_dat.bk_w2 = sxpar(w2_hdr, 'BKGRD')
     this_dat.bk_w3 = sxpar(w3_hdr, 'BKGRD')
     this_dat.bk_w4 = sxpar(w4_hdr, 'BKGRD')

     this_dat.i1 = i1[samp_x, samp_y]
     this_dat.i2 = i2[samp_x, samp_y]
     this_dat.e_i1 = rms_i1
     this_dat.e_i2 = rms_i2

     this_dat.bk_i1 = i1_bkgrd
     this_dat.bk_i2 = i2_bkgrd

     if found_st then begin
     this_dat.st = st[samp_x, samp_y]
     this_dat.e_st = rms_st
     endif

     if found_ns then begin
        this_dat.ns = ns[samp_x, samp_y]
        this_dat.e_ns = rms_ns
     endif     

     if n_elements(dat) eq 0 then begin
        dat = this_dat
     endif else begin
        dat = [dat, this_dat]
     endelse

  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SAVE DATABASE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  out_file = '../measurements/s4g_comp.fits'
  spawn, 'rm -rf '+out_file
  mwrfits, dat, out_file

end
