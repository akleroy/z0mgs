pro print_version2_call

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Z0MGS V1 RE-CALL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  z0mgs_tab = gal_data(tag='z0mgs', /full)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ALL GALAXIES WITHIN 50 MPC
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  lv_tab = gal_data(/all, /full)
  lv_tab = lv_tab[where(lv_tab.dist_mpc lt 50.)]
  
  keep = bytarr(n_elements(lv_tab))
  for ii = 0, n_elements(lv_tab)-1 do $
     keep[ii] = total(z0mgs_tab.pgc eq lv_tab[ii].pgc) eq 0
  lv_tab = lv_tab[where(keep)]

  z0mgs_tab = [lv_tab, z0mgs_tab]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CALIFA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  readcol $
     , '../orig_data/CALIFAList.csv' $
     , name, zstars, log_Mass, Re_kpc $
     , log_SFR_SF, lSFR, log_SFR_ssp, log_L_Ha_cor $
     , format='A,F,F,F,F,F,F,F'
  
  califa_tab = gal_data(strcompress(name,/rem),/full)
  califa_tab = califa_tab[where(califa_tab.pgc ne 0)]

  keep = bytarr(n_elements(califa_tab))
  for ii = 0, n_elements(califa_tab)-1 do $
     keep[ii] = total(z0mgs_tab.pgc eq califa_tab[ii].pgc) eq 0
  califa_tab = califa_tab[where(keep)]

  z0mgs_tab = [z0mgs_tab, califa_tab]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE THE CALL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  pgc_vec = strarr(n_elements(z0mgs_tab))
  for ii = 0, n_elements(z0mgs_tab)-1 do $
     pgc_vec[ii] = 'PGC'+str(z0mgs_tab[ii].pgc)
  ra_vec = z0mgs_tab.ra_deg
  dec_vec = z0mgs_tab.dec_deg
  size_vec = 6.*z0mgs_tab.r25_deg
  nan_ind = where(finite(size_vec) eq 0, nan_ct)
  if nan_ct gt 0 then size_vec[nan_ind] = 6.*(30./3600.)
  out_dir = 'z0mgs_v2/'
  
  size_thresh = 0.5
  mega_thresh = 2.0
  small_ind = where(size_vec le size_thresh, small_ct)
  large_ind = where((size_vec gt size_thresh) and $
                    (size_vec le mega_thresh), large_ct)
  mega_ind = where((size_vec gt mega_thresh), mega_ct)
  print, "Small, large, mega: ", small_ct, large_ct, mega_ct

  print_unwise_cutout_call $
     , galname = pgc_vec[small_ind] $
     , outdir = out_dir $
     , ra_deg = ra_vec[small_ind] $
     , dec_deg = dec_vec[small_ind] $
     , size_deg = size_vec[small_ind] $
     , outfile = '../measurements/z0mgs_v3_unwise_call_smallgals.txt'

  print_unwise_cutout_call $
     , galname = pgc_vec[large_ind] $
     , outdir = out_dir $
     , ra_deg = ra_vec[large_ind] $
     , dec_deg = dec_vec[large_ind] $
     , size_deg = size_vec[large_ind] $
     , outfile = '../measurements/z0mgs_v3_unwise_call_largegals.txt'

  stop
  
end
