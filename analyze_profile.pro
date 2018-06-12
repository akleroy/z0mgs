function analyze_profile $
   , infile=infile $
   , empty=empty

  nan = !values.f_nan
  phot_struct = $
     { $
     pgcname:'' $
     , band:'' $
     , sum_hybrid: nan $
     , sum_hybrid_bksub: nan $
     , r50_hybrid: nan $
     , r90_hybrid: nan $
     , sum_mean: nan $
     , sum_mean_bksub: nan $
     , r50_mean: nan $
     , r90_mean: nan $
     , sum_med: nan $
     , sum_med_bksub: nan $
     , r50_med: nan $
     , r90_med: nan $
     }
  if keyword_set(empty) then $
     return, phot_struct

  ap_in = 1.5
  ap_out = 3.0

  readcol, infile, format='X,A,A', numline = 10, key, val, /silent
  prof_target = strcompress(val[0], /rem)
  prof_band = strcompress(val[1], /rem)
  prof_file = strcompress(val[2], /rem)
  prof_pa = float(val[3])
  prof_incl = float(val[4])
  prof_ra = float(val[5])
  prof_dec = float(val[6])
  prof_srpix = float(val[7])^2
  prof_rbin = float(val[8])
  prof_r25as = float(val[9])
  if finite(prof_r25as) eq 0 then begin
     print, "No r25."
     prof_r25as = 60.0
  endif

  readcol, infile, comment = '#' $
           , prof_rmid, prof_flux, prof_mean, prof_med $
           , prof_counts, prof_mad, prof_std, prof_16th, prof_84th $
           , /silent

  prof_r25 = prof_rmid / prof_r25as
  
  bk_ind = where(prof_r25 ge ap_in and prof_r25 le ap_out, bk_ct)
  bkval = median(prof_med[bk_ind])
  
  prof_hybrid = prof_med
  prof_hybrid[0:2] = prof_mean[0:2]

  ap_ind = where(prof_r25 le ap_in, ap_ct)
;  print, bk_ct, ap_ct, n_elements(prof_r25), max(prof_r25,/nan), prof_r25as

  sum_med = total(prof_med[ap_ind]*prof_counts[ap_ind], /nan)*prof_srpix*1d6
  sum_mean = total(prof_mean[ap_ind]*prof_counts[ap_ind], /nan)*prof_srpix*1d6
  sum_hybrid =  total(prof_hybrid[ap_ind]*prof_counts[ap_ind], /nan)*prof_srpix*1d6

  sum_med_bksub = total((prof_med[ap_ind]-bkval)*prof_counts[ap_ind], /nan)*prof_srpix*1d6
  sum_mean_bksub = total((prof_mean[ap_ind]-bkval)*prof_counts[ap_ind], /nan)*prof_srpix*1d6
  sum_hybrid_bksub =  total((prof_hybrid[ap_ind]-bkval)*prof_counts[ap_ind], /nan)*prof_srpix*1d6

  cum_med = total(prof_med[ap_ind]*prof_counts[ap_ind], /cumul,/nan)*prof_srpix*1d6/sum_med
  r50_med = interpol(prof_rmid[ap_ind], cum_med, 0.5)
  r90_med = interpol(prof_rmid[ap_ind], cum_med, 0.9)

  cum_mean = total(prof_mean[ap_ind]*prof_counts[ap_ind], /cumul,/nan)*prof_srpix*1d6/sum_mean
  r50_mean = interpol(prof_rmid[ap_ind], cum_mean, 0.5)
  r90_mean = interpol(prof_rmid[ap_ind], cum_mean, 0.9)

  cum_hybrid = total(prof_hybrid[ap_ind]*prof_counts[ap_ind], /cumul,/nan)*prof_srpix*1d6/sum_hybrid
  r50_hybrid = interpol(prof_rmid[ap_ind], cum_hybrid, 0.5)
  r90_hybrid = interpol(prof_rmid[ap_ind], cum_hybrid, 0.9)

  phot_struct.pgcname = prof_target 
  phot_struct.band = prof_band
  phot_struct.sum_hybrid = sum_hybrid
  phot_struct.sum_hybrid_bksub = sum_hybrid_bksub
  phot_struct.r50_hybrid = r50_hybrid
  phot_struct.r90_hybrid = r90_hybrid
  phot_struct.sum_mean =sum_mean
  phot_struct.sum_mean_bksub = sum_mean_bksub
  phot_struct.r50_mean = r50_mean
  phot_struct.r90_mean = r90_mean
  phot_struct.sum_med = sum_med
  phot_struct.sum_med_bksub = sum_med_bksub
  phot_struct.r50_med = r50_med
  phot_struct.r90_med = r90_med

  stop

  return, phot_struct

end
