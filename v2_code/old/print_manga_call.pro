pro print_manga_call

  c = 2.99792458e5

  manga_tab = mrdfits('../orig_data/drpall-v3_1_1.fits',1,h)
  ra = manga_tab.objra
  dec = manga_tab.objdec
  z = manga_tab.nsa_z
  n = n_elements(manga_tab)
  ind = lonarr(n)-1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE THE CALL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  name_vec = strarr(n)
  for ii = 0, n-1 do $
     name_vec[ii] = 'MANGA'+str(manga_tab[ii].mangaid)
  ra_vec = ra
  dec_vec = dec
  size_vec = 6.*(30./3600.)+0.*ra_vec
  out_dir = 'unmanga/'

  print_unwise_cutout_call $
     , galname = name_vec $
     , outdir = out_dir $
     , ra_deg = ra_vec $
     , dec_deg = dec_vec $
     , size_deg = size_vec $
     , outfile = '../measurements/manga_unwise_call.txt'

  stop

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CROSS-MATCHING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  thresh = 20./3600.
  all_gals = gal_data(/all, /full)
  all_z = all_gals.vhel_kms/c

  hits = 0L
  misses = 0L
  default_z_mask = finite(all_gals.ra_deg)
  for ii = 0L, n-1 do begin
     counter, ii, n, 'galaxy'
     z_thresh = 0.005
     if z[ii] gt 0.0 then begin
        z_mask = abs(z[ii] - all_z) le z_thresh
     endif else begin
        z_mask = default_z_mask
     endelse
     dist = sphdist(all_gals.ra_deg, all_gals.dec_deg $
                    , ra[ii], dec[ii], /deg)+(z_mask eq 0)*999.
     mindist = min(dist,minind,/nan)
     if mindist le thresh then begin
        ind[ii] = minind
        hits += 1
     endif else begin
        misses += 1
     endelse
     print, hits, misses
  endfor

  print, total(ind ge 0)

  stop
  
end
