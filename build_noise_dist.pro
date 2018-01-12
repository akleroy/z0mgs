pro build_noise_dist

  restore, '../measurements/unwise_stats_with_dat.idl', /v
  atlas_dir = '../unwise/atlas/'

  xmin = -7.0
  xmax = 6.0
  xbinsize = 0.05
  xbins = (xmax-xmin)/(xbinsize)
  xmid = findgen(xbins)*xbinsize+0.5*xbinsize+xmin

  bmin = 0.
  bmax = 90.
  bbinsize = 10.
  bbins = (bmax-bmin)/bbinsize
  bmid = findgen(bbins)*bbinsize+0.5*bbinsize

  pos = fltarr(4, bbins, xbins)*0.0
  neg = fltarr(4, bbins, xbins)*0.0
  
  first = 1B
  n_gal = n_elements(dat)
  for ii = 0, n_gal-1 do begin
     counter, ii, n_gal, 'Reading galaxy'
     if ii mod 25 eq 0 then begin
        loadct, 33, /silent
        im = reform(alog10(pos[0,*,*]))
        for jj = 0, bbins-1 do $
           im[jj,*] = im[jj,*]/total(im[jj,*],/nan)
        disp, im, bmid, xmid
     endif

     w1_fname = atlas_dir+dat[ii].pgcname+'_w1_gauss15.fits'
     test = file_search(w1_fname, count=ct)
     if ct eq 0 then continue

     w1 = readfits(atlas_dir+dat[ii].pgcname+'_w1_gauss15.fits', w1_hdr, /silent)
     w2 = readfits(atlas_dir+dat[ii].pgcname+'_w2_gauss15.fits', w2_hdr, /silent)
     w3 = readfits(atlas_dir+dat[ii].pgcname+'_w3_gauss15.fits', w3_hdr, /silent)
     w4 = readfits(atlas_dir+dat[ii].pgcname+'_w4_gauss15.fits', w4_hdr, /silent)
     mask = readfits(atlas_dir+dat[ii].pgcname+'_mask.fits', mask_hdr, /silent)
     
     ind = where(mask ne 10 and mask ne 100 and $
                 finite(w1) and finite(w2) and finite(w3) and finite(w4) $
                 , ct)
     if ct eq 0 then continue

     w1_pos_hist = $
        bin_data(alog10(w1[ind]) > xmin $
                 , alog10(w1[ind])*0.0+1.0, /nan $
                 , xmin=xmin, xmax=xmax, binsize=xbinsize)
     w1_neg_hist = $
        bin_data(alog10(-1.0*w1[ind]) > xmin $
                 , alog10(-1.0*w1[ind])*0.0+1.0, /nan $
                 , xmin=xmin, xmax=xmax, binsize=xbinsize)

     w2_pos_hist = $
        bin_data(alog10(w2[ind]) > xmin $
                 , alog10(w2[ind])*0.0+1.0, /nan $
                 , xmin=xmin, xmax=xmax, binsize=xbinsize)
     w2_neg_hist = $
        bin_data(alog10(-1.0*w2[ind]) > xmin $
                 , alog10(-1.0*w2[ind])*0.0+1.0, /nan $
                 , xmin=xmin, xmax=xmax, binsize=xbinsize)

     w3_pos_hist = $
        bin_data(alog10(w3[ind]) > xmin $
                 , alog10(w3[ind])*0.0+1.0, /nan $
                 , xmin=xmin, xmax=xmax, binsize=xbinsize)
     w3_neg_hist = $
        bin_data(alog10(-1.0*w3[ind]) > xmin $
                 , alog10(-1.0*w3[ind])*0.0+1.0, /nan $
                 , xmin=xmin, xmax=xmax, binsize=xbinsize)

     w4_pos_hist = $
        bin_data(alog10(w4[ind]) > xmin $
                 , alog10(w4[ind])*0.0+1.0, /nan $
                 , xmin=xmin, xmax=xmax, binsize=xbinsize)
     w4_neg_hist = $
        bin_data(alog10(-1.0*w4[ind]) > xmin $
                 , alog10(-1.0*w4[ind])*0.0+1.0, /nan $
                 , xmin=xmin, xmax=xmax, binsize=xbinsize)
     
     b_bin = floor((abs(b[ii]) - bmin)/bbinsize)
     pos[0, b_bin, *] += w1_pos_hist.counts*1.0
     pos[1, b_bin, *] += w2_pos_hist.counts*1.0
     pos[2, b_bin, *] += w3_pos_hist.counts*1.0
     pos[3, b_bin, *] += w4_pos_hist.counts*1.0

     neg[0, b_bin, *] += w1_neg_hist.counts*1.0
     neg[1, b_bin, *] += w2_neg_hist.counts*1.0
     neg[2, b_bin, *] += w3_neg_hist.counts*1.0
     neg[3, b_bin, *] += w4_neg_hist.counts*1.0

  endfor
  
  save, file='../measurements/noise_dist.idl', pos, neg, xmid, bmid
  
end
