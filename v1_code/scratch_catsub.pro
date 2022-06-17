pro scratch_catsub

  wise1 = readfits('../unwise/atlas/PGC65001_w1_gauss15.fits', w1_hdr)
  wise2 = readfits('../unwise/atlas/PGC65001_w2_gauss15.fits', w2_hdr)
  wise3 = readfits('../unwise/atlas/PGC65001_w3_gauss15.fits', w3_hdr)
  wise4 = readfits('../unwise/atlas/PGC65001_w4_gauss15.fits', w4_hdr)
  table = z0mgs_read_ipac_table('../unwise/catalogs/PGC65001_allwise_cat.txt')

  make_axes, w1_hdr, ra=ra, da=da, ri=ri, di=di

  this_dat = gal_data('ngc6946')
  deproject, ri, di, gal=this_dat, rgrid=rgrid
  rad_mask = rgrid lt 2.*this_dat.r25_deg

  ind = where(table.w1mpro lt 14)
  adxy, w1_hdr, table.ra[ind], table.dec[ind], pt_x, pt_y

  pt_mask = finite(wise1)*0B
  pt_mask[pt_x, pt_y] = 1B
  pt_mask[where(rad_mask)] = 0B

  rms_wise1 = mad(wise1)
  rms_wise2 = mad(wise2)
  rms_wise3 = mad(wise3)
  rms_wise4 = mad(wise4)

  bright = (wise1 gt 5.*rms_wise1 or $
            wise2 gt 5.*rms_wise2 or $
            wise3 gt 5.*rms_wise3 or $
            wise4 gt 5.*rms_wise4)

  mask = grow_mask(pt_mask, constraint=bright)
  disp, wise1*(mask eq 0), max=1

  star_ind = $
     where(table.w1mpro lt 14 and $
           abs(table.w1mpro - table.w2mpro) lt 0.5 and $
;           (table.w1mpro - table.w3mpro) lt 0.5 and $
;           (table.w1mpro - table.w4mpro) lt 0.5 and $
           (table.ext_flg eq 0 or table.ext_flg eq 2 or table.ext_flg eq 4) $
           , star_ct)
  adxy, w1_hdr, table.ra[star_ind], table.dec[star_ind], star_x, star_y
  
  star_mask = finite(wise1)*0B
  sz = size(wise1)
  ximg = findgen(sz[1]) # (fltarr(sz[2])+1.)
  yimg = (fltarr(sz[1])+1.) # findgen(sz[2]) 
  for ii = 0, star_ct-1 do begin & $
     if (star_x[ii] lt 0 or star_y[ii] lt 0 or $
         star_x[ii] ge sz[1] or star_y[ii] ge sz[2]) then $
            continue & $
     peak_val = wise1[star_x[ii], star_y[ii]] & $
     rat = peak_val/rms_wise1 & $
     sig = (2.*alog(rat))^(0.5) & $
     this_rad = sig*(15./2.354/2.75) & $
     rimg = sqrt((ximg-star_x[ii])^2 + (yimg-star_y[ii])^2) & $
     star_mask[where(rimg lt this_rad)] = 1B & $
     endfor
  star_mask[star_x, star_y] = 1B
  
  disp, wise1*(mask eq 0)*(star_mask eq 0), max=1

end
