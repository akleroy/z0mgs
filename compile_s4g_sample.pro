pro compile_s4g_sample

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
  gdata = gdata[ind] 
  uniq_ind = uniq(gdata.pgc,sort(gdata.pgc))
  gdata = gdata[uniq_ind]
  n_gals = n_elements(gdata)

  nan = !values.f_nan
  empty = $
     { $
     pgc_name:'' $
     , ra_deg: nan $
     , dec_deg: nan $
     , fuv: nan $
     , fuv_rms: nan $
     , nuv: nan $
     , nuv_rms: nan $
     , wise1: nan $
     , wise1_rms: nan $
     , wise2: nan $
     , wise2_rms: nan $
     , wise3: nan $
     , wise3_rms: nan $
     , wise4: nan $
     , wise4_rms: nan $
     , irac1: nan $
     , irac1_rms: nan $
     , irac2: nan $
     , irac2_rms: nan $
     , stellar: nan $
     , stellar_rms: nan $
     , nonstellar: nan $
     , nonstellar_rms: nan $
     }

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD S4G SAMPLING FOR COMPARISON
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  loadct, 0
  !p.multi=0
  
  for ii = 0, n_gals-1 do begin
     
     counter, ii, n_gals, ' out of '
     
     pgcname = 'PGC' + str(gdata[ii].pgc)
     
     wise1_fname = atlas_dir+pgcname+'_w1.fits'
     if file_test(wise1_fname) eq 0 then $
        continue
     wise1 = readfits(wise1_fname, /silent, wise1_hdr)

     wise2_fname = atlas_dir+pgcname+'_w2.fits'
     if file_test(wise2_fname) eq 0 then $
        wise2 = !values.f_nan*wise1 $
     else $
        wise2 = readfits(wise2_fname,/silent, wise2_hdr)

     wise3_fname = atlas_dir+pgcname+'_w3.fits'
     if file_test(wise3_fname) eq 0 then $
        wise3 = !values.f_nan*wise1 $
     else $
        wise3 = readfits(wise3_fname, /silent, wise3_hdr)

     wise4_fname = atlas_dir+pgcname+'_w4.fits'
     if file_test(wise4_fname) eq 0 then $
        wise4 = !values.f_nan*wise1 $
     else $
        wise4 = readfits(wise4_fname, /silent, wise4_hdr)

     fuv_fname = atlas_dir+pgcname+'_fuv.fits'
     if file_test(fuv_fname) eq 0 then $
        fuv = !values.f_nan*wise1 $
     else $
        fuv = readfits(fuv_fname, /silent, fuv_hdr)

     nuv_fname = atlas_dir+pgcname+'_nuv.fits'
     if file_test(nuv_fname) eq 0 then $
        nuv = !values.f_nan*wise1 $
     else $
        nuv = readfits(nuv_fname, /silent, nuv_hdr)

     irac1_fname = out_dir + pgcname + '_irac1_bksub.fits'     
     if file_test(irac1_fname) eq 0 then $
        irac1 = !values.f_nan*wise1 $
     else $
        irac1 = readfits(irac1_fname, irac1_hdr, /silent)

     irac2_fname = out_dir + pgcname + '_irac2_bksub.fits'     
     if file_test(irac2_fname) eq 0 then $
        irac2 = !values.f_nan*wise1 $
     else $
        irac2 = readfits(irac2_fname, irac2_hdr, /silent)

     stellar_fname = out_dir + pgcname + '_stellar_bksub.fits'     
     if file_test(stellar_fname) eq 0 then $
        stellar = !values.f_nan*wise1 $
     else $
        stellar = readfits(stellar_fname, stellar_hdr, /silent)

     nonstellar_fname = out_dir + pgcname + '_nonstellar_bksub.fits'     
     if file_test(nonstellar_fname) eq 0 then $
        nonstellar = !values.f_nan*wise1 $
     else $
        nonstellar = readfits(nonstellar_fname, nonstellar_hdr, /silent)
     
     maskfile = unwise_dir+pgcname+'_mask.fits'
     test = file_search(maskfile, count=ct)
     if ct eq 0 then $
        continue
     mask = readfits(maskfile, /silent, mask_hdr)
     hastrom, mask, mask_hdr, wise1_hdr, interp=0
     
     samp_ind = where(mask eq 10 , samp_ct)
     if samp_ct eq 0 then continue
     
     out_ind = where(mask ne 10 and finite(irac1) and finite(irac2), out_ct)

     make_axes, wise1_hdr, ri=ri, di=di

     this_samp = replicate(empty, samp_ct)
     this_samp.ra_deg = ri[samp_ind]
     this_samp.dec_deg = di[samp_ind]
     this_samp.fuv = fuv[samp_ind]
     this_samp.fuv_rms = sxpar(fuv_hdr, 'MADALL')
     this_samp.nuv = nuv[samp_ind]
     this_samp.nuv_rms = sxpar(nuv_hdr, 'MADALL')
     this_samp.wise1 = wise1[samp_ind]
     this_samp.wise1_rms = sxpar(wise1_hdr, 'MADALL')
     this_samp.wise2 = wise2[samp_ind]
     this_samp.wise2_rms = sxpar(wise2_hdr, 'MADALL')
     this_samp.wise3 = wise3[samp_ind]
     this_samp.wise3_rms = sxpar(wise3_hdr, 'MADALL')
     this_samp.wise4 = wise4[samp_ind]
     this_samp.wise4_rms = sxpar(wise4_hdr, 'MADALL')
     this_samp.irac1 = irac1[samp_ind]
     if out_ct gt 10 then this_samp.irac1_rms = mad(irac1[out_ind])
     this_samp.irac2 = irac2[samp_ind]
     if out_ct gt 10 then this_samp.irac2_rms = mad(irac2[out_ind])
     this_samp.stellar = stellar[samp_ind]
     this_samp.nonstellar = nonstellar[samp_ind]
     
     if n_elements(samp) eq 0 then $
        samp = this_samp $
     else $
        samp = [samp, this_samp]

  endfor     

  save, file='../measurements/s4g_samples.idl' $
        , samp
  stop

end
