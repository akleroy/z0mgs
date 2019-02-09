pro compile_sample

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SETUP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  atlas_dir = '../delivery/'

  tab = mrdfits('../measurements/delivery_index_gauss15.fits', 1, h)
  n_gals = n_elements(tab)

  nan = !values.f_nan
  empty = $
     { $
     pgc:0L $
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
     }
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD S4G SAMPLING FOR COMPARISON
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  loadct, 0
  !p.multi=0

  sample = replicate(empty, 1e7)
  counter = 0L

  for ii = 0, n_gals-1 do begin
     
     if tab[ii].pgc eq 2557 then continue

     counter, ii, n_gals, ' out of '

     if tab[ii].has_wise1 eq 0 then continue
     
     pgc_name = strcompress(tab[ii].pgc_name, /rem)
     
     rgrid = readfits(atlas_dir+pgc_name+'_gauss15_rgrid.fits', rhdr, /silent)

     wise1_fname = atlas_dir+pgc_name+'_w1_gauss15.fits'
     wise1 = readfits(wise1_fname, wise1_hdr, /silent)
     wise1_rms = sxpar(wise1_hdr, 'RMS')

     if tab[ii].has_wise2 then begin
        wise2_fname = atlas_dir+pgc_name+'_w2_gauss15.fits'
        wise2 = readfits(wise2_fname, wise2_hdr, /silent)
        wise2_rms = sxpar(wise2_hdr, 'RMS')
     endif else begin
        wise2 = nan * wise1
        wise2_rms = nan
     endelse

     if tab[ii].has_wise3 then begin
        wise3_fname = atlas_dir+pgc_name+'_w3_gauss15.fits'
        wise3 = readfits(wise3_fname, wise3_hdr, /silent)
        wise3_rms = sxpar(wise3_hdr, 'RMS')
     endif else begin
        wise3 = nan * wise1
        wise3_rms = nan
     endelse

     if tab[ii].has_wise4 then begin
        wise4_fname = atlas_dir+pgc_name+'_w4_gauss15.fits'
        wise4 = readfits(wise4_fname, wise4_hdr, /silent)
        wise4_rms = sxpar(wise4_hdr, 'RMS')
     endif else begin
        wise4 = nan * wise1
        wise4_rms = nan
     endelse

     if tab[ii].has_nuv then begin
        nuv_fname = atlas_dir+pgc_name+'_nuv_gauss15.fits'
        nuv = readfits(nuv_fname, nuv_hdr, /silent)
        nuv_rms = sxpar(nuv_hdr, 'RMS')
     endif else begin
        nuv = nan * wise1
        nuv_rms = nan
     endelse

     if tab[ii].has_fuv then begin
        fuv_fname = atlas_dir+pgc_name+'_fuv_gauss15.fits'
        fuv = readfits(fuv_fname, fuv_hdr, /silent)
        fuv_rms = sxpar(fuv_hdr, 'RMS')
     endif else begin
        fuv = nan * wise1
        fuv_rms = nan
     endelse

     make_axes, rhdr, ri=ri, di=di
     samp_ind = where(rgrid lt sxpar(rhdr, 'FIDRAD'), samp_ct)

     sample[counter:(counter+samp_ct-1)].pgc = tab[ii].pgc
     sample[counter:(counter+samp_ct-1)].ra_deg = ri[samp_ind]
     sample[counter:(counter+samp_ct-1)].dec_deg = di[samp_ind]

     sample[counter:(counter+samp_ct-1)].fuv = fuv[samp_ind]
     sample[counter:(counter+samp_ct-1)].nuv = nuv[samp_ind]
     sample[counter:(counter+samp_ct-1)].wise1 = wise1[samp_ind]
     sample[counter:(counter+samp_ct-1)].wise2 = wise2[samp_ind]
     sample[counter:(counter+samp_ct-1)].wise3 = wise3[samp_ind]
     sample[counter:(counter+samp_ct-1)].wise4 = wise4[samp_ind]

     sample[counter:(counter+samp_ct-1)].fuv_rms = fuv_rms
     sample[counter:(counter+samp_ct-1)].nuv_rms = nuv_rms
     sample[counter:(counter+samp_ct-1)].wise1_rms = wise1_rms
     sample[counter:(counter+samp_ct-1)].wise2_rms = wise2_rms
     sample[counter:(counter+samp_ct-1)].wise3_rms = wise3_rms
     sample[counter:(counter+samp_ct-1)].wise4_rms = wise4_rms
     
     counter += samp_ct

  endfor     

  save, file='../measurements/z0mgs_samples.idl' $
        , sample

  stop

end
