pro add_photometry_to_index

  index = mrdfits('../measurements/delivery_index.fits',1,h)
  
  phot = mrdfits('../measurements/z0mgs_photometry.fits',1,h)
  n_pgc = n_elements(index)
  
  bands = ['FUV','NUV','WISE1','WISE2','WISE3','WISE4']
  n_bands = n_elements(bands)

  for ii = 0, n_pgc-1 do begin

     phot_ind = where(phot.pgc eq index[ii].pgc, phot_ct)
     if phot_ct eq 0 then continue
     
     phot_sub = phot[phot_ind]

     if index[ii].has_wise1 then begin
        ind = where(strcompress(phot_sub.band,/rem) eq 'WISE1', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                           
        index[ii].flux_wise1 = $
           phot_sub[ind].mean[1] + $
           (phot_sub[ind].med[2]-phot_sub[ind].med[1])
        index[ii].rms_flux_wise1 = phot_sub[ind].unc_stat[2]
        index[ii].std_flux_wise1 = phot_sub[ind].unc_conf[2]
     endif

     if index[ii].has_wise2 then begin
        ind = where(strcompress(phot_sub.band, /rem) eq 'WISE2', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                           
        index[ii].flux_wise2 = $
           phot_sub[ind].mean[1] + $
           (phot_sub[ind].med[2]-phot_sub[ind].med[1])
        index[ii].rms_flux_wise2 = phot_sub[ind].unc_stat[2]
        index[ii].std_flux_wise2 = phot_sub[ind].unc_conf[2]
     endif

     if index[ii].has_wise3 then begin
        ind = where(strcompress(phot_sub.band, /rem) eq 'WISE3', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                           
        index[ii].flux_wise3 = $
           phot_sub[ind].mean[1] + $
           (phot_sub[ind].med[2]-phot_sub[ind].med[1])
        index[ii].rms_flux_wise3 = phot_sub[ind].unc_stat[2]
        index[ii].std_flux_wise3 = phot_sub[ind].unc_conf[2]
     endif

     if index[ii].has_wise4 then begin
        ind = where(strcompress(phot_sub.band, /rem) eq 'WISE4', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                           
        index[ii].flux_wise4 = $
           phot_sub[ind].mean[1] + $
           (phot_sub[ind].med[2]-phot_sub[ind].med[1])
        index[ii].rms_flux_wise4 = phot_sub[ind].unc_stat[2]
        index[ii].std_flux_wise4 = phot_sub[ind].unc_conf[2]
     endif

     if index[ii].has_nuv then begin
        ind = where(strcompress(phot_sub.band, /rem) eq 'NUV', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                           
        index[ii].flux_nuv = $
           phot_sub[ind].mean[1] + $
           (phot_sub[ind].med[2]-phot_sub[ind].med[1])
        index[ii].rms_flux_nuv = phot_sub[ind].unc_stat[2]
        index[ii].std_flux_nuv = phot_sub[ind].unc_conf[2]
     endif

     if index[ii].has_fuv then begin
        ind = where(strcompress(phot_sub.band, /rem) eq 'FUV', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                           
        index[ii].flux_fuv = $
           phot_sub[ind].mean[1] + $
           (phot_sub[ind].med[2]-phot_sub[ind].med[1])
        index[ii].rms_flux_fuv = phot_sub[ind].unc_stat[2]
        index[ii].std_flux_fuv = phot_sub[ind].unc_conf[2]
     endif

  endfor

  mwrfits, index $
           , '../measurements/delivery_index.fits' $
           , /create

  stop

end
