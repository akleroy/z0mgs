pro add_photometry_to_index

  index15 = mrdfits('../measurements/delivery_index_gauss15.fits',1,h15)
  index7p5 = mrdfits('../measurements/delivery_index_gauss7p5.fits',1,h7p5)
  
  phot15 = mrdfits('../measurements/z0mgs_photometry_gauss15.fits',1,h)
  phot7p5 = mrdfits('../measurements/z0mgs_photometry_gauss7p5.fits',1,h)
  n_pgc = n_elements(index15)
  
  bands = ['FUV','NUV','WISE1','WISE2','WISE3','WISE4']
  n_bands = n_elements(bands)

  for ii = 0, n_pgc-1 do begin

     phot_ind = where(phot15.pgc eq index15[ii].pgc, phot_ct)
     if phot_ct eq 0 then continue
     
     phot15_sub = phot15[phot_ind]
     phot7p5_sub = phot7p5[phot_ind]

     if index7p5[ii].has_wise1 then begin
        ind = where(strcompress(phot15_sub.band,/rem) eq 'WISE1', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                           
        this_flux = phot7p5_sub[ind].mean[1] + $
                    (phot7p5_sub[ind].med[2]-phot7p5_sub[ind].med[1])
        this_noise = phot7p5_sub[ind].unc_stat[2]
        this_conf = phot7p5_sub[ind].unc_conf[2]
        this_outer = phot7p5_sub[ind].mean[5] - phot7p5_sub[ind].mean[4] 

        index15[ii].flux_wise1 = this_flux
        index15[ii].rms_flux_wise1 = this_noise
        index15[ii].std_flux_wise1 = this_conf
        ;index15[ii].outer_flux_wise1 = this_outer

        index7p5[ii].flux_wise1 = this_flux
        index7p5[ii].rms_flux_wise1 = this_noise
        index7p5[ii].std_flux_wise1 = this_conf
        ;index7p5[ii].outer_flux_wise1 = this_outer

     endif

     if index7p5[ii].has_wise2 then begin
        ind = where(strcompress(phot15_sub.band, /rem) eq 'WISE2', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                           

        this_flux = phot7p5_sub[ind].mean[1] + $
                    (phot7p5_sub[ind].med[2]-phot7p5_sub[ind].med[1])
        this_noise = phot7p5_sub[ind].unc_stat[2]
        this_conf = phot7p5_sub[ind].unc_conf[2]
        this_outer = phot7p5_sub[ind].mean[5] - phot7p5_sub[ind].mean[4] 

        index15[ii].flux_wise2 = this_flux
        index15[ii].rms_flux_wise2 = this_noise
        index15[ii].std_flux_wise2 = this_conf
        ;index15[ii].outer_flux_wise2 = this_outer

        index7p5[ii].flux_wise2 = this_flux
        index7p5[ii].rms_flux_wise2 = this_noise
        index7p5[ii].std_flux_wise2 = this_conf
        ;index7p5[ii].outer_flux_wise2 = this_outer

     endif

     if index7p5[ii].has_wise3 then begin

        ind = where(strcompress(phot15_sub.band, /rem) eq 'WISE3', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                           

        this_flux = phot7p5_sub[ind].mean[1] + $
                    (phot7p5_sub[ind].med[2]-phot7p5_sub[ind].med[1])
        this_noise = phot7p5_sub[ind].unc_stat[2]
        this_conf = phot7p5_sub[ind].unc_conf[2]
        this_outer = phot7p5_sub[ind].mean[5] - phot7p5_sub[ind].mean[4] 

        index15[ii].flux_wise3 = this_flux
        index15[ii].rms_flux_wise3 = this_noise
        index15[ii].std_flux_wise3 = this_conf
        ;index15[ii].outer_flux_wise3 = this_outer

        index7p5[ii].flux_wise3 = this_flux
        index7p5[ii].rms_flux_wise3 = this_noise
        index7p5[ii].std_flux_wise3 = this_conf
        ;index7p5[ii].outer_flux_wise3 = this_outer

     endif

     if index15[ii].has_wise4 then begin

        ind = where(strcompress(phot15_sub.band, /rem) eq 'WISE4', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                         

        this_flux = phot15_sub[ind].mean[1] + $
                    (phot15_sub[ind].med[2]-phot15_sub[ind].med[1])
        this_noise = phot15_sub[ind].unc_stat[2]
        this_conf = phot15_sub[ind].unc_conf[2]
        this_outer = phot15_sub[ind].mean[5] - phot15_sub[ind].mean[4] 

        index15[ii].flux_wise4 = this_flux
        index15[ii].rms_flux_wise4 = this_noise
        index15[ii].std_flux_wise4 = this_conf
        ;index15[ii].outer_flux_wise4 = this_outer

        index7p5[ii].flux_wise4 = this_flux
        index7p5[ii].rms_flux_wise4 = this_noise
        index7p5[ii].std_flux_wise4 = this_conf
        ;index7p5[ii].outer_flux_wise4 = this_outer

     endif

     if index7p5[ii].has_nuv then begin
        ind = where(strcompress(phot15_sub.band, /rem) eq 'NUV', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                           

        this_flux = phot7p5_sub[ind].mean[1] + $
                    (phot7p5_sub[ind].med[2]-phot7p5_sub[ind].med[1])
        this_noise = phot7p5_sub[ind].unc_stat[2]
        this_conf = phot7p5_sub[ind].unc_conf[2]
        this_outer = phot7p5_sub[ind].mean[5] - phot7p5_sub[ind].mean[4] 

        index15[ii].flux_nuv = this_flux
        index15[ii].rms_flux_nuv = this_noise
        index15[ii].std_flux_nuv = this_conf
        ;index15[ii].outer_flux_nuv = this_outer

        index7p5[ii].flux_nuv = this_flux
        index7p5[ii].rms_flux_nuv = this_noise
        index7p5[ii].std_flux_nuv = this_conf
        ;index7p5[ii].outer_flux_nuv = this_outer

     endif

     if index7p5[ii].has_fuv then begin

        ind = where(strcompress(phot15_sub.band, /rem) eq 'FUV', ct)
        if ct eq 0 then begin
           print, "I'm missing a measurement I should have."
           stop
        endif                           

        this_flux = phot7p5_sub[ind].mean[1] + $
                    (phot7p5_sub[ind].med[2]-phot7p5_sub[ind].med[1])
        this_noise = phot7p5_sub[ind].unc_stat[2]
        this_conf = phot7p5_sub[ind].unc_conf[2]
        this_outer = phot7p5_sub[ind].mean[5] - phot7p5_sub[ind].mean[4] 

        index15[ii].flux_fuv = this_flux
        index15[ii].rms_flux_fuv = this_noise
        index15[ii].std_flux_fuv = this_conf
        ;index15[ii].outer_flux_fuv = this_outer

        index7p5[ii].flux_fuv = this_flux
        index7p5[ii].rms_flux_fuv = this_noise
        index7p5[ii].std_flux_fuv = this_conf
        ;index7p5[ii].outer_flux_fuv = this_outer

     endif

  endfor

  mwrfits, index15 $
           , '../measurements/delivery_index_gauss15.fits' $
           , /create

  mwrfits, index7p5 $
           , '../measurements/delivery_index_gauss7p5.fits' $
           , /create

  stop

end
