pro add_mtol_to_index

  @constants.bat

  lsun_3p4 = 1.83d18
  nu_fuv = c/(154.d-9*1d2)
  nu_nuv = c/(231.d-9*1d2)
  nu_w1 = c/(3.4d-6*1d2)
  nu_w2 = c/(4.5d-6*1d2)
  nu_w3 = c/(12.d-6*1d2)
  nu_w4 = c/(22.d-6*1d2)
  
  index15 = mrdfits('../measurements/delivery_index_gauss15.fits',1,h15)
  index7p5 = mrdfits('../measurements/delivery_index_gauss7p5.fits',1,h7p5)

  for jj = 0, 1 do begin
     
     if jj eq 0 then this_index = index15
     if jj eq 1 then this_index = index7p5

     dat = gal_data(pgc=this_index.pgc)
     
     n_pgc = n_elements(this_index)

     has_w1 = finite(this_index.flux_wise1)
     has_w3 = finite(this_index.flux_wise3)
     has_nuv = finite(this_index.flux_nuv)
     
     nuv_lum = this_index.flux_nuv*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2
     w1_lum = this_index.flux_wise1*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2
     w3_lum = this_index.flux_wise3*1d-23*4.*!pi*(dat.dist_mpc*1d6*pc)^2    

     w3w1 = alog10(w3_lum/w1_lum)
     nuvw1 = alog10(nuv_lum/w1_lum) 
     lumw1 = alog10(w1_lum/lsun_3p4)

     have_all = where(has_w1 and has_w3 and has_nuv, all_ct)
     if all_ct gt 0 then begin
        mtol_all = $
           lookup_mtol(nuvw1=nuvw1[have_all] $
                       , w3w1=w3w1[have_all] $
                       , lumw1=lumw1[have_all] $
                       , unc=unc)
        this_index[have_all].mtol_w1 = mtol_all
        this_index[have_all].mtol_unc = unc
     endif

     just_wise = where(has_w1 and has_w3 and has_nuv eq 0, wise_ct)
     if wise_ct gt 0 then begin
        mtol_wise = $
           lookup_mtol(w3w1=w3w1[just_wise] $
                       , lumw1=lumw1[just_wise] $
                       , unc=unc) 
        this_index[just_wise].mtol_w1 = mtol_wise
        this_index[just_wise].mtol_unc = unc
     endif

     if jj eq 0 then index15 = this_index
     if jj eq 1 then index7p5 = this_index

  endfor
  
  mwrfits, index15 $
           , '../measurements/delivery_index_gauss15.fits' $
           , /create

  mwrfits, index7p5 $
           , '../measurements/delivery_index_gauss7p5.fits' $
           , /create  

  stop

end
