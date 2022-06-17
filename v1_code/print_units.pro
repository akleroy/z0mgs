pro print_units

  @constants.bat
  nu3p6 = c/3.6d-4
  nu3p4 = c/3.4d-4
  nu12 = c/12.d-4
  nu22 = c/22.d-4
  nu24 = c/24.d-4
  
  print, 'For 22um, 1MJy/sr is ', 1.*1d-23*1d6*nu22*4.*!pi*(1d6*pc)^2/(1d3)^2*10.^(-42.7), ' Msun/yr/kpc^2'
  print, 'For 24um, 1MJy/sr is ', 1.*1d-23*1d6*nu24*4.*!pi*(1d6*pc)^2/(1d3)^2*10.^(-42.7), ' Msun/yr/kpc^2'

  print, 'For 12um, 1MJy/sr is ', 1.*1d-23*1d6*nu12*4.*!pi*(1d6*pc)^2/(1d3)^2*10.^(-42.9), ' Msun/yr/kpc^2'

  nufuv = c/(154d-7)
  nunuv = c/(231d-7)
  print, 'For FUV and the KE12 conversion, 1 MJy/sr is ', 1.*1d-23*1d6*nufuv*4.*!pi*(1d6*pc)^2/(1d3)^2*10.^(-43.35), ' Msun/yr/kpc^2'
  print, 'For NUV and the KE12 conversion 1 MJy/sr is ', 1.*1d-23*1d6*nunuv*4.*!pi*(1d6*pc)^2/(1d3)^2*10.^(-43.17), ' Msun/yr/kpc^2'

  irac_zmag_jy = 280.9
  abs_sun_irac1_zmag = 3.24
  print, "For absolute ZMAG 3.24 mag:"
  lnu_sun_irac1 = 10.^(-3.24d/2.5)*280.9*4.*!pi*(10.d*pc)^2*1d-23
  print, "Specific luminosity of the Sun at 3.6 ... ", lnu_sun_irac1, " erg/s/Hz"
  print, "For which nu*Lnu is ... ", lnu_sun_irac1*nu3p6, " erg/s or ", lnu_sun_irac1*nu3p6/lsun, " Lsun"

  abs_sun_wise1_vega = 3.24
  abmag_sun_wise1 = abs_sun_wise1_vega + 2.699
  lnu_sun_wise1 = 10.^(-1.*abmag_sun_wise1/2.5)*3631.d*4.*!pi*(10.d*pc)^2*1d-23
  print, "Specific luminosity of the Sun at 3.4 ... ", lnu_sun_wise1, " erg/s/Hz"
  print, "For which nu*Lnu is ... ", lnu_sun_wise1*nu3p4, " erg/s or ", lnu_sun_wise1*nu3p4/lsun, " Lsun"
    
  print, "Then for a mass-to-light ratio of 0.5:"
  mtol = 0.5
  print, "For 3.6um, 1MJy/sr is ", 1.*1d-23*1d6*4.*!pi*(1d6*pc)^2/(1d6)^2/lnu_sun_irac1*mtol
  print, "For 3.4um, 1MJy/sr is ", 1.*1d-23*1d6*4.*!pi*(1d6*pc)^2/(1d6)^2/lnu_sun_wise1*mtol

;  dat = gal_data(/all)
;  lnu_w1 = dat.lum_w1/nu3p6
;  s4g_mtol = dat.s4g_mstar / (lnu_w1/lnu_sun_irac1)
;  vec = s4g_mtol[where(finite(s4g_mtol))]
;  vec = vec[sort(vec)]
;  n = n_elements(vec)
;  print, "S4G Querejeta mass-to-light ratio: "
;  print, vec[0.16*n], ' - ', vec[0.5*n], ' - ', vec[0.84*n] 

  stop

end
