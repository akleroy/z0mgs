pro v2_build_star_pred $
   , ra_ctr = ra_ctr $
   , dec_ctr = dec_ctr $
   , gaia_file = this_gaia_file $
   , ks_struct = ks_struct $
   , band=this_band $
   , infile=infile $
   , pred_flux_file=star_flux_file $
   , psf_file=psf_file $
   , pred_image_file=native_res_pred_file $
   , show=show $
   , pause=pause

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Initialize
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  nan = !values.f_nan
  
  if n_elements(this_pgc) eq 0 then begin
     this_pgc = -1
  endif
  
  if n_elements(infile) eq 0 then begin
     print, "infile required."
     return
  endif
  
  if file_test(infile) eq 0 then begin
     print, infile, " not found."
     return
  endif

  if n_elements(center_tol) eq 0 then begin
     center_tol = 3./3600.
  endif  
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Read the map and initialize the prediction
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  map = readfits(infile, hdr)
  pred_flux = finite(map)*0.0d
  sz = size(map)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Load the 2MASS bright stars  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if n_elements(ks_struct) eq 0 then begin
     
     print, "Loading 2MASS stars."
     restore, '../../measurements/2mass_stars.idl', /v
     ks_empty = { $
                name: '', $
                ra: nan, $
                dec: nan, $
                ks_mag: nan $
                }
     n_ks = n_elements(star_ra)
     ks_struct = replicate(ks_empty, n_ks)
     ks_struct.name = star_name
     ks_struct.ra = star_ra
     ks_struct.dec = star_dec
     ks_struct.ks_mag = star_km

     tab_2mass_file = $
        '../../measurements/tab_2mass_stars.fits'     
     if file_test(tab_2mass_file) eq 0 then begin
        print, "Writing 2MASS table"
        mwrfits, ks_struct, tab_2mass_file, /create
     endif
     
  endif  
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Load the GAIA file
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%    

  if file_test(this_gaia_file) then begin

     gaia_tab = mrdfits(this_gaia_file, 1, gaia_hdr)
     gaia_found = 1B
     
  endif else begin

     gaia_found = 0B
     
  endelse
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Get the conversions between magnitude and flux
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%    

  fid_res = 'gauss7p5'
  beam_area_sr = (7.5/3600.*!dtor/2.0)^2*!pi/alog(2)
  coef_intens_to_flux = beam_area_sr*1d6
  
  readcol $
     , '../../measurements/mag_to_intens.txt', comment='#' $
     , format='A,A,A,F' $
     , tab_mag, tab_band, tab_res, tab_coef $
     , /silent
  tab_mag = strcompress(tab_mag, /rem)
  tab_res = strcompress(tab_res, /rem)
  tab_band = strcompress(tab_band, /rem)
  
  gaia_ind = $
     where(tab_res eq fid_res and $
           tab_band eq this_band and $
           tab_mag eq 'gaiag', gaia_ct)
  if gaia_ct eq 0 then begin
     print, "No match for GAIA."
     stop
  endif
  gaia_coef = (tab_coef[gaia_ind])[0]*coef_intens_to_flux

  ks_ind = $
     where(tab_res eq fid_res and $
           tab_band eq this_band and $
           tab_mag eq 'ks', ks_ct)
  if ks_ct eq 0 then begin
     print, "No match for 2MASS."
     stop
  endif
  ks_coef = (tab_coef[ks_ind])[0]*coef_intens_to_flux

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Identify 2MASS stars in the image  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; PLACE ALL STARS IN PIXEL COORDS
  adxy, hdr $
        , ks_struct.ra, ks_struct.dec $
        , ks_x, ks_y

; IDENTIFY STARS IN THE IMAGE
  in_image = $
     (ks_x ge 0) and (ks_x lt sz[1]) and $
     (ks_y ge 0) and (ks_y lt sz[2])

; EXCLUDE STARS NEAR THE CENTER  
  away_from_center = $
     sqrt(((ks_struct.ra - ra_ctr)*cos(dec_ctr*!dtor))^2 + $
          ((ks_struct.dec - dec_ctr))^2) $
     gt center_tol

  ks_ind = where(in_image and away_from_center, n_ks_found)

  if n_ks_found gt 0 then begin

     x_to_blank = round(ks_x[ks_ind])
     y_to_blank = round(ks_y[ks_ind])

     this_pred_flux = 10.^(-1.0d*ks_struct[ks_ind].ks_mag/2.5)*ks_coef

     pred_flux[x_to_blank, y_to_blank] += this_pred_flux

  endif  
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Identify GAIA stars in the image  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(gaia_found) then begin

     gaia_s2n_cut = 3.5

     this_pred_flux = 10.^(-1.0d*gaia_tab.phot_g_mean_mag/2.5)*gaia_coef
     
     in_milky_way = $
        (gaia_tab.parallax ge gaia_s2n_cut*gaia_tab.parallax_error) or $
        (gaia_tab.pmra ge gaia_s2n_cut*gaia_tab.pmra_error) or $
        (gaia_tab.pmdec ge gaia_s2n_cut*gaia_tab.pmdec_error)
     
     away_from_center = $
        sqrt(((gaia_tab.ra - ra_ctr)*cos(dec_ctr*!dtor))^2 + $
             ((gaia_tab.dec - dec_ctr))^2) $
        gt center_tol
     
     adxy, hdr, gaia_tab.ra, gaia_tab.dec, gaia_x, gaia_y

     to_blank = $
        (gaia_x ge 0) and (gaia_x lt sz[1]) and $
        (gaia_y ge 0) and (gaia_y lt sz[2]) and $
        in_milky_way and $
        away_from_center and $
        finite(this_pred_flux)
     
     blank_ind = where(to_blank, n_found)
     
     if n_found gt 0 then begin

        pred_flux[round(gaia_x[blank_ind]), round(gaia_y[blank_ind])] += $
           this_pred_flux[blank_ind]
        
     endif     
     
  endif
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
; Write flux map to disk
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if n_elements(star_flux_file) ne 0 then begin

     flux_hdr = hdr
     sxdelpar, flux_hdr, 'BMAJ'
     sxdelpar, flux_hdr, 'BMIN'
     sxdelpar, flux_hdr, 'BPA'     
     sxaddpar, flux_hdr, 'BUNIT', 'JY/PIXEL'
     writefits, star_flux_file, pred_flux, flux_hdr
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Convolve to native resolution  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if file_test(psf_file) then begin

     xyad, hdr, 0, 0, ra0, dec0
     xyad, hdr, 0, 1, ra1, dec1     
     
     rad_per_pix = sphdist(ra0, dec0, ra1, dec1, /deg)*!dtor
     sr_per_pix = rad_per_pix^2     
     intens_map = pred_flux / sr_per_pix / 1e6

     native_hdr = hdr
     sxdelpar, native_hdr, 'BMAJ'
     sxdelpar, native_hdr, 'BMIN'
     sxdelpar, native_hdr, 'BPA'     
     sxaddpar, native_hdr, 'BUNIT', 'MJy/sr'

     scratch_file = '../../scratch/scratch_intens_map.fits'
     writefits, scratch_file, intens_map, native_hdr
     
     convolve_image $
        , kernel=psf_file $
        , image=scratch_file $
        , outfile=native_res_pred_file
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Show the results  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(show) then begin
     
     !p.multi=[0,2,2]

     intens_map = readfits(native_res_pred_file, intens_hdr)
     
     if total(finite(map)) lt 10 then $ 
        rms = 0.0 $
     else $
        rms = mad(map[where(finite(map))])
     
     loadct, 0
     disp, map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, intens_map, lev=[1.*rms], color=cgcolor('red'), /overplot
     
     loadct, 0
     disp, intens_map, max=5*rms, min=-5.*rms $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, intens_map, lev=[1.*rms], color=cgcolor('red'), /overplot
     
     loadct, 0
     disp, map-intens_map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, intens_map, lev=[1.*rms], color=cgcolor('red'), /overplot
     

     if keyword_set(pause) then begin
        ch = get_kbrd(1)
     endif
     
     !p.multi=0
     
  endif

  
end
