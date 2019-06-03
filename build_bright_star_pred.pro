pro build_bright_star_pred $
   , pgcname=pgcname $
   , galdata=this_data $
   , infile=infile $
   , band=band $
   , res=res $
   , outfile=outfile $
   , star_ra = star_ra $
   , star_dec = star_dec $
   , star_km = star_km $
   , n_found = n_found $
   , ra_found = ra_found $
   , dec_found = dec_found $
   , km_found = km_found $
   , nopad = nopad $
   , center_tol = center_tol $
   , show=show $
   , pause=pause

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS AND READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  nan = !values.f_nan

  gaia_dir = '../stars/gaia/'

  if n_elements(res) eq 0 then  begin
     res = 'gauss7p5'
  endif
  
  valid_res = ['gauss7p5', 'gauss15']
  if total(res eq valid_res) eq 0 then begin
     print, "Defaulting to 7.5as Gaussian"
     res = 'gauss7p5'
  endif

  valid_bands = ['w1','w2','w3','w4','nuv','fuv']
  if total(band eq valid_bands) eq 0 then begin
     print, "Defaulting to WISE1"
     band = 'w1'
  endif
  
  if n_elements(star_ra) eq 0 then begin
     print, "Loading 2MASS stars."
     restore, '../measurements/2mass_stars.idl', /v
  endif

  readcol $
     , 'mag_to_intens.txt', comment='#' $
     , format='A,A,A,F' $
     , tab_mag, tab_band, tab_res, tab_coef $
     , /silent
  tab_mag = strcompress(tab_mag, /rem)
  tab_res = strcompress(tab_res, /rem)
  tab_band = strcompress(tab_band, /rem)
  
  gaia_ind = $
     where(tab_res eq res and $
           tab_band eq band and $
           tab_mag eq 'gaiag', gaia_ct)
  if gaia_ct eq 0 then begin
     print, "No match for GAIA."
     stop
  endif
  gaia_coef = (tab_coef[gaia_ind])[0]

  ks_ind = $
     where(tab_res eq res and $
           tab_band eq band and $
           tab_mag eq 'ks', ks_ct)
  if ks_ct eq 0 then begin
     print, "No match for 2MASS."
     stop
  endif
  ks_coef = (tab_coef[ks_ind])[0]

  if file_test(infile) eq 0 then begin
     message, 'Target file not found and map and header not supplied.', /info
     return
  endif     
  
  if n_elements(center_tol) eq 0 then begin
     center_tol = 3./3600.
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  map = readfits(infile, hdr, /silent)  
  stars = finite(map)*0.0
  sz = size(stars)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IDENTIFY 2MASS STARS TO BLANK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; PLACE ALL STARS IN PIXEL COORDS
  adxy, hdr, star_ra, star_dec, star_x, star_y

; IDENTIFY STARS IN THE IMAGE
  in_image = $
     (star_x ge 0) and (star_x lt sz[1]) and $
     (star_y ge 0) and (star_y lt sz[2])
  
  away_from_center = $
     sqrt(((star_ra - this_data.ra_deg)*cos(this_data.dec_deg*!dtor))^2 + $
          ((star_dec - this_data.dec_deg))^2) $
     gt center_tol

  star_ind = where(in_image and away_from_center, n_found)

  if n_found gt 0 then begin

     x_to_blank = round(star_x[star_ind])
     y_to_blank = round(star_y[star_ind])
     pred_intens = 10.^(-1.0d*star_km[star_ind]/2.5)*ks_coef

     stars[x_to_blank, y_to_blank] += pred_intens

;    GRAB THE FOUND STARS, TO PASS BACK OUT
     
     ra_found = star_ra[star_ind]
     dec_found = star_dec[star_ind]
     km_found = star_km[star_ind]

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IDENTIFY GAIA STARS TO BLANK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  gaia_file = gaia_dir+pgcname+'_gaia.txt'
  readcol $
     , gaia_file $
     , delim=',' $
     , gaia_id, gaia_ra, gaia_dec, gaia_gmag $
     , gaia_parallax,gaia_parallax_error $
     , gaia_pmra,gaia_pmra_error $
     , gaia_pmdec,gaia_pmdec_error $
     , format='L,D,D,D,D,D,D,D,D,D' $
     , count=lines $
     , /silent    

  gaia_s2n_cut = 3.5

  if lines gt 0 then begin

     pred_intens = 10.^(-1.0d*gaia_gmag/2.5)*gaia_coef

     in_milky_way = $
        (gaia_parallax gt gaia_s2n_cut*gaia_parallax_error) or $
        (gaia_pmra gt gaia_s2n_cut*gaia_pmra_error) or $
        (gaia_pmdec gt gaia_s2n_cut*gaia_pmdec_error)
     
     away_from_center = $
        sqrt(((gaia_ra - this_data.ra_deg)*cos(this_data.dec_deg*!dtor))^2 + $
             ((gaia_dec - this_data.dec_deg))^2) $
        gt center_tol
     
     adxy, hdr, gaia_ra, gaia_dec, gaia_x, gaia_y

     to_blank = $
        (gaia_x ge 0) and (gaia_x lt sz[1]) and $
        (gaia_y ge 0) and (gaia_y lt sz[2]) and $
        in_milky_way and $
        away_from_center
     
     blank_ind = where(to_blank, n_found)
     
     if n_found gt 0 then begin

        stars[round(gaia_x[blank_ind]), round(gaia_y[blank_ind])] += $
           pred_intens[blank_ind]
        
     endif

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BLUR OUT THE MAP TO THE TARGET RESOLUTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if res eq 'gauss7p5' then $
     target_res_as  = 7.5
  if res eq 'gauss15' then $
     target_res_as  = 15.

  conv_with_gauss $
     , data=stars, hdr=hdr $
     , start=0.0*[1,1,0] $
     , target=target_res_as*[1.,1,0] $
     , out_data=smooth_stars $
     , /nonorm $
     , /quiet
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if keyword_set(show) then begin

     !p.multi=[0,2,2]

     if total(finite(map)) lt 10 then $ 
        rms =0.0 $
     else $
        rms = mad(map[where(finite(map))])

     loadct, 0
     disp, map, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, smooth_stars, lev=[3.*rms], color=cgcolor('red'), /overplot

     loadct, 0
     disp, smooth_stars, max=5*rms, min=-5.*rms $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, smooth_stars, lev=[3.*rms], color=cgcolor('red'), /overplot

     loadct, 0
     disp, map-smooth_stars, max=5*rms+median(map), min=-5.*rms+median(map) $
           , /sq, xstyle=1, ystyle=1, reserve=5
     contour, smooth_stars, lev=[3.*rms], color=cgcolor('red'), /overplot
     

     if keyword_set(pause) then begin
        ch = get_kbrd(1)
     endif

     !p.multi=0

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if n_elements(outfile) gt 0 then begin
     writefits, outfile, smooth_stars, hdr
  endif

end
