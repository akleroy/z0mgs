pro extract_gaia_training

;+
;
;Loop over our atlas at 7.5" resolution and extract a bunch of
;intensities associated with stars known from GAIA but outside the
;fiducial radius and not overlapping galaxies. Based on this build a
;vector of {g, fuv, nuv, w1, w2, w3} that can be used to reject or
;accept a "star-like" color on a galaxy.
;
;-  

  atlas_dir = '../delivery/'
  gaia_dir = '../stars/gaia/'

  tab = mrdfits('../measurements/delivery_index_gauss7p5.fits', 1, h)
  n_pgc = n_elements(tab)

  g_thresh = 19.0

  counter = 0L
  nstars = 1e7
  gmag = fltarr(nstars)*!values.f_nan
  w1 = fltarr(nstars)*!values.f_nan
  w2 = fltarr(nstars)*!values.f_nan
  w3 = fltarr(nstars)*!values.f_nan
  nuv = fltarr(nstars)*!values.f_nan
  fuv = fltarr(nstars)*!values.f_nan

  for ii = 0, n_pgc-1 do begin

     counter, counter, nstars, 'Star '

     if tab[ii].has_wise1 eq 0 then $
        continue

     pgc_name = strcompress(tab[ii].pgc_name, /rem)

     w1_im = readfits(atlas_dir+pgc_name+'_w1_gauss7p5.fits', w1_hdr, /silent)
     gal_mask = readfits(atlas_dir+pgc_name+'_gauss7p5_galaxies.fits', gal_hdr, /silent)
     rgrid = readfits(atlas_dir+pgc_name+'_gauss7p5_rgrid.fits', rad_hdr, /silent)

     gaia_file = gaia_dir+pgc_name+'_gaia.txt'
 
     readcol $
        , gaia_file $
        , delim=',' $
        , this_id, this_ra, this_dec, this_gmag $
        , format='L,D,D,D,X,X,X,X,X,X' $
        , count=lines $
        , /silent

     if lines eq 0 then $
        continue
     
     adxy, w1_hdr, this_ra, this_dec, this_x, this_y
    
     mask = rgrid gt sxpar(rad_hdr,'FIDRAD')*2.0 and $
            gal_mask eq 0

     ind = where(this_gmag lt g_thresh and $
                 mask[this_x,this_y] eq 1, ct)

     if ct eq 0 then $
        continue

     if counter+ct gt n_elements(w1) then $
        break

     x = this_x[ind]
     y = this_y[ind]
     
     lo = counter
     hi = (counter+ct-1)
     gmag[lo:hi] = this_gmag[ind]

     w1[lo:hi] = w1_im[x,y]

     if tab[ii].has_wise2 eq 1 then begin
        w2_im = readfits(atlas_dir+pgc_name+'_w2_gauss7p5.fits', w2_hdr, /silent)
        w2[lo:hi] = w2_im[x,y]
     endif

     if tab[ii].has_wise3 eq 1 then begin
        w3_im = readfits(atlas_dir+pgc_name+'_w3_gauss7p5.fits', w3_hdr, /silent)
        w3[lo:hi] = w3_im[x,y]
     endif

     if tab[ii].has_nuv eq 1 then begin
        nuv_im = readfits(atlas_dir+pgc_name+'_nuv_gauss7p5.fits', nuv_hdr, /silent)
        nuv[lo:hi] = nuv_im[x,y]
     endif

     if tab[ii].has_fuv eq 1 then begin
        fuv_im = readfits(atlas_dir+pgc_name+'_fuv_gauss7p5.fits', fuv_hdr, /silent)
        fuv[lo:hi] = fuv_im[x,y]
     endif

     counter += ct

  endfor

  save $
     , gmag, w1, w2, w3, w4, nuv, fuv $
     , file='../measurements/gaia_stars_gauss7p5.idl'

  stop

end
