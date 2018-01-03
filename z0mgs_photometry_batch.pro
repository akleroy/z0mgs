pro z0mgs_photometry_batch $
   , show=show


; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ISOPHOTAL FITTING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  if keyword_set(do_isophot) then begin

     if n_elements(tag) gt 0 then begin                
        all_data = gal_data(tag=tag)
     endif else begin
        all_data = gal_data(/all)
     endelse

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Isophotal fitting '

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        this_ind = where(all_data.pgc eq pgc_num[ii], ct)
        if ct eq 0 then continue
        this_dat = all_data[this_ind]

        band = 1
        infile = out_dir+pgc_name+'_w'+str(band)+'_gauss7p5.fits'
        z0mgs_isophot_fit $
           , outfile = '../measurements/'+pgc_name+'_isofit.txt' $
           , outimage = '../measurements/'+pgc_name+'_isofit.png' $
           , infile = infile $
           , gal_data = this_dat $
           , /show

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SLICES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

  if keyword_set(do_slice) then begin

     if n_elements(tag) gt 0 then begin                
        all_data = gal_data(tag=tag)
     endif else begin
        all_data = gal_data(/all)
     endelse

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Slice construction '

        pgc_name = pgc_list[ii]

        if n_elements(just) gt 0 then $
           if total(pgc_name eq just) eq 0 then $
              continue

        this_ind = where(all_data.pgc eq pgc_num[ii], ct)
        if ct eq 0 then continue
        this_dat = all_data[this_ind]

        band = 1
        infile = out_dir+pgc_name+'_w'+str(band)+'_gauss15.fits'
        z0mgs_image_slice $
           , outfile = '../measurements/'+pgc_name+'_slice.txt' $
           , outimage = '../measurements/'+pgc_name+'_slice.png' $
           , infile = infile $
           , gal_data = this_dat $
           , /show        

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PHOTOMETRY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GET THE LIST OF TARGETS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, 'survey_z0mgs.txt', format='A', z0mgs_name
  n_gals = n_elements(z0mgs_name)

  gals = gal_data(z0mgs_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  fuv_flux = fltarr(n_gals)*!values.f_nan
  nuv_flux = fltarr(n_gals)*!values.f_nan
  w1_flux = fltarr(n_gals)*!values.f_nan
  w2_flux = fltarr(n_gals)*!values.f_nan
  w3_flux = fltarr(n_gals)*!values.f_nan
  w4_flux = fltarr(n_gals)*!values.f_nan

  fuv_medflux = fltarr(n_gals)*!values.f_nan
  nuv_medflux = fltarr(n_gals)*!values.f_nan
  w1_medflux = fltarr(n_gals)*!values.f_nan
  w2_medflux = fltarr(n_gals)*!values.f_nan
  w3_medflux = fltarr(n_gals)*!values.f_nan
  w4_medflux = fltarr(n_gals)*!values.f_nan

  wise_dir = '../unwise/atlas/'
  galex_dir = '../galex/atlas/'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = 0, n_gals-1 do begin

     counter, ii, n_gals, 'Galaxy '

     for jj = 0, 5 do begin

        if jj eq 0 then begin
           map = readfits(wise_dir+z0mgs_name[ii]+'_w1_gauss15.fits' $
                          , hdr, /silent)
        endif
        if jj eq 1 then begin
           map = readfits(wise_dir+z0mgs_name[ii]+'_w2_gauss15.fits' $
                          , hdr, /silent)
        endif
        if jj eq 2 then begin
           map = readfits(wise_dir+z0mgs_name[ii]+'_w3_gauss15.fits' $
                          , hdr, /silent)
        endif
        if jj eq 3 then begin
           map = readfits(wise_dir+z0mgs_name[ii]+'_w4_gauss15.fits' $
                          , hdr, /silent)
        endif
        if jj eq 4 then begin
           map = readfits(galex_dir+z0mgs_name[ii]+'_nuv_gauss15.fits' $
                          , hdr, /silent)
        endif
        if jj eq 5 then begin
           map = readfits(galex_dir+z0mgs_name[ii]+'_fuv_gauss15.fits' $
                          , hdr, /silent)
        endif

        found = map[0] ne -1
        if found eq 0 then begin
           message, 'No file found for '+z0mgs_name[ii], /info
           continue
        endif
        
;       Convert to Jy/pix
        pix_sr = (sxpar(hdr,'CD1_1')*!dtor)^2
        if pix_sr eq 0 then $
           pix_sr = (sxpar(hdr,'CDELT1')*!dtor)^2
        map *= pix_sr*1d6
        
        incl = gals[ii].incl_deg
        if finite(incl) eq 0 then incl = 0.0
        if incl gt 80. then incl = 80.
        
        posang = gals[ii].posang_deg
        if finite(posang) eq 0 then begin
           posang = 0.0
           incl = 0.0
        endif

;       Work out galactocentric radius
        make_axes, hdr, ri=ri, di=di
        pos_vec = $
           [posang $
            , incl $
            ,gals[ii].ra_deg, gals[ii].dec_deg]
        if total(finite(pos_vec)) ne 4 then $
           continue
        deproject, ri, di, pos_vec $
                   , rgrid=rgrid

        r25 = gals[ii].r25_deg
        if finite(r25) eq 0 then r25 = 1./60.
        r25_grid = rgrid/r25

        mask = r25_grid lt 1.25
        if total(mask) eq 0 then stop
        
        if keyword_set(show) then begin
           if jj eq 0 then !p.multi = [0, 3, 2]
           disp, map, max=1d-4, min=0, title=z0mgs_name[ii]
           contour, mask, lev=[1], /overplot    
        endif 
        
        binsize = 7.5/3600./r25
        
        bins = bin_data(r25_grid, map $
                        , xmin=0.0, xmax=1.25, binsize=binsize)
        
        flux = total(map*mask, /nan)

        nbins = n_elements(bins)
        if nbins le 2 then begin
           medflux = flux
        endif else begin
           medflux = $
              total(bins[0:1].ymed*bins[0:1].counts) + $
              total(bins[1:(nbins-1)].ymed*bins[1:(nbins-1)].counts)
        endelse     

        if jj eq 0 then begin
           w1_flux[ii] = flux
           w1_medflux[ii] = medflux
        endif
        if jj eq 1 then begin
           w2_flux[ii] = flux
           w2_medflux[ii] = medflux
        endif
        if jj eq 2 then begin
           w3_flux[ii] = flux
           w3_medflux[ii] = medflux
        endif
        if jj eq 3 then begin
           w4_flux[ii] = flux
           w4_medflux[ii] = medflux
        endif
        if jj eq 4 then begin
           nuv_flux[ii] = flux
           nuv_medflux[ii] = medflux
        endif
        if jj eq 5 then begin
           fuv_flux[ii] = flux
           fuv_medflux[ii] = medflux
        endif

     endfor

  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  openw,1,'../measurements/z0mgs_photometry.txt'

  for ii = 0, n_gals-1 do begin

     line = z0mgs_name[ii]+' '
     line += string(w1_flux[ii])+' '
     line += string(w1_medflux[ii])+' '
     line += string(w2_flux[ii])+' '
     line += string(w2_medflux[ii])+' '
     line += string(w3_flux[ii])+' '
     line += string(w3_medflux[ii])+' '
     line += string(w4_flux[ii])+' '
     line += string(w4_medflux[ii])+' '
     line += string(nuv_flux[ii])+' '
     line += string(nuv_medflux[ii])+' '
     line += string(fuv_flux[ii])+' '
     line += string(fuv_medflux[ii])
     
     printf, 1, line

  endfor

  close, 1

end
