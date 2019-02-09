pro nearby_3color $
   , rebin=do_rebin $
   , mask=do_mask
  
  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  plot_dir = '../plots/'
  atlas_dir = '../delivery/'

  build_galaxy_list $
     , in_dir = in_dir $
     , tag=tag $
     , just=just $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = gal_data $
     , start = start_num $
     , stop = stop_num $
     , exclude = ['PGC17223']
  n_pgc = n_elements(pgc_list)
  sort_ind = sort((gal_data.gb_deg))
  gal_data = gal_data[sort_ind]

  cutoff = 10.*15./3600.
  ind = where(gal_data.r25_deg gt cutoff, ct)
  ngal = n_elements(ind)
  
  good_name = get_good_name(gal_data[ind].pgc)

  row = 0
  nrows = 10
  deltarow = 1.0/(nrows*1.0)

  column = 0
  ncolumns = 8
  deltacolumn = 1.0/(ncolumns*1.0)  
  
  page = 0

  for ii=0, ngal-1 do BEGIN

     if ii mod (nrows*ncolumns) eq 0 then begin

        if page gt 0 then begin
           ps, /xw  
           spawn, 'evince '+psfile+' &'
           spawn, 'convert -density 300x300 '+psfile+' '+pnfile
           !p.multi=0
        endif

        page += 1        
        psfile='../plots/test_3color_'+str(page)+'.eps'
        pnfile='../plots/test_3color_'+str(page)+'.png'
        ps, /def, file=psfile,/encapsulated,/portrait,/color,$
            xsize=8,ysize=10,/ps
        row = 0
        column = 0

     endif

     this_pgc = gal_data[ind[ii]].pgc
     print, 'PGC'+str(this_pgc)     
     w1_file = atlas_dir+'PGC'+str(this_pgc)+'_w1_gauss7p5.fits'
     w3_file = atlas_dir+'PGC'+str(this_pgc)+'_w3_gauss7p5.fits'
     nuv_file = atlas_dir+'PGC'+str(this_pgc)+'_nuv_gauss7p5.fits'
     if file_test(w1_file) eq 0 then begin
        print, "Did not find "+w1_file
        continue
     endif

     w1 = readfits(w1_file, w1_hdr)
     sz = size(w1, /dim)
     if n_elements(sz) eq 1 then begin
        print, "Header only for "+w1_file
        continue
     endif
     w3 = readfits(w3_file, w3_hdr)
          
     if file_test(nuv_file) then begin
        nuv = readfits(nuv_file,nuv_hdr)
        has_nuv = 1B
     endif else begin
        nuv = 0.0*w1+1d-10
        nuv_hdr = w1_hdr
        has_nuv = 0B
     endelse

     if keyword_set(do_mask) then begin
        galmask = atlas_dir+'PGC'+str(this_pgc)+'_gauss7p5_galaxies.fits'
        rgrid = readfits(atlas_dir+'PGC'+str(this_pgc)+'_gauss7p5_rgrid.fits', rhdr)

        w1_bright = atlas_dir+'PGC'+str(this_pgc)+'_w1_gauss7p5_bright_stars.fits'
        w3_bright = atlas_dir+'PGC'+str(this_pgc)+'_w3_gauss7p5_bright_stars.fits'
        if has_nuv then nuv_bright = atlas_dir+'PGC'+str(this_pgc)+'_nuv_gauss7p5_bright_stars.fits'

        w1_found = atlas_dir+'PGC'+str(this_pgc)+'_w1_gauss7p5_found_stars.fits'
        w3_found = atlas_dir+'PGC'+str(this_pgc)+'_w3_gauss7p5_found_stars.fits'
        if has_nuv then nuv_found = atlas_dir+'PGC'+str(this_pgc)+'_nuv_gauss7p5_found_stars.fits'        

        w1_mask = finite(w1) eq 0 
        w3_mask = finite(w3) eq 0
        if has_nuv then nuv_mask = finite(nuv) eq 0
        if file_test(galmask) then begin
           galmask = readfits(galmask)
           w1_mask = w1_mask or galmask
           w3_mask = w3_mask or galmask
           if has_nuv then nuv_mask = nuv_mask or galmask
        endif

        if file_test(w1_bright) then w1_mask = w1_mask or (readfits(w1_bright))
        if file_test(w3_bright) then w3_mask = w3_mask or (readfits(w3_bright))
        if has_nuv then $
           if file_test(nuv_bright) then nuv_mask = nuv_mask or (readfits(nuv_bright))

        if file_test(w1_found) then w1_mask = w1_mask or (readfits(w1_found))
        if file_test(w3_found) then w3_mask = w3_mask or (readfits(w3_found))
        if has_nuv then $
           if file_test(nuv_found) then nuv_mask = nuv_mask or (readfits(nuv_found))

        bad_w1 = where((w1_mask mod 10) eq 1 and rgrid gt sxpar(rhdr,'FIDRAD'), bad_ct)
        if bad_ct gt 0 then w1[bad_w1] = !values.f_nan
        bad_w3 = where((w3_mask mod 10) eq 1 and rgrid gt sxpar(rhdr,'FIDRAD'), bad_ct)
        if bad_ct gt 0 then w3[bad_w3] = !values.f_nan
        if has_nuv then begin
           bad_nuv = where((nuv_mask mod 10) eq 1 and rgrid gt sxpar(rhdr,'FIDRAD'), bad_ct)
           if bad_ct gt 0 then nuv[bad_nuv] = !values.f_nan
        endif

     endif

;     if (total(finite(nuv)) eq 0) or $
;        (total(size(nuv,/dimen) ne size(w1,/dimen)) ne 0) then begin
;     endif

     sz = size(w1,/dimen)
     xcen = sz[0]/2
     ycen = sz[1]/2

     delta = gal_data[ind[ii]].r25_deg*3600./2.75*1.25
     
     hextract,w1,w1_hdr,xcen-delta,xcen+delta,ycen-delta,ycen+delta
     hextract,w3,w3_hdr,xcen-delta,xcen+delta,ycen-delta,ycen+delta
     hextract,nuv,nuv_hdr,xcen-delta,xcen+delta,ycen-delta,ycen+delta
     
     sz = size(w1, /dimen)
     if keyword_set(do_rebin) then begin
        hrebin, w1, w1_hdr, outsize=[sz[0]*4, sz[0]*4]
        hrebin, w3, w3_hdr, outsize=[sz[0]*4, sz[0]*4]
        hrebin, nuv, nuv_hdr, outsize=[sz[0]*4, sz[0]*4]
     endif

     make_axes, w1_hdr, ra=ra, da=da, ri=ri, di=di

     loadct,0
     
     pos = [column*deltacolumn, 1.0-(row+1.0)*deltarow, $
            (column+1.0)*deltacolumn, 1.0-((row)*deltarow)]
     column += 1
     if column ge ncolumns then begin
        column = 0
        row += 1
     endif

     minval = [-2, -2, -3.0]
     maxval = [1.5, 1.5, -0.25]

     disp3,alog10(w3),alog10(w1),alog10(nuv), ra, da $
           , /sq,min=minval,max=maxval $ 
           , charthick=3,thick=3,xthick=3,ythick=3 $
           , xstyle=4, ystyle=4, position = pos $
           , /noerase
     pgcstring = 'PGC'+strcompress(str(gal_data[ind[ii]].pgc), /rem)
     name = strcompress(good_name[ii], /rem)
     al_legend, /top, /right, box=0, clear=0, lines=-99 $
                , [pgcstring,name] $
                , textcolor=cgcolor('white') $
                , charsize=0.5, charthick=1.0

     nuv_string = ''
     if has_nuv eq 0 then nuv_string='No GALEX'
     lat_string = '!8b!6 = '+strcompress(string(gal_data[ind[ii]].gb_deg,format='(3I)'),/rem)
     al_legend, /bottom, /left, box=0, clear=0, lines=-99 $
                , [nuv_string, lat_string] $
                , textcolor=cgcolor('white') $
                , charsize=0.5, charthick=1.0

  endfor

  ps, /xw  
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile
  !p.multi=0

  stop
end
