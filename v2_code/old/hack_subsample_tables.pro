pro hack_subsample_tables

  tdir = '../../../measurements/'

  nan = !values.f_nan
  empty = { $
          name:'', $
          pgc:-1L, $
          ctr_ra: nan, $
          ctr_dec: nan, $
          incl_deg: nan, $
          posang_deg: nan, $
          rgal_deg: nan, $
          imsize: nan, $
          unwise_ra: nan, $
          unwise_dec: nan $
  }
  
  for ii = 0, 4 do begin

     if ii eq 0 then begin
        infile = tdir+'unwise_v2_index_localgroup.fits'
        outfile = tdir+'unwise_v2_index_localgroup_v2.fits'
     endif

     if ii eq 1 then begin
        infile = tdir+'unwise_v2_index_localvolume.fits'
        outfile = tdir+'unwise_v2_index_localvolume_v2.fits'
     endif

     if ii eq 2 then begin
        infile = tdir+'unwise_v2_index_largeleda.fits'
        outfile = tdir+'unwise_v2_index_largeleda_v2.fits'
     endif

     if ii eq 3 then begin
        infile = tdir+'unwise_v2_index_smallleda.fits'
        outfile = tdir+'unwise_v2_index_smalleda_v2.fits'
     endif

     if ii eq 4 then begin
        infile = tdir+'unwise_v2_index_manga.fits'
        outfile = tdir+'unwise_v2_index_manga_v2.fits'
     endif

     ; Read the table
     tab = mrdfits(infile, 1)

     ; Get the galaxy data
     dat = gal_data(pgc=tab.pgc, /full)

     ; Get the WISE 1 file
     fname = '../'+tab.W1_FNAME

     n_gal = n_elements(tab)    
     new_tab = replicate(empty, n_gal)
     
     new_tab.name = tab.z0mgs_name
     new_tab.pgc = tab.pgc

     new_tab.incl_deg = dat.incl_deg
     new_tab.posang_deg = dat.posang_deg

     new_tab.ctr_ra = dat.ra_deg
     new_tab.ctr_dec = dat.dec_deg
     
     for jj = 0, n_gal-1 do begin
        h = headfits(strcompress(fname[jj],/rem))
        nx = sxpar(h, 'NAXIS1')
        ny = sxpar(h, 'NAXIS2')        
        cdeltx = sxpar(h, 'CD1_1')
        cdelty = sxpar(h, 'CD2_2')        
        new_tab[jj].imsize = (nx*cdeltx) > (ny*cdelty)
        mid_x = (nx / 2.0 - 0.5)
        mid_y = (ny / 2.0 - 0.5)
        xyad, h, mid_x, mid_y, ra, dec
        new_tab[jj].unwise_ra = ra
        new_tab[jj].unwise_dec = dec
     endfor

     new_tab.rgal_deg = dat.r25_deg
     
     mwrfits, new_tab, outfile, /create

     ;fasthist, new_tab.ctr_ra - new_tab.unwise_ra
     
     stop
     
  endfor
  
end
