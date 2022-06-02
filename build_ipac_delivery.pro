pro build_ipac_delivery $
   , rebuild=rebuild $
   , zapweights=zapweights

  tab = mrdfits('../measurements/delivery_index_gauss15.fits',1,h) 
  n_gals = n_elements(tab)

  if keyword_set(rebuild) then begin

     print, 'Making new directory'
     spawn, 'rm -rf ../ipac_delivery/'

     print, "Copying all data"
     spawn, 'cp -r ../delivery ../ipac_delivery'
     for ii = 0, n_gals-1  do begin

        counter, ii, n_gals, 'Galaxy'
        this_tab = tab[ii]
        this_name = strcompress(this_tab.pgc_name, /rem)

        print, "Removing rgrid files"
        spawn, 'rm -rf ../ipac_delivery/'+this_name+'_gauss15_rgrid.fits'
        spawn, 'rm -rf ../ipac_delivery/'+this_name+'_gauss7p5_rgrid.fits'
     endfor

  endif

  if keyword_set(zapweights) then begin

     jj = 0
     tot = total(tab.has_fuv eq 0)

     for ii = 0, n_gals-1 do begin

        this_tab = tab[ii]

        if this_tab.has_fuv eq 0 then begin
           jj += 1
           counter, jj, tot, 'Galaxy'
        
           wtname = '../ipac_delivery/'+$
                    strcompress(this_tab.pgc_name, /rem)+$
                    '_fuv_gauss15_weight.fits'
           spawn, 'rm -rf '+wtname

           wtname = '../ipac_delivery/'+$
                    strcompress(this_tab.pgc_name, /rem)+$
                    '_fuv_gauss7p5_weight.fits'
           spawn, 'rm -rf '+wtname
        endif

        if this_tab.has_nuv eq 0 then begin
           wtname = '../ipac_delivery/'+$
                    strcompress(this_tab.pgc_name, /rem)+$
                    '_nuv_gauss15_weight.fits'
           spawn, 'rm -rf '+wtname

           wtname = '../ipac_delivery/'+$
                    strcompress(this_tab.pgc_name, /rem)+$
                    '_nuv_gauss7p5_weight.fits'
           spawn, 'rm -rf '+wtname
        endif

     endfor

  endif

end
