pro go_get_2mass_lga $
   , checkfirst=checkfirst $
   , dryrun=dryrun

  readcol, 'lga_list_2mass.txt', format='A' $
           , lga_names

  lga_names = strcompress(lga_names, /rem)
  n_lga = n_elements(lga_names)
  
  root_url = 'https://irsa.ipac.caltech.edu/data/LGA/images/'

  outdir = '../2mass_lga/'
  cd, outdir

  for ii = 0, n_lga-1 do begin

     for kk = 0, 2 do begin
        if kk eq 0 then band = 'j'
        if kk eq 1 then band = 'h'
        if kk eq 2 then band = 'k'

        fname = lga_names[ii]+'_mosaic_'+band+'.fits'
        if keyword_set(checkfirst) then $
           if file_test(fname) eq 1 then $
              continue
        
        wget_call = 'wget '+root_url+lga_names[ii]+'/'+fname

        print, wget_call

        if keyword_set(dryrun) then $
           continue

        spawn, wget_call
     endfor

  endfor

end
