pro fetch_all_unwise_tiles $
   , start=start $
   , stop=stop $
   , out_dir=out_dir $
   , in_file=in_file

; Fetch the all GALEX tiles based off of CAS job output.

  if n_elements(out_dir) eq 0 then $
     out_dir = '/data/tycho/0/leroy.42/ellohess/data/unwise/tiles/data/'
  
  if n_elements(start) eq 0 then $
     start = 0L
  
  if n_elements(stop) eq 0 then $
     stop = 0L

; Read the input file

  tab = mrdfits('tiles.fits',1,hdr)

; Change diretory and set up the command

  full_ext_list = ["-w1-frames.fits" $
                   , "-w1-img-m.fits" $
                   , "-w1-img-u.fits" $
                   , "-w1-invvar-m.fits.gz" $
                   , "-w1-invvar-u.fits.gz" $
                   , "-w1-mask.tgz" $
                   , "-w1-n-m.fits.gz" $
                   , "-w1-n-u.fits.gz" $
                   , "-w1-std-m.fits.gz" $
                   , "-w1-std-u.fits.gz" $
                   , "-w2-frames.fits" $
                   , "-w2-img-m.fits" $
                   , "-w2-img-u.fits" $
                   , "-w2-invvar-m.fits.gz" $
                   , "-w2-invvar-u.fits.gz" $
                   , "-w2-mask.tgz" $
                   , "-w2-n-m.fits.gz" $
                   , "-w2-n-u.fits.gz" $
                   , "-w2-std-m.fits.gz" $
                   , "-w2-std-u.fits.gz" $
                   , "-w3-frames.fits" $
                   , "-w3-img-m.fits" $
                   , "-w3-img-u.fits" $ 
                   , "-w3-invvar-m.fits.gz" $
                   , "-w3-invvar-u.fits.gz" $
                   , "-w3-mask.tgz" $
                   , "-w3-n-m.fits.gz" $
                   , "-w3-n-u.fits.gz" $ 
                   , "-w3-std-m.fits.gz" $
                   , "-w3-std-u.fits.gz" $
                   , "-w4-frames.fits" $
                   , "-w4-img-m.fits" $
                   , "-w4-img-u.fits" $
                   , "-w4-invvar-m.fits.gz" $
                   , "-w4-invvar-u.fits.gz" $
                   , "-w4-mask.tgz" $
                   , "-w4-n-m.fits.gz" $
                   , "-w4-n-u.fits.gz" $
                   , "-w4-std-m.fits.gz" $
                   , "-w4-std-u.fits.gz"]

  ext_list = $
     ["-w1-img-m.fits" $
      , "-w2-img-m.fits" $
      , "-w3-img-m.fits" $
      , "-w4-img-m.fits"]

  cd, out_dir

  for ii = long(start), long(stop) do begin

     if ii gt n_elements(tab) then continue
   
     print, "-----------------------------------------"
     print, "I am pulling tile "+str(ii)+" and stopping at "+str(stop)
     print, "-----------------------------------------"

     fname = 'http://unwise.me/data/'
     fname += strmid(tab[ii].coadd_id,0,3)+'/'
     fname += tab[ii].coadd_id+'/unwise-'+tab[ii].coadd_id

     for kk = 0, n_elements(ext_list)-1 do begin
        this_fname = fname + ext_list[kk]
        line = 'wget -nd "'+this_fname+'"'        
        spawn, line
     endfor

  endfor

end
