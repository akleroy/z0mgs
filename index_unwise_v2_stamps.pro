pro index_unwise_v2_stamps   

; Put together basic data on our second UNWISE atlas run

  nan = !values.f_nan
  empty = $
     { $
     z0mgs_name:'' $
     , pgc:-1 $
     , subsample:'' $
     , ra_ctr: nan $
     , dec_ctr: nan $     
     , blc_ra:nan $
     , blc_dec:nan $
     , trc_ra:nan $
     , trc_dec:nan $
     }

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Loop over four datasets
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
  
  for vv = 0, 3 do begin

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Karachentsev catalog
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

     if vv eq 0 then begin
        
        dir = '../orig_data/unwise_v2/v2_karachentsev/'        
        flist = file_search(dir+'*/*/unwise-*-w1-img-m.fits', count=ct)
        subsample = 'localvolume'
        
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Large LEDA galaxies
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

     if vv eq 1 then begin
        
        dir = '../orig_data/unwise_v2/v2_largegals/'        
        flist = file_search(dir+'*/unwise-*-w1-img-m.fits', count=ct)  
        subsample = 'largeleda'
        
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Small LEDA galaxies
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

     if vv eq 2 then begin
        
        dir = '../orig_data/unwise_v2/v2_smallgals/'
        flist = file_search(dir+'*/*/unwise-*-w1-img-m.fits', count=ct)
        subsample = 'smallleda'
        
     endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MANGA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

     if vv eq 3 then begin
        
        dir = '../orig_data/unwise_v2/v2_unmanga/'
        flist = file_search(dir+'*/unwise-*-w1-img-m.fits', count=ct)
        subsample = 'manga'
        
     endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Loop over files
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
     
     dbase = replicate(empty, ct)

     for ii = 0, ct-1 do begin
        counter, ii, ct, 'Indexing galaxy '
        
        this_file = flist[ii]
        start_pos = strpos(this_file, 'unwise-')+strlen('unwise-')
        stop_pos = strlen(this_file) - strlen('-w1-img-m.fits')
        this_name = strmid(this_file,start_pos,stop_pos-start_pos)

        hdr = headfits(this_file)
        make_axes, hdr, ri=ri, di=di
        sz = size(ri)
        
        dbase[ii].z0mgs_name = this_name
        dbase[ii].subsample = subsample
        if (subsample eq 'smallleda') or (subsample eq 'largeleda') then begin
           dbase[ii].pgc = long(strmid(this_name,3,strlen(this_name)-3))
        endif
        
        dbase[ii].ra_ctr = ri[sz[1]/2,sz[2]/2]
        dbase[ii].dec_ctr = di[sz[1]/2,sz[2]/2]
        dbase[ii].blc_ra = ri[0,0]
        dbase[ii].blc_dec = di[0,0]
        dbase[ii].trc_ra = ri[sz[1]-1,sz[2]-1]
        dbase[ii].trc_dec = di[sz[1]-1,sz[2]-1]
        
     endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Write to disk
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
     
     mwrfits, dbase, '../measurements/unwise_v2_index_'+subsample+'.fits', /create

  endfor
  
  stop
  
end
