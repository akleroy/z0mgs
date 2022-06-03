pro index_unwise_v2_stamps   

; Put together basic data on our second UNWISE atlas run

  nan = !values.f_nan
  empty = $
     { $
     z0mgs_name:'' $
     , pgc:-1 $
     , subsample:'' $
     , w1_fname:'' $
     , ra_ctr: nan $
     , dec_ctr: nan $     
     , blc_ra:nan $
     , blc_dec:nan $
     , trc_ra:nan $
     , trc_dec:nan $
     , use_w1:1B $
     , use_w2:1B $
     , use_w3:1B $
     , use_w4:1B $
     }

  skip_galaxy = []
  skip_w1 = []
  skip_w2 = []
  skip_w3 = []
  skip_w4 = []    

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
        nx = sxpar(hdr,'NAXIS1')
        ny = sxpar(hdr,'NAXIS2')
        xyad, hdr, 0, 0, blc_ra, blc_dec
        xyad, hdr, nx-1, ny-1, trc_ra, trc_dec        
        xyad, hdr, nx/2., ny/2., ra_ctr, dec_ctr
        
        dbase[ii].z0mgs_name = this_name
        dbase[ii].subsample = subsample
        if (subsample eq 'smallleda') or (subsample eq 'largeleda') then begin
           dbase[ii].pgc = long(strmid(this_name,3,strlen(this_name)-3))
        endif
        dbase[ii].w1_fname = this_file
        dbase[ii].ra_ctr = ra_ctr
        dbase[ii].dec_ctr = dec_ctr
        dbase[ii].blc_ra = blc_ra
        dbase[ii].blc_dec = blc_dec
        dbase[ii].trc_ra = trc_ra
        dbase[ii].trc_dec = trc_ra

        if total(this_name eq skip_galaxy) gt 0 then begin
           dbase[ii].use_w1 = 0B
           dbase[ii].use_w2 = 0B
           dbase[ii].use_w3 = 0B
           dbase[ii].use_w4 = 0B
        endif       

        if total(this_name eq skip_w1) gt 0 then begin
           dbase[ii].use_w1 = 0B
        endif       

        if total(this_name eq skip_w2) gt 0 then begin
           dbase[ii].use_w2 = 0B
        endif       

        if total(this_name eq skip_w3) gt 0 then begin
           dbase[ii].use_w3 = 0B
        endif       
        
        if total(this_name eq skip_w4) gt 0 then begin
           dbase[ii].use_w4 = 0B
        endif       
        
     endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Write to disk
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
     
     mwrfits, dbase, '../measurements/unwise_v2_index_'+subsample+'.fits', /create

  endfor
  
  stop
  
end
