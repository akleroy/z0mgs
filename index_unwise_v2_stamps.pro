pro index_unwise_v2_stamps   

; Put together basic data on our second UNWISE atlas run

  nan = !values.f_nan
  empty = $
     { $
     'z0mgs_name':'' $
     , 'pgc':-1 $
     , 'tlc_ra':nan $
     , 'tlc_dec':nan $
     , 'brc_ra':nan $
     , 'brc_dec':nan $
     , 'has_w1':0B $
     , 'has_w2':0B $
     , 'has_w3':0B $
     , 'has_w4':0B $     
  }
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Karachentsev catalog
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  kdir = '../orig_data/unwise_v2/v2_karachentsev/'
  
  kflist = file_search(kdir+'*/*/unwise-*-w1-img-m.fits', count=kct)
  ; check - kflist_w4 = file_search(kdir+'*/*/unwise-*-w4-img-m.fits', count=kct_w4)

  kdbase = replicate(empty, kct)
  for ii = 0, kct-1 do begin
     this_file = kflist[ii]
     stop
  endfor

  mwrfits, kdbase, '../measurements/unwise_v2_index_karachentsev.fits', /create
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Large LEDA galaxies
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  ldir = '../orig_data/unwise_v2/v2_largegals/'

  lflist = file_search(ldir+'*/unwise-*-w1-img-m.fits', count=lct)  

  ldbase = replicate(empty, lct)
  for ii = 0, lct-1 do begin
     this_file = lflist[ii]
     stop
  endfor

  mwrfits, ldbase, '../measurements/unwise_v2_index_largegals.fits', /create
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; Small LEDA galaxies
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  sdir =  '../orig_data/unwise_v2/v2_smallgals/'

  sflist = file_search(sdir+'*/*/unwise-*-w1-img-m.fits', count=sct)  
  sdbase = replicate(empty, sct)
  for ii = 0, sct-1 do begin
     this_file = sflist[ii]
     stop
  endfor

  mwrfits, sdbase, '../measurements/unwise_v2_index_smallgals.fits', /create
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MANGA dir
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  mdir = '../orig_data/unwise_v2/v2_unmanga/'
  
  mflist = file_search(sdir+'*/unwise-*-w1-img-m.fits', count=mct)  
  mdbase = replicate(empty, mct)
  for ii = 0, mct-1 do begin
     this_file = mflist[ii]
     stop
  endfor
  
  mwrfits, mdbase, '../measurements/unwise_v2_index_manga.fits', /create

  stop
  
end
