pro show_big_targets

  
  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  out_dir = '../unwise/atlas/'

; BUILD A LIST OF PGC GALAXIES THAT WE HAVE RIGHT NOW
  flist = file_search(in_dir+'PGC*', count=file_ct)
  pgc_list = strarr(file_ct)
  pgc_num = lonarr(file_ct)
  for ii = 0, file_ct-1 do begin
     pgc_list[ii] = strmid(flist[ii],strlen(in_dir),strlen(flist[ii]))
     pgc_num[ii] = long(strmid(pgc_list[ii],3,strlen(pgc_list[ii])-3))
  endfor
  n_pgc = n_elements(pgc_list)

  all_data = gal_data(/all)

  out_dir = '../unwise/atlas/'

  for ii = 0, n_pgc-1 do begin

     pgc_name = pgc_list[ii]
    
     this_dat = all_data[where(all_data.pgc eq pgc_num[ii])]

     if this_dat.r25_deg lt 5./60. or finite(this_dat.r25_deg) eq 0 then continue

     !p.multi=[0,2,2]

     loadct, 33

     print, this_dat.name, this_dat.r25_deg

     for band = 1, 4 do begin
        
        infile = out_dir+pgc_name+'_w'+str(band)+'_bksub.fits'
        test = file_search(infile, count=ct)
        if ct eq 0 then $
           continue
        map = readfits(infile, hdr, /silent)
        rms = sxpar(hdr, 'NOISE')
        disp, map, /sq, min=-3.*rms, max=3.*rms
        
     endfor
     
     ch = get_kbrd(1)

     !p.multi=0

  endfor


end
