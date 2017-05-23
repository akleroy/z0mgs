pro compile_atlas $
   , units = do_units $
   , mask = do_mask $
   , bksub = do_bksub $
   , convol = do_convol $
   , show = show

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

; WRITE 

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY FILES AND CHANGE UNITS TO STAGE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_units) then begin
     
     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Unit conversion '

        pgc_name = pgc_list[ii]

        for band = 1, 4 do begin
           
           infile = in_dir+pgc_name+'/'+'unwise-'+pgc_name+'-w'+str(band)+'-img-m.fits'

           outfile = out_dir+pgc_name+'_w'+str(band)+'_mjysr.fits'

           convert_unwise_to_mjysr $
              , infile = infile $
              , outfile = outfile $
              , band = band
           
        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE A BASIC MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=[0,5,5]

  if keyword_set(do_mask) then begin

     all_data = gal_data(/all)

     for ii = 0, n_pgc-1 do begin

        ;counter, ii, n_pgc, 'Unit conversion '
        pgc_name = pgc_list[ii]
        print, ii*1.0/n_pgc, ':', pgc_name

        infile = out_dir+pgc_name+'_w3_mjysr.fits'       
        outfile = out_dir+pgc_name+'_mask.fits'

        test = file_search(infile, count=ct)
        if ct eq 0 then begin
           infile = out_dir+pgc_name+'_w1_mjysr.fits' 
           test = file_search(infile, count=ct)
           if ct eq 0 then begin
              infile = out_dir+pgc_name+'_w2_mjysr.fits' 
              test = file_search(infile, count=ct)
              if ct eq 0 then begin
                 infile = out_dir+pgc_name+'_w4_mjysr.fits' 
                 test = file_search(infile, count=ct)
                 if ct eq 0 then begin
                    continue
                 endif
              endif
           endif
        endif

;       Fracture this into a program        

        map = readfits(infile, hdr, /silent)
        if n_elements(map) eq 1 then $
           continue

        make_axes, hdr, ri=ri, di=di
        
        this_dat = all_data[where(all_data.pgc eq pgc_num[ii])]

        pa = this_dat.posang_deg
        incl = this_dat.incl_deg        

        if finite(incl) eq 0 then begin
           incl = 0.0           
        endif
        if incl gt 60 then begin
           incl = 60.
        endif

        if finite(pa) eq 0 then begin
           pa = 0.0
           incl = 0.0
        endif
        
        xctr = this_dat.ra_deg
        yctr = this_dat.dec_deg

        gal_vec = [pa, incl, xctr, yctr]
        deproject, ri, di, gal_vec, rgrid=rgrid

        fid_rad = this_dat.r25_deg        
        if fid_rad lt 10./3600. or finite(fid_rad) eq 0 then begin
           fid_rad = 10./3600.
        endif

        gal_ind = where(rgrid lt fid_rad*3., gal_ct)
        mask = finite(map)
        mask[gal_ind] = 10B
        bk_ind = where(rgrid gt fid_rad*3. and rgrid lt fid_rad*6., bk_ct)
        mask[bk_ind] = 100B

        if keyword_set(show) then begin
           disp, map, max=0.5, min=0, /sq, xstyle=5, ystyle=5
           contour, mask, lev=[2,11,101], /overplot, color=cgcolor('red')
        endif

        sxaddpar, hdr, 'BUNIT', 'MASK'
        writefits, outfile, mask, hdr

        ;ch = get_kbrd()

     endfor     

  endif

  !p.multi=0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RUN A BACKGROUND SUBTRACTION OUTSIDE THE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=[0,4,5]

  if keyword_set(do_bksub) then begin

     all_data = gal_data(/all)

     for ii = 0, n_pgc-1 do begin

        counter, ii, n_pgc, 'Background subtraction '

        pgc_name = pgc_list[ii]

        maskfile = out_dir+pgc_name+'_mask.fits'
        test = file_search(maskfile, count=ct)
        if ct eq 0 then $
           continue
        mask = readfits(maskfile, /silent, mask_hdr)
        
        for band = 1, 4 do begin
           
           infile = out_dir+pgc_name+'_w'+str(band)+'_mjysr.fits'
           test = file_search(infile, count=ct)
           if ct eq 0 then $
              continue
           map = readfits(infile, hdr, /silent)

           bkind = where(mask eq 100)
           bklev = median(map[bkind])
           rms = mad(map[bkind])
           map -= bklev
           sxaddpar, hdr, 'NOISE', rms
           sxaddpar, hdr, 'BKGRD', bklev

           outfile =  out_dir+pgc_name+'_w'+str(band)+'_bksub.fits'

           writefits, outfile, map, hdr

           if keyword_set(show) then begin
              disp, map, /sq, xstyle=5, ystyle=5, min=-3.*rms, max=3.*rms
           endif
           
        endfor
        
     endfor
          
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RUN THE CONVOLUTIONS TO MATCH BEAMS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=0

end
