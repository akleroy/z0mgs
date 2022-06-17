pro sort_tiles_galex $
   , mkdir=mkdir $
   , raw=raw $
   , copy=copy $
   , unzip=unzip

;+
; Sort the 150,000 GALEX tiles into smaller groups to be separated
; into separated directories.
;-

  hr_ra = lindgen(24)

  raw_dir = '../galex/all_tiles/'
  in_dir = '../galex/scratch_tiles/'
  out_dir = '../galex/sorted_tiles/'

  if keyword_set(mkdir) then begin
     for ii = 0, n_elements(hr_ra)-1 do begin
        command = 'mkdir '+out_dir+str(hr_ra[ii])+'h'
        spawn, command
     endfor
  endif

  if keyword_set(raw) then begin
;     flist = file_search(raw_dir+'*.fits.gz', count=count)
     flist = file_search(raw_dir+'*rrhr*.fits.gz', count=count)
     for jj = 0L, count-1 do begin
        counter, jj, count, 'Copying file '
        froot = strmid(flist[jj],strlen(in_dir))
        cpcmd = 'cp '+flist[jj]+' '+in_dir+'.'
        spawn, cpcmd
     endfor
  endif

  if keyword_set(unzip) then begin
     flist = file_search(in_dir+'*.fits.gz', count=count)
     for jj = 0L, count-1 do begin
        counter, jj, count, 'Unzipping file '
        command = 'gzip -d '+flist[jj]
        spawn, command
     endfor
  endif

  if keyword_set(copy) then begin
     flist = file_search(in_dir+'*-int.fits*', count=count)     
     for ii = 0L, count-1 do begin
        counter, ii, count, 'Copying file '
        this_hdr = headfits(flist[ii])
        dummy = sxpar(this_hdr, 'CRVAL1', count=crct)
        if crct eq 0 then begin
           message, "No CRVAL1 found - skipping tile "+flist[ii], /info
           continue
        endif

        nx = sxpar(this_hdr,'NAXIS1')
        ny = sxpar(this_hdr,'NAXIS2')
        
        xyad, this_hdr, nx/2, ny/2, mean_ra, mean_dec
                
        hour = long(floor(mean_ra/15.))
        hour_str = str(hour)+'h'

        froot = strmid(flist[ii],strlen(in_dir))

        dummy = file_search(out_dir+hour_str+'/'+froot, count=there)
        if there eq 0 then begin
           command = 'cp '+flist[ii]+' '+out_dir+hour_str+'/'+froot
           spawn, command
        endif

        flagroot = mg_streplace(froot,'-int.fits','-flags.fits')     
        originfile = mg_streplace(flist[ii],'-int.fits','-flags.fits')     

        dummy = file_search(out_dir+hour_str+'/'+flagroot, count=there)
        if there eq 0 then begin
           command = 'cp '+originfile+' '+out_dir+hour_str+'/'+flagroot
           spawn, command
        endif

        rrhrroot = mg_streplace(froot,'-int.fits','-rrhr.fits')     
        originfile = mg_streplace(flist[ii],'-int.fits','-rrhr.fits')     
        dummy = file_search(out_dir+hour_str+'/'+rrhrroot, count=there)
        if there eq 0 then begin
           command = 'cp '+originfile+' '+out_dir+hour_str+'/'+rrhrroot
           spawn, command
        endif

     endfor
  endif

end
