pro sort_tiles_unwise $
   , mkdir=mkdir $
   , copy=copy

;+
; Sort the UNWISE tiles into smaller groups to be separated
; into separated directories.
;-

  hr_ra = lindgen(24)

  in_dir = 'data/'
  out_dir = 'sorted_tiles/'

  if keyword_set(mkdir) then begin
     for ii = 0, n_elements(hr_ra)-1 do begin
        command = 'mkdir '+out_dir+str(hr_ra[ii])+'h'
        spawn, command
     endfor
  endif

  if keyword_set(copy) then begin
     flist = file_search(in_dir+'*.fits*', count=count)
     for ii = 0L, count-1 do begin
        counter, ii, count, 'Copying file '
        this_hdr = headfits(flist[ii])
        make_axes, this_hdr, ri=ri, di=di
        mean_ra = mean(ri)
        mean_dec = mean(di)
        hour = long(floor(mean_ra/15.))
        hour_str = str(hour)+'h'
        froot = strmid(flist[ii],strlen(in_dir))
        command = 'cp '+flist[ii]+' '+out_dir+hour_str+'/'+froot
        spawn, command
     endfor
  endif

end
