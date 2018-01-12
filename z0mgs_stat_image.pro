pro z0mgs_stat_image $
   , infile=infile $
   , outfile=outfile $
   , mask=maskfile $
   , reject=rejectfile $
   , weight=wtfile $
   , thresh=rej_thresh $
   , print=do_print $
   , galex=galex
  
  if n_elements(rej_thresh) eq 0 then $
     rej_thresh = 0.2

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  image = readfits(infile, hdr, /silent)
  mask = readfits(maskfile, mask_hdr, /silent)
  reject = readfits(rejectfile, reject_hdr, /silent)
  rej_mask = reject gt rej_thresh

  if keyword_set(galex) then begin
     weight = readfits(wtfile, wt_hdr, /silent)
     flat = image*sqrt(weight)
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; STAT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sz = size(image)
  quarters = mask*0B
  quarters[0:sz[1]/2,0:sz[2]/2] = 1
  quarters[0:sz[1]/2,sz[2]/2+1:*] = 2
  quarters[sz[1]/2+1:*,0:sz[2]/2] = 3
  quarters[sz[1]/2+1:*,sz[2]/2+1:*] = 4

  fin_image = finite(image)

  for ii = 0, 4 do begin
     if ii eq 0 then begin
        ind = where(mask ne 10 and mask ne 100 and $
                    fin_image, ct)
        ext = 'ALL'
     endif
     if ii eq 1 then begin
        ind = where(mask ne 10 and mask ne 100 and $
                    quarters eq 1 and fin_image, ct)
        ext = 'Q1'
     endif
     if ii eq 2 then begin
        ind = where(mask ne 10 and mask ne 100 and $
                    quarters eq 2 and fin_image, ct)
        ext = 'Q2'
     endif
     if ii eq 3 then begin
        ind = where(mask ne 10 and mask ne 100 and $
                    quarters eq 3 and fin_image, ct)
        ext = 'Q3'
     endif
     if ii eq 4 then begin
        ind = where(mask ne 10 and mask ne 100 and $
                    quarters eq 4 and fin_image, ct)
        ext = 'Q4'
     endif
     
     if ct lt 6 then begin
        sxaddpar, hdr, 'MED'+ext, -999.
        sxaddpar, hdr, 'MEAN'+ext, -999.
        sxaddpar, hdr, 'MAD'+ext, -999.
        sxaddpar, hdr, 'STD'+ext, -999.
        sxaddpar, hdr, 'MAX'+ext, -999.
        sxaddpar, hdr, 'REJ'+ext, -999.
        if keyword_set(galex) then begin
           sxaddpar, hdr, 'FLATMAD'+ext, -999.
           sxaddpar, hdr, 'FLATSTD'+ext, -999.           
        endif
        continue
     endif

     vec = image[ind]
     rej_frac = total(rej_mask[ind]*1.0 / (n_elements(ind)*1.0))
     if keyword_set(galex) then begin
        flat_vec = flat[ind]
     endif

     sxaddpar, hdr, 'STD'+ext, stddev(vec,/nan)
     sxaddpar, hdr, 'MAX'+ext, max(vec,/nan)
     sxaddpar, hdr, 'MED'+ext, median(vec)
     sxaddpar, hdr, 'MEAN'+ext, mean(vec,/nan)
     if keyword_set(galex) then begin
        sxaddpar, hdr, 'FLATSTD'+ext, stddev(flat_vec, /nan)
     endif

     ok_ind = where(rej_mask[ind] eq 0, ok_ct)
     if ok_ct gt 6 then begin
        sxaddpar, hdr, 'MAD'+ext, mad(vec[ok_ind])
        if keyword_set(galex) then begin
           sxaddpar, hdr, 'FLATMAD'+ext, mad(flat_vec[ok_ind])
        endif
     endif else begin
        sxaddpar, hdr, 'MED'+ext, -999.
        sxaddpar, hdr, 'MEAN'+ext, -999.
        sxaddpar, hdr, 'MAD'+ext, -999.
        sxaddpar, hdr, 'STD'+ext, -999.
        sxaddpar, hdr, 'MAX'+ext, -999.
        sxaddpar, hdr, 'REJ'+ext, -999.
        if keyword_set(galex) then begin
           sxaddpar, hdr, 'FLATMAD'+ext, -999.
           sxaddpar, hdr, 'FLATSTD'+ext, -999.           
        endif
        continue
     endelse
     sxaddpar, hdr, 'REJ'+ext, rej_frac

     if keyword_set(do_print) then begin
        print, ext, sxpar(hdr,'MED'+ext), sxpar(hdr,'MEAN'+ext) $
               , sxpar(hdr,'MAD'+ext), sxpar(hdr,'STD'+ext) $
               , sxpar(hdr,'REJ'+ext)
     endif

  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INTO HEADER AND WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  writefits, outfile, image, hdr

end
