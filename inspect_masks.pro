pro inspect_masks $
   , radius=do_radius $
   , galaxies=do_galaxies $
   , pause=pause $
   , band=band $
   , res=res $
   , start=start_num $
   , stop=stop_num

  if n_elements(band) eq 0 then band = 'w1'
  if n_elements(res) eq 0 then res = 'gauss7p5'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  unwise_dir = '../unwise/atlas/' 
  galex_dir = '../galex/atlas/'
  mask_dir = '../masks/'
  delivery_dir = '../delivery/'

  build_galaxy_list $
     , in_dir = in_dir $
     , tag=tag $
     , just=only $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = gal_data $
     , start = start_num $
     , stop = stop_num $
     , exclude = ['PGC17223','PGC89980','PGC917425']
  
  n_pgc = n_elements(pgc_list)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  counter = 0
  !p.multi=[0,3,2]

  if n_elements(start_num) eq 0 then start_num = 0
  if n_elements(stop_num) eq 0 then stop_num = n_pgc

  for ii = start_num, stop_num-1 do begin
     
     pgc_name = strcompress(pgc_list[ii], /rem)
     this_dat = gal_data[ii]
     
     if n_elements(just) gt 0 then $
        if total(pgc_name eq just) eq 0 then $
           continue

     image_file = delivery_dir+pgc_name+'_'+band+'_'+res+'.fits'
     rad_mask_file = delivery_dir+pgc_name+'_'+res+'_rgrid.fits'
     gal_mask_file = delivery_dir+pgc_name+'_'+res+'_galaxies.fits'
     
     image = readfits(image_file, image_hdr)
     loadct, 0
     disp, image, min=0., max=0.25

     if keyword_set(do_radius) then begin
        rad = readfits(rad_mask_file, rhdr)
        contour, rad, /overplot, color=cgcolor('green'), thick=3 $
                 , lev=sxpar(rhdr, 'FIDRAD')*[1,2]
     endif

     if keyword_set(do_galaxies) then begin
        gals = readfits(gal_mask_file, rhdr)
        contour, gals, /overplot, color=cgcolor('red'), thick=3 $
                 , lev=[1]
     endif

     if keyword_set(pause) then begin
        if counter eq 5 then begin
           print, 'Showing '+pgc_name+' or galaxy '+str(ii)+ ' - hit key to continue.'
           ch = get_kbrd(1)
           counter = 0
           !p.multi=[0,3,2]           
        endif else begin
           counter += 1
        endelse
     endif

  endfor

end
