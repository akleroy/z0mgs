pro make_sings_cutouts $
   , cutout=cutout $
   , copy=copy $
   , convolve=convolve $
   , align=align

  out_dir = '../cutouts/sings/'
  mips_dir = '/data/tycho/0/leroy.42/ellohess/data/mips/sings/'
  los_dir = '../../ellohess/code/index/'
  kdir = '/data/tycho/0/leroy.42/ellohess/kernels/Low_Resolution/'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE UNWISE CUTOUTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(cutout) then begin

     gals = gal_data(tag='SINGS')
     n_gals = n_elements(gals)
     for ii = 0, n_gals-1 do begin
        this_gal = gals[ii]
        
        extract_unwise_stamp $
           , band=1 $
           , ra_ctr=this_gal.ra_deg $
           , dec_ctr=this_gal.dec_deg $
           , size_deg=30*60./3600. $
           , image=out_image $
           , hdr=out_hdr $
           , /show
        
        writefits $
           , out_dir+strcompress(strupcase(this_gal.name),/rem) $
           +'_UNWISE1.FITS' $
           , out_image, out_hdr

        extract_unwise_stamp $
           , band=2 $
           , ra_ctr=this_gal.ra_deg $
           , dec_ctr=this_gal.dec_Deg $
           , size_deg=30*60./3600. $
           , image=out_image $
           , hdr=out_hdr $
           , /show
        
        writefits $
           , out_dir+strcompress(strupcase(this_gal.name),/rem) $
           +'_UNWISE2.FITS' $
           , out_image, out_hdr


        extract_unwise_stamp $
           , band=3 $
           , ra_ctr=this_gal.ra_deg $
           , dec_ctr=this_gal.dec_Deg $
           , size_deg=30*60./3600. $
           , image=out_image $
           , hdr=out_hdr $
           , /show
        
        writefits $
           , out_dir+strcompress(strupcase(this_gal.name),/rem) $
           +'_UNWISE3.FITS' $
           , out_image, out_hdr

        extract_unwise_stamp $
           , band=4 $
           , ra_ctr=this_gal.ra_deg $
           , dec_ctr=this_gal.dec_Deg $
           , size_deg=30*60./3600. $
           , image=out_image $
           , hdr=out_hdr $
           , /show
        
        writefits $
           , out_dir+strcompress(strupcase(this_gal.name),/rem) $
           +'_UNWISE4.FITS' $
           , out_image, out_hdr
        
     endfor
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY SINGS FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(copy) then begin
     
     readcol, los_dir + 'singsir_release.txt' $
              , comment='#', format='A,X,A,A' $
              , sname, sband, sfile
     sname = strupcase(strcompress(sname, /rem))
     sband = (strcompress(sband, /rem))
     gb_sname = sname
     for ii = 0, n_elements(sname)-1 do $
        gb_sname[ii] = $
        lookup_galbase_name(sname[ii], alias=alias_vec, name=name_vec)
     
     for ii = 0, n_elements(sname)-1 do begin

        if sband[ii] eq 'irac1' then begin
           outfile = out_dir + strcompress(gb_sname[ii], /rem) + $
                     '_IRAC1.FITS'
           spawn, 'cp '+los_dir+sfile[ii]+' '+outfile
        endif

        if sband[ii] eq 'irac4' then begin
           outfile = out_dir + strcompress(gb_sname[ii], /rem) + $
                     '_IRAC4.FITS'
           spawn, 'cp '+los_dir+sfile[ii]+' '+outfile
        endif

        if sband[ii] eq 'mips24' then begin
           outfile = out_dir + strcompress(gb_sname[ii], /rem) + $
                     '_MIPS24.FITS'
           spawn, 'cp '+los_dir+sfile[ii]+' '+outfile
        endif

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLVE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(convolve) then begin

     unwise1_kern = kdir+ $
                    'Kernel_LoRes_WISE_FRAME_3.4_to_Gauss_20.fits'

     unwise2_kern = kdir+ $
                    'Kernel_LoRes_WISE_FRAME_4.6_to_Gauss_20.fits'

     unwise3_kern = kdir+ $
                    'Kernel_LoRes_WISE_FRAME_11.6_to_Gauss_20.fits'
     
     unwise4_kern = kdir+ $
                    'Kernel_LoRes_WISE_FRAME_22.1_to_Gauss_20.fits'

     irac1_kern = kdir+ $
                    'Kernel_LoRes_IRAC_3.6_to_Gauss_20.fits'

     irac4_kern = kdir+ $
                  'Kernel_LoRes_IRAC_8.0_to_Gauss_20.fits'

     mips24_kern = kdir+ $
                   'Kernel_LoRes_MIPS_24_to_Gauss_20.fits'

     gals = gal_data(tag='SINGS')
     n_gals = n_elements(gals)
     for ii = 0, n_gals-1 do begin
        this_gal = gals[ii]
        print, "Convolving for "+this_gal.name

        out_dir = '../cutouts/sings/'

        for jj = 0, 3 do begin

           if jj eq 0 then begin
              in_unwise = out_dir + $
                          strupcase(strcompress(this_gal.name,/rem)) + $
                          '_UNWISE1.FITS'
              out_unwise = out_dir + $
                           strupcase(strcompress(this_gal.name,/rem)) + $
                           '_UNWISE1_GAUSS20.FITS'
              kern_unwise = unwise1_kern
           endif

           if jj eq 1 then begin
              in_unwise = out_dir + $
                          strupcase(strcompress(this_gal.name,/rem)) + $
                          '_UNWISE2.FITS'
              out_unwise = out_dir + $
                           strupcase(strcompress(this_gal.name,/rem)) + $
                           '_UNWISE2_GAUSS20.FITS'
              kern_unwise = unwise2_kern
           endif

           if jj eq 2 then begin
              in_unwise = out_dir + $
                          strupcase(strcompress(this_gal.name,/rem)) + $
                          '_UNWISE3.FITS'
              out_unwise = out_dir + $
                           strupcase(strcompress(this_gal.name,/rem)) + $
                           '_UNWISE3_GAUSS20.FITS'
              kern_unwise = unwise3_kern
           endif

           if jj eq 3 then begin
              in_unwise = out_dir + $
                          strupcase(strcompress(this_gal.name,/rem)) + $
                          '_UNWISE4.FITS'
              out_unwise = out_dir + $
                           strupcase(strcompress(this_gal.name,/rem)) + $
                           '_UNWISE4_GAUSS20.FITS'
              kern_unwise = unwise4_kern
           endif

           convolve_image $
              , kernel=kern_unwise $
              , image=in_unwise $
              , outfile=out_unwise

        endfor
           
        in_mips24 = out_dir + $
                     strupcase(strcompress(this_gal.name,/rem)) + $
                     '_MIPS24.FITS'

        out_mips24 = out_dir + $
                      strupcase(strcompress(this_gal.name, /rem)) + $
                      '_MIPS24_GAUSS20.FITS'

        convolve_image $
           , kernel=mips24_kern $
           , image=in_mips24 $
           , outfile=out_mips24
        
        in_irac1 = out_dir + $
                     strupcase(strcompress(this_gal.name,/rem)) + $
                     '_IRAC1.FITS'

        out_irac1 = out_dir + $
                      strupcase(strcompress(this_gal.name, /rem)) + $
                      '_IRAC1_GAUSS20.FITS'

        convolve_image $
           , kernel=irac1_kern $
           , image=in_irac1 $
           , outfile=out_irac1
        
        in_irac4 = out_dir + $
                     strupcase(strcompress(this_gal.name,/rem)) + $
                     '_IRAC4.FITS'

        out_irac4 = out_dir + $
                      strupcase(strcompress(this_gal.name, /rem)) + $
                      '_IRAC4_GAUSS20.FITS'

        convolve_image $
           , kernel=irac4_kern $
           , image=in_irac4 $
           , outfile=out_irac4
        

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ALIGN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(align) then begin

     gals = gal_data(tag='SINGS')
     n_gals = n_elements(gals)
     for ii = 0, n_gals-1 do begin
        this_gal = gals[ii]
     endfor

  endif  

end
