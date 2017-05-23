pro make_sfng_cutouts $
   , cutout=cutout $
   , convolve=convolve

  out_dir = '../cutouts/sfng/'
  los_dir = '../../ellohess/code/index/'
  kdir = '/data/tycho/0/leroy.42/ellohess/kernels/Low_Resolution/'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE UNWISE CUTOUTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  gal_list = $
     ['ngc0628' $
      , 'ngc1672' $
      , 'ngc3351' $
      , 'ngc3627' $
      , 'ngc4254' $
      , 'ngc4303' $
      , 'ngc4321' $
      , 'ngc4535' $
      , 'ngc5068' $
      , 'ngc6744' $
     ]

  if keyword_set(cutout) then begin

     gals = gal_data(gal_list)
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
           , dec_ctr=this_gal.dec_deg $
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
           , dec_ctr=this_gal.dec_deg $
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
; CONVOLVE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(convolve) then begin

     gals = gal_data(gal_list)

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

     gals = gal_data(gal_list)
     n_gals = n_elements(gals)
     for ii = 0, n_gals-1 do begin
        this_gal = gals[ii]
        print, "Convolving for "+this_gal.name

        out_dir = '../cutouts/sfng/'

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
        
     endfor

  endif

end
