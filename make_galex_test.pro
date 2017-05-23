pro make_galex_test

  glist = ['ngc2976', 'ngc0024', 'ugc05159']
  outdir = '../cutouts/galex_test/'

  for ii = 0, n_elements(glist)-1 do begin
     gal = glist[ii]

     s = gal_data(gal)

     for jj = 0, 1 do begin
        extract_galex_stamp $
           , fuv=(jj eq 0) $
           , ra_ctr=s.ra_deg $
           , dec_ctr=s.dec_deg $
           , size_deg=30./60. $
           , image=outim $
           , weight=weightim $
           , hdr=outhdr $
           , /show ;$
;           , /bksub

        if jj eq 0 then band = 'fuv' else band = 'nuv'
        writefits, outdir+gal+'_'+band+'_mjysr.fits', outim, outhdr
        sxaddpar, outhdr, 'BUNIT', 'WEIGHT'
        writefits, outdir+gal+'_'+band+'_weight.fits', weightim, outhdr
     endfor     

  endfor


  end
