pro nearby_3color $
   , rebin=do_rebin
  
  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  plot_dir = '../plots/'
  atlas_dir = '../delivery/'

  build_galaxy_list $
     , in_dir = in_dir $
     , tag=tag $
     , just=just $
     , pgc_list = pgc_list $
     , pgc_num = pgc_num $
     , dat = gal_data $
     , start = start_num $
     , stop = stop_num $
     , exclude = ['PGC17223']
  n_pgc = n_elements(pgc_list)
  sort_ind = sort((gal_data.gb_deg))
  gal_data = gal_data[sort_ind]

  cutoff = 10.*15./3600.
  ind = where(gal_data.r25_deg gt cutoff, ct)
  ngal = n_elements(ind)
  
  good_name = get_good_name(gal_data[ind].pgc)

  row = 0
  nrows = 10
  deltarow = 1.0/(nrows*1.0)

  column = 0
  ncolumns = 8
  deltacolumn = 1.0/(ncolumns*1.0)  
  
  page = 0

  for ii=0, ngal-1 do BEGIN

     if ii mod (nrows*ncolumns) eq 0 then begin

        if page gt 0 then begin
           ps, /xw  
           spawn, 'evince '+psfile+' &'
           spawn, 'convert -density 300x300 '+psfile+' '+pnfile
           !p.multi=0
        endif

        page += 1        
        psfile='../plots/test_3color_'+str(page)+'.eps'
        pnfile='../plots/test_3color_'+str(page)+'.png'
        ps, /def, file=psfile,/encapsulated,/portrait,/color,$
            xsize=8,ysize=10,/ps
        row = 0
        column = 0

     endif

     this_pgc = gal_data[ind[ii]].pgc
     if file_test(atlas_dir+'PGC'+str(this_pgc)+'_w1.fits') eq 0 then begin
        print, "Did not find PGC"+str(this_pgc)
        continue
     endif
     print, 'PGC'+str(this_pgc)
     
     w1 = readfits(atlas_dir+'PGC'+str(this_pgc)+'_w1.fits',w1hdr)
     w2 = readfits(atlas_dir+'PGC'+str(this_pgc)+'_w2.fits',w2hdr)
     w3 = readfits(atlas_dir+'PGC'+str(this_pgc)+'_w3.fits',w3hdr)
     w4 = readfits(atlas_dir+'PGC'+str(this_pgc)+'_w4.fits',w4hdr)
     fuv = readfits(atlas_dir+'PGC'+str(this_pgc)+'_fuv.fits',fuvhdr)
     nuv = readfits(atlas_dir+'PGC'+str(this_pgc)+'_nuv.fits',nuvhdr)

     has_nuv = 1B
     if (total(finite(nuv)) eq 0) or $
        (total(size(nuv,/dimen) ne size(w1,/dimen)) ne 0) then begin
        nuv = 0.0*w1+1d-10
        nuvhdr = w1hdr
        has_nuv = 0B
     endif

     sz = size(w1,/dimen)
     xcen = sz[0]/2
     ycen = sz[1]/2

     delta = gal_data[ind[ii]].r25_deg*3600./5.5*1.25

     hextract,w1,w1hdr,xcen-delta,xcen+delta,ycen-delta,ycen+delta
     hextract,w2,w2hdr,xcen-delta,xcen+delta,ycen-delta,ycen+delta
     hextract,w3,w3hdr,xcen-delta,xcen+delta,ycen-delta,ycen+delta
     hextract,w4,w4hdr,xcen-delta,xcen+delta,ycen-delta,ycen+delta
     hextract,nuv,nuvhdr,xcen-delta,xcen+delta,ycen-delta,ycen+delta

     sz = size(w1, /dimen)
     if keyword_set(do_rebin) then begin
        hrebin, w1, w1hdr, outsize=[sz[0]*4, sz[0]*4]
        hrebin, w2, w2hdr, outsize=[sz[0]*4, sz[0]*4]
        hrebin, w3, w3hdr, outsize=[sz[0]*4, sz[0]*4]
        hrebin, w4, w4hdr, outsize=[sz[0]*4, sz[0]*4]
        hrebin, nuv, nuvhdr, outsize=[sz[0]*4, sz[0]*4]
     endif

     make_axes, w1hdr, ra=ra, da=da, ri=ri, di=di

     loadct,0
     
     pos = [column*deltacolumn, 1.0-(row+1.0)*deltarow, $
            (column+1.0)*deltacolumn, 1.0-((row)*deltarow)]
     column += 1
     if column ge ncolumns then begin
        column = 0
        row += 1
     endif

     minval = [-2, -2, -3.0]
     maxval = [1.5, 1.5, -0.25]

     disp3,alog10(w3),alog10(w1),alog10(nuv),ra,da $
           , /sq,min=minval,max=maxval $ 
           , charthick=3,thick=3,xthick=3,ythick=3 $
           , xstyle=4, ystyle=4, position = pos $
           , /noerase
     pgcstring = 'PGC'+strcompress(str(gal_data[ind[ii]].pgc), /rem)
     name = strcompress(good_name[ii], /rem)
     al_legend, /top, /right, box=0, clear=0, lines=-99 $
                , [pgcstring,name] $
                , textcolor=cgcolor('white') $
                , charsize=0.5, charthick=1.0

     nuv_string = ''
     if has_nuv eq 0 then nuv_string='No GALEX'
     lat_string = '!8b!6 = '+strcompress(string(gal_data[ind[ii]].gb_deg,format='(3I)'),/rem)
     al_legend, /bottom, /left, box=0, clear=0, lines=-99 $
                , [nuv_string, lat_string] $
                , textcolor=cgcolor('white') $
                , charsize=0.5, charthick=1.0

  endfor

  ps, /xw  
  spawn, 'evince '+psfile+' &'
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile
  !p.multi=0

  stop
end
