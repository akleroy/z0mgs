pro build_unwise_mask $
   , pgc=pgc_name $
   , galdata=this_dat $
   , outfile=outfile $
   , show=show

  atlas_dir = '../unwise/atlas/'

  if n_elements(outfile) eq 0 then $
     outfile = atlas_dir+pgc_name+'_mask.fits'
  
  if n_elements(this_dat) eq 0 then begin
     this_dat = gal_data(pgc_name)
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE NATIVE RESOLUTION DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  map = readfits(atlas_dir+pgc_name+'_w1_mjysr.fits' , hdr, /silent)
  if map[0] eq -1 then found = 0 else found = 1

  if found eq 0 then begin

;    Based on my checks, no galaxy misses more than one band.
     map = readfits(atlas_dir+pgc_name+'_w2_mjysr.fits' , hdr, /silent)
     if map[0] eq -1 then found = 0 else found = 1

     message, 'No WISE image found in the atlas. Returning.', /info
     return
  endif


; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IDENTIFY THE REGION NEAR THE GALAXY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

; Work out coordinates
  make_axes, hdr, ri=ri, di=di

; Set inclination and position angle, default to face on lacking this
; information. Even for edge on systems, assume an inclination of ~60
; to avoid issues with the thickness of the disk.

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
  
; Work out the galactocentric radius from the orientation and center

  xctr = this_dat.ra_deg
  yctr = this_dat.dec_deg

  gal_vec = [pa, incl, xctr, yctr]
  deproject, ri, di, gal_vec, rgrid=rgrid

; Pick the radius to define the area around the galaxy. This will be
; avoided in a background subtraction.

  fid_rad = this_dat.r25_deg        
  if fid_rad lt 10./3600. or finite(fid_rad) eq 0 then begin
     fid_rad = 10./3600.
  endif

; Label < 3 r25 (galaxy) and 3-6 r25 (background)

  gal_ind = where(rgrid lt fid_rad*1., gal_ct)
  mask = finite(map)
  mask[gal_ind] = 10B
  bk_ind = where(rgrid ge fid_rad*1. and rgrid lt fid_rad*2., bk_ct)
  mask[bk_ind] = 100B
  bk_ind = where(rgrid ge fid_rad*2. and rgrid lt fid_rad*3., bk_ct)
  mask[bk_ind] = 110B
  bk_ind = where(rgrid ge fid_rad*3. and rgrid lt fid_rad*4., bk_ct)
  mask[bk_ind] = 120B
  bk_ind = where(rgrid ge fid_rad*4. and rgrid lt fid_rad*5., bk_ct)
  mask[bk_ind] = 130B
  bk_ind = where(rgrid ge fid_rad*5. and rgrid lt fid_rad*6., bk_ct)
  mask[bk_ind] = 140B

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW IF REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  if keyword_set(show) then begin

     disp, map, max=0.5, min=0, /sq, xstyle=5, ystyle=5

     contour, mask, lev=[2,11,101,111,121,131], /overplot, color=cgcolor('red')

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  sxaddpar, hdr, 'BUNIT', 'MASK'
  writefits, outfile, mask, hdr

end
