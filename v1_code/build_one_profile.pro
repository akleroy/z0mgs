function build_one_profile $
   , pgc_name = pgc_name $
   , gal_data = gal_data $
   , index = index $
   , data_dir = data_dir $
   , outfile=outfile $
   , bin_size=bin_size $
   , beam_size=beam_size $
   , vary_orient=vary_orient $
   , show=show $
   , pause=pause $
   , empty=empty

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINITIONS AND DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(bin_size) eq 0 then $
     bin_size = 7.5

  if n_elements(beam_size) eq 0 then $
     beam_size = 15.

  if n_elements(unc_posang) eq 0 then $
     unc_posang = 10.

  if n_elements(unc_incl) eq 0 then $
     unc_incl = 7.5

  max_bands = 6

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; EMPTY STRUCTURE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  nan = !values.f_nan
  ring = $
     { $
     pgc_name:'' $
     , ra_ctr: nan $
     , dec_ctr: nan $
     , no_posang: 0B $
     , posang_deg: nan $
     , delta_posang: 0.0 $
     , no_incl: 0B $
     , capped_incl: 0B $
     , incl_deg: nan $
     , delta_incl: 0.0 $
     , beam_arcsec: nan $
     , pix_arcsec: nan $
     , bin_arcsec: nan $
     , band: '' $
     , rms_image: nan $
     , std_image: nan $
     , ring_id: -1L $
     , ring_rad: nan $
     , ring_sum: nan $
     , ring_mean: nan $
     , ring_med: nan $     
     , ring_counts: 0L $     
     , ring_mad: nan $     
     , ring_std: nan $     
     , ring_perc16: nan $     
     , ring_perc84: nan $
     , ring_emean: nan $
     , ring_emed: nan $
     , label: '' $
     }

  if keyword_set(empty) then $
     return, ring
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; USE THE INDEX FILE TO READ THE RELEVANT DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  found_ref = 0B
  
  if index.has_fuv eq 0 and $
     index.has_nuv eq 0 and $
     index.has_wise1 eq 0 and $
     index.has_wise2 eq 0 and $
     index.has_wise3 eq 0 and $
     index.has_wise4 eq 0 then begin
     print, "No files expected. Returning."
     return, ring
  endif

  n_bands = 0

  if index.has_fuv then begin
     fuv_file = data_dir + pgc_name + '_fuv.fits'     
     test = file_search(fuv_file, count=fuv_ct)
     if fuv_ct eq 0 then begin
        print, "No FUV file found, but one expected. Stopping."
        stop
     endif
     fuv = readfits(fuv_file, fuv_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = fuv_hdr
     endif
     n_bands += 1
  endif

  if index.has_nuv then begin
     nuv_file = data_dir + pgc_name + '_nuv.fits'
     test = file_search(nuv_file, count=nuv_ct)
     if nuv_ct eq 0 then begin
        print, "No NUV file found, but one expected. Stopping."
        stop
     endif
     nuv = readfits(nuv_file, nuv_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = nuv_hdr
     endif
     n_bands += 1
  endif

  if index.has_wise1 then begin
     wise1_file = data_dir + pgc_name + '_w1.fits'
     test = file_search(wise1_file, count=wise1_ct)
     if wise1_ct eq 0 then begin
        print, "No WISE1 file found, but one expected. Stopping."
        stop
     endif
     wise1 = readfits(wise1_file, wise1_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = wise1_hdr
     endif
     n_bands += 1
  endif

  if index.has_wise2 then begin
     wise2_file = data_dir + pgc_name + '_w2.fits'
     test = file_search(wise2_file, count=wise2_ct)
     if wise2_ct eq 0 then begin
        print, "No WISE2 file found, but one expected. Stopping."
        stop
     endif
     wise2 = readfits(wise2_file, wise2_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = wise2_hdr
     endif
     n_bands += 1
  endif

  if index.has_wise3 then begin
     wise3_file = data_dir + pgc_name + '_w3.fits'
     test = file_search(wise3_file, count=wise3_ct)
     if wise3_ct eq 0 then begin
        print, "No WISE3 file found, but one expected. Stopping."
        stop
     endif
     wise3 = readfits(wise3_file, wise3_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = wise3_hdr
     endif
     n_bands += 1
  endif

  if index.has_wise4 then begin
     wise4_file = data_dir + pgc_name + '_w4.fits'
     test = file_search(wise4_file, count=wise4_ct)
     if wise4_ct eq 0 then begin
        print, "No WISE4 file found, but one expected. Stopping."
        stop
     endif
     wise4 = readfits(wise4_file, wise4_hdr)
     if found_ref eq 0 then begin
        found_ref = 1B
        ref_hdr = wise4_hdr
     endif
     n_bands += 1
  endif

  if found_ref eq 0 then begin
     print, "I should have a reference header at this point. Stopping."
     stop
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD ASTROMETRY AND RESAMPLE (HACK FOR PIXEL CLIPPING)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  nx = sxpar(ref_hdr, 'NAXIS1')
  ny = sxpar(ref_hdr, 'NAXIS2')

  binfac = 4

  hrebin, fltarr(nx, ny), ref_hdr, out = [binfac*nx, binfac*ny]

  if index.has_fuv then $
     hrebin, fuv, fuv_hdr, out = [binfac*nx, binfac*ny]

  if index.has_nuv then $
     hrebin, nuv, nuv_hdr, out = [binfac*nx, binfac*ny]

  if index.has_wise1 then $
     hrebin, wise1, wise1_hdr, out = [binfac*nx, binfac*ny]

  if index.has_wise2 then $
     hrebin, wise2, wise2_hdr, out = [binfac*nx, binfac*ny]

  if index.has_wise3 then $
     hrebin, wise3, wise3_hdr, out = [binfac*nx, binfac*ny]

  if index.has_wise4 then $
     hrebin, wise4, wise4_hdr, out = [binfac*nx, binfac*ny]

  make_axes, ref_hdr, ri=ri, di=di

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FILL IN HEADER
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  ring.pgc_name = pgc_name

  ring.ra_ctr = gal_data.ra_deg
  ring.dec_ctr = gal_data.dec_deg

  ring.posang_deg = gal_data.posang_deg
  ring.incl_deg = gal_data.incl_deg

  ring.posang_deg = gal_data.posang_deg
  ring.incl_deg = gal_data.incl_deg

  ring.beam_arcsec = beam_size
  ring.pix_arcsec = get_pixel_scale(ref_hdr)*3600.
  ring.bin_arcsec = bin_size

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SPECIFY BINS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  rmin = -0.5*bin_size
  delta = (max(di, /nan) - min(di, /nan))*3600./2.0
  rmax = delta
  nbins = long(ceil((rmax - rmin) / bin_size))

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVERT ASTROMETRY TO GALACTOCENTRIC RADIUS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  fid_posang = gal_data.posang_deg
  fid_incl = gal_data.incl_deg

  if finite(fid_posang) eq 0 or finite(fid_incl) eq 0. then begin
     fid_posang = 0.0
     fid_incl = 0.0
     ring.no_posang = 1B
     ring.no_incl = 1B
  endif

  if fid_incl gt 75. then begin
     fid_incl = 75.
     ring.capped_incl = 1B
  endif

  ring.incl_deg = fid_incl
  ring.posang_deg = fid_posang

  if keyword_set(vary_orient) then begin
     rings = replicate(ring, max_bands, nbins*5)
  endif else begin
     rings = replicate(ring, max_bands, nbins) 
  endelse

  for jj = 0, 4 do begin
        
     if keyword_set(vary_orient) eq 0 and $
        jj ne 0 then continue

     if jj eq 0 then begin
        this_ring = ring
        this_ring.label = 'BEST'
        this_ring.delta_posang = 0.0
        this_ring.delta_incl = 0.0
     endif

     if jj eq 1 then begin
        this_ring = ring
        this_ring.label = 'POSANGDOWN'
        this_ring.delta_posang = -1.*unc_posang
        this_ring.delta_incl = 0.0
     endif

     if jj eq 2 then begin
        this_ring = ring
        this_ring.label = 'POSANGUP'
        this_ring.delta_posang = +1.*unc_posang
        this_ring.delta_incl = 0.0
     endif

     if jj eq 3 then begin
        this_ring = ring
        this_ring.label = 'INCLDOWN'
        this_ring.delta_posang = 0.0
        this_ring.delta_incl = -1.*unc_incl
     endif

     if jj eq 4 then begin
        this_ring = ring
        this_ring.label = 'INCLUP'
        this_ring.delta_posang = 0.0
        this_ring.delta_incl = +1.*unc_incl
     endif

     adopted_incl = ((fid_incl + this_ring.delta_incl) < 85.) > 0.0
     adopted_posang = (fid_posang + this_ring.delta_posang)
        
     pos_vec = $
        [adopted_posang, adopted_incl $
         , this_ring.ra_ctr, this_ring.dec_ctr]
     deproject, ri, di, pos_vec, rgrid=rgrid, tgrid=tgrid
     rgrid *= 3600.
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER BANDS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
     for ii = 0, 5 do begin
        
        if ii eq 0 then begin
           if index.has_fuv eq 0 then $
              continue
           data = fuv
           band = 'FUV'
           rms = index.rms_fuv
           std = index.std_fuv
        endif
        
        if ii eq 1 then begin
           if index.has_nuv eq 0 then $
              continue
           data = nuv
           band = 'NUV'
           rms = index.rms_nuv
           std = index.std_nuv
        endif
        
        if ii eq 2 then begin
           if index.has_wise1 eq 0 then $
              continue
           data = wise1
           band = 'WISE1'
           rms = index.rms_wise1
           std = index.std_wise1
        endif
        
        if ii eq 3 then begin
           if index.has_wise2 eq 0 then $
              continue
           data = wise2
           band = 'WISE2'
           rms = index.rms_wise2
           std = index.std_wise2
        endif
        
        if ii eq 4 then begin
           if index.has_wise3 eq 0 then $
              continue
           data = wise3
           band = 'WISE3'
           rms = index.rms_wise3
           std = index.std_wise3
        endif
        
        if ii eq 5 then begin
           if index.has_wise4 eq 0 then $
              continue
           data = wise4
           band = 'WISE4'
           rms = index.rms_wise4
           std = index.std_wise4
        endif
        
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BIN THE IMAGE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
           
        bins = bin_data(rgrid, data, /nan $
                        , xmin=rmin, xmax=rmax, binsize=bin_size $
                        , perc=0.16)
        if n_elements(bins) ne nbins then begin
           print, "Number of bins inconsistent with estimate. Stopping."
           stop
        endif
        
        ring_id = lindgen(nbins)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SAVE TO STRUCTURE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        
        temp_rings = replicate(this_ring, nbins)
        temp_rings[*].band = band
        temp_rings[*].rms_image = rms
        temp_rings[*].std_image = std
        temp_rings[*].ring_id = ring_id
        temp_rings[*].ring_rad = bins.xmid
        temp_rings[*].ring_sum = bins.ysum
        temp_rings[*].ring_mean = bins.ymean
        temp_rings[*].ring_med = bins.ymed
        temp_rings[*].ring_counts = bins.counts
        temp_rings[*].ring_mad = bins.ymad
        temp_rings[*].ring_std = bins.ystd
        temp_rings[*].ring_perc16 = bins.ylo_perc
        temp_rings[*].ring_perc84 = bins.yhi_perc
        
        start_ind = jj*nbins
        stop_ind =  (jj+1)*nbins-1

        rings[ii,start_ind:stop_ind] = temp_rings
        
     endfor
     
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WORK OUT ERRORS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; POST-PROCESS THE MONTE CARLO SIDE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SAVE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(outfile) gt 0 then begin    
     mwrfits, rings, outfile, /create
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RETURN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  return, rings  

end
