pro z0mgs_photometry

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DIRECTORY AND BUILD GALAXY LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  in_dir = '../unwise/dlang_custom/z0mgs/PGC/'
  out_dir = '../measurements/'
;  sample_dir = '../measurements/samples/'
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

  bands = ['FUV','NUV','WISE1','WISE2','WISE3','WISE4']
  n_bands = n_elements(bands)
  
  spacing_deg = 7.5/3600.
  spacing_rad = spacing_deg*!dtor
  area_sr = 3.0d*sqrt(3.0d)/2.0d*(spacing_rad/2.d)^2
  
  limit_r25 = 2.0d
  min_limit = 30./3600.

  oversamp = 5.23

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE OUTPUT STRUCTURE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  nan = !values.f_nan
  empty = $
     {pgc:0L $
      , band: '' $
      , rms: nan $
      , std: nan $
      , ct: 0L $
      , unc_stat: nan $
      , unc_conf: nan $
      , limit_deg: nan $
      , posang_deg: nan $
      , incl_deg: nan $
      , med: nan $
      , med_plusrad: nan $
      , med_minusrad: nan $
      , mean: nan $
      , mean_plusrad: nan $
      , mean_minusrad: nan $
      , hybrid: nan $
      , hybrid_plusrad: nan $
      , hybrid_minusrad: nan $
     }
  phot = replicate(empty, n_pgc, n_bands)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = 0, n_pgc-1 do begin

     counter, ii, n_pgc, 'Photometry for galaxy '

     ;sample_file = $
     ;   sample_dir + $
     ;   strcompress('PGC'+str(gal_data[ii].pgc), /rem) + $
     ;   '_samples.fits'

     ;if file_test(sample_file) eq 0 then begin
     ;   print, "Did not find "+sample_file
     ;   ;stop
     ;   continue        
     ;endif

     ;tab = mrdfits(sample_file,1,h,/silent)
          
     this_pgc = gal_data[ii].pgc

     w1_file = atlas_dir+'PGC'+str(this_pgc)+'_w1.fits'
     has_w1 = file_test(w1_file)
     w2_file = atlas_dir+'PGC'+str(this_pgc)+'_w2.fits'
     has_w2 = file_test(w2_file)
     w3_file = atlas_dir+'PGC'+str(this_pgc)+'_w3.fits'
     has_w3 = file_test(w3_file)
     w4_file = atlas_dir+'PGC'+str(this_pgc)+'_w4.fits'
     has_w4 = file_test(w4_file)
     fuv_file = atlas_dir+'PGC'+str(this_pgc)+'_fuv.fits'
     has_fuv = file_test(fuv_file)
     nuv_file = atlas_dir+'PGC'+str(this_pgc)+'_nuv.fits'
     has_nuv = file_test(nuv_file)

     if has_w1 eq 0 then begin
        print, this_pgc, ' lacks WISE1. Should not.'
        continue
     endif
     w1 = readfits(w1_file,w1hdr, /silent)
     w2 = readfits(w2_file,w2hdr, /silent)
     w3 = readfits(w3_file,w3hdr, /silent)
     w4 = readfits(w4_file,w4hdr, /silent)

     area_sr = (sxpar(w1hdr, 'CD1_1')*!dtor)^2

     if has_nuv then begin
        nuv = readfits(nuv_file, nuvhdr, /silent)
     endif else begin
        nuv = w1*!values.f_nan
        nuvhdr = w1hdr
     endelse

     if has_fuv then begin
        fuv = readfits(fuv_file, fuvhdr, /silent)
     endif else begin
        fuv = w1*!values.f_nan
        fuvhdr = w1hdr
     endelse

     this_incl = gal_data[ii].incl_deg
     if finite(this_incl) eq 0 then $
        this_incl = 0.0
     this_incl = (this_incl - 10.0) > 0.0

     this_posang = gal_data[ii].posang_deg
     if finite(this_posang) eq 0 then begin
        this_incl = 0.0
        this_posang = 0.0
     endif

     galpos = [this_posang, this_incl, gal_data[ii].ra_deg, gal_data[ii].dec_deg]

     make_axes, w1hdr, ri=ri, di=di
     deproject, ri, di, galpos, rgrid=rad_deg
;     deproject, tab.ra_deg, tab.dec_deg, galpos, rgrid=rad_deg, /vec

     limit_deg = (limit_r25*gal_data[ii].r25_deg)
     if finite(limit_deg) eq 0 or limit_deg lt min_limit then $
        limit_deg = min_limit

     plus_limit = limit_deg*1.33
     minus_limit = limit_deg/1.33

     phot_ind = where(rad_deg le limit_deg, phot_ct)
     if phot_ct eq 0 then begin
        print, "No valid data for PGC "+str(gal_data[ii].pgc)
        continue
     endif

     plus_ind = where(rad_deg le plus_limit)
     minus_ind = where(rad_deg le minus_limit)

     for jj = 0, n_bands-1 do begin

        this_phot = empty
        this_phot.pgc = gal_data[ii].pgc
        this_phot.band = bands[jj]
        this_phot.limit_deg = limit_deg
        this_phot.posang_deg = this_posang
        this_phot.incl_deg = this_incl

        if jj eq 0 then begin
           y = fuv
        endif
        if jj eq 1 then begin
           y = nuv
        endif
        if jj eq 2 then begin
           y = w1
        endif
        if jj eq 3 then begin
           y = w2
        endif
        if jj eq 4 then begin
           y = w3
        endif
        if jj eq 5 then begin
           y = w4
        endif
        this_phot.ct = phot_ct

;       Adjust by "mega" factor
        y *= 1d6
        ;this_phot.rms *= 1d6
        ;this_phot.std *= 1d6

        ;this_phot.unc_stat = sqrt(this_phot.ct) * area_sr * this_phot.rms
        ;this_phot.unc_conf = sqrt(this_phot.ct) * area_sr * this_phot.std
        
        if total(finite(y[phot_ind])) eq 0 then $
           continue

        bins = $
           bin_data(rad_deg[phot_ind], y[phot_ind], /nan $
                    , xmin=-0.5*spacing_deg $
                    , xmax=limit_deg $
                    , binsize=spacing_deg)

        plus_bins = $
           bin_data(rad_deg[plus_ind], y[plus_ind], /nan $
                    , xmin=-0.5*spacing_deg $
                    , xmax=plus_limit $
                    , binsize=spacing_deg)

        minus_bins = $
           bin_data(rad_deg[minus_ind], y[minus_ind], /nan $
                    , xmin=-0.5*spacing_deg $
                    , xmax=minus_limit $
                    , binsize=spacing_deg)
        
        this_phot.mean = total(bins.ymean*area_sr*bins.counts, /nan)
        this_phot.med = total(bins.ymed*area_sr*bins.counts, /nan)
        if n_elements(bins) le 3 then begin
           this_phot.hybrid = this_phot.mean
        endif else begin
           this_phot.hybrid = $
              total(bins[0:2].ymean*area_sr*bins[0:2].counts, /nan) + $
              total(bins[3:*].ymed*area_sr*bins[3:*].counts, /nan)
        endelse

        this_phot.mean_plusrad = $
           total(plus_bins.ymean*area_sr*plus_bins.counts, /nan)
        this_phot.med_plusrad = $
           total(plus_bins.ymed*area_sr*plus_bins.counts, /nan)
        if n_elements(bins) le 3 then begin
           this_phot.hybrid_plusrad = this_phot.mean_plusrad
        endif else begin
           this_phot.hybrid_plusrad = $
              total(plus_bins[0:2].ymean*area_sr*plus_bins[0:2].counts, /nan) + $
              total(plus_bins[3:*].ymed*area_sr*plus_bins[3:*].counts, /nan)
        endelse

        this_phot.mean_minusrad = $
           total(minus_bins.ymean*area_sr*minus_bins.counts, /nan)
        this_phot.med_minusrad = $
           total(minus_bins.ymed*area_sr*minus_bins.counts, /nan)
        if n_elements(bins) le 3 then begin
           this_phot.hybrid_minusrad = this_phot.mean_minusrad
        endif else begin
           this_phot.hybrid_minusrad = $
              total(minus_bins[0:2].ymean*area_sr*minus_bins[0:2].counts, /nan) + $
              total(minus_bins[3:*].ymed*area_sr*minus_bins[3:*].counts, /nan)
        endelse

        phot[ii,jj] = this_phot

     endfor

  endfor
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  save, file='../measurements/z0mgs_photometry.idl', phot

  mwrfits, phot $
           , '../measurements/z0mgs_photometry.fits' $
           , /create

end
