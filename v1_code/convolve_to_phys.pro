pro convolve_to_phys $
   , just_res=just_res $
   , only_gal=only_gal $
   , skip_gal=skip_gal $
   , dry_run=dry_run

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINITIONS, DIRECTORIES, DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  tol = 0.1

  target_res_pc = [5e2, 1e3, 2e3]
  res_ext = ['500pc', '1kpc', '2kpc']
  out_dir = '../convolved/'+['500pc','1kpc','2kpc']+'/'
  n_res = n_elements(target_res_pc)

  bands = ['fuv','nuv','w1','w2','w3','w4']
  band_ext = ['_fuv_gauss7p5_interp.fits' $
              ,'_nuv_gauss7p5_interp.fits' $
              ,'_w1_gauss7p5_interp.fits' $
              ,'_w2_gauss7p5_interp.fits' $
              ,'_w3_gauss7p5_interp.fits' $
              ,'_w4_gauss15_interp.fits' $
             ]
  band_dir = [ $
             '/data/kant/0/leroy.42/allsky/postprocess/interpol/' $   ; fuv
             , '/data/kant/0/leroy.42/allsky/postprocess/interpol/' $ ; nuv
             , '/data/kant/0/leroy.42/allsky/postprocess/interpol/' $ ; w1
             , '/data/kant/0/leroy.42/allsky/postprocess/interpol/' $ ; w2
             , '/data/kant/0/leroy.42/allsky/postprocess/interpol/' $ ; w3
             , '/data/kant/0/leroy.42/allsky/postprocess/interpol/' $  ; w4
             ]
  band_res_as = [7.5, 7.5, 7.5, 7.5, 7.5, 15.]
  n_bands = n_elements(bands)


;             '/data/tycho/0/behrens.52/masked_interp/' $   ; fuv
;             , '/data/tycho/0/behrens.52/masked_interp/' $ ; nuv
;             , '/data/tycho/0/behrens.52/masked_interp/' $ ; w1
;             , '/data/tycho/0/behrens.52/masked_interp/' $ ; w2
;             , '/data/tycho/0/behrens.52/masked_interp/' $ ; w3
;             , '/data/kant/0/leroy.42/allsky/delivery/' $  ; w4

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD A LIST OF TARGETS AND DISTANCES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  index = mrdfits('../measurements/delivery_index_gauss15.fits',1,hi)
  
  check_missing = 0B
  if check_missing then begin
     dat = gal_data(tag='Z0MGS', /full)
     missing_pgc = bytarr(n_elements(index))
     for ii = 0L, n_elements(index)-1 do begin
        if total(index[ii].pgc eq dat.pgc) gt 0 then $
           continue
        missing_pgc[ii] = 1B
     endfor
     stop
  endif

  dat = gal_data(pgc=index.pgc)
  n_gal = n_elements(dat)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER RESOLUTION AND TARGET
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = 0L, n_res-1 do begin

     this_res = target_res_pc[ii]
     this_res_ext = res_ext[ii]
     this_out_dir = out_dir[ii]

     if n_elements(just_res) gt 0 then begin
        if (total(this_res_ext eq just_res) eq 0) and $
           (total(this_res eq just_res) eq 0) then begin
           continue
        endif
     endif

     for jj = 0L, n_gal-1 do begin
        
        this_gal = strcompress(dat[jj].pgcname,/rem)
        if this_gal eq '' then $
           continue

        if n_elements(skip_gal) gt 0 then begin
           if (total(this_gal eq skip_gal) gt 0) then begin
              continue
           endif
        endif

        if n_elements(only_gal) gt 0 then begin
           if (total(this_gal eq only_gal) eq 0) then begin
              continue
           endif
        endif

        this_dist = dat[jj].dist_mpc
        
        for kk = 0L, n_bands-1 do begin
           
           this_band = bands[kk]

           ;if this_band eq 'w4' then continue

           this_band_ext = band_ext[kk]
           this_band_dir = band_dir[kk]
           this_band_res_as = band_res_as[kk]

           this_res_pc = $
              this_band_res_as*!dtor/3600. * $
              this_dist * 1d6

           target_res_as = $
              this_res / (this_dist * 1d6) / !dtor * 3600.
           
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%        
           
           infile = $
              this_band_dir + $
              this_gal + $
              this_band_ext

           if file_test(infile) eq 0 then begin
              message $
                 , "Input file not found: "+infile $
                 , /info
              continue
           endif

           outfile = $
              this_out_dir + $
              this_gal + $
              '_'+this_band+'_'+ $
              this_res_ext + $
              '.fits'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DO THE CONVOLUTION OR COPYING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%        

; SKIP - RESOLUTION NOT POSSIBLE
           if this_res_pc gt this_res*(1.+tol) then begin
              message, "Resolution "+this_res_ext+ $
                       " not met for "+this_gal+" . Skipping.", /info
              continue
           endif

; COPY - RESOLUTION CLOSE ENOUGH
           if (this_res_pc le this_res*(1.+tol)) and $
              (this_res_pc ge this_res*(1.-tol)) $
           then begin

              print, "Current: ", this_res_pc, " target: ", this_res
              message, "Resolution "+this_res_ext+ $
                       " within tolerance for "+this_gal+" . Copying.", /info

              if keyword_set(dry_run) eq 0 then begin
                 cmd = '\cp '+infile+' '+outfile
                 spawn, cmd
              endif

           endif

; CONVOLVE - RESOLUTION NEEDS ADJUSTING

           if this_res_pc lt this_res*(1.-tol) then begin

              message, "Resolution "+this_res_ext+ $
                       " allows convolution for "+this_gal+ $
                       " . Convolving.", /info

              if keyword_set(dry_run) eq 0 then begin

                 conv_with_gauss $
                    , data=infile $
                    , target_beam=target_res_as*[1,1,0] $
                    , out_file=outfile

              endif

           endif

        endfor

     endfor

  endfor

end
