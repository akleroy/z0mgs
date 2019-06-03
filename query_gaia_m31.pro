pro query_gaia_m31 $
   , image=image $
   , ruwe=do_ruwe $
   , outfile=outfile $
   , send=do_send
  
  hdr = headfits(image)
  make_axes, hdr, ra=ra, da=da

  step_size = 1.0               ; degrees

  ramin = ra[n_elements(ra)-1]
  ramax = ra[0]
  if ramin gt ramax then begin
     print, "Not set up for a wrap in this case."
     stop
  endif
  nstep_ra = ceil((ramax-ramin)/step_size)

  decmin = da[0]
  decmax = da[n_elements(da)-1]
  nstep_dec = ceil((decmax-decmin)/step_size)
  
  first = 1B
  for ii = 0, nstep_ra-1 do begin
     for jj = 0, nstep_dec-1 do begin

        this_ramin = (ramin + (ii)*step_size)
        this_ramax = (ramin + (ii+1.)*step_size) < ramax

        this_decmin = (decmin + (jj)*step_size)
        this_decmax = (decmin + (jj+1.)*step_size) < decmax
        
        ra_low_string = strcompress(string(this_ramin,format='(F9.5)'), /rem)
        ra_high_string = strcompress(string(this_ramax,format='(F9.5)'), /rem)
        
        dec_low_string = strcompress(string(this_decmin,format='(F9.5)'), /rem)
        dec_high_string = strcompress(string(this_decmax,format='(F9.5)'), /rem)
        
        if n_elements(outfile) eq 0 then begin
           outfile = '../stars/gaia/dummy.txt'
        endif
        
        if keyword_set(first) then begin
           symbol = '>'
           first = 0
        endif else begin
           symbol = '>>'
        endelse

        if keyword_set(do_ruwe) eq 0 then begin
           call = 'curl "http://gea.esac.esa.int/tap-server/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=CSV&QUERY=' + $
                  'SELECT+source_id,ra,dec,phot_g_mean_mag,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error+FROM+gaiadr2.gaia_source+'+ $
                  'WHERE+ra+between+'+ra_low_string+'+and+'+ra_high_string+ $
                  '+and+dec+between+'+dec_low_string+'+and+'+dec_high_string+ $
                  '" '+symbol+' '+outfile
        endif else begin
           call = 'curl "http://dc.zah.uni-heidelberg.de/__system__/tap/run/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=CSV&QUERY=' + $
                  'SELECT+source_id,ra,dec,phot_g_mean_mag,ruwe+FROM+gaiadr2.gaia_source+'+ $
                  'WHERE+ra+between+'+ra_low_string+'+and+'+ra_high_string+ $
                  '+and+dec+between+'+dec_low_string+'+and+'+dec_high_string+ $
                  '" '+symbol+' '+outfile
        endelse
        
        if keyword_set(do_send) then begin
           spawn, call
        endif else begin
           print, call
        endelse

     endfor
  endfor

end
