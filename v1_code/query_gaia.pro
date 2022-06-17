pro query_gaia $
   , image=image $
   , ruwe=do_ruwe $
   , outfile=outfile $
   , send=do_send

  hdr = headfits(image)

; assume north/south

  make_axes, hdr, ra=ra, da=da

  ramin = ra[n_elements(ra)-1]
  ramax = ra[0]

  if ramin gt ramax then begin
     twopart_call = 1B 
     this_ramin = ramin
     this_ramax = 360.0
  endif else begin
     twopart_call = 0B
     this_ramin = ramin
     this_ramax = ramax
  endelse

  decmin = da[0]
  decmax = da[n_elements(da)-1]

; to strings
  ra_low_string = strcompress(string(this_ramin,format='(F9.5)'), /rem)
  ra_high_string = strcompress(string(this_ramax,format='(F9.5)'), /rem)

  dec_low_string = strcompress(string(decmin,format='(F9.5)'), /rem)
  dec_high_string = strcompress(string(decmax,format='(F9.5)'), /rem)

  if n_elements(outfile) eq 0 then begin
     outfile = '../stars/gaia/dummy.txt'
  endif

  if keyword_set(do_ruwe) eq 0 then begin
     call = 'curl "http://gea.esac.esa.int/tap-server/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=CSV&QUERY=' + $
            'SELECT+source_id,ra,dec,phot_g_mean_mag,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error+FROM+gaiadr2.gaia_source+'+ $
            'WHERE+ra+between+'+ra_low_string+'+and+'+ra_high_string+ $
            '+and+dec+between+'+dec_low_string+'+and+'+dec_high_string+ $
            '" > '+outfile
  endif else begin
     call = 'curl "http://dc.zah.uni-heidelberg.de/__system__/tap/run/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=CSV&QUERY=' + $
            'SELECT+source_id,ra,dec,phot_g_mean_mag,ruwe+FROM+gaiadr2.gaia_source+'+ $
            'WHERE+ra+between+'+ra_low_string+'+and+'+ra_high_string+ $
            '+and+dec+between+'+dec_low_string+'+and+'+dec_high_string+ $
            '" > '+outfile
  endelse

  if keyword_set(do_send) then begin
     spawn, call
  endif else begin
     print, call
  endelse

  if twopart_call then begin

     this_ramin = 0.0
     this_ramax = ramax

     ra_low_string = strcompress(string(this_ramin,format='(F9.5)'), /rem)
     ra_high_string = strcompress(string(this_ramax,format='(F9.5)'), /rem)
     
     if keyword_set(do_ruwe) eq 0 then begin
        call = 'curl "http://gea.esac.esa.int/tap-server/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=CSV&QUERY=' + $
               'SELECT+source_id,ra,dec,phot_g_mean_mag,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error+FROM+gaiadr2.gaia_source+'+ $
               'WHERE+ra+between+'+ra_low_string+'+and+'+ra_high_string+ $
               '+and+dec+between+'+dec_low_string+'+and+'+dec_high_string+ $
               '" >> '+outfile
     endif else begin
        call = 'curl "http://dc.zah.uni-heidelberg.de/__system__/tap/run/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=CSV&QUERY=' + $
               'SELECT+source_id,ra,dec,phot_g_mean_mag,ruwe+FROM+gaiadr2.gaia_source+'+ $
               'WHERE+ra+between+'+ra_low_string+'+and+'+ra_high_string+ $
               '+and+dec+between+'+dec_low_string+'+and+'+dec_high_string+ $
               '" >> '+outfile
     endelse
     
     if keyword_set(do_send) then begin
        spawn, call
     endif else begin
        print, call
     endelse

  endif

end
