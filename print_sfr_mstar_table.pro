pro print_sfr_mstar_table

  nan = !values.f_nan
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ INDEX
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  index15 = mrdfits('../measurements/delivery_index_gauss15.fits',1,h15)
  n_index = n_elements(index15)

  ngc_names = pgc_to_othername(index15.pgc, prefix='NGC')
  ugc_names = pgc_to_othername(index15.pgc, prefix='UGC')
  ic_names = pgc_to_othername(index15.pgc, prefix='IC')

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PRINT THE TEXT FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  outfile = '../measurements/z0mgs_fullsfrmstar.txt'
  openw, 1, outfile
  printf, 1, '# Column 1: PGC number'
  printf, 1, '# Column 2: NGC name [if available]'
  printf, 1, '# Column 3: UGC name [if available]'
  printf, 1, '# Column 4: IC name [if available]'
  printf, 1, '# Column 5: Adopted distance [Mpc]'
  printf, 1, '# Column 6: Uncertainty in log10 distance [dex]'
  printf, 1, '# Column 7: log10 stellar mass [Msun]'
  printf, 1, '# Column 8: Uncertainty (not including distance uncertainty) in stellar mass [dex]'
  printf, 1, '# Column 9: Adopted WISE1 mass-to-light ratio [Msun/Lsun]'
  printf, 1, '# Column 10: Method used to estimate WISE1 mass-to-light ratio'
  printf, 1, '# Column 11: log10 star formation rate [Msun/yr]'
  printf, 1, '# Column 12: Uncertainty (not including distance uncertainty) in star formation rate [dex]'
  printf, 1, '# Column 13: Method used to estimate star formation rate'
  printf, 1, '# Column 14: Offset from star forming main sequence [dex]'
  printf, 1, '# Column 16: Flags'

  for ii = 0, n_index-1 do begin

     line = ''
     line += string(index15[ii].pgc,format='(I8)')+', '
     
;       It seems like a good general rule is that we want the shortest
;       name that starts with NGC, IC, or UGC. Ask for some eyeballs
;       to check this.

     this_ngc = ngc_names[ii]
     if this_ngc ne '' then begin
        tokens = strsplit(this_ngc,';',/extract)
        len_tokens = strlen(tokens)
        dummy = min(len_tokens, minind)
        print_ngc = tokens[minind]
     endif else begin
        print_ngc = ''
     endelse
     line += string(print_ngc,format='(A10)')+', '

     this_ugc = ugc_names[ii]
     if this_ugc ne '' then begin
        tokens = strsplit(this_ugc,';',/extract)
        len_tokens = strlen(tokens)
        dummy = min(len_tokens, minind)
        print_ugc = tokens[minind]
     endif else begin
        print_ugc = ''
     endelse
     line += string(print_ugc,format='(A10)')+', '

     this_ic = ic_names[ii]
     if this_ic ne '' then begin
        tokens = strsplit(this_ic,';',/extract)
        len_tokens = strlen(tokens)
        dummy = min(len_tokens, minind)
        print_ic = tokens[minind]
     endif else begin
        print_ic = ''
     endelse
     line += string(print_ic,format='(A10)')+', '

     line += string(index15[ii].dist_mpc,format='(F5.1)')+', '
     line += string(index15[ii].e_dist_dex,format='(F5.2)')+', '

     if finite(index15[ii].logmass_color) then begin
        line += string(index15[ii].logmass_color,format='(F5.2)')+', '
        line += string(index15[ii].e_logmass_color,format='(F5.2)')+', '
        line += string(index15[ii].mtol_color,format='(F5.2)')+', '
     endif else begin
        line += string(index15[ii].logmass_color*nan,format='(F5.2)')+', '
        line += string(index15[ii].e_logmass_color*nan,format='(F5.2)')+', '
        line += string(index15[ii].mtol_color*nan,format='(F5.2)')+', '
     endelse
     line += string(index15[ii].method_mtol,format='(A10)')+', '     

     if finite(index15[ii].logsfr_fixed) then begin
        line += string(index15[ii].logsfr_fixed,format='(F5.2)')+', '
        line += string(index15[ii].e_logsfr_fixed,format='(F5.2)')+', '
     endif else begin
        line += string(index15[ii].logsfr_fixed*nan,format='(F5.2)')+', '
        line += string(index15[ii].e_logsfr_fixed*nan,format='(F5.2)')+', '
     endelse

     line += string(index15[ii].method_sfr,format='(A10)')+', '
     line += string(index15[ii].deltams,format='(F5.2)')+', '

     flags = ''
     if index15[ii].galaxy_overlap_flag then $
        flags+='G'
     if index15[ii].photometry_mismatch_flag then $
        flags+='P'
     if index15[ii].sat_effects_wise1 or $
        index15[ii].sat_effects_wise2 or $
        index15[ii].sat_effects_wise3 or $
        index15[ii].sat_effects_wise4 or $
        index15[ii].sat_effects_fuv or $
        index15[ii].sat_effects_nuv $
     then $
        flags+='A'
     if index15[ii].star_flag_wise1 or $
        index15[ii].star_flag_wise2 or $
        index15[ii].star_flag_wise3 or $
        index15[ii].star_flag_wise4 or $
        index15[ii].star_flag_fuv or $
        index15[ii].star_flag_nuv $
     then $
        flags+='S'
     line += string(flags,format='(A5)')

     printf, 1, line

  endfor

  close, 1
  
  stop
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PRINT THE LATEX STUB
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  outfile = '../measurements/z0mgs_sfrmstar.tex'
  openw, 1, outfile

  nlines = 0
  print_total = 10

  for ii = 0, n_index-1 do begin
     
     if nlines ge print_total then continue
     nlines += 1
     
     line = ''
     line += string(index15[ii].pgc,format='(I8)')+' & '
     
;       It seems like a good general rule is that we want the shortest
;       name that starts with NGC, IC, or UGC. Ask for some eyeballs
;       to check this.

     this_ngc = ngc_names[ii]
     if this_ngc ne '' then begin
        tokens = strsplit(this_ngc,';',/extract)
        len_tokens = strlen(tokens)
        dummy = min(len_tokens, minind)
        print_ngc = tokens[minind]
     endif else begin
        print_ngc = ''
     endelse
     line += string(print_ngc,format='(A10)')+' & '

     this_ugc = ugc_names[ii]
     if this_ugc ne '' then begin
        tokens = strsplit(this_ugc,';',/extract)
        len_tokens = strlen(tokens)
        dummy = min(len_tokens, minind)
        print_ugc = tokens[minind]
     endif else begin
        print_ugc = ''
     endelse
     line += string(print_ugc,format='(A10)')+' & '

     this_ic = ic_names[ii]
     if this_ic ne '' then begin
        tokens = strsplit(this_ic,';',/extract)
        len_tokens = strlen(tokens)
        dummy = min(len_tokens, minind)
        print_ic = tokens[minind]
     endif else begin
        print_ic = ''
     endelse
     line += string(print_ic,format='(A10)')+' & '

     line += '$'+string(index15[ii].dist_mpc,format='(F5.1)')+'$ & '
     line += '$'+string(index15[ii].e_dist_dex,format='(F5.2)')+'$ & '

     if finite(index15[ii].logmass_color) then begin
        line += '$'+string(index15[ii].logmass_color,format='(F5.2)')+'$ $\pm$ '
        line += '$'+string(index15[ii].e_logmass_color,format='(F5.2)')+'$ & '
        line += '$'+string(index15[ii].mtol_color,format='(F5.2)')+'$ & '
     endif else begin
        line += '$'+string(index15[ii].logmass_color*nan,format='(F5.2)')+'$ $\pm$ '
        line += '$'+string(index15[ii].e_logmass_color*nan,format='(F5.2)')+'$ & '
        line += '$'+string(index15[ii].mtol_color*nan,format='(F5.2)')+'$ & '
     endelse
     line += string(index15[ii].method_mtol,format='(A10)')+' & '     

     if finite(index15[ii].logsfr_fixed) then begin
        line += '$'+string(index15[ii].logsfr_fixed,format='(F5.2)')+'$ $\pm$ '
        line += '$'+string(index15[ii].e_logsfr_fixed,format='(F5.2)')+'$ & '
     endif else begin
        line += '$'+string(index15[ii].logsfr_fixed*nan,format='(F5.2)')+'$ $\pm$ '
        line += '$'+string(index15[ii].e_logsfr_fixed*nan,format='(F5.2)')+'$ & '
     endelse

     line += string(index15[ii].method_sfr,format='(A10)')+' & '
     line += '$'+string(index15[ii].deltams,format='(F5.2)')+'$ & '

     flags = ''
     if index15[ii].galaxy_overlap_flag then $
        flags+='G'
     if index15[ii].photometry_mismatch_flag then $
        flags+='P'
     if index15[ii].sat_effects_wise1 or $
        index15[ii].sat_effects_wise2 or $
        index15[ii].sat_effects_wise3 or $
        index15[ii].sat_effects_wise4 or $
        index15[ii].sat_effects_fuv or $
        index15[ii].sat_effects_nuv $
     then $
        flags+='A'
     if index15[ii].star_flag_wise1 or $
        index15[ii].star_flag_wise2 or $
        index15[ii].star_flag_wise3 or $
        index15[ii].star_flag_wise4 or $
        index15[ii].star_flag_fuv or $
        index15[ii].star_flag_nuv $
     then $
        flags+='S'
     line += string(flags,format='(A5)')

     line += ' \\'
     printf, 1, line

  endfor

  close, 1


end
