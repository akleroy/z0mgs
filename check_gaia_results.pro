pro check_gaia_results

  tab = mrdfits('../measurements/delivery_index_gauss15.fits', 1, h)
  n = n_elements(tab)
  pgcnamelist = strcompress(tab.pgc_name, /rem)
  gdir = '../stars/gaia/'
  for ii = 0, n-1 do begin
     counter, ii, n, 'Checking GAIA file'
     gfile = gdir+pgcnamelist[ii]+'_gaia.txt'
     if file_test(gfile) then begin
        nlines = file_lines(gfile)
     endif else begin
        nlines = 0
     endelse
     if nlines le 1 then begin
        if n_elements(problems) eq 0 then $
           problems = [gfile] $
        else $
           problems = [problems, gfile]
     endif
  endfor
  
  stop

; this file straddles the meridian
  image = '../unwise/atlas/PGC4_w1_gauss7p5_small.fits'
  hdr = headfits(image)
  query_gaia $
     , image=image, /send

; m31 is too big. break into parts.
  image = '../unwise/atlas/PGC2557_w1_gauss7p5_small.fits'
  hdr = headfits(image)
  query_gaia_m31 $
     , image=image, /send

; this file does not exist? has wise2 though. excluded.
  image = '../unwise/atlas/PGC89980_w1_gauss7p5_small.fits'
  hdr = headfits(image)
  query_gaia $
     , image=image

; this file straddles the meridian
  image = '../unwise/atlas/PGC166785_w1_gauss7p5_small.fits'
  hdr = headfits(image)
  query_gaia $
     , image=image
  
; this file does not exist? has wise2 though. excluded.
  image = '../unwise/atlas/PGC917425_w1_gauss7p5_small.fits'
  hdr = headfits(image)
  query_gaia $
     , image=image ;$
     ;, ruwe=do_ruwe $
     ;, outfile=outfile $
     ;, send=do_send

; ../stars/gaia/PGC4_gaia.txt 
; ../stars/gaia/PGC2557_gaia.txt 
; ../stars/gaia/PGC89980_gaia.txt 
; ../stars/gaia/PGC166785_gaia.txt 
; ../stars/gaia/PGC917425_gaia.txt
end
