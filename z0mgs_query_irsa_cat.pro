pro z0mgs_query_irsa_cat $
   , targetname_OR_coords $
   , catalog=catalog $
   , radius=radius $
   , radunits=radunits $
   , outfile=outfile $
   , query=url_q $
   , noread=skip_read $
   , table=table

;+
;
; *****************************************************************************
; This is a hacked version of the provided query_irsa_cat, which seems
; broken because IDLs web requests are either being rejected or the
; protocols are out of date. So instead, we just use this program to
; form a url and then spawn a wget command. The results go to a file
; that can be read with read_ipac_table.
; *****************************************************************************
;
;-

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFAULTS AND DEFINITIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

on_error,2
compile_opt idl2

if N_params() lt 1 then begin
  print,'Syntax - info = query_irsa_cat(targetname_or_coords,'
  print,'           [catalog=catalog,radius=radius,radunits=radunits,'
  print,'            outfile=outfile,/nowrite,/noread])'
endif

IF NOT(keyword_set(radius)) THEN $
   radius = 60
IF NOT(keyword_set(radunits)) THEN $
   radunits = 'arcsec'

if n_elements(outfile) eq 0 then begin
   outfile = "tmp_irsa.out"
endif

dummy = file_search(outfile, count=exists)
IF exists gt 0 THEN BEGIN
   spawn, 'rm -rf '+outfile
ENDIF 

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONSTRUCT THE PARTS OF THE QUERY STRING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

root = 'http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query'

; Catalog

IF keyword_set(catalog) THEN $
   catalog_name=catalog ELSE $
      catalog_name='fp_psc'

catstr='&catalog='+catalog_name

; Object

target = targetname_OR_coords

IF N_elements(target) EQ 2 THEN BEGIN 
   ra = double(target[0])
   dec = double(target[1])
   objstr = '&objstr='+strn(ra)+'+'+strn(dec)
ENDIF $
ELSE BEGIN 
   object = repstr(target,'+','%2B')
   object = repstr(strcompress(object),' ','+')
   objstr = '&objstr='+object
ENDELSE 

; No empty string
IF strlen(objstr) le 8 THEN BEGIN
  print, 'Empty object string not allowed.'
  return
ENDIF

; SEARCH SHAPE AND SIZE

spatial_str='Cone' 
spatial_param_name=['radius','radunits']
spatial_param_value_str = [strn(radius), radunits]

nspat = n_elements(spatial_param_name)

spatstr = '&spatial='+spatial_str
spatparstr = ''

FOR i = 0l, nspat-1 DO $
   spatparstr=spatparstr+'&'+spatial_param_name[i]+'='+spatial_param_value_str[i]
      
;;;; USE IPAC FORMAT

out_fmt = '?outfmt=1'

;;;; combine into query string

url_q = root+out_fmt+objstr+spatstr+spatparstr+catstr

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; EXECUTE THE QUERY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

command = "wget -O "+outfile+" '"+url_q+"'"
spawn, command

dummy = file_search(outfile, count=ct)
if ct ne 1 then begin
   message, "No match for WGET attempt. Proceeding.", /info        
   return
endif

;fail_string = "."
;if strpos(dummy, fail_string) ne -1 then begin
;   message, "Object not found in NED. Proceeding.", /info        
;   return, -1
;endif
;found = 1B

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if keyword_set(skip_read) eq 0 then begin

   table = z0mgs_read_ipac_table(outfile)

endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DONE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

return

END
