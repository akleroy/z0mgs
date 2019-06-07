pro shrink_figs $
   , in_dir = in_dir $
   , out_dir = out_dir $
   , resolution = resolution

  if n_elements(in_dir) eq 0 then $
     in_dir = '../plots/for_paper/'

  if n_elements(out_dir) eq 0 then $
     out_dir = '../plots/for_astroph/'
  
  if n_elements(resolution) eq 0 then $
     resolution = '150'

  fig_list_1 = file_search(in_dir+"*.eps")  
  fig_list_2 = file_search(in_dir+"*.png")  
  fig_list = [fig_list_1, fig_list_2]
  n_fig = n_elements(fig_list)

  for i = 0, n_fig-1 do begin
     last_slash = strpos(fig_list[i],'/',/reverse_search)
     fig_list[i] = strmid(fig_list[i],last_slash+1,strlen(fig_list[i])-(last_slash+1))
  endfor
     
  spawn, "rm -rf "+out_dir+"*.eps"
  spawn, "rm -rf "+out_dir+"*.jpg"

  for i = 0, n_fig-1 do begin
     print, "Shrinking "+fig_list[i]
     source = in_dir+fig_list[i]
     dest = out_dir+fig_list[i]
     spawn, 'gs -r'+resolution+' -dEPSCrop -dTextAlphaBits=4 -sDEVICE=jpeg -sOutputFile='+ $
            source+'.jpg -dBATCH -dNOPAUSE '+source
     spawn, 'jpeg2ps '+source+'.jpg > '+dest
     spawn, 'rm -rf '+source+'.jpg'
  endfor

;  spawn, "rm "+out_dir+"*.jpg"

end
