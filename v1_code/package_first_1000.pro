pro package_first_1000

  tab = mrdfits('../measurements/delivery_index_gauss7p5.fits', 1, h)

  spawn, 'rm -rf ../first_1000/'
  spawn, 'mkdir ../first_1000/'

  for ii = 0, 1000 do begin

     counter, ii, 1000, "Copying galaxy "

     this_pgc_name = strcompress(tab[ii].pgc_name, /rem)
     
     spawn, 'cp ../delivery/'+this_pgc_name+'_*.fits ../first_1000/.'

  endfor

end
