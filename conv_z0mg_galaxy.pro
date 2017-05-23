pro conv_z0mg_galaxy,$
	infile=infile,$
	start_psf=start_psf,$
	end_psf=end_psf,$
	outfile=outfile

;+
;  start_psf possibilities = 'w1', 'w2', 'w3', 'w4', 'nuv', 'fuv'
;  
;  end_psf possibilities = 'g7p5', 'g11', 'g15', 'g20', 'w2', 'nuv'
;
;  	* note that not all combinations are possible
;
;  Assumes infile and outfile are full paths
;
;	for now the uncertainty convolution isn't implemented but it can be
;-

	; get possible psf combinations and kernel paths
	readcol,'kernel_index.txt',psfin,psfout,kernelpath,format='A,A,A',/silent

	thiskernel = where(psfin eq start_psf and psfout eq end_psf,ct)
	if ct ne 1 then BEGIN
		print,'Not a valid combination of input and output psfs.'
		stop
	endif

	; convolution call
	z0mgconv,$
		file=infile,$
		kernel_file=kernelpath[thiskernel],$
		out_file=outfile

end	
