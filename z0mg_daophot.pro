pro z0mg_daophot,$
	file=file

	im = readfits(file,hdr)
	sky,im,skymode,skysig

	blah = strsplit(file,'.',/extract)

	in_hmin = 1.
	in_fwhm = 3.
	in_roundlim = [0.2,1.0]
	in_sharplim = [-1.0,1.0]

	find,im,x,y,flux,sharp,round,in_hmin,in_fwhm,in_sharplim,in_roundlim

	disp,im,/sq,max=10,min=0
	oplot,x,y,ps=1,color=getcolor('red')

	in_phpadu = 1
	in_apr = [5]
	in_skyrad = [7,9]
	aper,im,x,y,mag,errap,sky,skyerr,in_phpadu,in_apr,in_skyrad,/nan


	okstars = where(x gt 400 and y lt 200)

	in_ronois = 0.1
	in_phpadu = 1
	in_psfrad = 5
	in_fitrad = 3
	in_psfname = blah[0]+'_psf.fits'

	getpsf,im,x,y,mag,sky,in_ronois,in_phpadu,gauss,psf,okstars,$
		in_psfrad,in_fitrad,in_psfname

	in_rcrit = 7

	group,x,y,in_rcrit,ngroup
	nx = n_elements(x)

	nstar,im,indgen(nx),x,y,mag,sky,ngroup,in_phpadu,in_ronois,in_psfname,$
		magerr,iter,chisq,peak
	newim = im

	substar,newim,x,y,mag,indgen(nx),in_psfname

	writefits,blah[0]+'_substar.fits',newim,hdr

	xyad,hdr,x,y,ra,dec
	save,file=blah[0]+'_mags.sav',mag,ra,dec,sharp,round

end
