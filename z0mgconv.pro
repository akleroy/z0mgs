
;+
; NAME: 
; 	z0mgconv
;
; PURPOSE: 
; 	Convolves an image or a cube with a kernel, keeps track of 
; 	important details in headers.
;
; INPUTS:
; 	file: fits file containing image or cube to convolve
; 	kernel_file: fits file containing kernel
; 	out_file: name for output fits image or cube
;
; KEYWORDS:
; 	kfspec_format: deals properly with extensions for Herschel
; 		spectroscopy as done by KINGFISH
; 	start_fwhm/end_fwhm: these keyword let you account properly
; 		for the uncertainty convolution, by finding the number
; 		of independent resolution elements that were averaged
; 		together during the convolution
; 	
; OUTPUTS:
;
; EXAMPLE:
;
;
; To Do:
; 	- needs error checking on unc image convolution
; 	- figure out what is up with /no_pad in convolve
; 	- better NaN replacement, nearest neighbor?
; 	- possibility to rebin image instead of kernel? 
; 		maybe important if image is undersampled
;-

pro z0mgconv,$
	file=file,$
	unc_file=unc_file,$
	kernel_file=kernel_file,$
	out_file=out_file,$
	outunc_file=outunc_file,$
	silent=silent,$
	kfspec_format=kfspec_format,$
	gunc_format=gunc_format,$
	start_fwhm=start_fwhm,$
	end_fwhm=end_fwhm

	on_error,0

	; check that file, kernel_file and out_file have all been set
	if keyword_set(file) eq 0 or keyword_set(kernel_file) eq 0 or $
		keyword_set(out_file) eq 0 then $
			message,'Must set file, kernel_file and out_file!'

	; check that if unc_file is set that outunc_file is also set
	if keyword_set(unc_file) and (keyword_set(outunc_file) eq 0) then $
		message,'Need to name output uncertainty file.'
	
	; check that the input files actually point to files
	if file_test(file) eq 0 then message,'Input file not found!'
	if file_test(kernel_file) eq 0 then message,'Kernel file not found!'
	if keyword_set(unc_file) then BEGIN
		if file_test(unc_file) eq 0 then $
			message,'Uncertainty file not found!'
	endif

	; read in the files
	if keyword_set(kfspec_format) then BEGIN
		fhdr = mrdfits(file,0,fullhdr)
		im = mrdfits(file,'image',hdr)
		index = mrdfits(file,'ImageIndex',hdrind)
	endif else BEGIN
		fits_read,file,im,hdr
	endelse

	fits_read,kernel_file,kernel,kerhdr

	if keyword_set(unc_file) then BEGIN
		if keyword_set(kfspec_format) then BEGIN
			ufhdr = mrdfits(unc_file,0,ufullhdr)
			uncim = mrdfits(unc_file,'image',unchdr)
			uindex = mrdfits(unc_file,'ImageIndex',uhdrind)
		endif else BEGIN
			fits_read,unc_file,uncim,unchdr
		endelse

		if keyword_set(gunc_format) then BEGIN
			uncim = uncim[*,*,3]
		endif
	endif

	; figure out dimensions of the images
	imsize = size(im,/dimen)
	kersize = size(kernel,/dimen)

	; get astrometry info
	extast,hdr,ast_info
	
	if (n_elements(ast_info) eq 0) then message,'No astrometry in image!'

	getrot,hdr,angle,cdelt
	orig_cdelt = cdelt

	; here in karl's program there is derotation
	; for now I am going to ignore this
	
	; check if cdelts are equal, if not need more code
	if (abs(cdelt[0])-abs(cdelt[1]))/abs(cdelt[0]) gt 1e-3 then $
		message,'Need to match cdelts in image, code more!'

	; calculate image pixel scale
	image_scale = abs(cdelt)*3600.0

	if (n_elements(silent) EQ 0) then $
		print,'Image scale [arcsec] = ', image_scale

	; normalize the kernel
	kernel = kernel/total(kernel)
	
	; get the kernel pixel scale
	ker_scale = fxpar(kerhdr,'PIXSCALE',count=count)
	if count eq 0 then ker_scale = fxpar(kerhdr,'SECPIX',count=count)
	if count eq 0 then ker_scale = abs(fxpar(kerhdr,'CD1_1',count=count)*3600.0)
	if count eq 0 then ker_scale = abs(fxpar(kerhdr,'CDELT1',count=count)*3600.0)
	if count eq 0 then message,'No pixel scale found in kernel image.'

	if keyword_set(silent) eq 0 then $
		print,'Kernel scale [arcsec] = ',ker_scale

	; copy header to output hdr
	outhdr = hdr
	if keyword_set(unc_file) then outunchdr = unchdr

	; match pixel scales between image and kernel
	new_kernel = z0mgconv_matchpixscale(kernel,ker_scale[0],image_scale[0])

	; make sure new kernel is centered
	z0mgconv_ensure_psf_centered,new_kernel

	; make sure new kernel is normalized
	new_kernel = new_kernel/total(new_kernel)

	; locate NaNs and fill in 
	outim = im
	naninds = where(finite(outim) eq 0,nanct,complement=okinds)
	if nanct gt 0 then outim[naninds] = 0d
	
	if keyword_set(unc_file) then begin
		outuncim = uncim
		unaninds = where(finite(outuncim) eq 0,unanct,complement=uokinds)
		if unanct gt 0 then outuncim[unaninds] = 0d
	endif

	; do the convolution
	if n_elements(imsize) gt 2 then BEGIN

		for i=0,imsize[2]-1 do BEGIN
			slice = outim[*,*,i]
			outslice = convolve(double(slice),new_kernel,/no_pad)
			outim[*,*,i] = outslice

			if keyword_set(unc_file) then BEGIN
				uslice = (outuncim[*,*,i])^2.
				outuslice = convolve(double(uslice),new_kernel,/no_pad)
				outuncim[*,*,i] = sqrt(outuslice)*(start_fwhm/end_fwhm)
			endif
		endfor
	
	endif else BEGIN
		
		outim = convolve(double(outim),new_kernel,/no_pad)

		if keyword_set(unc_file) then BEGIN
			outuncim = sqrt(convolve(outuncim^2.,new_kernel,/no_pad))*(start_fwhm/end_fwhm)
		endif

	endelse

	; put the NaNs back in
	if nanct gt 0 then outim[naninds] = !values.f_nan
	if keyword_set(unc_file) then BEGIN
		if unanct gt 0 then outuncim[unaninds] = !values.f_nan
	endif

	; add comment to hdr
	sxaddpar,outhdr,'COMMENT','Convolved with z0mgconv.'
	sxaddpar,outunchdr,'COMMENT','Convolved with z0mgconv.'

	; write output file
	if keyword_set(kfspec_format) then BEGIN
		mwrfits,0,out_file,fullhdr,/create
		mwrfits,outim,out_file,outhdr
		mwrfits,index,out_file,hdrind
	endif else BEGIN
		writefits,out_file,outim,outhdr
	endelse

	if keyword_set(unc_file) then BEGIN
		if keyword_set(kfspec_format) then BEGIN
			mwrfits,0,outunc_file,ufullhdr,/create
			mwrfits,outuncim,outunc_file,outunchdr
			mwrfits,uindex,outunc_file,uhdrind
		endif else BEGIN
			writefits,outunc_file,outuncim,outunchdr
		endelse
	endif

end

