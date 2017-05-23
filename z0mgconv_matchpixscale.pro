
;+
; NAME:
; 	z0mgconv_matchpixscale
;
; PURPOSE:
; 	rebins an input image to have the requested pixel scale
;
; INPUTS:
; 	inpsf: input image to rebin
; 	indx: input pixel scale
; 	refdx: requested pixel scale
;
; KEYWORDS:
;	normalize: normalize the output image to have total of 1
;		useful for psf images
;
; OUTPUTS:
; 	rebinned image
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Fixed up from previous code by Karin (8 Aug 2013)
; 	Attempting to be responsible & verison control with Git.
;
; To Do:
; 	- needs a lot of checking
;	- does this handle increase and decrease of scale?
;
;-

function z0mgconv_matchpixscale,inpsf,indx,refdx,$
	normalize=normalize

	; find size of input image
	insize = size(inpsf,/dimen)
	
	; calculate what size the rebinned image has
	refsize = (insize*indx)/refdx

	; make sure the rebinned image will have odd numbers of pixels
	if floor(refsize[0]) mod 2 eq 0 then $
		refsize=floor(refsize)+1 else refsize=floor(refsize)

	; pick out center pixel
	maxlocin = insize/2
	maxlocref = refsize/2

	; create x and y index vectors
	xvecin = dindgen(insize[0])*indx
	xvecin = xvecin - xvecin[maxlocin[0]]
	yvecin = dindgen(insize[1])*indx
	yvecin = yvecin - yvecin[maxlocin[1]]
	xarrin = rebin(xvecin,insize[0],insize[1])
	yarrin = rebin(reform(yvecin,1,insize[1]),insize[0],insize[1])

	; do the same for the requested scale
	xvecref = dindgen(refsize[0])*refdx
	xvecref = xvecref - xvecref[maxlocref[0]]
	yvecref = dindgen(refsize[1])*refdx
	yvecref = yvecref - yvecref[maxlocref[1]]
	xarrref = rebin(xvecref,refsize[0],refsize[1])
	yarrref = rebin(reform(yvecref,1,refsize[1]),refsize[0],refsize[1])

	; call rebinning program
	newpsf = z0mgconv_regrid(inpsf,xarrin,yarrin,xarrref,yarrref)

	; renormalize
	if keyword_set(normalize) then newpsf = newpsf/total(newpsf)

	return,newpsf
end
