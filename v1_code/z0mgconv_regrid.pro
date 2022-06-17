
;+
; NAME: 
; 	z0mgconv_regrid
;
; PURPOSE:
; 	regrids an input image to match a given xy grid
;
; INPUTS:
; 	arr: array to regrid
; 	xarr, yarr: arrays of x and y locations
; 	xout, yout: arrays for output x and y locations
;
; KEYWORDS:
;
; OUTPUTS:
; 	regridded array
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Modified from previous code by Karin (8 Aug 2013)
; 	Attempting to be responsible & version control with GIT
;
; To Do:
; 	- not clear that this works with increasing and decreasing
; 		pixel scale
; 	- add different interpolation options?
;-

function z0mgconv_regrid,arr,xarr,yarr,xout,yout

	; get pix scale
	dxin = xarr[1,0] - xarr[0,0]
	dyin = yarr[0,1] - yarr[0,0]
	
	size_out = size(xout,/dimen)
	xout_0 = xout[0,0] & yout_0 = yout[0,0]
	xout_max = xout[size_out[0]-1,size_out[1]-1] 
	yout_max = yout[size_out[0]-1,size_out[1]-1]

	; get x fractional index of xout_0
	close_xout_0 = value_locate(reform(xarr[*,0]),xout_0)
	if close_xout_0 eq -1 then close_xout_0 = 0
	fi_xout_0 = close_xout_0 + ((xout_0 - xarr[close_xout_0,0])/dxin)

	; get x fractional index of xout_max
	close_xout_max = value_locate(reform(xarr[*,0]),xout_max)
	fi_xout_max = close_xout_max + ((xout_max - xarr[close_xout_max,0])/dxin)

	; create vector of x fractional indices
;	fi_xout = scale_vector(dindgen(size_out[0]),fi_xout_0,fi_xout_max)
	fi_xout = cgScaleVector(dindgen(size_out[0]),fi_xout_0,fi_xout_max)

	; get y fractional index of yout_0
	close_yout_0 = value_locate(reform(yarr[0,*]),yout_0)
	if close_yout_0 eq -1 then close_yout_0 = 0
	fi_yout_0 = close_yout_0 + ((yout_0 - yarr[0,close_yout_0])/dyin)
	
	; get y fractional index of yout_max
	close_yout_max = value_locate(reform(yarr[0,*]),yout_max)
	fi_yout_max = close_yout_max + ((yout_max - yarr[0,close_yout_max])/dyin)
	
	; create vector of y fractional indices
;	fi_yout = scale_vector(dindgen(size_out[1]),fi_yout_0,fi_yout_max)
	fi_yout = cgScaleVector(dindgen(size_out[1]),fi_yout_0,fi_yout_max)

	fi_xarr = rebin(fi_xout,size_out[0],size_out[1])
	fi_yarr = rebin(reform(fi_yout,1,size_out[1]),size_out[0],size_out[1])

	; do the interpolation
	newarr = interpolate(arr,fi_xarr,fi_yarr,missing=0d)

	return,newarr

end
	
