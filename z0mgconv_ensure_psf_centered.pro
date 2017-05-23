
;+
; NAME:
; 	z0mgconv_ensure_psf_centered
;
; PURPOSE:
; 	center kernels and psfs before using to convolve
;
; INPUTS: 
; 	image: the psf or kernel in question
;
; KEYWORDS:
;
; OUTPUTS: 
; 	centered image
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Stolen from K. Gordon by Karin (9 Aug 2013)
; 	Attempting to be responsible & version control with GIT
;
; To Do:
;
;-

pro z0mgconv_ensure_psf_centered,image

	size_image = size(image)
	; check that image is square if not make it square
	; by adding rows or columns of zeroes
	if (size_image[1] NE size_image[2]) then begin
	    if (size_image[1] GT size_image[2]) then begin
	        new_image = replicate(0.0,size_image[1],size_image[1])
	    endif else begin
	        new_image = replicate(0.0,size_image[2],size_image[2])
	    endelse
	    new_image[0:size_image[1]-1,0:size_image[2]-1] = image
	    image = new_image
	    size_image = size(image)
	endif

	; make sure image has an odd number of elements
	if ((size_image[1] mod 2) EQ 0) then begin
	    new_image = replicate(0.0,size_image[1]+1,size_image[2]+1)
	    new_image[0:size_image[1]-1,0:size_image[2]-1] = image
	    image = new_image
	    size_image = size(image)
	endif

	; now find where the maxium in the image is at
	indxs = where(finite(image) EQ 1,n_indxs)
	if (n_indxs GT 0) then begin
	    indxs = where(image EQ max(image[indxs]))
	endif else begin
	    print,'no finite elements to psf image'
	    stop
	endelse
	max_y = indxs[0]/size_image[1]
	max_x = indxs[0] mod size_image[1]
	print,'maximum at (' + strtrim(string(max_x),2) + $
		',' + strtrim(string(max_y),2) + ')'
	print,'image size is (' + strtrim(string(size_image[1]),2) + $ 
		',' + strtrim(string(size_image[2]),2) + ')'

	; determine the needed shifts
	shift_x = 0
	shift_y = 0
	if (max_x NE size_image[1]/2) then shift_x = size_image[1]/2 - max_x
	if (max_y NE size_image[2]/2) then shift_y = size_image[2]/2 - max_y

	; make the shift if nonzero
	if ((shift_x NE 0) OR (shift_y NE 0)) then begin
	    print,'shifting center by (' + strtrim(string(shift_x),2) + $
			',' + strtrim(string(shift_y),2) + ')'
	    new_image = replicate(0.0,size_image[1],size_image[2])
	    if (shift_x GT 0) then begin
	        x1 = 0
	        x2 = size_image[1] - abs(shift_x) - 1
	        nx1 = abs(shift_x)
	        nx2 = size_image[1] - 1
	    endif else begin
	        nx1 = 0
	        nx2 = size_image[1] - abs(shift_x) - 1
	        x1 = abs(shift_x)
	        x2 = size_image[1] - 1
	    endelse
	    if (shift_y GT 0) then begin
	        y1 = 0
	        y2 = size_image[2] - abs(shift_y) - 1
	        ny1 = abs(shift_y)
	        ny2 = size_image[2] - 1
	    endif else begin
	        ny1 = 0
	        ny2 = size_image[2] - abs(shift_y) - 1
		    y1 = abs(shift_y)
        	y2 = size_image[2] - 1
    	endelse

	    new_image[nx1:nx2,ny1:ny2] = image[x1:x2,y1:y2]
	
	    image = new_image
	endif

end
