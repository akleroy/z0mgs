; This does EVERYTHING - your probably do NOT want to do this.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD TARGET LIST
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Rebuild our galaxy list from the ground up (the exclusion is the
; LMC, which breaks things on image size).

build_galaxy_list, /rebuild $
                   , exclude = ['PGC17223']

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GALEX ATLAS CONSTRUCTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Make the GALEX cutouts
build_galex_atlas, /cutouts

; Make r25-based masks
build_galex_atlas, /mask

; Record the integration times and look for missing images. Record
; this in the header.
build_galex_atlas, /inventory

; Make the 15" version of the atlas
build_galex_atlas, /convol

; Subtract a background
build_galex_atlas, /bksub

; Compile statistics for the GALEX images and save them in the header
build_galex_atlas, /stats

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; UNWISE ATLAS CONSTRUCTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Convert the unwise images to have the right units
build_unwise_atlas, /units

; Build r25 masks on the unwise astrometry
build_unwise_atlas, /mask

; Carry out an iterative plane fit to remove a local background
build_unwise_atlas, /bksub

; Convolve to have a 15" Gaussian PSF
build_unwise_atlas, /convol

; Carry out a second background subtraction after the convolution
build_unwise_atlas, /second

; Measure statistics for the final images and record them in the header
build_unwise_atlas, /stats

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COMPILE A RELEASE AND REDUCE PIXEL SCALE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Wipe the previous delivery
build_delivery, /reset

; Copy and down-sample the GALEX data
build_delivery, /galex

; Copy and down-sample the UNWISE data
build_delivery, /wise

; Index the atlas into a large table
build_delivery, /index

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD THE INTENSITY DATABASE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Sample the maps and build the intensity table.
build_intensity_table




