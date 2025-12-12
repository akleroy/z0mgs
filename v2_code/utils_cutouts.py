# Routines to extract cutouts for specific surveys. contains a lot of 

def extract_galex_stamp(
        band = 'fuv',
        ra_ctr = 0.0,
        dec_ctr = 0.0,
        size_deg = 0.01,
        index_file = '../../working_data/galex/index/galex_tile_index.fits',
        index_tab = None,
        use_int_files = True,
        outfile_image = None,
        outfile_weight = None,
        show = False,
        overwrite = True):
    """
    This builds one GALEX image from the individual calibrated tiles. 

    This is slow.
    """
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Make a target header
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    center_coord = SkyCoord(ra=ra_ctr*u.deg, dec=dec_ctr*u.deg, frame='icrs')    
    pix_scale = np.array([1.5/3600.,1.5/3600.])
    nx = int(np.ceil(size_deg / pix_scale[0]))
    ny = int(np.ceil(size_deg / pix_scale[1]))
    print("... pixel scale, nx, ny: ", pix_scale, nx, ny)
    
    target_hdr = make_simple_header(
        center_coord, pix_scale, nx=nx, ny=ny, return_header=True)    
    target_hdr['BUNIT'] = 'MJy/sr'
    target_wcs = wcs.WCS(target_hdr)
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Find contributing tiles
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    overlap_tab = find_index_overlap(
        index_file = index_file,
        center_coord = center_coord,
        image_extent = size_deg,
        selection_dict = {'filter':band},
    )    
    n_overlap = len(overlap_tab)
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Initialize output
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    weight_image = np.zeros((ny, nx))
    sum_image = np.zeros((ny, nx))

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Loop over tiles
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    
    for ii, this_overlap_row in enumerate(overlap_tab):
        
        print("... processing tile ", ii, " of ", n_overlap)

        # Can use background subtracted or integrated images, we
        # prefer to do the background subtraction ourself.
        
        if use_int_files:
            this_image_fname = this_overlap_row['fname'].strip()
        else:
            this_image_fname = this_overlap_row['bgsub_fname'].strip()            

        # Relative response file name
        this_rrhr_fname = this_overlap_row['rrhr_fname'].strip()

        # Flag file name
        this_flag_fname = this_overlap_row['flag_fname'].strip()

        # Read in the image
        image_hdu = fits.open(this_image_fname)
        image = image_hdu[0].data
        
        # Check if the image is empty (in which case we skip this
        # file). Could refine the index to do this.
        
        if np.sum(np.isfinite(image)*(image != 0)) == 0:
            print("... ... no finite or non-zero values. Proceeding.")
            continue

        tile_hdr = image_hdu[0].header
        if 'CDELT1' not in tile_hdr:
            print("... ... no astrometry found. Proceeding.")

        # Read in the response
        rrhr_hdu = fits.open(this_rrhr_fname)
        rrhr = rrhr_hdu[0].data
        rrhr_hdr = rrhr_hdu[0].header

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        # Read in and apply flags
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        print("... ... reading and applying flags.")
        
        # Read in the flags
        flag_hdu = fits.open(this_flag_fname)

        # Align the flags to the image using nearest neighbor
        # alignment to preserve precise pixel values.
        tile_wcs = wcs.WCS(tile_hdr)
        aligned_flags, footprint = reproject_interp(
            flag_hdu[0], tile_wcs,
            order='nearest-neighbor')
        aligned_flags[np.where(footprint == 0)] = 1024

        flag_mask = (((aligned_flags % 256) % 128) % 64) >= 32

        image[flag_mask] = np.nan
        rrhr[flag_mask] = np.nan
        
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        # Convert units
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        print("... ... unit conversion.")
                
        # Counts per second to Jy transformation. Based on
        # http://galexgi.gsfc.nasa.gov/docs/galex/FAQ/counts_background.html
        #
        # for GALEX FUV:
        # mAB = -2.5 log10 CPS + 18.82
        # for GALEX NUV:
        # mAB = -2.5 log10 CPS + 20.08
        #
        # then mAB = -2.5 log10 (f_nu / 3630.8)

        if band == 'fuv':
            cpstojy = 0.000107647
            target_hdr['CPSTOJY'] = (cpstojy, 'FUV value')
        else:
            cpstojy = 3.37289e-05
            target_hdr['CPSTOJY'] = (cpstojy, 'NUV value')

        image = cpstojy*image/1E6            
            
        # Assume square and will be decimal degrees
        pix_scale_deg = proj_plane_pixel_scales(tile_wcs)[0]
        pix_sr = (np.pi/180.*pix_scale_deg)**2
        image /= pix_sr

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        # Align the image and response to the target header
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

        print("... ... align and accumulate.")
        
        # Image
        aligned_image, footprint = reproject_interp(
            (image, tile_hdr), target_wcs,
            order='bilinear')        
        #aligned_image, footprint = reproject_adaptive(
        #    (image, tile_hdr), target_wcs,
        #    bad_value_mode='ignore')
        aligned_image[np.where(footprint == 0)] = np.nan

        # RRHR
        aligned_rrhr, footprint = reproject_interp(
            (rrhr, rrhr_hdr), target_wcs,
            order='bilinear')
        #aligned_rrhr, footprint = reproject_adaptive(
        #    (rrhr, tile_hdr), target_wcs,
        #    bad_value_mode='ignore')
        aligned_rrhr[np.where(footprint == 0)] = 0.0
        
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        # Accumulate
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        im_ind = \
            np.where((aligned_rrhr > 0.0)*np.isfinite(aligned_rrhr)* \
                     np.isfinite(aligned_image))

        sum_image[im_ind] = sum_image[im_ind]+ \
            (aligned_image[im_ind]*aligned_rrhr[im_ind])
        weight_image[im_ind] = weight_image[im_ind]+ \
            aligned_rrhr[im_ind]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Construct final image
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    final_image = sum_image/weight_image

    # Construct HDUs and write if requested
    
    final_image_hdu = fits.PrimaryHDU(data=final_image, header=target_hdr)
    if outfile_image is not None:
        final_image_hdu.writeto(outfile_image, overwrite=overwrite)

    response_hdr = target_hdr
    response_hdr['BUNIT'] = 'WEIGHT'
    final_weight_hdu = fits.PrimaryHDU(data=weight_image, header=response_hdr)
    if outfile_weight is not None:
        final_weight_hdu.writeto(outfile_weight, overwrite=overwrite)

    return(final_image_hdu, final_weight_hdu)
