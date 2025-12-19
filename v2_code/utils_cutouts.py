# Routines to extract cutouts for specific surveys. contains a lot of
# specifics related to those data sets.

import os

from astropy.table import Table
from astropy.io import fits
from astropy import wcs
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.stats import mad_std

from scipy.ndimage import binary_dilation

from reproject import reproject_interp, reproject_adaptive

from utils_z0mgs_images import *

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# GALEX
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def extract_galex_stamp(
        band = 'fuv',
        ctr_ra = 0.0,
        ctr_dec = 0.0,
        size_deg = 0.01,
        index_file = '../../working_data/galex/index/galex_tile_index.fits',
        index_tab = None,
        use_int_files = True,
        outfile_image = None,
        outfile_weight = None,
        overwrite = True):
    """
    This builds one GALEX image from the individual calibrated tiles.

    This is slow.
    """
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Make a target header
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    center_coord = SkyCoord(ra=ctr_ra*u.deg, dec=ctr_dec*u.deg, frame='icrs')    
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

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# UNWISE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def make_invvar_mask(
        invvar_fname = None,
        invvar_hdu = None,
        invvar_image = None,
        thresh = 5.0):
    """
    Construct a mask based on the inverse variance image.
    """

    if invvar_image is None:

        if invvar_hdu is None:
            invvar_hdu = fits.open(invvar_fname)[0]

            invvar_image = invvar_hdu.data

    rms = np.nanstd(invvar_image)
    med_val = np.nanmedian(invvar_image)
    resid = (invvar_image - med_val)/rms
    mask = binary_dilation(resid >= thresh, iterations=5)
            
    return(mask)

def unwise_counts_to_mjysr(
        band = 'w1',
        pix_scale_as = 2.75):
    """
    Convert counts to MJy/sr.
    """

    vega_to_ab = {
        'w1':2.683,
        'w2':3.319,
        'w3':5.242,
        'w4':6.604,
    }

    norm_mag_vega = 22.5

    countspix_to_jypix = \
        10.**(-1.*(norm_mag_vega+vega_to_ab[band])/2.5)*3631.

    jypix_to_mjysr = 1./1E6/(pix_scale_as/3600.*np.pi/180.)**2

    countspix_to_mjysr = countspix_to_jypix*jypix_to_mjysr
    
    return(countspix_to_mjysr)

def extract_unwise_stamp(
        band = 'w1',
        ctr_ra = 0.0,
        ctr_dec = 0.0,
        method = 'copy',
        tol_deg = 10./3600.,
        size_deg = 0.01,
        survey = 'unwise',
        index_dir = '../../working_data/unwise/index/',
        index_file = None,
        index_tab = None,
        outfile_image = None,
        outfile_mask = None,
        overwrite = True,
):
    """
    """

    valid_methods = ['copy','reproject']
    if method not in valid_methods:
        print("Invalid method: ", method)
        print("... allowed methods ", valid_methods)
        return(None)

    valid_surveys = ['unwise', 'unwise_custom', 'allwise', 'neowise']
    if survey not in valid_surveys:
        print("Invalid survey: ", survey)
        print("... allowed surveys ", valid_surveys)
        return(None)

    if survey == 'unwise':
        index_file = index_dir + 'unwise_custom_index.fits'

    if survey == 'unwise_custom':
        index_file = index_dir + 'unwise_custom_index.fits'

    if survey == 'allwise':
        index_file = index_dir + 'unwise_allwise_index.fits'

    if survey == 'neowise':
        index_file = index_dir + 'unwise_neowise_index.fits'
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Copying case
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if method == 'copy':

        print("Copying closest matching tile...")

        center_coord = SkyCoord(ra=ctr_ra*u.deg, dec=ctr_dec*u.deg, frame='icrs')
        
        overlap_tab, separations = \
            find_index_overlap(
                index_file = index_file,
                center_coord = center_coord,
                image_extent = size_deg,
                selection_dict = {'filter':band},
                force_tolerance = tol_deg,
                return_separations = True,
            )    
        n_overlap = len(overlap_tab)

        print("... found ", n_overlap, " tiles")

        if n_overlap > 1:

            min_ind = np.argmin(separations)
            overlap_tab = overlap_tab[min_ind]
            separations = separations[min_ind]

        print("Closest tile center distance: ", separations)
            
        this_fname = overlap_tab[0]['fname'].strip()
        this_hdu = fits.open(this_fname)[0]
        this_hdr = this_hdu.header
        this_image = this_hdu.data

        print("... converting units.")
        this_wcs = wcs.WCS(this_hdr)
        pix_scale_deg = proj_plane_pixel_scales(this_wcs)[0]
        pix_scale_as = pix_scale_deg*3600.
        print("... ... pixel scale [as]: ", pix_scale_as)
        
        conv_fac = \
            unwise_counts_to_mjysr(
                band = band,
                pix_scale_as = pix_scale_as)

        this_image *= conv_fac
        this_hdr['BUNIT'] = 'MJy/sr'
        
        print("... making an inverse variance mask.")
        this_invvar_fname = this_fname.replace(
            '-img-m.fits','-invvar-m.fits')
        if os.path.isfile(this_invvar_fname) == False:
            this_invvar_fname = this_fname.replace(
                '-img-m.fits','-invvar-m.fits.gz')
        if os.path.isfile(this_invvar_fname) == False:
            print("Cannot find invvar file")
                
        invvar_mask = make_invvar_mask(
            invvar_fname=this_invvar_fname)
        
        # ................................................
        # Write to disk
        # ................................................
        
        new_image_hdu = fits.PrimaryHDU(data=this_image, header=this_hdr)
        
        if outfile_image is not None:
            new_image_hdu.writeto(outfile_image, overwrite=overwrite)

        mask_hdr = this_hdr
        mask_hdr['BUNIT'] = 'Mask'
        invvar_mask_hdu = fits.PrimaryHDU(data=invvar_mask.astype(int), header=mask_hdr)

        if outfile_mask is not None:        
            invvar_mask_hdu.writeto(outfile_mask, overwrite=overwrite)

        return(new_image_hdu, invvar_mask_hdu)
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Reprojection case
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if method == 'reproject':

        # ................................................
        # Define header
        # ................................................        
        
        center_coord = SkyCoord(ra=ctr_ra*u.deg, dec=ctr_dec*u.deg, frame='icrs')    
        pix_scale = np.array([2.75/3600.,2.75/3600.])
        nx = int(np.ceil(size_deg / pix_scale[0]))
        ny = int(np.ceil(size_deg / pix_scale[1]))
        print("... pixel scale, nx, ny: ", pix_scale, nx, ny)
    
        target_hdr = make_simple_header(
            center_coord, pix_scale, nx=nx, ny=ny, return_header=True)    
        target_hdr['BUNIT'] = 'MJy/sr'
        target_wcs = wcs.WCS(target_hdr)
    
        # ................................................
        # Find overlap
        # ................................................

        print("Finding all overlapping tiles within ", size_deg*60., " arcmin ...")
        
        overlap_tab = \
            find_index_overlap(
                index_file = index_file,
                center_coord = center_coord,
                image_extent = size_deg,
                selection_dict = {'filter':band},
            )
        print(overlap_tab)
        n_overlap = len(overlap_tab)

        print("... found ", n_overlap, " tiles")
        
        # ................................................
        # Initialize output
        # ................................................

        weight_image = np.zeros((ny, nx))
        sum_image = np.zeros((ny, nx))
        mask_image = np.zeros((ny, nx))
        
        # ................................................
        # Loop over tiles, convert, mask then reproject
        # ................................................        

        for ii, this_overlap_row in enumerate(overlap_tab):

            print("... processing frame ", ii, " of ", n_overlap)

            # Read and convert units
            this_fname = this_overlap_row['fname'].strip()
            this_hdu = fits.open(this_fname)[0]
            this_hdr = this_hdu.header
            this_image = this_hdu.data
            
            print("... converting units.")
            this_wcs = wcs.WCS(this_hdr)
            pix_scale_deg = proj_plane_pixel_scales(this_wcs)[0]
            pix_scale_as = pix_scale_deg*3600.
            print("... ... pixel scale [as]: ", pix_scale_as)
        
            conv_fac = \
                unwise_counts_to_mjysr(
                    band = band,
                    pix_scale_as = pix_scale_as)

            this_image *= conv_fac
            this_hdr['BUNIT'] = 'MJy/sr'

            # Make a mask based on the inverse variance
            
            print("... making an inverse variance mask.")
            
            this_invvar_fname = this_fname.replace(
                '-img-m.fits','-invvar-m.fits')
            if os.path.isfile(this_invvar_fname) == False:
                this_invvar_fname = this_fname.replace(
                    '-img-m.fits','-invvar-m.fits.gz')
            if os.path.isfile(this_invvar_fname) == False:
                print("Cannot find invvar file")
            
            invvar_mask = make_invvar_mask(
                invvar_fname=this_invvar_fname)

            print("... reproject and accumulate.")
            
            # Reproject image
            aligned_image, footprint_image = reproject_interp(
                (this_image, this_hdr), target_wcs,
                order='bilinear')
            aligned_image[footprint_image==0] = 0.0
            
            print("Sum - ", np.sum(footprint_image), np.sum(aligned_image))
            
            # Reproject mask
            aligned_mask, footprint_mask = reproject_interp(
                (invvar_mask, this_hdr), target_wcs,
                order='nearest-neighbor')
            aligned_mask[footprint_mask==0] = 0.0
            
            # Accumulate
            sum_image += aligned_image
            weight_image += footprint_image*1.0
            mask_image += aligned_mask
            
        # ................................................
        # Create full image and write to disk
        # ................................................        

        full_image = sum_image / weight_image
        full_mask = (mask_image >= 1.0)
        
        full_image_hdu = fits.PrimaryHDU(data=full_image, header=target_hdr)
        
        if outfile_image is not None:
            full_image_hdu.writeto(outfile_image, overwrite=overwrite)

        mask_hdr = target_hdr
        mask_hdr['BUNIT'] = 'Mask'
        
        full_mask_hdu = fits.PrimaryHDU(data=full_mask.astype(int),
                                          header=mask_hdr)
        
        if outfile_mask is not None:        
            full_mask_hdu.writeto(outfile_mask, overwrite=overwrite)
            
        return(full_image_hdu, full_mask_hdu)
    
    print("Should not reach here.")
    return(None)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# SDSS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def extract_sdss_stamp(
        band = 'g',
        ctr_ra = 0.0,
        ctr_dec = 0.0,
        size_deg = 0.01,
        index_file = '../../working_data/sdss/index/sdss_tile_index.fits',
        index_tab = None,
        use_int_files = True,
        outfile_image = None,
        outfile_weight = None,
        overwrite = True):

    pass
