# Module to accept a table that's the result of a GAIA query and
# extract cutouts at the location of each star.

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Imports
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import numpy as np
from reproject import reproject_interp
from astropy.table import Table
import astropy.io.fits as fits
import astropy.wcs as wcs
from astropy.utils.console import ProgressBar

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Cutout machinery
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def build_centered_header(
        ra_ctr=None, dec_ctr=None,
        pix_scale=None, half_width_pix=None,
        return_header=False):
    """Build a fiducial 2-d header centered at the specified (ra, dec)
    location.
    """
    full_width_pix = 2*half_width_pix+1
    
    hdu = fits.PrimaryHDU()

    hdu.header['NAXIS'] = 2
    hdu.header['NAXIS1'] = full_width_pix
    hdu.header['NAXIS2'] = full_width_pix
    hdu.header['CTYPE1'] = 'RA---TAN'
    hdu.header['CRVAL1'] = ra_ctr
    hdu.header['CRPIX1'] = half_width_pix+1
    hdu.header['CDELT1'] = -1.0*pix_scale
    hdu.header['CTYPE2'] = 'DEC--TAN'
    hdu.header['CRVAL2'] = dec_ctr
    hdu.header['CRPIX2'] = half_width_pix+1
    hdu.header['CDELT2'] = pix_scale
    hdu.header['EQUINOX'] = 2000.

    if return_header:
        return(hdu.header)
    else:
        return(hdu)

def extract_one_gaia_cutout(
        hdu_image, hdu_target,
        order='bilinear', missing=np.nan, return_hdu=True):
    """Extract a single cutout from an image using (by default) bilinear
    interpolation. If requested return a new HDU.
    """
    
    reprojected_image, footprint = reproject_interp(
        hdu_image, hdu_target.header, order=order
    )
    reprojected_image[footprint == 0] = missing

    if return_hdu == False:
        return(reprojected_image)
    
    new_header = hdu_target.header.copy()
    for this_key in ['BMAJ','BMIN','BPA','BUNIT']:
        try:
            new_header[this_key] = hdu_image.header[this_key]
        except KeyError:
            pass
    new_hdu = fits.PrimaryHDU(reprojected_image, new_header)

    return(new_hdu)

def extract_gaia_stack(
        image_fname,
        gaia_table_fname=None,
        out_fname=None, overwrite=True,
        oversamp_fac=2.0, half_width_pix=20,
        order='bilinear'):
    """Extract a stack of gaia cutouts an place them into a three-d stack
    in a relative astrometry.
    """

    # Open the target image
    hdu_image = fits.open(image_fname)[0]

    # Extract the pixel scale
    this_wcs = wcs.WCS(hdu_image)    
    pix_scale = wcs.utils.proj_plane_pixel_scales(this_wcs)
    target_pix_scale = pix_scale[0]/oversamp_fac
    print("Pixel scale [as] of image, target: ",
          pix_scale[0]*3600.,target_pix_scale*3600.)

    # Open the GAIA table (or any other table with ra, dec)
    gaia_table = Table.read(gaia_table_fname)
    n_rows = len(gaia_table)

    # Initialize the stack
    stack_hdu = build_centered_header(
        ra_ctr=0.0, dec_ctr=0.0,
        pix_scale=target_pix_scale,
        half_width_pix=half_width_pix)
    stack_hdu.header['NAXIS'] = 3
    stack_hdu.header['NAXIS3'] = n_rows
    stack_hdu.header['CRVAL3'] = 0
    stack_hdu.header['CRPIX3'] = 1
    stack_hdu.header['CDELT3'] = 1
    stack_hdu.header['CTYPE3'] = 'PLANE'
    full_width_pix = 2*half_width_pix+1
    stack = np.zeros((n_rows,
                      full_width_pix,
                      full_width_pix))*np.nan

    # Loop over all cutout coordinate locations
    ii = 0
    for this_row in ProgressBar(gaia_table):

        # Build the target astrometry for this object
        this_ra = float(this_row['ra'])
        this_dec = float(this_row['dec'])
        target_hdu = build_centered_header(
            ra_ctr=this_ra, dec_ctr=this_dec,
            pix_scale=target_pix_scale,
            half_width_pix=half_width_pix)

        # Extract the cutout, here we don't need a new HDU
        extracted_image = extract_one_gaia_cutout(
            hdu_image, target_hdu,
            order=order, return_hdu=False)
        stack[ii,:,:] = extracted_image
        ii += 1

    # Form an HDU
    stack_hdu = fits.PrimaryHDU(stack, stack_hdu.header)
    
    # Write if requested then return
    if out_fname is not None:
        stack_hdu.writeto(out_fname, overwrite=overwrite)
    
    return(stack_hdu)
        
