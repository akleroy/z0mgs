#!/usr/bin/env python3

import os, glob, sys
import numpy as np
import warnings
from enum import IntFlag, auto

import urllib.request
import urllib.error

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.table import Table, QTable
from astropy import units as u
from astropy.stats import sigma_clipped_stats
import astropy.constants as const

# Astroquery for IRSA access
from astroquery.ipac.irsa import Irsa

# Reproject for image alignment
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs

from astropy.utils.console import ProgressBar

#import warnings
#warnings.filterwarnings('ignore', category=UserWarning)
#import astropy.utils.exceptions
#warnings.simplefilter('ignore', category=astropy.utils.exceptions.AstropyWarning)                     

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to query data from IRSA
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def check_spherex_collections():
    """Query IRSA's collection list to see which collections hold spherex
    data. Use this to inform the image search.
    """
    
    collections = Irsa.list_collections(servicetype='SIA', filter='spherex')

    print(f"Found {len(collections)} SPHEREX collections:")
    for collection in collections['collection']:
        print(f"  - {collection}")

    return()


def search_spherex_images(
        target: str = None,
        coordinates: SkyCoord = None,
        radius: u.Quantity = 10*u.arcmin,
        collection: str = 'spherex_qr2',
        verbose: bool = True) -> Table:
    
    """
    From Adam Ginsburg (pared down a lot)

    Search for SPHEREX images in IRSA archive.

    Parameters
    ----------
        
    target : str, optional
    Target name (e.g., 'M31', 'NGC 1234')
    
    coordinates : SkyCoord, optional
    Target coordinates

    radius : Quantity
    Search radius

    collection : str (default 'spherex_qr2')
    Collection holding the images, will need to be managed    

    Returns
    -------
    Table
    
    Table of available SPHEREX images
    """
    
    if verbose:
        print(f"Searching for SPHEREX images...")

    # Make sure we have target coordinates
    if coordinates is None and target is not None:
        coordinates = SkyCoord.from_name(target)
    elif coordinates is None:
        raise ValueError("Either target name or coordinates must be provided")

    if verbose:
        print(f"Target coordinates: {coordinates.to_string('hmsdms')}")
        print(f"Search radius: {radius}")

    # Search for images in the specified collection
    try:
        images = Irsa.query_sia(pos=(coordinates, radius), collection=collection)
    except Exception as e:
        if verbose:
            print(f"SIA query failed: {e}")

    if verbose:
        print(f"Found {len(images)} images")
        if len(images) > 0:
            print("Available columns (first 10): ", images.colnames[:10])  # Show first 10 columns

    return(images)


def download_images(
        images: Table,
        max_images: int = None,
        outdir: str = '',
        incremental: bool = True,
        verbose: bool = True) -> list[str]:
    """Download SPHEREX images from a table produced by the IRSA query above.

    Modified from version by Adam Ginsburg.
    
    Parameters
    ----------
    
    images : Table
             Table of images from search_spherex_images
        
    max_images : int, optional
                 Maximum number of images to download

    outdir : str, optional
             Location where the download directs

    incremental : bool, default True
                  Check if file is present before download

    verbose : bool, default True
              Turn on extra print statements
    
    Returns
    -------
    
    list[str] : List of downloaded file paths

    """

    # If an empty table, return
    if len(images) == 0:
        if verbose:
            print("No images to download")
        return([])

    # If we're capped in files to download impose that
    if max_images is not None:
        images = images[:max_images]

    # Initialize record
    downloaded_files = []

    # Loop and download
    if verbose:
        print(f"Downloading {len(images)} images...")
        if len(images) > 10:
            print("  (This may take a while for large datasets...)")

    # Loop over the list of images
    for ii, this_row in enumerate(images):

        # Determine the access URL        
        try:
            this_url = None
            # ... try the different potential column names
            for url_col in ['access_url', 'cloud_access', 'download_url', 'url']:
                if url_col in this_row.colnames and this_row[url_col]:
                    this_url = str(this_row[url_col])
                    break

            # ... catch the error case
            if this_url is None:
                if verbose:
                    print(f"  Skipping image {ii+1}: No access URL found")
                continue

            # Generate filename

            # AKL: Note that you cannot just use the obsid because two
            # images with different filters are obtained with the same
            # obsid - parse the full file name.

            # Observation ID + detector should be unique, but
            # processing time in the ID stamp is also potentially of
            # interest.
           
            #obs_id = this_row.get('obs_id', f'spherex_{ii:04d}')
            
            obs_fname = this_url.split('/')[-1]            
            this_fname = outdir + obs_fname

            # Check if the filename is present already and either delete
            # the file or proceed with a notification.
            if os.path.isfile(this_fname):
                if incremental:
                    if verbose:
                        print(f"  Image {ii+1}/{len(images)}: {this_fname} already exists.")
                        print(f"  Incremental is set to TRUE.")
                        downloaded_files.append(str(this_fname))
                    continue
                else:
                    if verbose:
                        print(f"  Image {ii+1}/{len(images)}: {this_fname} already exists.")
                        print(f"  Incremental is set to FALSE. Will proceed and overwrite.")
                    os.system('rm -rf '+this_fname)

            # Download the file
            if verbose:
                print(f"  Downloading image {ii+1}/{len(images)}: {this_fname}")

            # Download using urllib first, then open with astropy
            try:                
                # Download to a temporary location first
                temp_fname = this_fname + '.tmp'
                urllib.request.urlretrieve(this_url, temp_fname)
            
                # Verify it's a valid FITS file by opening it
                with fits.open(temp_fname) as this_hdu:                
                    # Save as final file
                    this_hdu.writeto(this_fname, overwrite=True)

                # Remove temporary file
                os.system('rm -rf '+temp_fname)

            except urllib.error.URLError as e:
                if verbose:
                    print(f"    URL error: {e}")
                continue            
            except Exception as fits_error:            
                # Clean up temp file if it exists
                temp_fname = this_fname + '.tmp'
                if os.path.isfile(temp_fname):
                    os.system('rm -rf '+temp_fname)
                    # Will crash the program?
                raise fits_error

            # If we made it here we're successful
            downloaded_files.append(str(this_fname))
        
        except Exception as e:
            if verbose:
                print(f"  Failed to download image {ii+1}: {e}")
            continue

    if verbose:
        print(f"Successfully downloaded {len(downloaded_files)} images")

    return(downloaded_files)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to support build a cube
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def make_wavelength_image(
        hdu_list = None,
        use_hdu = None,
):
    """Use astropy WCS to construct images of wavelength and bandwidth
    across a SPHEREx image.

    Feed in the HDUlist for an image and the HDU holding the image of
    interest. Program returns 

    lam, bw

    The central wavelength and bandwidth of the image in microns.

    Based on IRSA tutorial.
    """

    # Testing shows that the various image extensions do produce the
    # same output images.
    
    this_header = hdu_list[use_hdu].header.copy()
    this_shape = hdu_list[use_hdu].data.shape

    # Remove SIP coefficients. Mostly this just avoids an annoying
    # error message, since setting spectral_wcs.sip = None below turns
    # off their application anyways.
    
    keywords_to_remove = ['A_?_?', 'B_?_?', 'AP_?_?', 'BP_?_?',
                          '*_ORDER']
    for this_pattern in keywords_to_remove:
        for this_keyword in this_header.copy():
            if glob.fnmatch.fnmatch(this_keyword, this_pattern):
                del this_header[this_keyword]
    
    # Key call. Feed the header associated with the desired image,
    # also pass the parent HDU list to handle distortions/wavelength
    # lookup, and use the "W"avelength transform.

    # This all relies on astropy being wired under the hood to use the
    # WCS-WAVE table in the final extension.
    
    spectral_wcs = WCS(this_header, fobj=hdu_list, key="W")

    # Turn off spatial distortion terms
    spectral_wcs.sip = None

    # Create a grid of pixel coords
    yy, xx = np.indices(this_shape)
    
    # Evaluate the WCS to get the wavelength and bandwidth
    lam, bw = spectral_wcs.pixel_to_world(xx, yy)

    # Get the units (not used right now)
    lam_unit = spectral_wcs.wcs.cunit[0]
    bw_unit = spectral_wcs.wcs.cunit[1]

    # Return images of lambda and bandwidth
    return((lam, bw))

def make_cube_header(
        center_coord,
        pix_scale = 6. / 3600.,
        extent = None, extent_x = None, extent_y = None,
        nx = None, ny = None,
        lam_min = 0.7, lam_max = 5.2, lam_step = 0.02,
        return_header=False):
    """Make a #D FITS header centered on the coordinate of interest with a
    user-specififed pixel scale and extent and wavelength axis.

    
    Parameters
    ----------

    center_coord : `~astropy.coordinates.SkyCoord` object or
        array-like Sky coordinates of the image center. If array-like
        then (ra, dec) in decimal degrees assumes.

    pix_scale : required. Size in decimal degrees of a pixel. Can be
        an array in which case it is pixel scale along x and y (e.g.,
        as returned by proj_pixel_scales).

    extent : the angular extent of the image in decimal degrees for both x and y.

    extent_x : the angular extent of the image along the x coordinate

    extent_y : the angular extent of the image along the y coordinate

    nx : the number of x pixels (not needed with extent_x and pix_scale)

    ny : the number of y pixels (not needed with extent_y and pix_scale)

    lam_min, lam_max, lam_step : minimum, maximum, and channel width
        for wavelength axis in microns

    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Figure out the center
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if isinstance(center_coord, str):
        coordinates = SkyCoord.from_name(center_coord)
        ra_ctr = coordinates.ra.degree
        dec_ctr = coordinates.dec.degree
    elif isinstance(center_coord, SkyCoord):
        ra_ctr = center_coord.ra.degree
        dec_ctr = center_coord.dec.degree
    else:
        ra_ctr, dec_ctr = center_coord
        if hasattr(ra_ctr, 'unit'):
            ra_ctr = ra_ctr.to(u.deg).value
            dec_ctr = dec_ctr.to(u.deg).value

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Figure out the pixel scale
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    if pix_scale is None:
        print("Pixel scale not specified. Returning.")
        return()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Figure out image extent
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if extent is not None:
        extent_x = extent
        extent_y = extent
    
    if (nx is not None) and (ny is not None):
        if isinstance(pix_scale, np.ndarray):
            extent_x = pix_scale[0] * nx
            extent_y = pix_scale[1] * ny
        else:
            extent_x = pix_scale * nx
            extent_y = pix_scale * ny
    elif (extent_x is not None) and (extent_y is not None):
        if isinstance(pix_scale, np.ndarray):
            nx = int(np.ceil(extent_x*0.5 / pix_scale[0]) * 2 + 1)
            ny = int(np.ceil(extent_y*0.5 / pix_scale[1]) * 2 + 1)
        else:
            nx = int(np.ceil(extent_x*0.5 / pix_scale) * 2 + 1)
            ny = int(np.ceil(extent_y*0.5 / pix_scale) * 2 + 1)            
    else:
        print("Extent not specified. Returning.")
        return()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Make the wavelength axis
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    lam_array = np.arange(lam_min, lam_max + lam_step, lam_step)
    nz = len(lam_array)
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Make the FITS header
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    
    hdu = fits.PrimaryHDU()    
    hdu.header = fits.Header()
    hdu.header['NAXIS'] = 3

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Spatial information
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    hdu.header['NAXIS1'] = nx
    hdu.header['NAXIS2'] = ny
    
    hdu.header['CTYPE1'] = 'RA---TAN'
    hdu.header['CRVAL1'] = ra_ctr
    hdu.header['CRPIX1'] = np.float16((nx / 2) * 1 - 0.5)

    hdu.header['CTYPE2'] = 'DEC--TAN'
    hdu.header['CRVAL2'] = dec_ctr
    hdu.header['CRPIX2'] = np.float16((ny / 2) * 1 - 0.5)
    
    if isinstance(pix_scale, np.ndarray):    
        hdu.header['CDELT1'] = -1.0 * pix_scale[0]
        hdu.header['CDELT2'] = 1.0 * pix_scale[1]
    else:
        hdu.header['CDELT1'] = -1.0 * pix_scale
        hdu.header['CDELT2'] = 1.0 * pix_scale
            
    hdu.header['EQUINOX'] = 2000.0
    hdu.header['RADESYS'] = 'FK5'
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Spectral information
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    hdu.header['NAXIS3'] = nz
    hdu.header['CTYPE3'] = 'WAVE'
    hdu.header['CUNIT3'] = 'um'
    hdu.header['CRPIX3'] = 1
    hdu.header['CRVAL3'] = lam_array[0]
    hdu.header['CDELT3'] = lam_step

    # Return header or HDU    
    if return_header:
        return(hdu.header)
    else:
        return(hdu)

# Potentially useful but interface with numpy is not straightforward.

class SpherexFlag(IntFlag):
    TRANSIENT = auto()
    OVERFLOW = auto()
    SUR_ERROR = auto()
    BLANK_1 = auto()
    PHANTOM = auto()
    REFERENCE = auto()
    NONFUNC = auto()
    DICHROIC = auto()
    BLANK_2 = auto()
    MISSING_DATA = auto()
    HOT = auto()
    COLD = auto()    
    FULLSAMPLE = auto()
    BLANK_3 = auto()
    PHANMISS = auto()
    NONLINEAR = auto()
    BLANK_4 = auto()    
    PERSIST = auto()
    BLANK_5 = auto()        
    OUTLIER = auto()
    BLANK_6 = auto()        
    SOURCE = auto()

def spherex_flag_dict():
    """
    """

    # From the header, flag definitions
    
    #HIERARCH MP_TRANSIENT = 0 / Transient detected during SUR                       
    #HIERARCH MP_OVERFLOW = 1 / Overflow detected during SUR                         
    #HIERARCH MP_SUR_ERROR = 2 / Error in onboard processing                         
    #HIERARCH MP_PHANTOM = 4 / Phantom pixel                                         
    #HIERARCH MP_REFERENCE = 5 / Reference pixel                                     
    #HIERARCH MP_NONFUNC = 6 / Permanently unusable                                  
    #HIERARCH MP_DICHROIC = 7 / Low efficiency due to dichroic                       
    #HIERARCH MP_MISSING_DATA = 9 / Onboard data lost                                
    #MP_HOT  =                   10 / Hot pixel                                      
    #MP_COLD =                   11 / Anomalously low signal                         
    #HIERARCH MP_FULLSAMPLE = 12 / Pixel full sample history is available            
    #HIERARCH MP_PHANMISS = 14 / Phantom correction was not applied                  
    #HIERARCH MP_NONLINEAR = 15 / Nonlinearity correction cannot be applied reliably 
    #HIERARCH MP_PERSIST = 17 / Persistent charge above threshold                    
    #HIERARCH MP_OUTLIER = 19 / Pixel flagged by Detect Outliers                     
    #HIERARCH MP_SOURCE = 21 / Pixel mapped to a known source    
    
    this_dict = {
        'TRANSIENT' : 0,
        'OVERFLOW' : 1,
        'SUR_ERROR' : 2,
        'PHANTOM' : 4,
        'REFERENCE' : 5,
        'NONFUNC' : 6,
        'DICHROIC' : 7,
        'MISSING_DATA' : 9,
        'HOT' : 10,
        'COLD' : 11,
        'FULLSAMPLE' : 12,
        'PHANMISS' : 14,
        'NONLINEAR' : 15,
        'PERSIST' : 17,
        'OUTLIER' : 19,
        'SOURCE' : 21,        
    }
    
    return(this_dict)
    
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routine to actually build a cube
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def extract_spherex_sed(
        target_coord,        
        image_list = [],
        outfile = None,
        overwrite=True):
    """
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Get the coordinates
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    if isinstance(target_coord, str):
        coordinates = SkyCoord.from_name(target_coord)
        ra_deg = coordinates.ra.degree
        dec_deg = coordinates.dec.degree
    elif isinstance(target_coord, SkyCoord):
        ra_deg = target_coord.ra.degree
        dec_deg = target_coord.dec.degree
    else:
        ra_deg, dec_deg = target_coord
        if hasattr(ra_deg, 'unit'):
            ra_deg = ra_deg.to(u.deg).value
            dec_deg = dec_deg.to(u.deg).value

    target_coord = SkyCoord(
        ra = ra_deg * u.deg, dec=dec_deg* u.deg,
        frame = 'icrs')
                                        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Loop over the image list
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    val_ra = np.zeros(len(image_list))*np.nan
    zodi_ra = np.zeros(len(image_list))*np.nan
    lam_ra = np.zeros(len(image_list))*np.nan
    bw_ra = np.zeros(len(image_list))*np.nan
    file_list = []
    
    counter = 0
    
    for this_fname in ProgressBar(image_list):        
        
        this_hdu_list = fits.open(this_fname)
        hdu_image = this_hdu_list['IMAGE']
        hdu_zodi = this_hdu_list['ZODI']        
        image_header = hdu_image.header
        this_wcs = WCS(image_header)
        
        lam, bw = make_wavelength_image(
            hdu_list = this_hdu_list,
            use_hdu = 'IMAGE',
        )

        y_pix, x_pix = this_wcs.world_to_array_index(target_coord)
        this_shape = hdu_image.data.shape
        if (y_pix < 0) | (y_pix >= this_shape[0]) | \
           (x_pix < 0) | (x_pix >= this_shape[1]):
            continue
        
        this_val = hdu_image.data[y_pix, x_pix]
        this_zodi = hdu_zodi.data[y_pix, x_pix]        
        this_lam = lam.data[y_pix, x_pix]
        this_bw = bw.data[y_pix, x_pix]

        val_ra[counter] = this_val
        zodi_ra[counter] = this_zodi
        lam_ra[counter] = this_lam
        bw_ra[counter] = this_bw
        file_list.append(this_fname.split('/')[-1])
        
        counter += 1

    flist_ra = np.array(file_list)
        
    keep_ind = np.where(np.isfinite(val_ra))
    val_ra = val_ra[keep_ind]
    zodi_ra = zodi_ra[keep_ind]
    lam_ra = lam_ra[keep_ind]
    bw_ra = bw_ra[keep_ind]
    flist_ra = flist_ra[keep_ind]
        
    tab = QTable([lam_ra*u.um, bw_ra*u.um, val_ra*u.MJy/u.sr, zodi_ra*u.MJy/u.sr,
                  flist_ra],
                 names=['lam','bw','val','zodi','fname'])

    if outfile is not None:
        tab.write(outfile, overwrite=overwrite, format='ascii.ecsv')
    
    return(tab)

def make_mask_from_flags(
        flag_image,
        flags_to_use = ['SUR_ERROR','NONFUNC','MISSING_DATA',
                        'HOT','COLD','NONLINEAR','PERSIST'],        
):
    """
    """

    # From the explanatory supplement

    # Suggested flags for background estimation masking:

    # OVERFLOW
    # SUR_ERROR
    # NONFUNC
    # MISSING_DATA
    # HOT
    # COLD
    # NONLINEAR
    # PERSIST
    # OUTLIER
    # SOURCE
    # TRANSIENT

    # Suggested flags for on source photometry:

    # SUR_ERROR
    # NONFUNC
    # MISSING_DATA
    # HOT
    # COLD
    # NONLINEAR
    # PERSIST

    use_flag_ind = []
    flag_dict = spherex_flag_dict()
    for this_flag in flags_to_use:
        try:
            this_ind = flag_dict[this_flag.upper()]
        except KeyError:
            print("Flag unrecognized: ", this_flag)
            continue
        use_flag_ind.append(this_ind)
    
    n_flags = 22

    # Make an array of powers of 2 to cover each relevant bit
    powers_of_2 = 1 << np.arange(n_flags)

    # Copy the flag image to AND against each flag
    flag_cube = np.repeat(flag_image[:,:,np.newaxis], n_flags, axis=2)
    
    # Use bitwise AND to hash against each flag
    mask_cube = (flag_cube & powers_of_2) != 0
    
    # Initialize the image mask
    mask = np.zeros_like(flag_image, dtype=bool)

    # Loop over and accumulate the requested flags
    for ii in np.arange(n_flags):

        if ii not in use_flag_ind:
            continue

        mask = mask | (mask_cube[:,:,ii])
        
    return(mask)

def grid_spherex_cube(
        target_hdu = None,
        image_list = [],
        sub_zodi = True,
        outfile = None,
        overwrite=True):
    """TBD #1: handle the wavelengths better. They're offset right now and
    don't pay attention to the bandwidth.

    TBD #2: Make a not-a-cube spectrum at each location holding the
    SED and wavelength.

    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Initialize the output
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    target_header = target_hdu.header
    nx = target_header['NAXIS1']
    ny = target_header['NAXIS2']
    nz = target_header['NAXIS3']

    target_header_2d = target_header.copy()
    target_header_2d['NAXIS'] = 2
    del target_header_2d['NAXIS3']
    del target_header_2d['CRVAL3']
    del target_header_2d['CDELT3']
    del target_header_2d['CRPIX3']
    del target_header_2d['CTYPE3']
    del target_header_2d['CUNIT3']
    
    lam_step = target_header['CDELT3']    
    lam_array = (np.arange(nz)-(target_header['CRPIX3']-1))* \
        lam_step + target_header['CRVAL3']
    
    sum_cube = np.zeros((nz,ny,nx),dtype=np.float32)
    weight_cube = np.zeros((nz,ny,nx),dtype=np.float32)
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Loop over the image list
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for this_fname in ProgressBar(image_list):        
        
        this_hdu_list = fits.open(this_fname)
        hdu_image = this_hdu_list['IMAGE']
        hdu_flags = this_hdu_list['FLAGS']
        hdu_zodi = this_hdu_list['ZODI']        
        image_header = hdu_image.header
        
        lam, bw = make_wavelength_image(
            hdu_list = this_hdu_list,
            use_hdu = 'IMAGE',
        )
        
        hdu_lam = fits.PrimaryHDU(lam, image_header)
        hdu_bw = fits.PrimaryHDU(bw, image_header)

        # Implement flags and subtract ZODI if requested

        mask = make_mask_from_flags(
            hdu_flags.data,
            flags_to_use = ['SUR_ERROR','NONFUNC','MISSING_DATA',
                        'HOT','COLD','NONLINEAR','PERSIST']
        )

        if sub_zodi:
            masked_data = hdu_image.data - hdu_zodi.data
        else:
            masked_data = hdu_image.data
            
        masked_data[mask] = np.nan
        hdu_masked_image = fits.PrimaryHDU(masked_data, image_header)
        
        # This is pretty annoyingly inefficient to repeat
        missing = np.nan
        
        reprojected_image, footprint_image = \
            reproject_interp(hdu_masked_image, target_header_2d, order='bilinear')
        reprojected_image[footprint_image == 0] = missing

        reprojected_lam, footprint_lam = \
            reproject_interp(hdu_lam, target_header_2d, order='bilinear')
        reprojected_lam[footprint_lam == 0] = missing

        reprojected_bw, footprint_bw = \
            reproject_interp(hdu_bw, target_header_2d, order='bilinear')
        reprojected_bw[footprint_bw == 0] = missing        

        min_lam = np.nanmin(reprojected_lam - reprojected_bw - lam_step)
        max_lam = np.nanmax(reprojected_lam + reprojected_bw + lam_step)

        overlap_pix = np.sum(footprint_image)
        pix_in_cube = 0.0
        
        for zz, this_lam in enumerate(lam_array):

            # Skip this channel if there's no overlap with the current image
            
            if this_lam < min_lam:
                continue

            if this_lam > max_lam:
                continue            

            # Compare each pixel to the center of the current model
            
            delta = np.abs(this_lam - reprojected_lam)
            weight = delta / lam_step
            
            y_ind, x_ind = \
                np.where(delta <= (lam_step))

            pix_in_cube += len(y_ind)
            
            z_ind = np.zeros_like(y_ind,dtype=int)+zz
            
            sum_cube[z_ind, y_ind, x_ind] = \
                sum_cube[z_ind, y_ind, x_ind] + \
                (reprojected_image[y_ind, x_ind]*weight[y_ind, x_ind])

            weight_cube[z_ind, y_ind, x_ind] = \
                weight_cube[z_ind, y_ind, x_ind] + weight[y_ind, x_ind]
            
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Output and return
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    cube = sum_cube / weight_cube
    cube[np.where(weight_cube == 0.0)] = np.nan
    
    new_hdu = fits.PrimaryHDU(cube, target_header)
    if outfile is not None:
        new_hdu.writeto(outfile, overwrite=overwrite)
    
    return(new_hdu)
