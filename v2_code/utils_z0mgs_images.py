# Utilities related to z0MGS atlas construction. These routines focus
# on image manipulation and the actual process of atlas construction.

# IMAGE BUILDING
# - routines to make new headers
# - routine to identify which memebers of an index matter to a new image
# - mosaicking routines

# STARS
# - routines to query Gaia and create stacks of images of stars
# - routines to predict stellar fluxes from Gaia
# - routines to make an image of stellar fluxes for various bands

# MASKS
# - make a mask of galaxies based on a file

# BACKGROUNDS
# - routines to fit linear and planar backgrounds

# Basics
import os
import numpy as np
import math

# Astropy stuff
from astropy.convolution import convolve_fft
from astropy.coordinates import SkyCoord
import astropy.io.fits as fits
from astropy.table import Table, vstack
import astropy.units as u
from astropy.utils.console import ProgressBar
import astropy.wcs as wcs
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.stats import mad_std

# This can be deprecated once we pull Gaia down
from astroquery.gaia import Gaia

# Convolution and reprojection
from reproject import reproject_interp, reproject_exact

from spectral_cube import SpectralCube, LazyMask, Projection
from radio_beam import Beam

from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import curve_fit

# Visualization
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.visualization import simple_norm, LogStretch, PercentileInterval

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Define a new image header
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def make_simple_header(center_coord, pix_scale,
                       extent_x = None, extent_y = None,
                       nx = None, ny = None,
                       return_header=False):
    """Make a simple 2D FITS header centered on the coordinate of interest
    with a user-specififed pixel scale and extent.

    
    Parameters
    ----------

    center_coord : `~astropy.coordinates.SkyCoord` object or
        array-like Sky coordinates of the image center. If array-like
        then (ra, dec) in decimal degrees assumes.

    pix_scale : required. Size in decimal degrees of a pixel. Can be
        an array in which case it is pixel scale along x and y (e.g.,
        as returned by proj_pixel_scales).

    extent_x : the angular extent of the image along the x coordinate

    extent_y : the angular extent of the image along the y coordinate

    nx : the number of x pixels (not needed with extent_x and pix_scale)

    ny : the number of y pixels (not needed with extent_y and pix_scale)

    """
    
    # Figure out the center, working with types and units
    if isinstance(center_coord, SkyCoord):
        ra_ctr = center_coord.ra.degree
        dec_ctr = center_coord.dec.degree
    else:
        ra_ctr, dec_ctr = center_coord
        if hasattr(ra_ctr, 'unit'):
            ra_ctr = ra_ctr.to(u.deg).value
            dec_ctr = dec_ctr.to(u.deg).value

    if pix_scale is None:
        print("Pixel scale not specified. Returning.")
        return()
            
    # Figure out extent
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

    hdu = fits.PrimaryHDU()
    
    hdu.header = fits.Header()
    hdu.header['NAXIS'] = 2
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

    if return_header:
        return(hdu.header)
    else:
        return(hdu)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Manage index-image overlaps
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def find_index_overlap(
        index_file = None,
        index_tab = None,
        index_coords = None,
        index_extent = None,
        center_coord = None,
        image_extent = None,
        force_tolerance = None,
        selection_dict = {},
        return_separations = False,
):
    """Given a tabular index of tiles and a coord + extent for a new
    image, find all tiles in the index that could contribute to the
    new tile.

    Parameters
    ----------

    index_file : 

    index_tab : 

    index_coords: 

    skip_galaxy : List of galaxies to omit. Default empty.
    
    center_coord : 

    image_extent : 

    force_tolerance :

    selection dict : 

    return_separations : 

    """
    
    if index_tab is None:
        index_tab = (Table.read(index_file, format='fits'))

    if index_coords is None:
        index_coords = SkyCoord(
            ra=np.array(index_tab['ctr_ra'])*u.deg,
            dec=np.array(index_tab['ctr_dec'])*u.deg,
            frame='icrs')

    # Center to corner extent of the image
    if index_extent is None:
        index_extent = np.sqrt(
            (0.5*index_tab['nx']*index_tab['pix_scale_x'])**2 + 
            (0.5*index_tab['ny']*index_tab['pix_scale_y'])**2)
        
    if not isinstance(center_coord, SkyCoord):
        ra_ctr, dec_ctr = center_coord
        if hasattr(ra_ctr, 'unit'):
            ra_ctr = ra_ctr.to(u.deg).value
            dec_ctr = dec_ctr.to(u.deg).value
        center_coord = SkyCoord(ra=ra_ctr*u.deg, dec=dec_ctr*u.deg, frame='icrs')

    # Find tiles within half of the image extent + tile extent (i.e.,
    # that overlap) or a user-set tolerance if supplied
    
    separations = np.array(index_coords.separation(center_coord))
    if force_tolerance is None:
        tolerance = np.array(index_extent + 0.5*image_extent)
    else:
        tolerance = force_tolerance
    tiles_overlap = separations < tolerance
        
    for this_key, this_value in selection_dict.items():
        tiles_overlap *= (index_tab[this_key] == this_value)

    overlap_tab = index_tab[tiles_overlap]

    if return_separations:
        return(overlap_tab, separations)
    else:
        return(overlap_tab)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Query GAIA for use identifying foreground stars
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# This builds an SQL query to Gaia that can be used to find all stars
# in the target footprint.

# TBD improvements: the square is a bit inelegant, could probably
# switch to a cone query and downselect to the square after the fact
# or just use a cone.

def query_gaia(
        ra_min=None,
        ra_max=None,
        dec_min=None,
        dec_max=None,
        outfile=None,
        add_fields=[],
        omit_fields=[],
        quiet=False,
        dry_run=True,
        skip_if_present=False):
    """Query GAIA DR3 with parameters of interest to extragalactic images
    in a rectangular footprint. This executes a single query using SQL
    and astroquery and writes the result as a FITS table to outfile.
    """

    if outfile is not None:
        if (os.path.isfile(outfile) == True) and \
           (skip_if_present == True):
            return
            
    tab = None
    delim=' '

    # Construct the base SQL query
    source = 'FROM'+delim+'gaiadr3.gaia_source'+delim+'AS'+delim+'src'+delim+'JOIN'+delim+'gaiadr3.astrophysical_parameters'+delim+'AS'+delim+'ap'+delim+'ON'+delim+'src.source_id'+delim+'='+delim+'ap.source_id'+delim+''

    # Set the list of fields to query
    fields = [
        'src.source_id',
        'src.ra',
        'src.dec',
        'ap.classprob_dsc_combmod_quasar',
        'ap.classprob_dsc_combmod_star',
        'ap.classprob_dsc_combmod_galaxy',
        'ap.classprob_dsc_combmod_binarystar',
        'ap.classprob_dsc_combmod_whitedwarf',
        'ap.teff_gspphot',
        'ap.teff_gspspec',
        'src.phot_g_mean_flux',
        'src.phot_g_mean_flux_error',
        'src.phot_g_mean_mag',
        'src.phot_bp_mean_mag',
        'src.phot_rp_mean_mag',                
        'src.parallax',
        'src.parallax_error',
        'src.pmra',
        'src.pmra_error',
        'src.pmdec',
        'src.pmdec_error',
        'src.ruwe',
    ]

    for field_to_add in add_fields:
        if fields.count(field_to_add) == 0:
            fields.append(field_to_add)
    for field_to_omit in omit_fields:
        if fields.count(field_to_omit) == 1:
            fields.remove(field_to_omit)

    # Add the fields to the SQL query
    all_fields = 'SELECT'+delim+''
    first=True
    for this_field in fields:
        if first == False:
            all_fields += ','
        all_fields += this_field.strip()
        first = False
    
    # Deal with longitude wrap (probably could do this better with
    # SkyCoords - come back and fix it)
    if ra_min > ra_max:
        twopart_call = True
        this_ramin = ra_min
        this_ramax = 360.0
    else:
        twopart_call = False
        this_ramin = ra_min
        this_ramax = ra_max

    # Convert RA and Dec to strings
    ra_low_string = (f'{this_ramin:.5f}').strip()
    ra_high_string = (f'{this_ramax:.5f}').strip()

    dec_low_string = (f'{dec_min:.5f}').strip()
    dec_high_string = (f'{dec_max:.5f}').strip()

    # Now convert to a selection to add to the SQL query
    selection = 'WHERE'+delim+'src.ra'+delim+'BETWEEN'+ \
        delim+''+ra_low_string+''+delim+'AND'+delim+''+ra_high_string+ \
        ''+delim+'AND'+delim+'src.dec'+delim+'BETWEEN'+ \
        delim+''+dec_low_string+''+delim+'AND'+delim+''+dec_high_string
        
    query = all_fields + ' ' + source + ' ' + selection

    # Execute the query
    if dry_run == False:
        if quiet == False:
            print(query)
        if outfile is not None:
            if (os.path.isfile(outfile) == False) or \
               (skip_if_present == False):            
                job = Gaia.launch_job_async(query)
                tab = job.get_results()
                print(tab)
                if twopart_call == False:            
                    tab.write(outfile, format='fits', overwrite=True)
    else:
        if quiet == False:
            print(query)

    if twopart_call == False:
        return(query)

    # From here only worry about the case where we wrapped around the
    # meridian and need to make a second call and then stitch two
    # separate tables back together.

    # ... build the second call
    
    first_part = query
    
    this_ramin = 0.0
    this_ramax = ra_max

    # Convert to strings
    ra_low_string = (f'{this_ramin:.5f}').strip()
    ra_high_string = (f'{this_ramax:.5f}').strip()

    selection = 'WHERE'+delim+'src.ra'+delim+'BETWEEN'+delim+''+ra_low_string+''+delim+'AND'+delim+''+ra_high_string+ \
        ''+delim+'AND'+delim+'src.dec'+delim+'BETWEEN'+delim+''+dec_low_string+''+delim+'AND'+delim+''+dec_high_string

    query = all_fields + ' ' + source + ' ' + selection

    # ... execute it
    
    if dry_run == False:
        if quiet==False:
            print(query)
        if outfile is not None:            
            if (os.path.isfile(outfile) == False) or \
               (skip_if_present == False):                        
                job = Gaia.launch_job_async(query)
                tab2 = job.get_results()

                # ... combine the two tables
                tab = vstack([tab,tab2])
                tab.write(outfile, format='fits', overwrite=True)
    else:
        if quiet==False:
            print(query)

    return([first_part,query]) 

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Make cutouts around GAIA sources
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Make small cutouts around a set of Gaia sources. This is in
# principle useful to test the PSF in an image or to model how to
# convert stellar fluxes to the atlas bands.

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
    stack_hdu = make_simple_header(
        center_coord=(0.0*u.deg, 0.0*u.deg),
        pix_scale=target_pix_scale,
        nx=half_width_pix*2+1, ny=half_width_pix*2+1)
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

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Related to predicting a stellar flux from Gaia or 2MASS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Predict a stellar flux at WISE or GALEX from Gaia or 2MASS Ks
# catalog magnitudes. For a table and a header make these predictions
# and construct a Jy/pixel image.

def pred_star_flux_from_gaia(
        g_mag = None,
        bp = None,
        rp = None,
        parallax = None,
):
    """Convert from Gaia G magnitude to a predicted flux in the GALEX and
    WISE bands. BP, RP, and parallax not currently used but could be.

    Parameters
    ----------

    g_mag : 

    bp : 

    rp : 

    parallax :
    
    """
    
    # Hard-code the conversions from z0mgs paper but in flux units
    # (i.e., predict Jy)
    gauss7p5_sr = (7.5/3600.*np.pi/180./2.0)**2*np.pi/np.log(2)
    gauss7p5_intens_to_flux = gauss7p5_sr*1e6

    gauss15_sr = (15/3600.*np.pi/180./2.0)**2*np.pi/np.log(2)
    gauss15_intens_to_flux = gauss15_sr*1e6

    # coefficient in flux = 10.**(-1.0*mag/2.5)*coefficient
    mag_to_jy = {
        'fuv': 27.07*gauss7p5_intens_to_flux,
        'nuv': 13060.2*gauss7p5_intens_to_flux,
        'w1': 829914.*gauss7p5_intens_to_flux,
        'w2': 450745.*gauss7p5_intens_to_flux,
        'w3': 77349.7*gauss7p5_intens_to_flux,
        'w4': 8398.*gauss15_intens_to_flux,
    }

    flux_dict = {}
    for this_band in mag_to_jy.keys():
        flux_dict[this_band]= mag_to_jy[this_band]* \
            10.**(-1.0*g_mag/2.5)

    return(flux_dict)
    
def pred_star_flux_from_ks(
        ks_mag = None
):
    """Convert from 2MASS Ks magnitude to a predicted flux in the GALEX and
    WISE bands.    
    """
    
    # Hard-code the conversions from z0mgs paper but in flux units
    # (i.e., predict Jy)    
    gauss7p5_sr = (7.5/3600.*np.pi/180./2.0)**2*np.pi/np.log(2)
    gauss7p5_intens_to_flux = gauss7p5_sr*1e6

    gauss15_sr = (15/3600.*np.pi/180./2.0)**2*np.pi/np.log(2)
    gauss15_intens_to_flux = gauss15_sr*1e6

    # coefficient in flux = 10.^(-1.0d*mag/2.5)*coefficient
    mag_to_jy = {
        'fuv': 1.35*gauss7p5_intens_to_flux,
        'nuv': 202.0*gauss7p5_intens_to_flux,
        'w1': 164900.*gauss7p5_intens_to_flux,
        'w2': 91800.*gauss7p5_intens_to_flux,
        'w3': 16800.*gauss7p5_intens_to_flux,
        'w4': 1550.*gauss15_intens_to_flux,
    }

    flux_dict = {}
    for this_band in mag_to_jy.keys():
        flux_dict[this_band]= mag_to_jy[this_band]* \
            10.**(-1.0*ks_mag/2.5)
        
    return(flux_dict)

def build_star_flux_image(
        template_file = None,
        template_hdu = None,
        outfile = None,
        band = 'w1',
        gaia_file = None,
        gaia_tab = None,
        ks_file = '../../measurements/tab_2mass_stars.fits',
        ks_tab = None,
        gaia_s2n_cut = 3.5,
        center_coord = None,
        center_tol = 3.0*u.arcsec,
        overwrite=True,
):
    """Build a flux image in Jy/pixel containing predicted fluxes of stars
    known from Gaia and 2MASS in the specified band.

    # TBD Interp is untested.

    # TBD The models could be improved!

    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Make an image matched to the template
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if template_hdu is None:

        template_hdu = fits.open(template_file)[0]

    # Initialize output. Copy the header, make a map, and get the WCS
    # for the new output map.
    
    out_hdr = template_hdu.header

    # Remove beam information and update units
    if 'BMAJ' in out_hdr:
        out_hdr.remove('BMAJ')
    if 'BMIN' in out_hdr:
        out_hdr.remove('BMIN')
    if 'BPA' in out_hdr:
        out_hdr.remove('BPA')
    out_hdr['BUNIT'] = 'Jy/pixel'
    out_hdr['BAND'] = band

    # Initialize an empty map and WCS object
    out_map = np.zeros_like(template_hdu.data)
    out_wcs = wcs.WCS(out_hdr)
    
    ny, nx = out_map.shape
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # If supplied add Gaia-based stars to the image
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    # Read the Gaia file

    if gaia_tab is None and gaia_file is not None:

        gaia_tab = Table.read(gaia_file, format='fits')

    if gaia_tab is not None:

        # Identify the Gaia sources within the image
        gaia_ra = gaia_tab['ra']
        gaia_dec =  gaia_tab['dec']
        gaia_coords = SkyCoord(ra=gaia_ra, dec=gaia_dec, frame='icrs')
        
        gaia_pix = out_wcs.world_to_pixel(gaia_coords)
        gaia_x = gaia_pix[0]
        gaia_y = gaia_pix[1]

        in_image = (gaia_x >= 0) & (gaia_x <= (nx-1)) \
            & (gaia_y >= 0) & (gaia_y <= (ny-1))

        gaia_tab = gaia_tab[in_image]
        gaia_coords = gaia_coords[in_image]
        gaia_x = gaia_x[in_image]
        gaia_y = gaia_y[in_image]
        
        # Select the Gaia stars likely to be foreground stars
        foreground = \
            (gaia_tab['parallax'] >= gaia_s2n_cut*gaia_tab['parallax_error']) & \
            (gaia_tab['pmra'] >= gaia_s2n_cut*gaia_tab['pmra_error']) & \
            (gaia_tab['pmdec'] >= gaia_s2n_cut*gaia_tab['pmdec_error'])

        gaia_tab = gaia_tab[foreground]
        gaia_coords = gaia_coords[foreground]
        gaia_x = gaia_x[foreground]
        gaia_y = gaia_y[foreground]        

        # Avoid the galaxy center if requested        
        if center_coord is not None:
        
            away_from_center = \
                (center_coord.separation(gaia_coords) > center_tol)

            gaia_tab = gaia_tab[away_from_center]
            gaia_coords = gaia_coords[away_from_center]
            gaia_x = gaia_x[away_from_center]
            gaia_y = gaia_y[away_from_center]
                    
        # Predict the flux for each 
        gaia_fluxes = pred_star_flux_from_gaia(
            g_mag = gaia_tab['phot_g_mean_mag'],
            bp = gaia_tab['phot_bp_mean_mag'],
            rp = gaia_tab['phot_rp_mean_mag'],
            parallax = gaia_tab['parallax']
        )
        
        # Add the flux of each Gaia star to the image. Must be a
        # better way to do this.
        for ii, this_flux in enumerate(gaia_fluxes[band]):
            # TBD - interpolate linearly between four pix
            method = 'nearest'
            if method == 'nearest':
                this_x = int(round(gaia_x[ii]))
                this_y = int(round(gaia_y[ii]))
                out_map[this_y, this_x] += this_flux
            if method == 'interp':
                low_x = int(math.floor(gaia_x[ii]))
                high_x = int(math.ceil(gaia_x[ii]))
                frac_x = (gaia_x[ii] - low_x)/(1.0*high_x - 1.0*low_x)
                low_y = int(math.floor(gaia_y[ii]))
                high_y = ceil(math.ceil(gaia_y[ii]))                
                frac_y = (gaia_y[ii] - low_y)/(1.0*high_y - 1.0*low_y)
                out_map[low_y, low_x] = (1.0-frac_x)*(1.0-frac_y)*this_flux
                out_map[low_y, high_x] = frac_x*(1.0-frac_y)*this_flux
                out_map[high_y, low_x] = (1.0-frac_x)*frac_y*this_flux
                out_map[high_y, high_x] = frac_x*frac_y*this_flux               
                
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # If supplied add the 2MASS based stars to the image
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # Read the 2MASS file or table
    if ks_tab is None and ks_file is not None:

        ks_tab = Table.read(ks_file, format='fits')
        
        # Identify the Gaia sources within the image
        ks_ra = ks_tab['RA'] * u.deg
        ks_dec =  ks_tab['DEC'] * u.deg
        ks_coords = SkyCoord(ra=ks_ra, dec=ks_dec, frame='icrs')
        
        ks_pix = out_wcs.world_to_pixel(ks_coords)
        ks_x = ks_pix[0]
        ks_y = ks_pix[1]

        in_image = (ks_x >= 0) & (ks_x < nx) \
            & (ks_y >= 0) & (ks_y < ny)

        ks_tab = ks_tab[in_image]
        ks_coords = ks_coords[in_image]
        ks_x = ks_x[in_image]
        ks_y = ks_y[in_image]
        
        # Avoid the galaxy center if requested
        if center_coord is not None:
        
            away_from_center = \
                (center_coord.separation(ks_coords) > center_tol)

            ks_tab = ks_tab[away_from_center]
            ks_coords = ks_coords[away_from_center]
            ks_x = ks_x[away_from_center]
            ks_y = ks_y[away_from_center]
                    
        # Predict the flux for each 
        ks_fluxes = pred_star_flux_from_ks(
            ks_mag = ks_tab['KS_MAG'],
        )
        
        # Add the flux of each 2MASS star to the image. Must be a
        # better way to do this.
        for ii, this_flux in enumerate(ks_fluxes[band]):
            method = 'nearest'
            if method == 'nearest':            
                this_x = int(round(ks_x[ii]))
                this_y = int(round(ks_y[ii]))
                out_map[this_y, this_x] += this_flux        
            if method == 'interp':
                low_x = int(math.floor(ks_x[ii]))
                high_x = int(math.ceil(ks_x[ii]))
                frac_x = (ks_x[ii] - low_x)/(1.0*high_x - 1.0*low_x)
                low_y = int(math.floor(ks_y[ii]))
                high_y = ceil(math.ceil(ks_y[ii]))                
                frac_y = (ks_y[ii] - low_y)/(1.0*high_y - 1.0*low_y)
                out_map[low_y, low_x] = (1.0-frac_x)*(1.0-frac_y)*this_flux
                out_map[low_y, high_x] = frac_x*(1.0-frac_y)*this_flux
                out_map[high_y, low_x] = (1.0-frac_x)*frac_y*this_flux
                out_map[high_y, high_x] = frac_x*frac_y*this_flux               

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Return and write to disk if requested
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    out_hdu = fits.PrimaryHDU(out_map, out_hdr)
    
    # Write to disk
    if outfile is not None:
        out_hdu.writeto(outfile, overwrite=overwrite)

    return(out_hdu)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Build masks
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def build_galaxy_mask(
        template_file = None,
        template_hdu = None,
        outfile = None,
        gal_tab = None,
        gal_tab_file = '/home/leroy.42/idl/galbase/gal_data/gal_base.fits',
        this_pgc = None,        
        pgc_to_skip = None,
        field_for_rad = 'R25_DEG',
        min_rad = 7.5/3600.*u.deg,
        rad_fac_to_blank = 1.0,
        use_orient = True,
        max_incl = 70.*u.deg,
        show = False,
        pause = False,
        overwrite = True,
):
    """Use an external galaxy database table to build a mask of galaxies
    (other than the target) in a provided field of view.

    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Read the table of galaxies
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    if gal_tab is None:

        gal_tab = Table.read(gal_tab_file, format='fits')

    gal_ra = gal_tab['RA_DEG'] * u.deg
    gal_dec = gal_tab['DEC_DEG'] * u.deg
    gal_coords = SkyCoord(ra=gal_ra, dec=gal_dec, frame='icrs')
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Make an image matched to the template
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if template_hdu is None:

        template_hdu = fits.open(template_file)[0]

    # Initialize output. Copy the header, make a map, and get the WCS
    # for the new output map.
    
    out_hdr = template_hdu.header

    # Remove beam information and update units
    if 'BMAJ' in out_hdr:
        out_hdr.remove('BMAJ')
    if 'BMIN' in out_hdr:
        out_hdr.remove('BMIN')
    if 'BPA' in out_hdr:
        out_hdr.remove('BPA')
    out_hdr['BUNIT'] = 'Mask'

    # Initialize an empty map and WCS object
    out_mask = np.zeros_like(template_hdu.data, dtype=int)
    out_wcs = wcs.WCS(out_hdr)
    
    ny, nx = out_mask.shape

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Identify galaxies that may overlap the image
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # Calculate coordinates of galaxies in the table
    gal_pix = out_wcs.world_to_pixel(gal_coords)
    gal_x = gal_pix[0]
    gal_y = gal_pix[1]

    # TBD Note the padding - build this out later. Depends on the size
    # of the galaxy in pixels. So calculate pixel size and divide
    # galaxy size by that.
    pix_scale = wcs.utils.proj_plane_pixel_scales(out_wcs)[0]
    pad_pix = np.array(gal_tab[field_for_rad] / pix_scale)
    print(pad_pix)
    
    # Galaxies of interest
    in_image = (gal_x >= (0-pad_pix)) & \
        (gal_x < (nx+pad_pix)) & \
        (gal_y >= (0-pad_pix)) & \
        (gal_y < (ny+pad_pix))
    
    # Downsample to galaxies of interest
    this_gal_tab = gal_tab[in_image]

    print("Building a galaxy mask")
    print("... found galaxies in image: ", len(this_gal_tab))
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Add those galaxies to the mask
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for this_gal in this_gal_tab:

        # Skip the target galaxy and any user flagged galaxies
        
        if this_pgc is not None:
            if this_gal['PGC'] == this_pgc:
                continue

        if pgc_to_skip is not None:
            if this_gal['PGC'] in pgc_to_skip:
                continue

        # Find the footprint of the galaxy in the image
            
        this_ra = this_gal['RA_DEG'] * u.deg
        this_dec = this_gal['DEC_DEG'] * u.deg
        this_coords = SkyCoord(ra=this_ra, dec=this_dec, frame='icrs')
        
        if use_orient:
            this_incl = this_gal['INCL_DEG']*u.deg
            this_pa = this_gal['POSANG_DEG']*u.deg
        else:
            this_incl = 0.*u.deg
            this_pa = 0.*u.deg
                        
        if np.isfinite(this_incl) == False | np.isfinite(this_pa) == False:
            this_incl = 0.*u.deg
            this_pa = 0.*u.deg

        if max_incl is not None:
            if this_incl > max_incl:
                max_incl = max_incl
                
        radius_deg, projang_deg = \
            deproject(
                center_coord=this_coords, incl=this_incl, pa=this_pa,
                template_header = out_hdr, return_offset = False)
                
        # Create the mask for this galaxy

        this_rad_deg = (this_gal[field_for_rad]*u.deg).to_value(u.deg)

        this_mask = radius_deg <= (this_rad_deg*rad_fac_to_blank)
                
        # Add to the overall mask

        out_mask[this_mask] = 1
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Return and write to disk if requested
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    out_hdu = fits.PrimaryHDU(out_mask, out_hdr)
    
    # Write to disk
    if outfile is not None:
        out_hdu.writeto(outfile, overwrite=overwrite)

    return(out_hdu)

def build_star_mask(
        image_file = None,
        image_hdu = None,
        star_file = None,
        star_hdu = None,
        outfile = None,
        clip_level = None,
        rms_fac = 3.0,
        rms_value = None,
        show = False,
        pause = False,
        overwrite = True,
):
    """Accept an input prediction file and create a mask based on a
    threshold.

    TBD could be simplified, removing the noise image.
    """
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Read the image and star prediction
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if image_hdu is None:        
        image_hdu = fits.open(image_file)[0]

    if star_hdu is None:        
        star_hdu = fits.open(star_file)[0]
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Determine the clip level
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if clip_level is None:

        if rms_value is None:

            rms_value = mad_std(image_hdu.data)

        clip_level = rms_value * rms_fac

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Make the mask
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    star_mask = (star_hdu.data >= clip_level)*1.
    star_mask_hdr = star_hdu.header

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Return and write
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    star_mask_hdu = fits.PrimaryHDU(data=star_mask, header=star_mask_hdr)
    
    # Write to disk
    if outfile is not None:
        star_mask_hdu.writeto(outfile, overwrite=overwrite)

    return(star_mask_hdu)

def stack_masks(
        mask_fnames = None,
        mask_hdus = None,        
        outfile = None,
        overwrite = True,
        # TBD - enable nearest neighbor reprojection
):
    """
    Join together maps that are already.

    """

    # Loop over file names

    have_first_mask = False
    
    for this_fname in mask_fnames:
        
        this_hdu = fits.open(this_fname)[0]
        this_mask = this_hdu.data
        if not have_first_mask:
            mask = this_hdu.data
            mask_hdr = this_hdu.header
            have_first_mask = True
        else:
            mask = np.logical_or(mask, this_mask)

    # Loop over HDUs
        
    for this_other_hdu in mask_hdus:

        this_mask = this_other_hdu.data
        if not have_first_mask:
            mask = this_hdu.data
            mask_hdr = this_hdu.header
            have_first_mask = True
        else:
            mask = np.logical_or(mask, this_mask)            

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Outfile
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    mask_hdu = fits.PrimaryHDU(mask.astype(int), header = mask_hdr)
    
    # Write to disk
    if outfile is not None:
        mask_hdu.writeto(outfile, overwrite=overwrite)

    return(mask_hdu)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Related to deprojection
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def deproject(
        center_coord=None,
        incl=0*u.deg,
        pa=0*u.deg,
        template_header=None,
        template_wcs=None,
        template_naxis=None,
        template_ra=None,
        template_dec=None,
        return_offset=False,
        verbose=False):
    """Calculate deprojected radii and projected angles in a disk.

    This function deals with projected images of astronomical objects
    with an intrinsic disk geometry. Given sky coordinates of the disk
    center, disk inclination and position angle, this function
    calculates deprojected radii and projected angles based on

    (1) a FITS header (`header`), or

    (2) a WCS object with specified axis sizes (`wcs` + `naxis`), or
    
    (3) RA and DEC coodinates (`ra` + `dec`).
    
    Both deprojected radii and projected angles are defined relative
    to the center in the inclined disk frame. For (1) and (2), the
    outputs are 2D images; for (3), the outputs are arrays with shapes
    matching the broadcasted shape of `ra` and `dec`.

    Parameters
    ----------
    center_coord : `~astropy.coordinates.SkyCoord` object or array-like
        Sky coordinates of the disk center
    incl : `~astropy.units.Quantity` object or number, optional
        Inclination angle of the disk (0 degree means face-on)
        Default is 0 degree.
    pa : `~astropy.units.Quantity` object or number, optional
        Position angle of the disk (red/receding side, North->East)
        Default is 0 degree.
    header : `~astropy.io.fits.Header` object, optional
        FITS header specifying the WCS and size of the output 2D maps
    wcs : `~astropy.wcs.WCS` object, optional
        WCS of the output 2D maps
    naxis : array-like (with two elements), optional
        Size of the output 2D maps
    ra : array-like, optional
        RA coordinate of the sky locations of interest
    dec : array-like, optional
        DEC coordinate of the sky locations of interest
    return_offset : bool, optional
        Whether to return the angular offset coordinates together with
        deprojected radii and angles. Default is to not return.

    Returns
    -------
    deprojected coordinates : list of arrays
        If `return_offset` is set to True, the returned arrays include
        deprojected radii, projected angles, as well as angular offset
        coordinates along East-West and North-South direction;
        otherwise only the former two arrays will be returned.

    Notes
    -----
    This is the Python version of an IDL function `deproject` included in the `cpropstoo` package. See URL below:

    https://github.com/akleroy/cpropstoo/blob/master/cubes/deproject.pro

    Convention on the in-plane position angle w.r.t. the receding node may be flipped.

    Python routine from Jiayi Sun. Modified for compatibility with
    z0mgs namespace.

    """

    if isinstance(center_coord, SkyCoord):
        x0_deg = center_coord.ra.degree
        y0_deg = center_coord.dec.degree
    else:
        x0_deg, y0_deg = center_coord
        if hasattr(x0_deg, 'unit'):
            x0_deg = x0_deg.to(u.deg).value
            y0_deg = y0_deg.to(u.deg).value
    if hasattr(incl, 'unit'):
        incl_deg = incl.to(u.deg).value
    else:
        incl_deg = incl
    if hasattr(pa, 'unit'):
        pa_deg = pa.to(u.deg).value
    else:
        pa_deg = pa

    if template_header is not None:
        wcs_cel = wcs.WCS(template_header).celestial
        naxis1 = template_header['NAXIS1']
        naxis2 = template_header['NAXIS2']
        # create ra and dec grids
        ix = np.arange(naxis1)
        iy = np.arange(naxis2).reshape(-1, 1)
        ra_deg, dec_deg = wcs_cel.wcs_pix2world(ix, iy, 0)
    elif (template_wcs is not None) and \
         (template_naxis is not None):
        wcs_cel = template_wcs.celestial
        naxis1, naxis2 = template_naxis
        # create ra and dec grids
        ix = np.arange(naxis1)
        iy = np.arange(naxis2).reshape(-1, 1)
        ra_deg, dec_deg = wcs_cel.wcs_pix2world(ix, iy, 0)
    else:
        if template_ra.ndim == 1:
            ra_deg, dec_deg = \
                np.broadcast_arrays(template_ra, template_dec)
        else:
            ra_deg, dec_deg = template_ra, template_dec
            if verbose:
                print("ra ndim != 1")
        if hasattr(template_ra, 'unit'):
            ra_deg = template_ra.to(u.deg).value
            dec_deg = template_dec.to(u.deg).value
    
    
    #else:
        #ra_deg, dec_deg = np.broadcast_arrays(ra, dec)
        #if hasattr(ra_deg, 'unit'):
            #ra_deg = ra_deg.to(u.deg).value
            #dec_deg = dec_deg.to(u.deg).value

    # recast the ra and dec arrays in term of the center coordinates
    # arrays are now in degrees from the center
    dx_deg = (ra_deg - x0_deg) * np.cos(np.deg2rad(y0_deg))
    dy_deg = dec_deg - y0_deg

    # rotation angle (rotate x-axis up to the major axis)
    rotangle = np.pi/2 - np.deg2rad(pa_deg)

    # create deprojected coordinate grids
    deprojdx_deg = (dx_deg * np.cos(rotangle) +
                    dy_deg * np.sin(rotangle))
    deprojdy_deg = (dy_deg * np.cos(rotangle) -
                    dx_deg * np.sin(rotangle))
    deprojdy_deg /= np.cos(np.deg2rad(incl_deg))

    # make map of deprojected distance from the center
    radius_deg = np.sqrt(deprojdx_deg**2 + deprojdy_deg**2)

    # make map of angle w.r.t. position angle
    projang_deg = np.rad2deg(np.arctan2(deprojdy_deg, deprojdx_deg))

    if return_offset:
        return radius_deg, projang_deg, deprojdx_deg , deprojdy_deg
    else:
        return radius_deg, projang_deg
    

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Related to convolution
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def get_pixscale(hdu):
    """From PJPIPE. Helper function used in convolve. Get pixel scale from
    header. Checks HDU header and returns a pixel scale

    Args:

        hdu: hdu to get pixel scale for

    """

    PIXEL_SCALE_NAMES = ["XPIXSIZE", "CDELT1", "CD1_1", "PIXELSCL"]

    for pixel_keyword in PIXEL_SCALE_NAMES:
        try:
            try:
                pix_scale = np.abs(float(hdu.header[pixel_keyword]))
            except ValueError:
                continue
            if pixel_keyword in ["CDELT1", "CD1_1"]:
                pix_scale = wcs.WCS(hdu.header).proj_plane_pixel_scales()[0].value * 3600
                # pix_scale *= 3600
            return pix_scale
        except KeyError:
            pass

    raise Warning("No pixel scale found")

def convolve_image_with_kernel(
        image_file=None,
        image_hdu=None,
        kernel_file=None,
        kernel_hdu=None,
        outfile=None,
        blank_zeros=False,
        force_jwst_syntax=False,
        overwrite=True,
):
    """Convolves input image with an input kernel, and writes to
    disk. Moderately edited from PJPIPE (Tom Williams) version to
    allow more flexible handling of extensions and remove the
    reprojection.

    Args:

    image_hdu: HDU to be convolved
    image_file: Path to image file (if HDU not supplied)

    kernel_hdu: HDU of kernel
    kernel_file: Path to kernel for convolution (if HDU not supplied)

    outfile: Path to output file (optional)

    blank_zeros: If True, then all zero values will be set to NaNs. Defaults to True
    force_jwst_syntax: Force use of SCI and ERR extensions, else try to be smart
    overwrite: overwrite the output file if relevant

    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Read the kernel
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    
    if kernel_hdu is None:
        kernel_hdu = fits.open(kernel_file)

    # Get the pixel scale
    kernel_pix_scale = get_pixscale(kernel_hdu[0])
    
    # Note the shape and grid of the kernel as input
    kernel_data = kernel_hdu[0].data
    kernel_hdu_length = kernel_hdu[0].data.shape[0]
    original_central_pixel = (kernel_hdu_length - 1) / 2

    # 1-d vector of angular offset values
    original_grid = \
        (np.arange(kernel_hdu_length) - original_central_pixel) * \
        kernel_pix_scale

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Read the image
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if image_hdu is None:
        image_hdu = fits.open(image_file)
    
    # Detect whether the JWST SCI/ERR extension syntax is
    # present. If not, then default to just convolving the Primary
    # HDU in the file.
        
    if force_jwst_syntax:
        sci_ext = 'SCI'
    else:
        hdu_dict = image_hdu.info(False)
        ext_list = []
        for this_hdu in hdu_dict:
            ext_list.append(this_hdu[1])
        if 'SCI' in ext_list:
            sci_ext = 'SCI'
        else:
            # Most common other case
            sci_ext = 'PRIMARY'
            
    if 'ERR' in ext_list:
        use_err = True
    else:
        use_err = False

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Handle blanks in the image
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    if blank_zeros:
        
        # make sure that all zero values were set to NaNs, which
        # astropy convolution handles with interpolation
        
        image_hdu[sci_ext].data[(image_hdu[sci_ext].data == 0)] = np.nan
        if use_err:
            image_hdu["ERR"].data[(image_hdu[sci_ext].data == 0)] = np.nan

    #print(np.nansum(image_hdu[sci_ext].data))
    #plt.imshow(image_hdu[sci_ext].data, origin='lower')
    #plt.imshow(kernel_data, origin='lower')
    #plt.show()
            
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Regrid the kernel for use with convolve
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    image_pix_scale = get_pixscale(image_hdu[sci_ext])

    # Calculate kernel size after interpolating to the image pixel
    # scale. Because sometimes there's a little pixel scale
    # rounding error, subtract a little bit off the optimum size
    # (from Tom Williams).

    interpolate_kernel_size = (
        np.floor(kernel_hdu_length * kernel_pix_scale / image_pix_scale) - 2
    )

    # Ensure the kernel has a central pixel
    
    if interpolate_kernel_size % 2 == 0:
        interpolate_kernel_size -= 1

    # Define a new coordinate grid onto which to project the kernel
    # but using the pixel scale of the image

    new_central_pixel = (interpolate_kernel_size - 1) / 2
    new_grid = (
        np.arange(interpolate_kernel_size) - new_central_pixel
    ) * image_pix_scale
    x_coords_new, y_coords_new = np.meshgrid(new_grid, new_grid)

    # Do the reprojection from the original kernel grid onto the new
    # grid with pixel scale matched to the image

    grid_interpolated = RegularGridInterpolator(
        (original_grid, original_grid),
        kernel_data,
        bounds_error=False,
        fill_value=0.0,
    )
    kernel_interp = grid_interpolated(
        (x_coords_new.flatten(), y_coords_new.flatten())
    )
    kernel_interp = kernel_interp.reshape(x_coords_new.shape)

    # Ensure the interpolated kernel is normalized to a sum of 1
    
    kernel_interp = kernel_interp / np.nansum(kernel_interp)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Convolve
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    # Now with the kernel centered and matched in pixel scale to the
    # input image use the FFT convolution routine from astropy to
    # convolve.

    conv_im = convolve_fft(
        image_hdu[sci_ext].data,
        kernel_interp,
        allow_huge=True,
        boundary='fill',
        fill_value=0.0,
        # This is debateable either way but I think for most z0mgs
        # cases this mimics previous work and is what we want.
        nan_treatment='fill',
        preserve_nan=True,
        psf_pad=True,
    )
    
    # If an error map is present, convolve errors (with kernel**2,
    # do not normalize it).  This, however, doesn't account for
    # covariance between pixels
    
    if use_err:
        conv_err = np.sqrt(
            convolve_fft(
                image_hdu["ERR"].data ** 2,
                kernel_interp ** 2,
                allow_huge=True,
                normalize_kernel=False,
                boundary='fill',
                fill_value=0.0,
                # TBD - The fill value here is no good but no good options
                nan_treatment='fill',
                preserve_nan=True,
                psf_pad=True,
            )
        )
        
    image_hdu[sci_ext].data = conv_im
    if use_err:
        image_hdu["ERR"].data = conv_err
        
    if outfile is not None:
        image_hdu.writeto(outfile, overwrite=overwrite)

    return(image_hdu)

def convolve_image_with_gauss(
        image_file=None,
        image_hdu=None,
        starting_res=None,
        target_res=None,
        nan_treatment='fill',
        dtype=np.float32,
        outfile=None,
        overwrite=True
):
    """
    Docs
    """

    
    print("Convolving with Gaussian")
            
    if image_hdu is None:
        print("... file: ", image_file)
        image_hdu = fits.open(image_file)

    proj = Projection.from_hdu(image_hdu)

    if starting_res is not None:
        old_proj = proj
        starting_beam = Beam(starting_res)
        proj = proj.with_beam(starting_beam)
        print("... starting from ", starting_beam)
    
    proj.allow_huge_operations = True

    target_beam = Beam(major=target_res)
    print("... target ", target_beam)
    
    convolved_proj = proj.convolve_to(
        target_beam,
        nan_treatment=nan_treatment, fill_value=0.0,
        allow_huge=True)

    out_hdu = fits.PrimaryHDU(
        np.array(convolved_proj.filled_data[:], dtype=dtype),
        header=convolved_proj.header)
        
    if outfile is not None:
        out_hdu.writeto(outfile, overwrite=overwrite)

    return(out_hdu)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Related to basic image manipulation
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def jypix_to_mjysr(
        image_fname=None,
        image_hdu=None,
        hdu_to_use=0,
        outfile=None,
        overwrite=True):
    """
    Docs
    """
    
    if image_hdu is None:
        image_hdu = fits.open(image_fname)[hdu_to_use]

    # Could logic check on input units. Now assume Jy/beam

    # Calculate pixel area        
    this_wcs = wcs.WCS(image_hdu.header)
    pix_scale_deg = proj_plane_pixel_scales(this_wcs)
    pix_sr = np.abs(pix_scale_deg[0]*pix_scale_deg[1])* \
        (np.pi/180.)**2

    # Rescale from Jy/pix to MJy/sr
    rescaled_image = image_hdu.data/1E6/pix_sr

    # Update header
    updated_hdr = image_hdu.header
    updated_hdr['BUNIT'] = 'MJy/sr'

    # Create new HDU
    out_hdu = fits.PrimaryHDU(data=rescaled_image, header=updated_hdr)    
    if outfile is not None:
        out_hdu.writeto(outfile, overwrite=overwrite)
    return(out_hdu)
    
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Related to background subtraction
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def plane_func(data, a, b, c):
    """
    Function for plane fitting optimization.
    """
    x, y, = data
    return(a * x + b * y + c)
    
def fit_z0mgs_background(
        image_fname=None,
        image_hdu=None,
        mask_fnames=[],
        mask_hdus=[],
        weight_fname=None,
        weight_hdu=None,
        rad_fname = None,
        rad_hdu = None,
        fid_rad = 15./3600.*u.deg,
        outfile = None,
        outfile_bkgrd = None,        
        overwrite = True,
        methods=['itermed'],
        clip_thresh=3.0,
        niter=5,
        ):
    """
    Fit a background to an image.
    """

    # Read the image
    if image_hdu is None:
        image_hdu = fits.open(image_fname)[0]
    image = image_hdu.data
    image_hdr = image_hdu.header
        
    # Read weights
    if weight_hdu is None:
        if weight_fname is None:
            weight_image = np.isfinite(image_hdu.data)*1.0
        else:
            weight_hdu = fits.open(weight_fname)[0]
            weight_image = weight_hdu.data
        
    # Initialize the background
    bkgrd = np.zeros_like(image)
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Masking
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    have_masks = False
    if (len(mask_fnames) > 0) | (len(mask_hdus) > 0):
        have_masks = True
        
    # Read and stack supplied list of masks
    if have_masks:
        mask_hdu = stack_masks(mask_fnames = mask_fnames, mask_hdus = mask_hdus)
        mask = mask_hdu.data

    # Read the aperture/coordinate definition and mask the galaxy
    if rad_hdu is None:
        if rad_fname is not None:
            rad_hdu = fits.open(rad_fname)[0]
            rad_image = rad_hdu.data*u.deg
            rad_mask = rad_image <= fid_rad
            mask = np.logical_or(mask, rad_mask).astype(int)

    # Mask out regions with no weight
    mask = np.logical_or(mask, (weight_image == 0)).astype(int)
            
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Estimate the noise
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # Scale by normalized weights (if provided) to get a local noise
    # map (especially important for GALEX). Assumption is that weight
    # is integration time-like so that noise \propto weight^(-0.5)

    unmasked_ind = np.where(mask == 0)
    med_weight_unmasked = np.nanmedian(weight_image[unmasked_ind])
    noiselike = 1./np.sqrt(weight_image/med_weight_unmasked)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Initialize background
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    bkgrd = np.zeros_like(image)
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Iteratively subtract a median
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    if 'itermed' in methods:

        print("Using iterative median for background estimation.")
                
        # Initialize an aperture and a rejected pixels mask
        aperture = (mask == 0)
        rejected = np.zeros_like(aperture,dtype=bool)

        converged = False
        prev_med = np.nan
        prev_rms = np.nan
        prev_num_outliers = 0
        
        for ii in range(niter):

            if converged:
                continue
            
            print("... iteration ", ii+1, " of ", niter)
            
            bkgrd_ind = np.where(aperture*(rejected == False))
            n_bkgrd_pix = len(bkgrd_ind[0])
            if n_bkgrd_pix == 0:
                print("No background pixels.")
                continue
            print("... background pixels: ", n_bkgrd_pix) 
            
            data_vec = image[bkgrd_ind]
            noiselike_vec = noiselike[bkgrd_ind]

            # Calculate the median and noise
            rms = mad_std(data_vec, ignore_nan=True)
            med_val = np.nanmedian(data_vec)

            # Find the outliers
            outliers = np.abs(data_vec - med_val) > \
                (clip_thresh * rms / noiselike_vec)

            # Mask them in the aperture
            num_outliers = np.sum(outliers)
            if num_outliers > 0:
                rejected[bkgrd_ind] = \
                    np.logical_or(rejected[bkgrd_ind], outliers)

            converged = False
            if (med_val == prev_med) & (rms == prev_rms) & \
               (prev_num_outliers == num_outliers):
                converged = True
            else:
                prev_med = med_val            
                prev_rms = rms
                prev_num_outliers = num_outliers
                
            print("... rejected outliers: ", num_outliers)
            print("... median value, rms: ", med_val, rms)
            print("... converged: ", converged)

        # Incorporate rejected pixels into the mask
        mask = np.logical_or(mask, rejected).astype(int)

        # Record the background
        bkgrd = bkgrd + med_val
        image = image - med_val
            
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Plane fit
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    if 'planefit' in methods:

        print("Using iterative plane fit for background estimation.")

        # Coordinate images
        ny, nx = image.shape
        x_axis = np.arange(nx)
        y_axis = np.arange(ny)
        x_img, y_img = np.meshgrid(x_axis, y_axis)

        # Empty plane to start
        plane_img = np.zeros_like(image)

        # Initialize mask and rejected pixels
        aperture = (mask == 0)
        rejected = np.zeros_like(aperture,dtype=bool)

        converged = False
        prev_a = np.nan
        prev_b = np.nan
        prev_c = np.nan        
        prev_num_outliers = 0
 
        for ii in range(niter):

            if converged:
                continue
            
            print("... iteration ", ii+1, " of ", niter)

            bkgrd_ind = np.where(aperture*(rejected == False))
            n_bkgrd_pix = len(bkgrd_ind[0])
            if n_bkgrd_pix == 0:
                print("No background pixels.")
                continue
            print("... background pixels: ", n_bkgrd_pix) 

            # Vectorize the data and points
            data_vec = image[bkgrd_ind]
            x_vec = x_img[bkgrd_ind]
            y_vec = y_img[bkgrd_ind]

            # Fit
            xy_stack = np.vstack([x_vec, y_vec])
            fit_params, fit_covariance = \
                curve_fit(plane_func, xy_stack, data_vec)
            a, b, c = fit_params
      
            # Whole background
            plane_img = a * x_img + b * y_img + c
            resid = image - plane_img
            rms = mad_std(resid, ignore_nan=True)
                  
            # Find outliers
            outliers = (np.abs(resid) > (clip_thresh*rms/noiselike))

            # Mask outliers
            num_outliers = np.sum(outliers[bkgrd_ind])
            if num_outliers > 0:
                  rejected[np.where(outliers)] = 1

            converged = False
            if (prev_a == a) & (prev_b == b) & (prev_c == c) \
               & (prev_num_outliers == num_outliers):
                converged = True
            else:
                prev_a = a
                prev_b = b
                prev_c = c
                prev_num_outliers = num_outliers

            print("... fitted plane equation: z = ", a, " * x + ",
                  b, " * y + ", c)
            print("... rejected outliers ", num_outliers)
            print("... converged: ", converged)
        
        # Incorporate rejected pixels into the mask
        mask = np.logical_or(mask, rejected).astype(int)

        image = image - plane_img
        bkgrd = bkgrd + plane_img

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Mode/histogram based calculation
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    if 'mode' in methods:

        print("Using mode of histogram for background estimation.")
        
        aperture = (mask == 0)
        bkgrd_ind = np.where(aperture)
        data_vec = image[bkgrd_ind]

        rms = mad_std(data_vec, ignore_nan=True)
        med_val = np.nanmedian(data_vec)
        
        binsize = 0.05 * rms
        bin_centers = np.arange(-50,51,1.)*binsize + med_val
        bin_edges = np.concatenate([bin_centers-0.5*binsize,
                                    np.array([bin_centers[-1]+0.5*binsize])])
        
        hist_counts, bin_edges_out = np.histogram(data_vec, bins=bin_edges)
        max_bin_ind = np.argmax(hist_counts)
        max_bin_val = bin_centers[max_bin_ind]

        print("... binsize: ", binsize)
        print("... max_bin_x, counts: ", max_bin_val, hist_counts[max_bin_ind])

        bkgrd = bkgrd + max_bin_val
        image = image - max_bin_val
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Median radial profile
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    # TBD

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Write to disk if desired
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    bksub_hdu = fits.PrimaryHDU(image, image_hdr)
    bkgrd_hdu = fits.PrimaryHDU(bkgrd, image_hdr)    
    
    # Write to disk
    if outfile is not None:
        bksub_hdu.writeto(outfile, overwrite=overwrite)

    if outfile_bkgrd is not None:
        bkgrd_hdu.writeto(outfile_bkgrd, overwrite=overwrite)
                  
    return(bksub_hdu,bkgrd_hdu)
        
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Related to visualization
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def show_z0mgs_image(
        image_fname = None,
        image_hdu = None,
        image_hdu_to_use = 0,
        rms = None,
        mask_hdu = None,
        mask_fname = None,
        mask_hdu_to_use = 0,
        mask_levels = None,
        show = False,
        outfile = None,
        title = 'Z0MGS Image',
        value_string = 'image',
        ):
    """Plot a z0mgs image on some standard stretches optionally with a
    mask contour.
    """

    if image_hdu is None:
        image_hdu = fits.open(image_fname)[image_hdu_to_use]

    image_data = image_hdu.data
    image_wcs = wcs.WCS(image_hdu.header)
    if np.sum(np.isfinite(image_data)) == 0:
        print("Image empty: ", image_fname)
        image_data = np.zeros_like(image_data)
    
    mask_data = None
    if mask_hdu is None:
        if mask_fname is not None:
            mask_hdu = fits.open(mask_fname)[mask_hdu_to_use]
    if mask_hdu is not None:
        mask_data = mask_hdu.data
    if mask_levels is None:
        mask_levels = [1.0]
        
    fig = plt.figure(figsize=(14, 8))

    this_norm = simple_norm(image_data, 'log', percent=99.5)
    
    ax1 = fig.add_subplot(1, 2, 1, projection=image_wcs)
    this_im = ax1.imshow(image_data, origin='lower', cmap='Greys', norm=this_norm)
    fig.colorbar(this_im, ax=ax1, fraction=0.05, pad=0.05, label='log10 '+value_string)
    if mask_data is not None:
        ax1.contour(mask_data*1.0, levels=mask_levels, colors='red', alpha=0.7)
    ax1.set_title(title + ' (log scale)')
    ax1.coords[0].set_axislabel('R.A.')
    ax1.coords[1].set_axislabel('Dec.')

    if rms is None:
        ind_for_mad = np.where((image_data != 0.0) & np.isfinite(image_data))
        rms = mad_std(image_data[ind_for_mad])
    meanval = np.nanmean(image_data)
        
    ax2 = fig.add_subplot(1, 2, 2, projection=image_wcs)
    this_im = ax2.imshow(image_data, origin='lower', cmap='Greys',
                         vmin=-5.*rms+meanval, vmax=+10.*rms+meanval)
    fig.colorbar(this_im, ax=ax2, fraction=0.05, pad=0.05, label=value_string)
    if mask_data is not None:
        ax2.contour(mask_data*1.0, levels=mask_levels, colors='red', alpha=0.7)
    ax2.set_title(title + ' (linear)')
    ax2.coords[0].set_axislabel('R.A.')
    ax2.coords[1].set_axislabel('Dec.')

    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
    if show:
        plt.show()

    plt.close()
        
    return()

# TBD - show image, radial profile, high stretch, histogram.

def show_z0mgs_background(
):
    pass
