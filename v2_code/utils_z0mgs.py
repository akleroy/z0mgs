# Utilities related to z0MGS atlas construction

# What's in here:

# - routines to create tables of galaxies to loop over
# - routines to query Gaia and create stacks of images of stars
# - routines to predict stellar fluxes from Gaia
# - routines to make an image of stellar fluxes for various bands

# TBD
# - routines to fit linear and planar backgrounds

# Imports
import os
import numpy as np

from astropy.convolution import convolve_fft
from astropy.coordinates import SkyCoord
import astropy.io.fits as fits
from astropy.table import Table, vstack
import astropy.units as u
from astropy.utils.console import ProgressBar
import astropy.wcs as wcs
from astropy.wcs.utils import proj_plane_pixel_scales

from astroquery.gaia import Gaia

from reproject import reproject_interp

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Construct tables of galaxies to be processed
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# This manages tables with the relevant format to create a table that
# the atlas construction pipelines will loop over. It assumes that the
# information about the galaxies comes from somewhere external.

def build_target_table(
        subsamples=None,
        table_dir='../../measurements/',
        just_galaxy=None,
        skip_galaxy=None,
        start_galaxy=None,
        stop_galaxy=None,
):
    """
    Build and return a table of targets for use in z0mgs optical, UV,
    or other atlas construction.
    """

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Definitions
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # These are tables where we keep data on the z0mgs subsamples.

    subsample_tables = {
        'localgroup':table_dir+'unwise_v2_index_localgroup.fits',
        'localvolume':table_dir+'unwise_v2_index_localvolume.fits',
        'largeleda':table_dir+'unwise_v2_index_largeleda.fits',
        'smallleda':table_dir+'unwise_v2_index_smallleda.fits',
        'manga':table_dir+'unwise_v2_index_manga.fits',
    }
    subsample_list = subsample_tables.keys()

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Selections
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    
    # Define which subsamples we are working with
    
    if subsamples is None:

        subsamples = ['all']

    if not isinstance(subsamples, list):
        subsamples = [subsamples]

    if subsamples == ['all']:

        subsamples = subsample_list

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Load the tables
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    targets_tab = None
    for this_subsample in subsamples:

        if this_subsample not in subsample_list:
            print("Invalid subsample: ", this_subsample)
        
        this_tab = Table.read(subsample_tables[this_subsample],
                              format='fits')
        this_tab['SUBSAMPLE'] = this_subsample

        if targets_tab is None:
            targets_tab = this_tab
        else:
            targets_tab = vstack(targets_tab, this_tab)
    
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Downselect
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # Initialize a mask
    targets_tab['COUNTER'] = 0
    targets_tab['USE_THIS_ROW'] = True
    
    # Loop on just/skip
    for ii, this_row in enumerate(targets_tab):

        this_name = this_row['Z0MGS_NAME'].strip()

        if just_galaxy is not None:            
            if this_name not in just_galaxy:
                this_row['USE_THIS_ROW'] = False
                
        if skip_galaxy is not None:
            if this_name in skip_galaxy:
                this_row['USE_THIS_ROW'] = False

    # Down select to manual selection
    targets_tab = targets_tab[targets_tab['USE_THIS_ROW']]
    
    # Loop on just/skip
    for ii, this_row in enumerate(targets_tab):

        targets_tab['COUNTER'] = ii
        
        if start_galaxy is not None:
            if ii < start_galaxy:
                this_row['USE_THIS_ROW'] = False

        if stop_galaxy is not None:
            if ii > stop_galaxy:
                this_row['USE_THIS_ROW'] = False

    # Down select to manual selection
    targets_tab =  targets_tab[targets_tab['USE_THIS_ROW']]
    
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Return
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    return(targets_tab)
    
def build_tab_for_one_target(
        name = 'GALNAME',
        pgc = -1,
        subsample = '',
        w1_fname = '',
        ra_ctr = 12.0,
        dec_ctr = 30.0,
        extent_arcmin = 1.0,
):
    """Routine to make a table for just one target based on coords, to be
    fed into atlas creation. Helper routine to allow the pipeline to
    be run one target at a time.
    """
    
    ctr_coords = SkyCoord(ra=ra_ctr*u.deg, dec=dec_ctr*u.deg, frame='icrs')
    offset_deg = extent_arcmin / 60. * np.sqrt(2) * u.deg
    trc_pa = -45.*u.deg
    blc_pa = 135.*u.deg
    trc_coords = ctr_coords.directional_offset_by(trc_pa, offset_deg)
    blc_coords = ctr_coords.directional_offset_by(blc_pa, offset_deg)
    
    gal_dict = \
        [{'Z0MGS_NAME':name,
          'PGC':pgc,
          'SUBSAMPLE':subsample,
          'W1_FNAME':w1_fname,
          'RA_CTR':ctr_coords.ra.value,
          'DEC_CTR':ctr_coords.dec.value,
          'BLC_RA':blc_coords.ra.value,
          'BLC_DEC':blc_coords.dec.value,
          'TRC_RA':trc_coords.ra.value,
          'TRC_DEC':trc_coords.dec.value,
          'USE_FUV':1,
          'USE_NUV':1,
          'USE_W1':1,
          'USE_W2':1,
          'USE_W3':1,
          'USE_W4':1}]

    gal_tab = Table(gal_dict)

    return(gal_tab)

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
# Make cutouts around GAIA sources (or any sources)
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Make small cutouts around a set of Gaia sources. This is in
# principle useful to test the PSF in an image or to model how to
# convert stellar fluxes to the atlas bands.

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

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Related to predicting a stellar flux from Gaia
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
    """
    
    # Hard-code the conversions from z0mgs paper but in flux units
    # (i.e., predict Jy)
    gauss7p5_sr = (7.5/3600.*np.pi/180./2.0)**2*np.pi/np.log(2)
    gauss7p5_intens_to_flux = gauss7p5_sr*1e6

    gauss15_sr = (15/3600.*np.pi/180./2.0)**2*np.pi/np.log(2)
    gauss15_intens_to_flux = gauss15_sr*1e6

    # coefficient in flux = 10.**(-1.0*mag/2.5)*coefficient
    mag_to_jy = {
        'fuv': 27.07*gauss_7p5_intens_to_flux,
        'nuv': 13060.2*gauss_7p5_intens_to_flux,
        'w1': 829914.*gauss_7p5_intens_to_flux,
        'w2': 450745.*gauss_7p5_intens_to_flux,
        'w3': 77349.7*gauss_7p5_intens_to_flux,
        'w4': 8398.*gauss_15_intens_to_flux,
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
        'fuv': 1.35*gauss_7p5_intens_to_flux,
        'nuv': 202.0*gauss_7p5_intens_to_flux,
        'w1': 164900.*gauss_7p5_intens_to_flux,
        'w2': 91800.*gauss_7p5_intens_to_flux,
        'w3': 16800.*gauss_7p5_intens_to_flux,
        'w4': 1550.*gauss_15_intens_to_flux,
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
):
    """
    Build an image in Jy with flux.
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
        gaia_ra = gaia_tab['ra'] * u.deg
        gaia_dec =  gaia_tab['dec'] * u.deg
        gaia_coords = SkyCoord(ra=gaia_ra, dec=gaia_dec, frame='icrs')
        
        gaia_pix = out_wcs.world_to_pixel(gaia_coords)
        gaia_x = gaia_pix[0]
        gaia_y = gaia_pix[1]

        in_image = (gaia_x >= 0) & (gaia_x < nx) \
            & (gaia_y >= 0) & (gaia_y < ny)

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
            g_mag = gaia_tab['g_mag'],
            bp = gaia_tab['bp'],
            rp = gaia_tab['rp'],
            parallax = gaia_tab['parallax']
        )
        
        # Add the flux of each Gaia star to the image. Must be a
        # better way to do this.
        for ii, this_flux in enumerate(gaia_fluxes):
            this_x = gaia_x[ii]
            this_y = gaia_y[ii]
            out_map[this_y, this_x] += this_flux
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # If supplied add the 2MASS based stars to the image
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # Read the 2MASS file or table
    if ks_tab is None and ks_file is not None:

        ks_tab = Table.read(ks_file, format='fits')

        # Identify the Gaia sources within the image
        ks_ra = ks_tab['ra'] * u.deg
        ks_dec =  ks_tab['dec'] * u.deg
        ks_coords = SkyCoord(ra=ks_ra, dec=ks_dec, frame='icrs')
        
        ks_pix = out_wcs.world_to_pixel(ks_coords)
        ks_x = pixel_coords[0]
        ks_y = pixel_coords[1]

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
            ks_mag = ks_tab['ks_mag'],
        )
        
        # Add the flux of each Gaia star to the image. Must be a
        # better way to do this.
        for ii, this_flux in enumerate(ks_fluxes):
            this_x = ks_x[ii]
            this_y = ks_y[ii]
            out_map[this_y, this_x] += this_flux        

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Return and write to disk if requested
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    out_hdu = fits.PrimaryHDU(out_map, out_hdr)
    
    # Write to disk
    if outfile is not None:
        out_hdu.writeto(outfile)

    return(out_hdu)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Build galaxy masks
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
    """Use the galaxy database table to build a mask of galaxies (other
    than the target).
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

        print(this_gal['INCL_DEG'], this_gal['POSANG_DEG'])
        
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
    
    Both deprojected radii and projected angles are defined relative to the center in the inclined disk frame. For (1) and (2), the outputs are 2D images; for (3), the outputs are arrays with shapes matching the broadcasted shape of `ra` and `dec`.

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

def z0mgs_kernel_name(
        from_band = None,
        to_band = None):
    """
    Return kernel name for use in convolution.
    """

    # TBD - add file checking and more flexibility on directories
    
    psf_dir = '/data/bell-tycho/leroy.42/ellohess/kernels/PSF_FITS_Files/'

    psf_dict = {}
    psf_dict['fuv'] = 'PSF_Corrected_GALEX_FUV_added_wing.fits'
    psf_dict['nuv'] = 'PSF_Corrected_GALEX_NUV_added_wing.fits'
    psf_dict['w1'] = 'PSF_Corrected_WISE_ATLAS_3.4_added_wing.fits'
    psf_dict['w2'] = 'PSF_Corrected_WISE_ATLAS_4.6_added_wing.fits'
    psf_dict['w3'] = 'PSF_Corrected_WISE_ATLAS_11.6_added_wing.fits'
    psf_dict['w4'] = 'PSF_Corrected_WISE_ATLAS_22.1_added_wing.fits'

    if to_band == 'native' or to_band == 'psf':
        return(psf_dir + psf_dict)

    kernel_dir = '../../kernels/'

    kernel_bands = {
        'fuv': 'GALEX_FUV',
        'nuv': 'GALEX_NUV',
        'w1': 'WISE_FRAME_3.4',
        'w2': 'WISE_FRAME_4.6',
        'w3': 'WISE_FRAME_11.6',
        'w4': 'WISE_FRAME_22.1',
        'gauss7p5': 'Gauss_07.5',
        'gauss11': 'Gauss_11',
        'gauss15': 'Gauss_15',
        'gauss20': 'Gauss_20',
    }
    
    kernel_name = kernel_dir + 'Kernel_LoRes_' + \
        kernel_bands[from_band] + '_to_' + \
        kernel_bands[to_band] + '.fits'
    
    return(kernel_name)

def convolve_image_with_kernel(
        image_file=None,
        image_hdu=None,
        kernel_file=None,
        kernel_hdu=None,
        outfile=None,
        blank_zeros=True,
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
        preserve_nan=True,
        fill_value=np.nan,
    )
    
    # If an error map is present, convolve errors (with kernel**2,
    # do not normalize it).  This, however, doesn't account for
    # covariance between pixels
    
    if use_err:
        conv_err = np.sqrt(
            convolve_fft(
                image_hdu["ERR"].data ** 2,
                kernel_interp ** 2,
                preserve_nan=True,
                allow_huge=True,
                normalize_kernel=False,
            )
        )
        
    image_hdu[sci_ext].data = conv_im
    if use_err:
        image_hdu["ERR"].data = conv_err
        
    if outfile is not None:
        image_hdu.writeto(file_out, overwrite=overwrite)

    return(image_hdu)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Related to background subtraction
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def subtract_z0mgs_background():

    pass
