# Utilities for the nuts-and-bolts specifics of z0mgs atlas
# construction. Creates tables to loop over, fetches PSFs, etc.. This
# is where most specifics related to the atlas goes.

# What's in here:

# This set of routines handles specific directories and tables.

# FILE, TABLE, and GALAXY management
# - routines to create tables of galaxies to loop over
# - routines to manage the directory structure
# - routines to get files for PSFs or kernels

# Basics
import os
import numpy as np

# Astropy stuff
import astropy.io.fits as fits
from astropy.table import Table, vstack
import astropy.units as u
from astropy.utils.console import ProgressBar
import astropy.wcs as wcs
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.stats import mad_std

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Construct tables of galaxies to be processed
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# This manages tables and subsample definitions to create a table that
# the atlas construction pipelines will loop over. It assumes that the
# information about the galaxies comes from somewhere external.

# TBD - refactor this to be simpler in terms of what's in the tables
# of galaxies of interest. Each subsample gets a list and center
# coords and a size. The size can have use-defined aspects.

def build_target_table(
        subsamples=['all'],
        table_dir='../../measurements/',
        just_galaxy=None,
        skip_galaxy=None,
        start_galaxy=None,
        stop_galaxy=None,
):
    """Build and return a table of targets for use in z0mgs optical, UV,
    or other atlas construction.

    Parameters
    ----------

    subsamples : Default ['all']. List of subsamples (strings) to
    include in the output table.

    table_dir : Directory where the tables defining the subsamples are
    found.

    just_galaxy : List of galaxies to include. Default empty (i.e., include all).

    skip_galaxy : List of galaxies to omit. Default empty.
    
    start_galaxy : Index of galaxy to start at. Useful to break apart
    large calls into smaller ones.

    stop_galaxy : Index of galaxy to stop at. Useful to break apart
    large calls into smaller ones.

    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=            
    # Definitions
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    # These are tables where we keep the definition of the z0mgs
    # subsamples. Could revise this later.

    subsample_tables = {
        'localgroup':table_dir+'unwise_v2_index_localgroup.fits',
        'localvolume':table_dir+'unwise_v2_index_localvolume.fits',
        'largeleda':table_dir+'unwise_v2_index_largeleda.fits',
        'smallleda':table_dir+'unwise_v2_index_smallleda.fits',
        'manga':table_dir+'unwise_v2_index_manga.fits',
    }
    subsample_list = subsample_tables.keys()

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
    # Selections
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    # Define which subsamples we are working with
    
    if subsamples is None:

        subsamples = ['all']

    if not isinstance(subsamples, list):
        subsamples = [subsamples]

    if subsamples == ['all']:

        subsamples = subsample_list

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Load the tables
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

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

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Downselect    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

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
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Return
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
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
    """Routine to make a table for just one target based on
    coords. Contains the necessary information to be fed into atlas
    creation, but otherwise minimal. 

    Helper routine to allow the pipeline to be run one new target at a
    time. To run on just one galaxy already in the subsample tables
    instead use build_target_table with a 'just' parameter.

    Parameter definitions (TBD)

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

# Make/check directory structure (convenience function)

def make_z0mgs_directories(
        root_dir='../../working_data/',
        surveys=['galex','unwise','sdss']):
    """Routine to make the z0mgs directory structure expected by the other
    programs. Includes a hard-coded list of subdirectories and
    samples.
    """
    
    subsample_list = [
        'largeleda','smallleda',
        'localgroup','localvolume',
        'manga','other']

    stages_list = [
        'staged', 'convolved',
        'final', 'index', 'masks',
        'coords', 'bkgrd', 'star_stacks',
    ]

    for this_survey in surveys:

        if not os.path.isdir(root_dir+this_survey):
            os.system('mkdir '+root_dir+this_survey)

        for this_stage in stages_list:
            this_dir = root_dir + this_survey + '/' + this_stage
            if not os.path.isdir(this_dir):
                os.system('mkdir '+this_dir)

            for this_subsample in subsample_list:
                this_dir = root_dir + this_survey + '/' + \
                    this_stage + '/' + this_subsample
                if not os.path.isdir(this_dir):
                    os.system('mkdir '+this_dir)

    return()

def z0mgs_psf_name(
        band = None,
):
    """
    Return a file path to a PSF file for the specified band.

    TBD Error trapping.
    """
    
    psf_dir = '/data/bell-tycho/leroy.42/ellohess/kernels/PSF_FITS_Files/'

    psf_dict = {}
    psf_dict['fuv'] = 'PSF_Corrected_GALEX_FUV_added_wing.fits'
    psf_dict['nuv'] = 'PSF_Corrected_GALEX_NUV_added_wing.fits'
    psf_dict['w1'] = 'PSF_Corrected_WISE_ATLAS_3.4_added_wing.fits'
    psf_dict['w2'] = 'PSF_Corrected_WISE_ATLAS_4.6_added_wing.fits'
    psf_dict['w3'] = 'PSF_Corrected_WISE_ATLAS_11.6_added_wing.fits'
    psf_dict['w4'] = 'PSF_Corrected_WISE_ATLAS_22.1_added_wing.fits'

    if band in psf_dict:
        psf_name = psf_dir+psf_dict[band]
        return(psf_name)
    else:
        return(None)

def z0mgs_kernel_name(
        from_res = None,
        to_res = None):
    """Return kernel name for use in convolution from the band or kernel
    'from_res' to the band or kernel 'to_res.'
    """

    # TBD - add file checking and more flexibility on directories
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
        kernel_bands[from_res] + '_to_' + \
        kernel_bands[to_res] + '.fits'
    
    return(kernel_name)
