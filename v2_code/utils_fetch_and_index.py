# Routines to fetch images and data from across the web to provide the
# raw material for the z0MGS atlas, run checks, and index these files
# for use in atlas constuction.

# What's in here:

# FETCH: wget the galex, unWISE, and allWISE data, rsync SDSS
# imaging. Check that we have all the files for each survey. For
# unwise this will implement a wget to retrieve any missing files.

# INDEX: create tables listing files and indexes of these files for
# each survey. These indexes are appropriate to create mosaics for
# individual galaxies.

# Indexing takes several hours (less than a day) for GALEX and the
# unWISE all-sky surveys. It is less than an hour for the unwise
# custom call. It takes O(1 week) for SDSS.

# This all includes some directory structure that's local to each
# system (and you need a place to stash O(20 TB) of stuff). It's not
# expected that this needs to be run in many places.

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Imports
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import os, glob, time
import numpy as np

from astropy.table import Table, Column, vstack
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import astropy.wcs as wcs
from astropy.utils.console import ProgressBar
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord, ICRS
from astropy.coordinates.representation import CartesianRepresentation, SphericalRepresentation

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Fetch various surveys from their archives
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def unwise_fetch(
        version = 'neo9',
        incremental = False,
        dry_run = True
        ):
    """Use WGET to transfer the all-sky versions of unwise NEOWISE 9
    (deepest current UNWISE for bands 1 and 2) and ALLWISE (has W3 and
    W4 but also has the median filter on). You will want to modify
    this to your own local system (recommend scratch space since these
    data are available on the web).
    """
    
# Alternative to unwise.me:
# https://portal.nersc.gov/project/cosmo/data/unwise/neo9/unwise-coadds/fulldepth/
    
    if version == 'neo9':
        prefix = '/data/fft_scratch/leroy.42/allsky/unwise_neo9/'
        target = 'https://unwise.me/data/neo9/unwise-coadds/fulldepth/'
        list_file = '../support_files/unwise_neowise_find.lst'
        extra_nonsense = 'data/neo9/unwise-coadds/fulldepth/'
    if version == 'allwise':
        prefix = '/data/fft_scratch/leroy.42/allsky/unwise_allwise/'
        target = 'https://unwise.me/data/allwise/unwise-coadds/fulldepth/'
        list_file = '../support_files/unwise_allwise_find.lst'
        extra_nonsense = 'data/allwise/unwise-coadds/fulldepth/'

    if incremental:
        tab = ascii.read(list_file, data_start=0, names=['fname'])
        fname_root = prefix + extra_nonsense
        tab['found'] = False
        print("Identifying missing files.")
        for this_row in ProgressBar(tab):
            this_fname = fname_root + this_row['fname'].strip()
            if os.path.isfile(this_fname):
                this_row['found'] = True
        needed_tab = tab[tab['found']==False]
        for this_row in needed_tab:
            out_file = (fname_root + this_row['fname'].strip())
            out_dir = out_file[:(-1*len((out_file.split('/'))[-1]))]            
            wget_call = 'wget -P '+ out_dir + ' '+target+this_row['fname'].strip()
            os.system(wget_call)
    else:
        print("This is a big call, you should run this by hand. I will print it but not run it.")
        this_call = 'wget -nH -r -P '+prefix+' '+target
        print(this_call)
        #if dry_run == False:
        #    os.system(this_call)

    return()

def sdss_fetch(
        cd_to_dir = '/data/fourier/leroy.42/allsky/sdss_dr12/',
        out_dir = 'frames/',
        dry_run = True,
        ):
    """
    Execute an rsync call to pull SDSS DR 12 images.
    """

    os.system(cd_to_dir)

    # See:
    # ../support_files/sdss_rsync_call.sh
    
    dry_run_call = 'rsync -avz --dry-run --include="*/" --include="frame-*.fits.bz2" --exclude="*" rsync://dtn.sdss.org/dr9/boss/photoObj/frames/ '+out_dir
    full_call =    'rsync -avz --include="*/" --include="frame-*.fits.bz2" --exclude="*" rsync://dtn.sdss.org/dr9/boss/photoObj/frames/ '+out_dir

    if dry_run:
        this_call = dry_run_call
    else:
        this_call = full_call

    print(this_call)
    if dry_run == False:
        print("This is a big call, you should run this by hand. I will print it but not run it.")        
        #os.system(this_call)

    return()

def galex_fetch():

    # See ../support_files/wget_galex_list.txt

    # TBD implement incremental pull / file checking
    
    pass

def gaia_fetch():

    # This can be done via globus from Flatiron, otherwise write a
    # call here to wget the ESA DR3.
    
    pass

def gaia_convert_to_fits():
    """Gaia serves their tables as ecsv but these are slow for astropy to
    read. Convert to FITS tables which read quickly (could also
    presumably write a sharper CSV reader). This reduces read time by
    about 30x.
    """
    
    working_dir = '/data/fft_scratch/leroy.42/allsky/gaia/'
    selection = '*.csv'

    flist = glob.glob(working_dir+selection)
    for this_fname in ProgressBar(flist):
        this_tab = Table.read(this_fname, format='ascii.ecsv')
        outfile = this_fname.replace('.csv','.fits')
        print("Will write to: ", outfile)
        if os.path.isfile(outfile):
            print("... it exists already.")
            continue
        this_tab.write(outfile, format='fits')

    return()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Index the fetched surveys
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def compile_list_of_images(
        survey=None,
        root_dir=None,
        selection=None,
        tab_dir=None,
        tab_file=None):
    """Compile a list of image files for one of our surveys. Creates a
    table that can then be used in indexing.

    """

    # Specify input and output for surveys
    
    if survey is not None:

        if survey == 'sdss':
            print("")
            print("Compiling list of all files in SDSS DR12 frames directory.")
            print("")
            
            root_dir = '/data/fourier/leroy.42/allsky/sdss_dr12/frames/301/'

            # Search for only g band files, the other bands cover the
            # same footprints and we can swap out names.

            # selection = '*/*/*.fits.bz2'
            selection = '*/*/frame-g*.fits.bz2'
            tab_dir = '../../working_data/sdss/index/'
            tab_file = 'sdss_frame_list.fits'

        if survey == 'galex':
            print("")
            print("Compiling list of all files in GALEX directory.")
            print("")
            
            root_dir = '/data/fourier/leroy.42/allsky/all_galex_tiles/'
            selection = '*d/*-int.fits.gz'
            tab_dir = '../../working_data/galex/index/'
            tab_file = 'galex_tile_list.fits'

        if survey == 'unwise_custom':
            print("")
            print("Compiling list of all files in UNWISE custom directory.")
            print("(this is a recursive call)")
            print("")

            root_dir = '../../orig_data/unwise_v2/v2_karachentsev/'
            selection = '*/*/unwise-*-w*-img-m.fits'
            tab_dir = '../../working_data/unwise/index/'
            tab_file = 'unwise_custom_localvolume_list.fits'
            print("... localvolume")
            tab_localvolume = compile_list_of_images(
                root_dir=root_dir, selection=selection,
                tab_dir=tab_dir, tab_file=tab_file)

            root_dir = '../../orig_data/unwise_v2/v2_largegals/'
            selection = '*/unwise-*-w*-img-m.fits'
            tab_dir = '../../working_data/unwise/index/'
            tab_file = 'unwise_custom_largeleda_list.fits'
            print("... largeleda")
            tab_largeleda = compile_list_of_images(
                root_dir=root_dir, selection=selection,
                tab_dir=tab_dir, tab_file=tab_file)            

            root_dir = '../../orig_data/unwise_v2/v2_smallgals/'
            selection = '*/*/unwise-*-w*-img-m.fits'
            tab_dir = '../../working_data/unwise/index/'
            tab_file = 'unwise_custom_smallleda_list.fits'
            print("... smallleda")
            tab_smallleda = compile_list_of_images(
                root_dir=root_dir, selection=selection,
                tab_dir=tab_dir, tab_file=tab_file)            

            root_dir = '../../orig_data/unwise_v2/v2_unmanga/'
            selection = '*/unwise-*-w*-img-m.fits'
            tab_dir = '../../working_data/unwise/index/'
            tab_file = 'unwise_custom_manga_list.fits'
            print("... manga")
            tab_manga = compile_list_of_images(
                root_dir=root_dir, selection=selection,
                tab_dir=tab_dir, tab_file=tab_file)            

            root_dir = '../../orig_data/unwise_v2/v2_localgroup/'
            selection = '*/*/unwise-*-w*-img-m.fits'
            tab_dir = '../../working_data/unwise/index/'
            tab_file = 'unwise_custom_localgroup_list.fits'
            print("... localgroup")
            tab_localgroup = compile_list_of_images(
                root_dir=root_dir, selection=selection,
                tab_dir=tab_dir, tab_file=tab_file)            

            tab_largeleda['subsample'] = ' ' * 15
            tab_largeleda['subsample'] = 'largeleda'

            tab_smallleda['subsample'] = ' ' * 15
            tab_smallleda['subsample'] = 'smallleda'

            tab_localvolume['subsample'] = ' ' * 15
            tab_localvolume['subsample'] = 'localvolume'

            tab_localgroup['subsample'] = ' ' * 15
            tab_localgroup['subsample'] = 'localgroup'

            tab_manga['subsample'] = ' ' * 15
            tab_manga['subsample'] = 'manga'
        
            tab_fname = '../../working_data/unwise/index/unwise_custom_list.fits'
            print("I wrote the combined stack to: ", tab_fname)
            tab = vstack([tab_largeleda, tab_smallleda, tab_localvolume, tab_localgroup, tab_manga])
            tab.write(tab_fname, format='fits', overwrite=True)
            
            return()
            
        if survey == 'unwise_allwise':
            print("")
            print("Compiling list of all files in UNWISE allwise directory.")
            print("")

            root_dir = '/data/fft_scratch/leroy.42/allsky/unwise_allwise/data/allwise/unwise-coadds/fulldepth/'
            selection = '*/*/*-w*-img-m.fits'
            tab_dir = '../../working_data/unwise/index/'
            tab_file = 'unwise_allwise_list.fits'
            
        if survey == 'unwise_neowise':
            print("")
            print("Compiling list of all files in NEOWISE neowise directory.")
            print("")

            root_dir = '/data/fft_scratch/leroy.42/allsky/unwise_neo9/data/neo9/unwise-coadds/fulldepth/'
            selection = '*/*/*-w*-img-m.fits'
            tab_dir = '../../working_data/unwise/index/'
            tab_file = 'unwise_neowise_list.fits'

        if survey == 'gaia':
            print("")
            print("Compiling list of all files in GAIA directory.")
            print("")

            root_dir = '/data/fft_scratch/leroy.42/allsky/gaia/'
            #selection = '*.csv'
            selection = '*.fits'
            tab_dir = '../../working_data/gaia/index/'
            tab_file = 'gaia_list.fits'
            
    # Find the relevant files
    
    flist = glob.glob(root_dir+selection)

    start_time = time.time()
    
    n_files = len(flist)
    id = np.arange(0,n_files,dtype=np.int32)

    tab = Table([id,np.array(flist)],names=('id','fname'))
    print("... I found ", len(tab), " images.")
    
    table_fname = tab_dir+tab_file
    print("... Writing them to ", table_fname)
    tab.write(table_fname, format='fits', overwrite=True)

    stop_time = time.time()

    print("This took : ", stop_time-start_time, " seconds.")

    return(tab)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Index the properties of the files
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def index_image_list(
        survey = None,
        list_table = None,
        outfile = None,
        do_all_frames=False,
        incremental=False):
    """
    Read a table  of FITS files and build these into an index.
    """

    if survey is not None:

        if survey == 'sdss':
            print("")
            print("Indexing SDSS tiles")
            print("")
            index_dir = '../../working_data/sdss/index/'
            list_table = index_dir + 'sdss_frame_list.fits'
            outfile = index_dir+'sdss_frame_index.fits'

        if survey == 'galex':
            print("")
            print("Indexing GALEX tiles")
            print("")
            index_dir = '../../working_data/galex/index/'
            list_table = index_dir + 'galex_tile_list.fits'
            outfile = index_dir+'galex_tile_index.fits'            

        if survey == 'unwise_custom':
            print("")
            print("Indexing custom UNWISE tiles")
            print("")
            index_dir = '../../working_data/unwise/index/'
            list_table = index_dir + 'unwise_custom_list.fits'
            outfile = index_dir+'unwise_custom_index.fits'

        if survey == 'unwise_allwise':
            print("")
            print("Indexing UNWISE allwise tiles.")
            print("")
            
            index_dir = '../../working_data/unwise/index/'
            list_table = index_dir + 'unwise_allwise_list.fits'
            outfile = index_dir+'unwise_allwise_index.fits'

        if survey == 'unwise_neowise':
            print("")
            print("Indexing UNWISE neowise tiles.")
            print("")

            index_dir = '../../working_data/unwise/index/'
            list_table = index_dir + 'unwise_neowise_list.fits'
            outfile = index_dir+'unwise_neowise_index.fits'

        if survey == 'gaia':
            print("")
            print("Indexing GAIA tables.")
            print("")

            index_dir = '../../working_data/gaia/index/'
            list_table = 'gaia_list.fits'
            outfile = 'gaia_index.fits'

            index_gaia_tables(
                list_table=index_dir+list_table,
                outfile=index_dir+outfile)

            return()
            
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Loop over table
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    start_time = time.time()

    # Read table
    tab = Table.read(list_table, format='fits')
    n_files = len(tab)
    
    # Initialize output
    for this_field in ['ctr_ra','ctr_dec']:
        tab[this_field] = np.nan
    for this_field in ['nx','ny']:
        tab[this_field] = int(0)
    for this_field in ['pix_scale_x','pix_scale_y']:
        tab[this_field] = np.nan

    tab['filter'] = ' ' * 100
        
    if survey == 'sdss':
        tab['u'] = False
        tab['g'] = False
        tab['r'] = False
        tab['i'] = False
        tab['z'] = False
        tab['exptime'] = np.nan
        tab['filled'] = False
        increment_to_write = 1000

    if survey.count('unwise') > 0:
        pass
        
    if survey == 'galex':
        tab['rrhr_fname'] = ' ' * 150
        tab['flag_fname'] = ' ' * 150
        tab['bgsub_fname'] = ' ' * 150

    # If doing an sdss incremental call read the table file from disk
    # to pick up where we left off.
    if (survey == 'sdss') & incremental:
        print("Reading the incremental file from disk.")
        tab = Table.read(outfile, format='fits')
        
    counter = 0

    for this_row in ProgressBar(tab):

        # If we're in incremental mode, skip finished rows
        if (survey == 'sdss') & incremental:
            if this_row['filled']:
                continue            
            
        # Read, extract header and wcs from file
        this_fname = this_row['fname']
        this_hdulist = fits.open(this_fname)
        this_header = this_hdulist[0].header
        w = wcs.WCS(this_header)

        if survey.count('unwise') > 0:
            if this_fname.count('w1') > 0:
                this_row['filter'] = 'w1'
            if this_fname.count('w2') > 0:
                this_row['filter'] = 'w2'
            if this_fname.count('w3') > 0:
                this_row['filter'] = 'w3'
            if this_fname.count('w4') > 0:
                this_row['filter'] = 'w4'

        if survey == 'galex':
            this_rrhr_fname = this_fname.replace(
                '-int.fits.gz', '-rrhr.fits.gz')
            if os.path.isfile(this_rrhr_fname):
                this_row['rrhr_fname'] = this_rrhr_fname
            
            this_flag_fname = this_fname.replace(
                '-int.fits.gz', '-flags.fits.gz')
            if os.path.isfile(this_flag_fname):
                this_row['flag_fname'] = this_flag_fname
                              
            this_bgsub_fname = this_fname.replace(
                '-int.fits.gz', '-intbgsub.fits.gz')
            if os.path.isfile(this_bgsub_fname):
                this_row['bgsub_fname'] = this_bgsub_fname
            
            if this_fname.count('-nd-') > 0:
                this_row['filter'] = 'nuv'
            if this_fname.count('-fd-') > 0:
                this_row['filter'] = 'fuv'
            
        # Figure out and record corner and center coordinates
        nx = this_header['NAXIS1']
        ny = this_header['NAXIS2']
        tab['nx'] = nx
        tab['ny'] = ny
            
        # Convert to world coordinates
        ctr_world =  w.wcs_pix2world(nx/2.,ny/2.,0)
        this_row['ctr_ra'] = ctr_world[0]
        this_row['ctr_dec'] = ctr_world[1]

        # Pixel scale
        pix_scale = wcs.utils.proj_plane_pixel_scales(w)
        this_row['pix_scale_x'] = pix_scale[0]
        this_row['pix_scale_y'] = pix_scale[1]

        #for this_field in this_row.colnames:
        #    print(this_field, this_row[this_field])

        # Survey specific items
        if survey == 'sdss':            
            this_row['exptime'] = float(this_header['EXPTIME'])
            for this_filter in ['u','g','r','i','z']:
                filter_fname = this_fname.replace('frame-g','frame-'+this_filter)
                if os.path.isfile(filter_fname):
                    this_row[this_filter] = True
                else:
                    print("Missing ", filter_fname)

            this_row['filled'] = True
            
        counter += 1

        # For the big SDSS job write our progress
        if (survey == 'sdss'):
            if counter % increment_to_write == 0:
                print("")
                print("Writing incremental table.")
                print("")
                tab.write(outfile, format='fits', overwrite=True)
                
    tab.write(outfile, format='fits', overwrite=True)

    stop_time = time.time()

    print("This took : ", stop_time-start_time, " seconds.")
    
    return(tab)

def index_gaia_tables(
        list_table = None,
        outfile = None,
):
    """
    Index the gaia tables.
    """

    start_time = time.time()

    # Read table
    index = Table.read(list_table, format='fits')
    n_files = len(index)
    
    # Initialize output
    for this_field in ['ctr_ra','ctr_dec','extent']:
        index[this_field] = np.nan

    for this_row in ProgressBar(index):

        this_tab = Table.read(this_row['fname'], format='ascii.ecsv')

        fid_dist = 10.*u.mpc
        eq_coords = SkyCoord(ra=this_tab['ra'], dec=this_tab['dec'], distance=fid_dist, frame='icrs')
        cart_coords = eq_coords.cartesian

        mean_x = np.mean(cart_coords.x)
        mean_y = np.mean(cart_coords.y)
        mean_z = np.mean(cart_coords.z)

        cart_rep = CartesianRepresentation(mean_x, mean_y, mean_z)
        sphere_rep = cart_rep.represent_as(SphericalRepresentation)
        ctr_coord = SkyCoord(sphere_rep).transform_to(ICRS)
        
        #ctr_coord = (SkyCoord(x=mean_x, y=mean_y, z=mean_z, representation_type='cartesian')).transform_to('icrs')
        this_row['ctr_ra'] = ctr_coord[0]
        this_row['ctr_dec'] = ctr_coord[1]

        separations = eq_coords.separation(ctr_coord)
        max_extent = np.max(separations)

        this_row['extent'] = max_extent
        print(this_row)
        
    index.write(outfile, format='fits', overwrite=True)
        
    stop_time = time.time()
    
    print("This took : ", stop_time-start_time, " seconds.")
    
    return(index)

    
