# Routines to fetch images and data from across the web to provide the
# raw material for the atlas, run checks, and index these files for
# use in atlas constuction.

# What's in here:

# - wget the galex, unwise, and allwise data, rsync SDSS imaging
# - create tables listing files and indexes of these files

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Imports
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import os, glob, time
import numpy as np
from astropy.table import Table
import astropy.io.fits as fits
import astropy.wcs as wcs
from astropy.utils.console import ProgressBar
from astropy.wcs.utils import proj_plane_pixel_scales

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
    W4 but also has the median filter on).

    TBD adding file checking against the provided file list.
    """
    
# Alternative:
# https://portal.nersc.gov/project/cosmo/data/unwise/neo9/unwise-coadds/fulldepth/
    
    if version == 'neo9':
        prefix = '/data/fft_scratch/leroy.42/allsky/unwise_neo9/'
        target = 'https://unwise.me/data/neo9/unwise-coadds/fulldepth/'
        list_file = 'unwise_neowise_find.lst'

    if version == 'allwise':
        prefix = '/data/fft_scratch/leroy.42/allsky/unwise_allwise/'
        target = 'https://unwise.me/data/neo9/unwise-coadds/fulldepth/'
        list_file = 'unwise_allwise_find.lst'

    if incremental:
        tab = ascii.read(list_file, data_start=0, names=['fname'])
        # TBD
    else:
        this_call = 'wget -nH -r -P '+prefix+' '+target
        print(this_call)
        if dry_run == False:
            os.system(this_call)

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
    
    dry_run_call = 'rsync -avz --dry-run --include="*/" --include="frame-*.fits.bz2" --exclude="*" rsync://dtn.sdss.org/dr9/boss/photoObj/frames/ '+out_dir
    full_call =    'rsync -avz --include="*/" --include="frame-*.fits.bz2" --exclude="*" rsync://dtn.sdss.org/dr9/boss/photoObj/frames/ '+out_dir

    if dry_run:
        this_call = dry_run_call
    else:
        this_call = full_call

    print(this_call)
    if dry_run == False:
        os.system(this_call)

    return()

def galex_fetch():

    # wget_galex_list.txt
    
    pass

def gaia_fetch():
    pass

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Index the fetched surveys
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def compile_list_of_images(
        survey=None,
        root_dir=None,
        selection=None,
        tab_dir=None,
        tab_file=None):
    """
    Compile a list of image files.
    """

    # Specify input and output for surveys
    
    if survey is not None:

        if survey == 'sdss':
            print("")
            print("Compiling list of all files in SDSS DR12 frames directory.")
            print("")
            
            root_dir = '/data/fourier/leroy.42/allsky/sdss_dr12/frames/'
            selection = '*/*/*/*.fits.bz2'
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
            compile_list_of_images(
                root_dir=root_dir, selection=selection,
                tab_dir=tab_dir, tab_file=tab_file)

            root_dir = '../../orig_data/unwise_v2/v2_largegals/'
            selection = '*/*/unwise-*-w*-img-m.fits'
            tab_dir = '../../working_data/unwise/index/'
            tab_file = 'unwise_custom_largeleda_list.fits'
            print("... largeleda")
            compile_list_of_images(
                root_dir=root_dir, selection=selection,
                tab_dir=tab_dir, tab_file=tab_file)            

            root_dir = '../../orig_data/unwise_v2/v2_smallgals/'
            selection = '*/*/unwise-*-w*-img-m.fits'
            tab_dir = '../../working_data/unwise/index/'
            tab_file = 'unwise_custom_smallleda_list.fits'
            print("... smallleda")
            compile_list_of_images(
                root_dir=root_dir, selection=selection,
                tab_dir=tab_dir, tab_file=tab_file)            

            root_dir = '../../orig_data/unwise_v2/v2_unmanga/'
            selection = '*/*/unwise-*-w*-img-m.fits'
            tab_dir = '../../working_data/unwise/index/'
            tab_file = 'unwise_custom_manga_list.fits'
            print("... manga")
            compile_list_of_images(
                root_dir=root_dir, selection=selection,
                tab_dir=tab_dir, tab_file=tab_file)            

            root_dir = '../../orig_data/unwise_v2/v2_localgroup/'
            selection = '*/*/unwise-*-w*-img-m.fits'
            tab_dir = '../../working_data/unwise/index/'
            tab_file = 'unwise_custom_localgroup_list.fits'
            print("... localgroup")
            compile_list_of_images(
                root_dir=root_dir, selection=selection,
                tab_dir=tab_dir, tab_file=tab_file)            
            
            return()
            
        if survey == 'unwise_allwise':
            print("")
            print("Compiling list of all files in UNWISE allwise directory.")
            print("")

            # TBD
            
        if survey == 'unwise_neowise':
            print("")
            print("Compiling list of all files in NEOWISE neowise directory.")
            print("")

            # TBD

        if survey == 'gaia':
            print("")
            print("Compiling list of all files in GAIA directory.")
            print("")

            # TBD
            
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
# Loop over the files
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def index_image_list(
        survey = None,
        list_table = None,
        outfile = None,
        do_all_frames=False):
    """
    Read a list of FITS files and build these into an index.
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

        if survey == 'unwise_custom' and (list_table == None):
            print("")
            print("Indexing custom UNWISE tiles")
            print("(this is a recursive call)")
            print("")
            index_dir = '../../working_data/galex/index/'
            for subsample in ['localvolume','localgroup',
                              'largeleda', 'smallleda', 'manga']:
                list_table = index_dir + 'unwise_custom_' + \
                    subsample+'_list.fits'
                outfile = index_dir+'unwise_custom_'+subsample+'_index.fits'
                index_image_list(
                    survey='unwise_custom',
                    list_table = list_table, outfile=outfile,
                    start_id = start_id, stop_id = stop_id,
                    do_all = do_all)
                return()

        if survey == 'unwise_allwise':
            print("")
            print("Indexing UNWISE allwise tiles.")
            print("")

            pass

        if survey == 'unwise_neowise':
            print("")
            print("Indexing UNWISE neowise tiles.")
            print("")

            pass

        if survey == 'gaia':
            print("")
            print("Indexing GAIA tables.")
            print("")

            pass

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
    if this_field in ['pix_scale_x','pix_scale_y']:
        tab[this_field] = np.nan

    tab['filter'] = ''
        
    if survey == 'sdss':
        tab['exptime'] = np.nan

    if survey.count('unwise') > 0:
        pass
        
    if survey == 'galex':
        tab['rrhr_fname'] = ''
        tab['flag_fname'] = ''
        tab['bgsub_fname'] = ''
        
    counter = 0

    for this_row in ProgressBar(tab):

        # Read, extract header and wcs from file
        this_fname = this_file['fname']
        this_hdulist = fits.open(this_fname)
        this_header = this_hdulist[0].header
        w = wcs.WCS(this_header)

        # Survey specific items
        if survey == 'sdss':            
            this_row['filter'] = this_header['FILTER'].strip()
            this_row['exptime'] = float(this_header['EXPTIME'])

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
        tab['nx'] = this_header['NAXIS1']
        tab['ny'] = this_header['NAXIS2']
            
        # Convert to world coordinates
        ctr_pix = [[(int(nx/2),int(ny/2))]]        
        ctr_world =  w.wcs_pix2world(ctr_pix,0)
        this_row['ctr_ra'] = ctr_world[0]
        this_row['ctr_dec'] = ctr_world[1]

        # Pixel scale
        pix_scale = wcs.utils.proj_plane_pixel_scales(w)
        this_row['pix_scale_x'] = pix_scale[0]
        this_row['pix_scale_y'] = pix_scale[1]        
        
    tab.write(outfile, format='fits', overwrite=True)

    stop_time = time.time()

    print("This took : ", stop_time-start_time, " seconds.")
    
    return(filled_tab)
        
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Run the compilation
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

do_fetch = False
do_check = False
do_flist = True
do_index = False

do_unwise = False
do_sdss = False
do_galex = True
do_gaia = False

if do_fetch:
    if do_unwise:
        unwise_fetch(dry_run=True, incremental=False, version='neo9')
        unwise_fetch(dry_run=True, incremental=False, version='allwise')

    if do_galex:
        galex_fetch(dry_run=True, incremental=False)

    if do_sdss:
        sdss_fetch(dry_run=True)

    if do_gaia:
        gaia_fetch(dry_run=True)
        
if do_check:
    pass
        
if index:
    if do_sdss:
        if do_flist:
            test = compile_list_of_images(survey='sdss')
        if do_index:
            test = index_image_list(survey='sdss')

    if do_galex:
        if do_flist:
            test = compile_list_of_images(survey='galex')
        if do_index:
            test = index_image_list(survey='galex')

    if do_unwise:
        if do_flist:
            test = compile_list_of_images(survey='unwise_custom')
            #test = compile_list_of_images(survey='unwise_allwise')
            #test = compile_list_of_images(survey='unwise_neowise')
        if do_index:
            test = index_image_list(survey='unwise_custom')
            #test = index_image_list(survey='unwise_allwise')
            #test = index_image_list(survey='unwise_neowise')

    if do_gaia:
        if do_flist:
            test = compile_list_of_images(survey='gaia')
        if do_index:
            test = index_image_list(survey='gaia')
            
#test = compile_list_of_images(survey='unwise_allwise')
#test = compile_list_of_images(survey='unwise_neowise')
#test = compile_list_of_images(survey='gaia')

#test = index_sdss_frames()
#test = index_sdss_frames(do_all=True, identifier='')
