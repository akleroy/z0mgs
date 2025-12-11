# GALEX atlas creation

# What's in here:

# - infrastructure to loop over subsamples
# - a pipeline function to reduce one galaxy
# - galex-specific functions

# TBD:

# - image construction itself seems to be the bottleneck

# - galaxy mask creation can be sped by doing it only once.

# - can save reading the tiles index over and over but this is not
# - likely to be a large speedup.

# - do better on incremental stuff

import os

from astropy.table import Table
from astropy.io import fits
from astropy import wcs
from astropy.wcs.utils import proj_plane_pixel_scales

from reproject import reproject_interp, reproject_adaptive

from utils_tabs_and_dirs import *
from utils_z0mgs_images import *

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Atlas construction loop
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def galex_build_atlas(
        tasks='all',
        subsamples=None,
        just_galaxy=None,
        skip_galaxy=None,
        start_galaxy=None,
        stop_galaxy=None,
        root_dir='../../working_data/galex/',
        table_dir='../../measurements/',
        bands=['fuv','nuv'],
        pause=False,
        incremental=False,
        show = False,
        overwrite = True):
    """Loop to construct the GALEX atlas. Manages construction of the
    list of targets and juggling directories then calls the pipeline
    for a single galaxy.
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Define tasks
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if not isinstance(tasks, list):
        tasks = [tasks]
    
    working_dirs = {
        'staged':root_dir+'staged/',
        'bkgrd':root_dir+'bkgrd/',        
        'masks':root_dir+'masks/',
        'coords':root_dir+'coords/',        
        'convolved':root_dir+'convolved/',
        'final':root_dir+'final/',
        'gaia':root_dir+'../gaia/',
    }

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Define targets
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
        
    target_table = build_target_table(
        subsamples=subsamples, table_dir=table_dir,
        just_galaxy=just_galaxy, skip_galaxy=skip_galaxy,
        start_galaxy=start_galaxy, stop_galaxy=stop_galaxy)
    n_targets = len(target_table)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=            
    # Define Gaussian resolutions
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=            

    res_dict = {
        'gauss7p5':7.5,
        'gauss15':15.,
        'gauss20':20.,
        }

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Loop over targets
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        
    for ii, this_target_row in enumerate(target_table):

        print("")
        print("Processing galaxy "+str(ii)+" of "+str(n_targets)+ \
              " "+this_target_row["Z0MGS_NAME"])
        print("")

        this_working_dirs = working_dirs.copy()
        if this_target_row['SUBSAMPLE'] is not None:
            for this_key in this_working_dirs.keys():
                this_dir = this_working_dirs[this_key]
                this_working_dirs[this_key] = this_dir + \
                    this_target_row['SUBSAMPLE']+'/'        

        galex_process_one_galaxy(
            target = this_target_row,
            working_dirs = this_working_dirs,
            res_dict=res_dict,
            tasks = tasks,
            bands = bands,
            pause = pause,
            incremental=incremental,
            show = show,
            overwrite = overwrite)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Single galaxy pipeline
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        
def galex_process_one_galaxy(
        target=None,
        tasks=['all'],
        bands=['fuv','nuv'],
        working_dirs='./',
        res_dict={},
        pause=False,
        incremental=False,
        show = False,
        overwrite = True):
    """
    Build the GALEX z0MGS products for one galaxy.
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Manage the task list and directories
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 

    if not isinstance(tasks, list):
        tasks = [tasks]
    
    if tasks == ['all']:
        tasks = [ 
            'stage',
            'plot_stage',
            'convolve',
            'plot_convolve',
            'coord_mask',
            'plot_coord_mask',             
            'galaxy_mask',
            'plot_galaxy_mask',
            'star_pred',
            'star_mask',
            'plot_star_mask',
            'bkgrd',
            'plot_bkgrd',
        ]

    if not isinstance(working_dirs, dict):
        working_dirs = {
            'staged':working_dirs,
            'bkgrd':working_dirs,
            'masks':working_dirs,
            'coords':working_dirs,            
            'convolved':working_dirs,
            'final':working_dirs,
            'gaia':working_dirs,
        }

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Manage the target
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    if not isinstance(target, dict):
        target = dict(target)
    
    this_name = target['Z0MGS_NAME'].strip()

    # TBD - revise this infrastructure to a set of field tags
    # consistent with the tile index syntax.
                
    ra_ctr = target['RA_CTR']
    dec_ctr = target['DEC_CTR']
    size_deg = target['TRC_DEC'] - target['BLC_DEC']
    
    if target['PGC'] == 2557:
        print("... special case of M31")
        size_deg = 2.0

    posang_deg = target[]
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Stage the images
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    # Identify which tiles in the full stack of GALEX tiles overlap
    # this target and use these to construct an intensity and a weight
    # image for the galaxy.

    # TBD - This step is slow. Explore options. Including possibly
    # unzipping GALEX?
    
    if 'stage' in tasks:
        
        for this_band in bands:

            outfile_image = working_dirs['staged']+ \
                this_name+'_'+this_band+'_mjysr.fits'

            outfile_weight = working_dirs['staged']+ \
                this_name+'_'+this_band+'_weight.fits'

            skip = False
            if incremental:
                image_present = os.system.isfile(outfile_image)
                weight_present = os.system.isfile(outfile_weight)
                skip = image_present * weight_present

            print("Staging: ", this_name)
            print("... RA, Dec center: ", ra_ctr, dec_ctr)
            print("... size in deg: ", size_deg)
            print("... to: ", outfile_image)

            index_file = \
                '../../working_data/galex/index/galex_tile_index.fits'
                
            if not skip:
               image, hdr = \
                   extract_galex_stamp(
                       band=this_band,
                       ra_ctr=ra_ctr,
                       dec_ctr=dec_ctr,
                       size_deg=size_deg,
                       index_file = index_file,
                       use_int_files = True,
                       outfile_image = outfile_image,
                       outfile_weight = outfile_weight,
                       show = show,
                       pause = pause,
                       overwrite = True)

    # Plot the images produced in staging. We will make fancier plots
    # later but this provides a quick look.
               
    if 'plot_stage' in tasks:

        for this_band in bands:
              
            staged_image_file = working_dirs['staged']+ \
                this_name+'_'+this_band+'_mjysr.fits'
            
            staged_weight_file = working_dirs['staged']+ \
                this_name+'_'+this_band+'_weight.fits'

            if this_band == 'fuv':
                this_rms = 10.**(-3.5)

            if this_band == 'nuv':
                this_rms = 10.**(-3.0)

            print("Plots for: ", this_name, ' ', this_band)
            
            show_z0mgs_image(
                image_fname = staged_image_file,
                show = False,
                outfile = staged_image_file.replace('.fits','.png'),
                title = this_name+' '+this_band,
                value_string = this_band+' [MJy/sr]',
                rms = this_rms
            )

            show_z0mgs_image(
                image_fname = staged_weight_file,
                show = False,
                outfile = staged_weight_file.replace('.fits','.png'),
                title = this_name+' '+this_band,
                value_string = this_band+' [relative response]',
                rms = this_rms
            )


    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
    # Make supporting coordinate images
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    # Create images of galactocentric radius and offset from the major
    # and minor axes.
    
    if 'coord_mask' in tasks:

        skip = False
            
        staged_image_file = working_dirs['staged']+ \
            this_name+'_nuv_mjysr.fits'

        if os.path.isfile(staged_image_file) == False:
            staged_image_file = working_dirs['staged']+ \
                this_name+'_fuv_mjysr.fits'
            if os.path.isfile(staged_image_file) == False:
                skip = True

        rad_file = working_dirs['coords']+ \
            this_name+'_rgal.fits'

        posang_file = working_dirs['coords']+ \
            this_name+'_posang.fits'

        major_file = working_dirs['coords']+ \
            this_name+'_major.fits'

        minor_file = working_dirs['coords']+ \
            this_name+'_minor.fits'        

        if incremental and (skip == False):

            skip = True
            if not os.path.isfile(rad_file):
                skip = False
            if not os.path.isfile(posang_file):
                skip = False
            if not os.path.isfile(major_file):
                skip = False
            if not os.path.isfile(minor_file):
                skip = False
                
        if not skip:

            template_hdu = fits.open(staged_image_file)[0]
            template_hdr = template_hdu.header
            template_hdr['BUNIT'] = 'Deg'
            
            radius_deg, projang_deg, major_deg, minor_deg = \
                deproject(
                    center_coord=None,
                    incl=0*u.deg,
                    pa=0*u.deg,
                    template_header=template_hdr,
                    return_offset=True,
                    verbose=False)

            rad_hdu = PrimaryHDU(radius_deg, 
        
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Do the convolutions
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    # Convolve to some specified Gaussian resolutions.

    # TBD - In principle this could come after star masking and could
    # blank them out of the convolution. This needs to be done before
    # background subtraction. Add ability to accept a mask.
    
    if 'convolve' in tasks:

        print("Convolving images ...")

        for this_band in bands:

            print("... ... ", this_band)
            
            for target_res in res_dict.keys():
                print("... to "+target_res)
                
                kern_to_this_res = z0mgs_kernel_name(
                    from_res=this_band, to_res=target_res)
                
                staged_image_file = working_dirs['staged']+ \
                    this_name+'_'+this_band+'_mjysr.fits'
                
                convolved_image_file = working_dirs['convolved']+ \
                    this_name+'_'+this_band+'_mjysr_'+target_res+'.fits'

                # TBD - Add masks.
                
                convolve_image_with_kernel(
                    image_file=staged_image_file,
                    outfile=convolved_image_file,
                    kernel_file=kern_to_this_res,
                    blank_zeros=False,
                    overwrite=True
                )

    # Make images of the convolved data.
                
    if 'plot_convolve' in tasks:

        for this_band in bands:

            print("... ... ", this_band)
            
            for this_res in res_dict.keys():
                print("... res "+this_res)
                                
                convolved_image_file = working_dirs['convolved']+ \
                    this_name+'_'+this_band+'_mjysr_'+this_res+'.fits'

                show_z0mgs_image(
                    image_fname = convolved_image_file,
                    show = False,
                    outfile = convolved_image_file.replace('.fits','.png'),
                    title = this_name+' '+this_band+' '+this_res,
                    value_string = this_band+' [MJy/sr]',
                )

    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Make galaxy masks
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # Make a mask of known other galaxies in the field. Leverages
    # existing galaxy database. Could farm these tuning parameters out
    # to a config file. My current impression is that they are a bit
    # too small for some of the bigger galaxies.
    
    if 'galaxy_mask' in tasks:

        gal_tab_file = '/home/leroy.42/idl/galbase/gal_data/gal_base.fits'
        pgc_to_skip = []

        # Find a template image

        skip = False
            
        staged_image_file = working_dirs['staged']+ \
            this_name+'_nuv_mjysr.fits'

        if os.path.isfile(staged_image_file) == False:
            staged_image_file = working_dirs['staged']+ \
                this_name+'_fuv_mjysr.fits'
            if os.path.isfile(staged_image_file) == False:
                skip = True
        
        galaxy_mask_file = working_dirs['masks']+ \
            this_name+'_galmask.fits'

        if incremental:
            mask_present = os.system.isfile(galaxy_mask_file)
            skip = mask_present

        # Build the galaxy mask
        if not skip:

            print("Making a galaxy mask for for: ", this_name)
            
            build_galaxy_mask(
                template_file = staged_image_file,
                outfile = galaxy_mask_file,
                gal_tab_file = gal_tab_file,
                this_pgc = target['PGC'],
                pgc_to_skip = pgc_to_skip,
                field_for_rad = 'R25_DEG',
                min_rad = 7.5/3600.*u.deg,
                rad_fac_to_blank = 1.0,
                use_orient = True,
                max_incl = 70.*u.deg,
                show = show,
                pause = pause,
                overwrite = overwrite)

            
        if 'plot_galaxy_mask' in tasks:

            galaxy_mask_image = \
                galaxy_mask_file.replace('.fits','.png')
            
            if incremental:
                image_present = os.system.isfile(galaxy_mask_image)
                skip = image_present

            this_rms = 10.**(-3.5)
            if not skip:                
                show_z0mgs_image(
                    image_fname = staged_image_file,
                    mask_fname = galaxy_mask_file,
                    mask_levels = [0.99],
                    show = False,
                    outfile = galaxy_mask_image,
                    title = this_name+' with galaxy mask',
                    value_string = this_band+' [MJy/sr]',
                    rms = this_rms
                )
                
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Make star predictions
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    # Make a prediction in flux/pixel units of the flux of foreground
    # star point sources in the image and convolve this with a few
    # relevant PSFs to make images. Leverages existing GAIA query or
    # can make the query if the file is missing.

    # TBD Order here is debatable relative to convolution especially.
    
    if 'star_pred' in tasks:

        print("Predicting stellar fluxes.")
        skip_query_if_present = True

        ks_file = '../../measurements/tab_2mass_stars.fits'
            
        for this_band in bands:

            staged_image_file = working_dirs['staged']+ \
                this_name+'_'+this_band+'_mjysr.fits'
            
            gaia_file = working_dirs['gaia']+ \
                this_name+'_gaia_dr3.fits'

            star_flux_file = working_dirs['masks']+\
                this_name+'_'+this_band+'_starflux.fits'

            star_intens_file = working_dirs['masks']+\
                this_name+'_'+this_band+'_starintens.fits'
            
            # Query GAIA if needed but skip if present            

            # TBD patched
            gaia_file = working_dirs['gaia'].replace('test_data','working_data') + \
                this_name+'_gaia_dr3.fits'

            print("... ensuring Gaia catalog present")

            # TBD - replace this with a disk based catalog query and
            # construction. Reading those tables is a time hit, so it
            # may be worth avoiding doing this per band.
            
            #query_gaia(
            #    ra_min = target['BLC_RA'],
            #    ra_max = target['TRC_RA'],
            #    dec_min = target['RLC_RA'],
            #    dec_max = target['TRC_RA'],
            #    outfile=gaia_file,
            #    skip_if_present=skip_query_if_present)

            print("... building a stellar flux image")
                        
            # Build an image of stellar flux
            build_star_flux_image(
                template_file = staged_image_file,
                outfile = star_flux_file,
                band = this_band,
                gaia_file = gaia_file,
                ks_file = ks_file,
                gaia_s2n_cut = 3.5,
                center_coord = None,
                center_tol = 3.0*u.arcsec,
            )

            print("... convolving the stellar flux image to intensity")
            
            star_intens_hdu = fits.HDUList(
                [jypix_to_mjysr(image_fname=star_flux_file)])

            pix_to_native_kernel = \
                z0mgs_psf_name(band=this_band)

            print("... ... to native resolution")
            
            convolve_image_with_kernel(                
                image_hdu=star_intens_hdu,
                outfile=star_intens_file,
                kernel_file=pix_to_native_kernel,
                blank_zeros=False,
            )

            for res_tag, res_in_arcsec in res_dict.items():
                print("... ... to resolution ", res_tag)                
                convolve_image_with_gauss(
                    image_hdu=star_intens_hdu,
                    starting_res = 0. * u.arcsec,
                    target_res = res_in_arcsec * u.arcsec,
                    outfile = star_intens_file.replace(
                        '.fits','_'+res_tag+'.fits'),
                    overwrite=True,
                )
        
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Make star masks
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%        

    # Use the existing stellar prediction images to create masks that
    # blank pixels associated with foreground stars.

    # TBD this requires knowledge of a specific threshold. For GALEX
    # this works best with a convolved image, which makes it desirable
    # to add a step above aggregative statistics from the convolved
    # maps. Here I have hard-coded things.
    
    if 'star_mask' in tasks:

        for this_band in bands:
            
            staged_image_file = working_dirs['staged']+ \
                this_name+'_'+this_band+'_mjysr.fits'
                
            star_intens_file = working_dirs['masks']+\
                this_name+'_'+this_band+'_starintens.fits'
                
            star_mask_file = working_dirs['masks']+\
                this_name+'_'+this_band+'_mask.fits'

            if this_band == 'fuv':
                this_rms = 10.**(-3.5)
            
            if this_band == 'nuv':
                this_rms = 10.**(-3.0)
                
            build_star_mask(
                image_file = staged_image_file,
                star_file = star_intens_file,
                outfile = star_mask_file,
                rms_value = this_rms,
                rms_fac = 3.0,
                overwrite = True,
            )
        
    if 'plot_star_mask' in tasks:

        for this_band in bands:
            
            star_intens_file = working_dirs['masks']+\
                this_name+'_'+this_band+'_starintens.fits'
        
            # TBD

   
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Fit and subtract a background
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    
    if 'bkgrd' in tasks:

        fit_z0mgs_background(
            image_fname=None,
            image_hdu=None,
            mask_fnames=[],
            mask_hdus=[],
            weight_fname=None,
            weight_hdu=None,
            rad_fname = None,
            rad_hdu = None,
            fid_rad = 15./3600.,
            outfile = None,
            methods=['itermed'],
            clip_thresh=3.0,
            niter=5,
        )        
        
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    # Extract a subimage
    # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        
    if 'extract' in tasks:
        
        pass

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Individual pipeline steps
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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
        pause = False,
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
