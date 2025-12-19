# Loop for atlas creation

# What's in here:

# - infrastructure to loop over subsamples
# - a pipeline function to reduce one galaxy

# TBD:

# - image construction itself seems to be the bottleneck

# - can save reading the tiles index over and over but this is not
# - likely to be a large speedup.

# - do better on incremental stuff

# - the order of ops for correctness could be improved to:

# stage ->
# convolve + mask ->
# bkfit ->
# mask + bksub (+ interpolate?) ->
# convolve

# right now convolving with filled 0s and then subtracting leaves edge
# artifact where we fill 0s outside the image. Could also patch this
# (imperfectly) by filling with the image median.

import os

from astropy.table import Table
from astropy.io import fits
from astropy import wcs
from astropy.wcs.utils import proj_plane_pixel_scales

from reproject import reproject_interp, reproject_adaptive

from utils_tabs_and_dirs import *
from utils_z0mgs_images import *
from utils_cutouts import *

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Atlas construction loop
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def z0mgs_build_atlas(
        survey='galex',
        tasks='all',
        subsamples=None,
        just_galaxy=None,
        skip_galaxy=None,
        start_galaxy=None,
        stop_galaxy=None,
        root_dir='../../working_data/',
        table_dir='../../measurements/',
        bands=None,
        unwise_method='reproject',
        incremental=False,
        overwrite = True):
    """Loop to construct one of the z0mgs atlases. Manages construction of
    the list of targets and juggling directories then calls the
    pipeline for a single galaxy.
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Define tasks and directories
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    if isinstance(survey,list):
        survey = survey[0]
    valid_surveys = ['galex','unwise','unwise_custom','allwise','neowise','sdss']
    if survey not in valid_surveys:
        print("Invalid survey: ", survey)
        print("... valid options: ", valid_surveys)
        return

    if survey == 'galex':
        valid_bands = ['fuv','nuv']
        if bands is None:
            band = valid_bands
        survey_dir = 'galex/'

    if (survey == 'unwise') | (survey == 'unwise_custom') \
       | (survey=='allwise'):
        valid_bands = ['w1','w2','w3','w4']
        if bands is None:
            band = valid_bands
        survey_dir = 'unwise/'

    if survey == 'neowise':
        valid_bands = ['w1','w2','w3','w4']
        if bands is None:
            band = valid_bands
        survey_dir = 'unwise/'
            
    if survey == 'sdss':
        valid_bands = ['u','g','r','i','z']
        if bands is None:
            band = valid_bands
        survey_dir = 'sdss/'
    
    if not isinstance(tasks, list):
        tasks = [tasks]
    
    working_dirs = {
        'staged':root_dir+survey_dir+'staged/',
        'bkgrd':root_dir+survey_dir+'bkgrd/',        
        'masks':root_dir+survey_dir+'masks/',
        'coords':root_dir+survey_dir+'coords/',        
        'convolved':root_dir+survey_dir+'convolved/',
        'final':root_dir+survey_dir+'final/',
        'gaia':root_dir+'gaia/',
    }

    # Make sure the directories are present
    make_z0mgs_directories(
        root_dir=root_dir, surveys=[survey])    

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Define targets
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    # Get the list of galaxies. Does not depend on the parent survey.
    
    target_table = build_target_table(
        subsamples=subsamples, table_dir=table_dir,
        just_galaxy=just_galaxy, skip_galaxy=skip_galaxy,
        start_galaxy=start_galaxy, stop_galaxy=stop_galaxy)
    n_targets = len(target_table)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=            
    # Define target Gaussian resolutions
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
              " "+this_target_row["NAME"])
        print("")

        this_working_dirs = working_dirs.copy()
        if this_target_row['SUBSAMPLE'] is not None:
            for this_key in this_working_dirs.keys():
                this_dir = this_working_dirs[this_key]
                this_working_dirs[this_key] = this_dir + \
                    this_target_row['SUBSAMPLE']+'/'        

        # Make sure that we have a position, orientation, etc.
        target_row = clean_up_target_row(this_target_row)
                
        z0mgs_process_one_galaxy(
            target = this_target_row,
            tasks = tasks,
            survey = survey,
            bands = bands,
            working_dirs = this_working_dirs,
            res_dict=res_dict,
            unwise_method=unwise_method,
            incremental=incremental,
            overwrite = overwrite)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Single galaxy pipeline
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
        
def z0mgs_process_one_galaxy(
        target=None,
        tasks=['all'],
        survey='galex',
        bands=None,
        working_dirs='./',
        res_dict= {
            'gauss7p5':7.5,
            'gauss15':15.,
            'gauss20':20.,
        },
        unwise_method='reproject',
        incremental=False,
        overwrite = True):
    """
    Build the GALEX z0MGS products for one galaxy.
    """

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Manage the survey, bands, tasks, and directories
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 

    # Check survey validity
    valid_surveys = ['galex','unwise','allwise','neowise','sdss']
    if survey not in valid_surveys:
        print("Invalid survey: ", survey)
        print("... valid options: ", valid_surveys)
        return

    # Key data for survey

    an_unwise_survey = \
        ((survey == 'unwise') | (survey == 'neowise') \
         | (survey == 'allwise'))
    
    # ... by default fit the background to the native res native res
    res_ext_for_bkgrd = ''
    preconv_res_dict = {}
    
    if survey == 'galex':        
        valid_bands = ['fuv','nuv']
        fid_rms = 5E-3
        # ... for galex smooth before background
        res_ext_for_bkgrd = '_gauss20'
        preconv_res_dict = {
            'gauss20':20.,
        }
    if survey == 'sdss':
        valid_bands = ['u','g','r','i','z']
        fid_rms = 1E-3
    if survey == 'unwise':
        valid_bands = ['w1','w2','w3','w4']
        fid_rms = 1E-2
    if survey == 'neowise':
        valid_bands = ['w1','w2']
        fid_rms = 1E-2
    if survey == 'allwise':
        valid_bands = ['w1','w2','w3','w4']
        fid_rms = 1E-2
        
    if bands is None:
        bands = valid_bands
    else:
        for this_band in bands:
            if this_band not in valid_bands:
                print("Invalid band for ", survey, ": ", this_band)
                print("... valid bands are ", valid_bands)
                
    # Check task list
    if not isinstance(tasks, list):
        tasks = [tasks]

    valid_tasks = [ 
        'stage',
        'plot_stage',
        'convolve_for_bkgrd',
        'coord_mask',
        'galaxy_mask',
        'star_pred',
        'star_mask',
        'bkgrd',
        'plot_bkgrd',
        'convolve',
        'plot_results',
    ]
        
    if tasks == ['all']:
        tasks = valid_tasks

    use_tasks = []
    for this_task in tasks:
        if this_task not in valid_tasks:
            print("Invalid task ", this_task)
            print("Returning ...")
            return()
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Manage the target
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    if not isinstance(target, dict):
        target = dict(target)

    # Extract values from the target table
    this_name = target['NAME'].strip()
    this_pgc = target['PGC']                
    ctr_ra = target['CTR_RA']
    ctr_dec = target['CTR_DEC']
    center_coord = SkyCoord(ra=ctr_ra*u.deg, dec=ctr_dec*u.deg, frame='icrs')
    pa = target['POSANG_DEG']*u.deg
    incl = target['INCL_DEG']*u.deg
    rgal = target['RGAL_DEG']*u.deg
    size_deg = target['IMSIZE']

    # TBD - overrides file
    if target['PGC'] == 2557:
        print("... special case of M31")
        size_deg = 2.0

    print("Processing: ", this_name)
    print("... R.A., Dec. center: ", ctr_ra, ctr_dec)
    print("... size in deg: ", size_deg)
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Stage the images
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    # Call the cutout construction routine. This varies by survey.

    # TBD - This step is slow for galex. Explore options.
    
    if 'stage' in tasks:        
        
        for this_band in bands:
            
            outfile_image = working_dirs['staged']+ \
                this_name+'_'+this_band+'_mjysr.fits'

            skip = False
            if incremental:
                image_present = os.system.isfile(outfile_image)
                skip = image_present
                if survey == 'galex':
                    weight_present = os.system.isfile(outfile_weight)
                    skip = skip * weight_present

            if survey == 'galex' and not skip:

                outfile_weight = working_dirs['staged']+ \
                    this_name+'_'+this_band+'_weight.fits'
                
                extract_galex_stamp(
                    band=this_band,
                    ctr_ra=ctr_ra,
                    ctr_dec=ctr_dec,
                    size_deg=size_deg,
                    use_int_files = True,
                    outfile_image = outfile_image,
                    outfile_weight = outfile_weight,
                    overwrite = True)

            if an_unwise_survey and not skip:

                outfile_mask = working_dirs['staged']+ \
                    this_name+'_'+this_band+'_mask.fits'
                
                extract_unwise_stamp(
                    band=this_band,
                    ctr_ra=ctr_ra,
                    ctr_dec=ctr_dec,
                    size_deg=size_deg,
                    outfile_image = outfile_image,
                    outfile_mask = outfile_mask,
                    method=unwise_method,
                    survey=survey,
                    overwrite = True)

            if survey == 'sdss' and not skip:

                print("SDSS staging not implemented yet.")
                
                pass            

    # Plot the images produced in staging. We will make fancier plots
    # later but this provides a quick look.
               
    if 'plot_stage' in tasks:

        for this_band in bands:
              
            staged_image_file = working_dirs['staged']+ \
                this_name+'_'+this_band+'_mjysr.fits'

            print("Plots for: ", this_name, ' ', this_band)
            
            show_z0mgs_image(
                image_fname = staged_image_file,
                show = False,
                outfile = staged_image_file.replace('.fits','.png'),
                title = this_name+' '+this_band,
                value_string = this_band+' [MJy/sr]',
                rms = fid_rms,
            )

            if survey == 'galex':
                
                staged_weight_file = working_dirs['staged']+ \
                    this_name+'_'+this_band+'_weight.fits'
            
                show_z0mgs_image(
                    image_fname = staged_weight_file,
                    outfile = staged_weight_file.replace('.fits','.png'),
                    title = this_name+' '+this_band,
                    value_string = this_band+' [relative response]',
                    rms = fid_rms
                )

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Convolve (if needed) before background subtraction
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 

    # Convolve to a specified Gaussian resolution (if relevant) before
    # background subtraction. The final convolutions will come after
    # background subtraction below.
    
    if 'convolve' in tasks:

        print("Convolving images for background subtraction ...")

        for this_band in bands:

            print("... ... ", this_band)
            
            for target_res in preconv_res_dict.keys():
                print("... to "+target_res)
                
                kern_to_this_res = z0mgs_kernel_name(
                    from_res=this_band, to_res=target_res)
                if os.path.isfile(kern_to_this_res) == False:
                    print("... ... kernel not found: ", kern_to_this_res)
                    print("... ... continuing without it.")
                    print("... ... create it if it should be there.")
                    continue
                
                staged_image_file = working_dirs['staged']+ \
                    this_name+'_'+this_band+'_mjysr.fits'
                
                convolved_image_file = working_dirs['staged']+ \
                    this_name+'_'+this_band+'_mjysr_'+target_res+'.fits'
                
                convolve_image_with_kernel(
                    image_file=staged_image_file,
                    outfile=convolved_image_file,
                    kernel_file=kern_to_this_res,
                    blank_zeros=False,
                    overwrite=True
                )
                
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
    # Make supporting coordinate images
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    # Create images of galactocentric radius and offset from the major
    # and minor axes.
    
    if 'coord_mask' in tasks:

        print("Constructing coordinate images ... ")
        
        # Find a template
        skip = True
        template_file = None
        for this_band in bands:
            if template_file is not None:
                continue
            try_file = working_dirs['staged']+ \
                this_name+'_'+this_band+'_mjysr.fits'
            if os.path.isfile(try_file):
                template_file = try_file
                skip = False

        # Define coord files
        rad_file = working_dirs['coords']+ \
            this_name+'_rgal.fits'

        theta_file = working_dirs['coords']+ \
            this_name+'_theta.fits'

        major_file = working_dirs['coords']+ \
            this_name+'_major.fits'

        minor_file = working_dirs['coords']+ \
            this_name+'_minor.fits'        

        # Proceed if any files are missing
        if incremental & (skip == False):
            skip = True
            for this_file in [rad_file, posang_file,
                              major_file, minor_file]:
                if not os.path.isfile(this_file):
                    skip = False
                
        if not skip:

            template_hdu = fits.open(template_file)[0]
            template_hdr = template_hdu.header
            template_hdr['BUNIT'] = 'Deg'

            radius_deg, theta_deg, major_deg, minor_deg = \
                deproject(
                    center_coord = center_coord, incl=incl,  pa=pa,
                    template_header=template_hdr, return_offset=True)

            rad_hdu = fits.PrimaryHDU(data=radius_deg, header=template_hdr)
            rad_hdu.writeto(rad_file, overwrite=True)
            
            theta_hdu = fits.PrimaryHDU(data=theta_deg, header=template_hdr)
            theta_hdu.writeto(theta_file, overwrite=True)

            major_hdu = fits.PrimaryHDU(data=major_deg, header=template_hdr)
            major_hdu.writeto(major_file, overwrite=True)
            
            minor_hdu = fits.PrimaryHDU(data=minor_deg, header=template_hdr)
            minor_hdu.writeto(minor_file, overwrite=True)        

        else:

            print("... ... skipping")

            
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Make galaxy masks
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    # Make a mask of known other galaxies in the field. Leverages
    # existing galaxy database. Could farm these tuning parameters out
    # to a config file. My current impression is that they are a bit
    # too small for some of the bigger galaxies.
    
    if 'galaxy_mask' in tasks:

        gal_tab_file = '/home/leroy.42/idl/galbase/gal_data/gal_base.fits'
        pgc_to_skip = []

        # Find a template
        skip = True
        template_file = None
        for this_band in bands:
            if template_file is not None:
                continue
            try_file = working_dirs['staged']+ \
                this_name+'_'+this_band+'_mjysr.fits'
            if os.path.isfile(try_file):
                template_file = try_file
                skip = False
        
        galaxy_mask_file = working_dirs['masks']+ \
            this_name+'_galmask.fits'

        # Skip if present if incremental mode requested
        if incremental & (skip == False):
            mask_present = os.system.isfile(galaxy_mask_file)
            skip = mask_present

        # Build the galaxy mask
        if not skip:

            print("Making a galaxy mask for for: ", this_name)
            
            build_galaxy_mask(
                template_file = template_file,
                outfile = galaxy_mask_file,
                gal_tab_file = gal_tab_file,
                this_pgc = target['PGC'],
                pgc_to_skip = pgc_to_skip,
                field_for_rad = 'R25_DEG',
                min_rad = 7.5/3600.*u.deg,
                rad_fac_to_blank = 1.0,
                use_orient = True,
                max_incl = 70.*u.deg,
                overwrite = overwrite)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Make star predictions    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
    # Make a prediction in flux/pixel units of the flux of foreground
    # star point sources in the image and convolve this with a few
    # relevant PSFs to make images. Leverages existing GAIA query or
    # can make the query if the file is missing.
    
    if 'star_pred' in tasks:

        print("Predicting stellar fluxes.")
        skip_query_if_present = True

        ks_file = '../../measurements/tab_2mass_stars.fits'
            
        for this_band in bands:

            template_file = working_dirs['staged']+ \
                this_name+'_'+this_band+'_mjysr.fits'
            
            gaia_file = working_dirs['gaia']+ \
                this_name+'_gaia_dr3.fits'

            star_flux_file = working_dirs['masks']+\
                this_name+'_'+this_band+'_starflux.fits'

            star_intens_file = working_dirs['masks']+\
                this_name+'_'+this_band+'_starintens.fits'
            
            # Query GAIA if needed but skip if present            

            # TBD patched
            gaia_file = working_dirs['gaia'] + this_name + '_gaia_dr3.fits'
            if gaia_file.count('test_data') > 0:
                gaia_file = gaia_file.replace('test_data','working_data') 

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
                template_file = template_file,
                outfile = star_flux_file,
                band = this_band,
                gaia_file = gaia_file,
                ks_file = ks_file,
                gaia_s2n_cut = 3.5,
                center_coord = center_coord,
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
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    # Make star masks
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

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
                
            build_star_mask(
                image_file = staged_image_file,
                star_file = star_intens_file,
                outfile = star_mask_file,
                rms_value = fid_rms,
                rms_fac = 3.0,
                overwrite = True,
            )
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
    # Fit and subtract a background
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
    
    if 'bkgrd' in tasks:

        print("Fitting backgrounds ...")

        for this_band in bands:

            print("... ... ", this_band)

            image_file = working_dirs['staged']+ \
                this_name+'_'+this_band+'_mjysr'+ \
                res_ext_for_bkgrd+'.fits'
                            
            rad_file = working_dirs['coords']+ \
                this_name+'_rgal.fits'

            bkgrd_file = working_dirs['bkgrd']+ \
                this_name+'_'+this_band+'_mjysr_bkgrd.fits'

            if survey == 'galex':
                weight_file = working_dirs['staged']+ \
                    this_name+'_'+this_band+'_weight.fits'
            else:
                weight_file = None

            galaxy_mask_file = working_dirs['masks']+ \
                this_name+'_galmask.fits'
            
            star_mask_file = working_dirs['masks']+\
                this_name+'_'+this_band+'_mask.fits'
                
            mask_fname_list = [
                galaxy_mask_file,
                star_mask_file,
            ]
            
            if survey == 'galex':
                bkgrd_methods = ['itermed','mode']

            if survey == 'sdss':
                bkgrd_methods = ['itermed','mode']        

            if an_unwise_survey and (this_band in ['w1','w2']):
                bkgrd_methods = ['itermed','mode']
                
            if an_unwise_survey and (this_band in ['w3','w4']):
                bkgrd_methods = ['itermed','mode','planefit']
            
            fit_z0mgs_background(
                image_fname = image_file,
                mask_fnames = mask_fname_list,
                weight_fname = weight_file,
                rad_fname = rad_file,
                fid_rad = rgal,
                outfile_bkgrd = bkgrd_file,                
                methods = bkgrd_methods,
                clip_thresh = 3.0,
                niter = 5,
            )

            print("... applying background subtraction.")

            # Read the background
            bkgrd_hdu = fits.open(bkgrd_file)
            bkgrd = bkgrd_hdu[0].data

            print("... ... to native")

            # Background subtract the native res image
            staged_image_fname = working_dirs['staged']+ \
                this_name+'_'+this_band+'_mjysr.fits'
            staged_image_hdu = fits.open(staged_image_fname)[0]
            staged_image = staged_image_hdu.data
            staged_image_hdr = staged_image_hdu.header
            
            bksub_image_fname = working_dirs['bkgrd']+ \
                this_name+'_'+this_band+'_mjysr_bksub.fits'
            bksub_image_hdu = fits.PrimaryHDU(
                data=staged_image - bkgrd, header=staged_image_hdr)
            bksub_image_hdu.writeto(bksub_image_fname, overwrite=True)

    if 'plot_bkgrd' in tasks:        

        for this_band in bands:

            print("Plotting background subtraction.")
            print("... to be implemented.")
        
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Do the convolutions
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        

    # Convolve the final products to to some specified Gaussian
    # resolutions.

    # TBD - masking
    
    if 'convolve' in tasks:

        print("Convolving images ...")

        for this_band in bands:

            print("... ... ", this_band)
            
            for target_res in res_dict.keys():
                print("... to "+target_res)
                
                kern_to_this_res = z0mgs_kernel_name(
                    from_res=this_band, to_res=target_res)
                
                if os.path.isfile(kern_to_this_res) == False:
                    print("... ... kernel not found: ", kern_to_this_res)
                    print("... ... continuing without it.")
                    print("... ... create it if it should be there.")
                    continue
                
                staged_image_file = working_dirs['bkgrd']+ \
                    this_name+'_'+this_band+'_mjysr_bksub.fits'
                
                convolved_image_file = working_dirs['bkgrd']+ \
                    this_name+'_'+this_band+'_mjysr_bksub_'+target_res+'.fits'
                
                convolve_image_with_kernel(
                    image_file=staged_image_file,
                    outfile=convolved_image_file,
                    kernel_file=kern_to_this_res,
                    blank_zeros=False,
                    overwrite=True
                )
                
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
    # Plot everything together
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    

    # Plot the final images.
    
    if 'plot_results' in tasks:
        
        pass

