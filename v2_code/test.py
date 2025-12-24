from z0mgs_build_atlas import *

z0mgs_build_atlas(
    subsamples='largeleda',
    #tasks = ['stage'],
    #tasks = 'all',
    tasks = [
        #'stage',
        #'plot_stage',
        #'convolve_for_bkgrd',
        #'coord_mask',
        #'galaxy_mask',
        #'star_pred',
        #'star_mask',
        #'bkgrd',
        #'plot_bkgrd',
        'convolve',
        'plot_results',
    ],
    #bands=['w3','w4',],
    bands = ['u','g','r','i','z'],
    root_dir='../../test_data/',
    #survey='galex',
    #survey='unwise_custom',
    survey='sdss',
    just_galaxy='PGC47404',
)
    
#from utils_gaia_cutouts import *

#image_fname = \
#    '../../working_data/unwise/staged/largeleda/PGC47404_w1_mjysr.fits'
#gaia_table_fname = \
#    '../../working_data/gaia/largeleda/PGC47404_gaia_dr3.fits'
#stack = extract_gaia_stack(
#    image_fname, gaia_table_fname,
#    oversamp_fac=2.0, half_width_pix=20)
#stack.writeto('../../scratch/test_gaia.fits',overwrite=True)
