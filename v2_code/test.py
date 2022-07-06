from utils_gaia_cutouts import *

image_fname = \
    '../../working_data/unwise/staged/largeleda/PGC47404_w1_mjysr.fits'
gaia_table_fname = \
    '../../working_data/gaia/largeleda/PGC47404_gaia_dr3.fits'
stack = extract_gaia_stack(
    image_fname, gaia_table_fname,
    oversamp_fac=2.0, half_width_pix=20)
stack.writeto('../../scratch/test_gaia.fits',overwrite=True)
