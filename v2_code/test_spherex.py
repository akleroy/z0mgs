import os, glob
from utils_spherex import *

root_dir = '../../test_data/spherex/'
this_gal = 'ngc5194'
#this_gal = 'ngc0253'
gal_dir = root_dir + this_gal + '/'
outdir = gal_dir + 'raw/'

re_download = False
make_sed = False
grid_cube = True

if re_download:
    
    os.system('mkdir '+gal_dir)
    os.system('mkdir '+outdir)

    image_tab = \
        search_spherex_images(
            target = this_gal,
            radius = 30*u.arcmin,
            collection = 'spherex_qr2',
            verbose = True)

    downloaded_images = \
        download_images(
            image_tab,
            outdir = outdir,
            incremental = True,
            verbose = True)

if make_sed:

    im_list = glob.glob(outdir+'*.fits')
    
    extract_spherex_sed(
        target_coord = this_gal,        
        image_list = im_list,
        outfile = gal_dir + this_gal + '_spherex_ctrsed.ecsv',
        overwrite=True)
    
if grid_cube:
    
    cube_hdu = make_cube_header(
        center_coord = this_gal,
        pix_scale = 6. / 3600.,
        extent = 30. / 60., 
        lam_min = 0.75, lam_max = 5.2, lam_step = 0.1,
        return_header=False)

    im_list = glob.glob(outdir+'*.fits')
    #im_list = glob.glob(outdir+'level2_2025W21_1B_0519_4D4_spx_l2b-v20-2025-248.fits')
    grid_spherex_cube(
        target_hdu = cube_hdu,
        image_list = im_list,
        outfile = gal_dir + this_gal+'_spherex_cube.fits',
        overwrite = True)
