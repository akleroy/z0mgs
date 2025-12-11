# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Imports
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import os, sys, glob
import numpy as np

import matplotlib.pyplot as plt

import astropy.io.fits as fits
from astropy.table import Table
import astropy.table
from astropy.utils.console import ProgressBar

from utils_gaia_cutouts import *

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Control flow and definitions
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Master tables for each sample

tdir = '../../measurements/'
tab_dict = {
    'localgroup':tdir+'unwise_v2_index_localgroup.fits',
    'largeleda':tdir+'unwise_v2_index_largeleda.fits',
    'smallleda':tdir+'unwise_v2_index_smallleda.fits',
    'localvolume':tdir+'unwise_v2_index_localvolume.fits',
    'manga':tdir+'unwise_v2_index_manga.fits',
}

# Skip (set here to the SMC by hand)
skip_pgc = [3085]

# Option to only consider a short list
just_pgc = []

# Get list of samples to run from command line ...
sample = []
for ii, this_arg in enumerate(sys.argv):
    if this_arg in tab_dict.keys():
        sample.append(this_arg)

# ... or default to large LEDA sample
if len(sample) == 0:
    sample = ['largeleda']

print("Building out PSF estimates tables for sample: ", sample)

# Run all bands by default
bands = ['fuv','nuv','w1','w2','w3','w4']
bands = ['w4']

# Location of images of galaxies in each band
working_dir = {
    'fuv':'../../working_data/galex/',
    'nuv':'../../working_data/galex/',
    'w1':'../../working_data/unwise/',
    'w2':'../../working_data/unwise/',
    'w3':'../../working_data/unwise/',
    'w4':'../../working_data/unwise/',
}

# Location of the GAIA query files
gaia_root = '../../working_data/gaia/'

# Minimum pixels to compute a background value
min_pixels_for_bkgrd = 100

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Execute loop
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Loop over the samples (large LEDA, small LEDA, local group, etc.)

for this_sample in sample:

    # Read the table for this sample ...
    tab = Table.read(tab_dict[this_sample])

    # ... then loop over individual galaxies
    for this_row in ProgressBar(tab):

        print('')
        print('Filling in table ', this_row['Z0MGS_NAME'].strip())

        # Check against list to skip
        if this_row['PGC'] in skip_pgc:
            print("Skipping ", this_row['Z0MGS_NAME'].strip())
            continue

        # Check against list of targets
        if len(just_pgc) > 0:
            if this_row['PGC'] not in just_pgc:
                continue

        # Define table name
        gaia_table_name = gaia_root+ this_sample+'/'+ \
            this_row['Z0MGS_NAME'].strip() + '_gaia_dr3_filled.fits'

        # Verify file existence
        if os.path.isfile(gaia_table_name) == False:
            print("Gaia table not found, skipping: ", gaia_table_name)
            continue

        # Read the GAIA table
        tab = Table.read(gaia_table_name)

        # Select the stars to consider
        tab['phot_g_mean_mag'].fill_value = np.nan
        tab['phot_bp_mean_mag'].fill_value = np.nan
        tab['phot_rp_mean_mag'].fill_value = np.nan
        
        g = np.array(tab['phot_g_mean_mag'])
        b = np.array(tab['phot_bp_mean_mag'])
        r = np.array(tab['phot_rp_mean_mag'])
        bmr = b-r                
        
        tab['parallax_error'].fill_value = np.nan
        tab['parallax'].fill_value = np.nan
        tab['ruwe'].fill_value = np.nan
        tab['pmra'].fill_value = np.nan
        tab['pmra_error'].fill_value = np.nan
        tab['pmdec'].fill_value = np.nan
        tab['pmdec_error'].fill_value = np.nan
        tab['classprob_dsc_combmod_star'].fill_value = np.nan
        
        e_mu = np.array(tab['parallax_error'].filled())
        mu = np.array(tab['parallax'].filled())
        dm = 5.*np.log10(1./(np.pi*1e-3))-5.0
        
        e_pmra = np.array(tab['pmra_error'].filled())
        pmra = np.array(tab['pmra'].filled())
        e_pmdec = np.array(tab['pmdec_error'].filled())
        pmdec = np.array(tab['pmdec'].filled())
        has_proper_motion = np.logical_or(pmra > 3.0*e_pmra, pmdec > 3.0*e_pmdec)
        
        #prob_star = np.array(tab['classprob_dsc_combmod_star'].filled())
        prob_star = np.array(tab['classprob_dsc_combmod_star'])
        ruwe = np.array(tab['ruwe'].filled())
        
        use_star = np.logical_or(mu > 3.*e_mu,has_proper_motion)*(prob_star > 0.9)*(ruwe < 1.2)*(g < 20.)
        print("I selected ... ", np.sum(use_star), " stars to use.")        
        
        # Loop over bands being considered
        for this_band in bands:

            print("... building psf for ", this_band)            

            if this_band == 'w1':
                vmin=-0.001
                vmax=0.01
                log_levs = [-2.9,-2.7,-2.5,-2.3,-2.,-1,-0.3]

            if this_band == 'w2':
                vmin=-0.001
                vmax=0.01
                log_levs = [-2.9,-2.7,-2.5,-2.3,-2.,-1,-0.3]

            if this_band == 'w3':
                vmin=-0.01
                vmax=0.1
                log_levs = [-2.9,-2.7,-2.5,-2.3,-2.,-1,-0.3]

            if this_band == 'w4':
                vmin=-0.01
                vmax=0.1
                log_levs = [-2.9,-2.7,-2.5,-2.3,-2.,-1,-0.3]
                
            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            # Get the stellar cutouts
            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            
            # Identify file holding the cutouts around stars for this band ...
            this_stack_dir = \
                working_dir[this_band]+'star_stacks/'+this_sample+'/'
            
            stack_file_name = this_stack_dir + \
                this_row['Z0MGS_NAME'].strip()+ \
                '_'+this_band+ \
                '_starstack.fits'

            # Check file existence
            if os.path.isfile(stack_file_name) == False:
                print("Stack not found, skipping: ", stack_file_name)
                continue
            
            # Open the stack
            hdu_image = fits.open(stack_file_name)[0]
            
            # Read the image stack and select the subset of planes we are using.
            sub_tab = tab[use_star]
            sub_stack = hdu_image.data[use_star,:,:]
            
            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            # Background-subtract and scale the image
            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            nz, ny, nx = sub_stack.shape

            peak_vec = np.array(sub_tab[this_band]).reshape(nz,1,1)
            bkgrd_vec = np.array(sub_tab[this_band+'_bkgrd_mean']).reshape(nz,1,1)
            
            norm_vec = peak_vec-bkgrd_vec
            weight_vec = norm_vec**2
            
            print("Cube stack shape: ", sub_stack.shape)
            print("Factor vector shape: ", norm_vec.shape)

            bksub_stack = sub_stack - bkgrd_vec
            norm_stack = bksub_stack / norm_vec

            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            # Collapse
            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            # with sums
            weighted_stack = norm_stack * weight_vec
            collapsed_image = np.nansum(weighted_stack, axis=0)
            collapsed_weight = np.nansum(weight_vec)
            model_psf = collapsed_image / collapsed_weight

            # with a sort
            weight_stack = np.ones_like(norm_stack)*weight_vec
            ai = np.argsort(norm_stack,axis=0)
            sorted_norm = np.take_along_axis(norm_stack,ai,axis=0)
            sorted_weights = np.take_along_axis(weight_stack,ai,axis=0)

            print("... interpolating")
            model_psf = np.zeros((ny,nx))*np.nan
            for ii in np.arange(ny):
                for jj in np.arange(nx):
                    weight_vec = sorted_weights[:,ii,jj]
                    value_vec = sorted_norm[:,ii,jj]
                    keep = np.isfinite(weight_vec*value_vec)
                    weight_vec = weight_vec[keep]
                    value_vec = value_vec[keep]
                    cum_weight_vec = np.cumsum(weight_vec)/np.sum(weight_vec)
                    #print(cum_weight_vec.shape)
                    #print(value_vec.shape)
                    model_psf[ii,jj] = np.interp(0.5,cum_weight_vec,value_vec)
            
            plt.clf()
            plt.imshow(model_psf, vmin=vmin, vmax=vmax, origin='lower', cmap='viridis')
            plt.contour(np.log10(model_psf), levels=log_levs, 
                        origin='lower', colors='red',linestyle='solid')
            plt.show()
            
