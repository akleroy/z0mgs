# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Imports
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import os, sys
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

print("Building out GAIA tables for sample: ", sample)

# Run all bands by default
bands = ['fuv','nuv','w1','w2','w3','w4']

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

        # Define input and output names
        gaia_table_name = gaia_root+ this_sample+'/'+ \
            this_row['Z0MGS_NAME'].strip() + '_gaia_dr3.fits'
        
        out_table_name = gaia_root+ this_sample+'/'+ \
            this_row['Z0MGS_NAME'].strip() + '_gaia_dr3_filled.fits'

        # Verify file existence
        if os.path.isfile(gaia_table_name) == False:
            print("Gaia table not found, skipping: ", gaia_table_name)
            continue

        # Read the GAIA table
        tab = Table.read(gaia_table_name)

        # Loop over bands being considered
        for this_band in bands:

            print("... filling in ", this_band)
            
            # Initialize fields
            tab[this_band] = np.nan
            tab[this_band+'_bksub'] = np.nan
            tab[this_band+'_bkgrd'] = np.nan

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
            
            # Read the middle pixel
            nz, ny, nx = hdu_image.shape
            mid_y = int(np.floor(ny/2))
            mid_x = int(np.floor(nx/2))
            tab[this_band] = np.array(hdu_image.data[:,mid_y,mid_x])

            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            # Construct a mask that identifies the background
            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            
            # ... width of edge mask in pixels
            edge = 5

            edge_mask = np.zeros((ny,nx),dtype=int)
            edge_mask[0:edge,:] = 1
            edge_mask[:,0:edge] = 1
            edge_mask[-1*(edge+1):,:] = 1
            edge_mask[:,-1*(edge+1):] = 1

            plt.clf()
            plt.imshow(edge_mask, origin='lower')
            plt.show()
            
            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            # Subtract the background from each plane of the image
            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            # ... read the stack in, this copy will be our background subtracted version
            star_stack = hdu_image.data
            bkgrd_mean = np.nan*np.zeros(nz)
            bkgrd_median = np.nan*np.zeros(nz)

            # ... loop over planes
            for zz in range(nz):
                plane = star_stack[zz,:,:]
                
                # ... if the plane is empty, skip it
                if np.sum(np.isfinite(plane)) == 0:                    
                    continue

                # ... also check that we get a background
                if np.sum(np.isfinite(plane[edge_mask])) < min_pixels_for_bkgrd:                    
                    #print("Not enough pixels for background calculation.")
                    continue

                # Estimate the median within the mask
                this_bkgrd_median = np.nanmedian(plane[edge_mask])                
                this_bkgrd_mean = np.nanmean(plane[edge_mask])                
                
                bkgrd_median[zz] = this_bkgrd_median
                bkgrd_mean[zz] = this_bkgrd_mean
                
            # ... save into the file            
            tab[this_band+'_bkgrd_mean'] = \
                np.array(bkgrd_mean)
            tab[this_band+'_bkgrd_median'] = \
                np.array(bkgrd_median)
                        
        # Write the updated table to disk as a new file
        tab.write(out_table_name, format='fits', overwrite=True)
