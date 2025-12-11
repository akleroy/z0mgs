import os, sys
import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
import astropy.table
from utils_gaia_cutouts import *
from astropy.utils.console import ProgressBar

# Master tables
tdir = '../../measurements/'
tab_dict = {
    'localgroup':tdir+'unwise_v2_index_localgroup.fits',
    'largeleda':tdir+'unwise_v2_index_largeleda.fits',
    'smallleda':tdir+'unwise_v2_index_smallleda.fits',
    'localvolume':tdir+'unwise_v2_index_localvolume.fits',
    'manga':tdir+'unwise_v2_index_manga.fits',
}

# Skip the SMC
skip_pgc = [3085]

sample = []
for ii, this_arg in enumerate(sys.argv):
    if this_arg in tab_dict.keys():
        sample.append(this_arg)

if len(sample) == 0:
    sample = ['largeleda']

print("Fitting stars for.")
    
bands = ['fuv','nuv','w1','w2','w3','w4']

working_dir = {
    'fuv':'../../working_data/galex/',
    'nuv':'../../working_data/galex/',
    'w1':'../../working_data/unwise/',
    'w2':'../../working_data/unwise/',
    'w3':'../../working_data/unwise/',
    'w4':'../../working_data/unwise/',
}

gaia_root = '../../working_data/gaia/'

# Loop over target samples
for this_sample in sample:

    tab = Table.read(tab_dict[this_sample])
    
    for this_row in ProgressBar(tab):

        print('')
        print('Filling in table ', this_row['Z0MGS_NAME'].strip())
        
        if this_row['PGC'] in skip_pgc:
            print("Skipping ", this_row['Z0MGS_NAME'].strip())
            continue

        gaia_table_name = gaia_root+ this_sample+'/'+ \
            this_row['Z0MGS_NAME'].strip() + '_gaia_dr3.fits'
        
        out_table_name = gaia_root+ this_sample+'/'+ \
            this_row['Z0MGS_NAME'].strip() + '_gaia_dr3_filled.fits'
        
        if os.path.isfile(gaia_table_name) == False:
            print("Gaia table not found, skipping: ", gaia_table_name)
            continue

        tab = Table.read(gaia_table_name)
        
        for this_band in bands:

            print("... filling in ", this_band)
            
            # Initialize fields
            tab[this_band] = np.nan
            tab[this_band+'_BKSUB'] = np.nan
            
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
            
            # Do a rough background subtraction and save that value
            edge = 5

            edge_mask = np.zeros((ny,nx),dtype=int)
            edge_mask[0:edge-1,:] = 1
            edge_mask[:,0:edge-1] = 1
            edge_mask[-1*(edge+1):,:] = 1
            edge_mask[:,-1*(edge+1):] = 1

            bksub_stack = hdu_image.data
            for zz in range(nz):
                plane = bksub_stack[zz,:,:]
                if np.sum(np.isfinite(plane)) == 0:
                    bksub_stack[zz,:,:] = plane
                    continue
                bkgrd = np.nanmedian(plane[edge_mask])
                bksub_stack[zz,:,:] = plane-bkgrd
            
            tab[this_band+'_BKSUB'] = \
                np.array(bksub_stack[:,mid_y,mid_x])
                        
        # Write the updated table to disk as a new file
        tab.write(out_table_name, format='fits', overwrite=True)
