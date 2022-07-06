import os, sys
from astropy.table import Table
import astropy.table
from utils_gaia_cutouts import *
from astropy.utils.console import ProgressBar

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

bands = ['fuv','nuv','w1','w2','w3','w4']

working_dir = {
    'fuv':'../../working_data/galex/',
    'nuv':'../../working_data/galex/',
    'w1':'../../working_data/unwise/',
    'w2':'../../working_data/unwise/',
    'w3':'../../working_data/unwise/',
    'w4':'../../working_data/unwise/',
}

print("Gaia cutouts for sample: ", sample)

gaia_root = '../../working_data/gaia/'


for this_sample in sample:

    tab = Table.read(tab_dict[this_sample])
    
    for this_row in ProgressBar(tab):

        print('')
        print('Querying ', this_row['Z0MGS_NAME'].strip())
        
        if this_row['PGC'] in skip_pgc:
            print("Skipping ", this_row['Z0MGS_NAME'].strip())
            continue

        gaia_table_name = gaia_root+ this_sample+'/'+ \
            this_row['Z0MGS_NAME'].strip() + '_gaia_dr3.fits'
        
        if os.path.isfile(gaia_table_name) == False:
            print("Gaia table not found, skipping: ", gaia_table_name)
            continue
        
        for this_band in bands:
        
            this_staged_dir = \
                working_dir[this_band]+'staged/'+this_sample+'/'
            this_stack_dir = \
                working_dir[this_band]+'star_stacks/'+this_sample+'/'
            
            image_file_name = this_staged_dir+ \
                this_row['Z0MGS_NAME'].strip() + '_' + \
                this_band + '_mjysr.fits'
            
            out_file_name = this_stack_dir + \
                this_row['Z0MGS_NAME'].strip()+ \
                '_'+this_band+ \
                '_starstack.fits'

            # Check file existence

            if os.path.isfile(image_file_name) == False:
                print("Image file not found, skipping: ", image_file_name)
                continue
            
            # Obtain and write the stack
            extract_gaia_stack(
                image_file_name,
                gaia_table_fname=gaia_table_name,
                out_fname=out_file_name, overwrite=True,
                oversamp_fac=2.0, half_width_pix=20,
                order='bilinear')
            
