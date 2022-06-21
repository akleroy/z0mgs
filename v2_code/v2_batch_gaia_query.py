from astropy.table import Table
import astropy.table
from utils_query_gaia import *
from astropy.utils.console import ProgressBar

sample = ['largeleda']

tdir = '../../measurements/'

tab_dict = {
    'localgroup':tdir+'unwise_v2_index_localgroup.fits',
    'largeleda':tdir+'unwise_v2_index_largeleda.fits',
    'smallleda':tdir+'unwise_v2_index_smallleda.fits',
    'localvolume':tdir+'unwise_v2_index_localvolume.fits',
    'manga':tdir+'unwise_v2_index_manga.fits',
}

out_root = '../../working_data/gaia/'

for this_sample in sample:

    out_dir = out_root + this_sample + '/'
    tab = Table.read(tab_dict[this_sample])

    for this_row in ProgressBar(tab):

        out_file_name = out_dir+this_row['Z0MGS_NAME'].strip()+'_gaia_dr3.fits'
        
        query_gaia(
            ra_min=float(this_row['TRC_RA']),
            ra_max=float(this_row['BLC_RA']),
            dec_min=float(this_row['BLC_DEC']),
            dec_max=float(this_row['TRC_DEC']),
            outfile=out_file_name,
            add_fields=[],
            omit_fields=[],
            quiet=True,
            dry_run=False)
    
