import os
from astropy.table import vstack
from astroquery.gaia import Gaia

# def query_gaia_for_image()

# def query_gaia_inchunks()

def query_gaia(
        ra_min=None,
        ra_max=None,
        dec_min=None,
        dec_max=None,
        outfile=None,
        add_fields=[],
        omit_fields=[],
        quiet=False,
        dry_run=True,
        skip_if_present=False):
    """
    Query GAIA DR3 with parameters of interest to extragalactic images in a rectangular footprint.
    """

    tab = None
    delim=' '

    source = 'FROM'+delim+'gaiadr3.gaia_source'+delim+'AS'+delim+'src'+delim+'JOIN'+delim+'gaiadr3.astrophysical_parameters'+delim+'AS'+delim+'ap'+delim+'ON'+delim+'src.source_id'+delim+'='+delim+'ap.source_id'+delim+''
    fields = [
        'src.source_id',
        'src.ra',
        'src.dec',
        'ap.classprob_dsc_combmod_quasar',
        'ap.classprob_dsc_combmod_star',
        'ap.classprob_dsc_combmod_galaxy',
        'ap.classprob_dsc_combmod_binarystar',
        'ap.classprob_dsc_combmod_whitedwarf',
        'ap.teff_gspphot',
        'ap.teff_gspspec',
        'src.phot_g_mean_mag',
        'src.phot_bp_mean_mag',
        'src.phot_rp_mean_mag',        
        'src.parallax',
        'src.parallax_error',
        'src.pmra',
        'src.pmra_error',
        'src.pmdec',
        'src.pmdec_error',
        'src.ruwe',
    ]
    for field_to_add in add_fields:
        if fields.count(field_to_add) == 0:
            fields.append(field_to_add)
    for field_to_omit in omit_fields:
        if fields.count(field_to_omit) == 1:
            fields.remove(field_to_omit)
    
    all_fields = 'SELECT'+delim+''
    first=True
    for this_field in fields:
        if first == False:
            all_fields += ','
        all_fields += this_field.strip()
        first = False
    
    # Deal with longitude wrap
    if ra_min > ra_max:
        twopart_call = True
        this_ramin = ra_min
        this_ramax = 360.0
    else:
        twopart_call = False
        this_ramin = ra_min
        this_ramax = ra_max

    # Convert to strings
    ra_low_string = (f'{this_ramin:.5f}').strip()
    ra_high_string = (f'{this_ramax:.5f}').strip()

    dec_low_string = (f'{dec_min:.5f}').strip()
    dec_high_string = (f'{dec_max:.5f}').strip()

    selection = 'WHERE'+delim+'src.ra'+delim+'BETWEEN'+delim+''+ra_low_string+''+delim+'AND'+delim+''+ra_high_string+ \
        ''+delim+'AND'+delim+'src.dec'+delim+'BETWEEN'+delim+''+dec_low_string+''+delim+'AND'+delim+''+dec_high_string
        
    query = all_fields + ' ' + source + ' ' + selection

    if dry_run == False:
        if quiet == False:
            print(query)
        if outfile is not None:
            if (os.path.isfile(outfile) == False) or \
               (skip_if_present == False):            
                job = Gaia.launch_job_async(query)
                tab = job.get_results()
                if twopart_call == False:            
                    tab.write(outfile, format='fits', overwrite=True)
    else:
        if quiet == False:
            print(query)

    if twopart_call == False:
        return(query)

    # From here only worry about the case where we wrapped around the meridian
  
    first_part = query
    
    this_ramin = 0.0
    this_ramax = ra_max

    # Convert to strings
    ra_low_string = (f'{this_ramin:.5f}').strip()
    ra_high_string = (f'{this_ramax:.5f}').strip()

    selection = 'WHERE'+delim+'src.ra'+delim+'BETWEEN'+delim+''+ra_low_string+''+delim+'AND'+delim+''+ra_high_string+ \
        ''+delim+'AND'+delim+'src.dec'+delim+'BETWEEN'+delim+''+dec_low_string+''+delim+'AND'+delim+''+dec_high_string

    query = all_fields + ' ' + source + ' ' + selection
    
    if dry_run == False:
        if quiet==False:
            print(query)
        if outfile is not None:            
            if (os.path.isfile(outfile) == False) or \
               (skip_if_present == False):                        
                job = Gaia.launch_job_async(query)
                tab2 = job.get_results()
                tab = vstack([tab,tab2])
                tab.write(outfile, format='fits', overwrite=True)
    else:
        if quiet==False:
            print(query)

    return([first_part,query]) 
        

