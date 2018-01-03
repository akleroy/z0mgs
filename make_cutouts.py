import numpy as np
import gal_data
import extract_stamp
import warnings
import argparse
from pdb import set_trace


_GALDATA_DIR = '/n/home00/lewis.1590/research/galbase/gal_data/'

def get_args():
    """ Get command line arguments """
    parser = argparse.ArgumentParser(description='Create cutouts of a given size around each galaxy center.')
    parser.add_argument('--imtype', default='int', help='input images to use. Default: int')
    parser.add_argument('--wttype', default='rrhr', help='images to use for the weighting. Default: rrhr')
    parser.add_argument('--size', default=30, help='cutout size in arcminutes. Default: 30.')
    parser.add_argument('--desired_pix_scale', default=1.5, help='desired pixel scale of output image. Default: 1.5 (GALEX)')
    parser.add_argument('--band', default=None, help='waveband. Default: None (does all)')
    parser.add_argument('--model_bg', default=True, help='model the background to match all images as best as possible. Default: True')
    parser.add_argument('--weight_ims', default=True, help='weight the input images by the desired weight images. Default: True')
    parser.add_argument('--convert_mjysr', action='store_true', help='set to convert images to MJy/sr. Default: False')
    parser.add_argument('--galaxy_list', default=None, nargs='+', help='Galaxy name if doing a single cutout or list of names. Default: None')
    parser.add_argument('--all_galaxies', action='store_true', help='run all galaxies in database. Default: False; include flag to store_true')
    parser.add_argument('--tag', default=None, help='tag to select galaxies, i.e., SINGS, HERACLES, etc. Default: None')
    parser.add_argument('--inds', nargs=2, type=int, help='index the all galaxy array')
    return parser.parse_args()


def main(**kwargs):
    """ Create cutouts using all data and combining images using Montage
    
    Parameters
    ----------
    imtype : str
        input image type to use from galex (Default: int)
    wttype : str
        input weights image type to use from galex (Default: rrhr)
    size : float
        cutout size in arcminutes (Default: 30.0)
    desired_pix_scale : float
        Desired pixel scale of output image. Default is currently set to GALEX pixel scale (Default: 1.5)
    band : str
        the band in which the cutout is made, either a single band or a list (Default: fuv)
    model_bg : bool
        model the background in Montage
    weight_ims : bool
        weight the input images with the weights images
    convert_mjysr : bool
        convert input images from counts/sec to MJy/sr
    galaxy_list : list
        list of one or more galaxies for which to make cutouts. Do not set if you want to make cutouts for all galaxies (Default: None)
    all_galaxies : bool
        Make cutouts for all galaxies in the galbase
    tag : str
        A tag to select a subset of galaxies; i.e., SINGS, HERACLES, etc. (Default: None)
    inds : int
        List of two ints to index the galaxy array from [int1:int2]

    Example:
    This code can be run from the command line or imported and run within a separate program.
    The following example creates cutouts of the SINGS sample that are 30x30 arcminutes in the FUV with 
    modeled backgrounds, images weighted by the exposure time, and converted to MJy/sr

    Usage:
    %run make_cutouts.py --size 30 --band fuv --convert_mjysr --tag SINGS
    (the model_bg and weight_ims flags do not need to be explicitly set as they are True by default)

    or

    import make_cutouts
    make_cutouts.main(size=30, band='fuv', convert_mjysr=True, tag='SINGS')
    """

    warnings.filterwarnings('ignore')
    wband = kwargs['band']
    if kwargs['inds']:
        kwargs['all_galaxies'] = True

    #get data from galbase
    gals = gal_data.gal_data(names=kwargs['galaxy_list'], data=None, all=kwargs['all_galaxies'], 
                             galdata_dir=_GALDATA_DIR, tag=kwargs['tag']) 

    if kwargs['inds']:
        ind_start, ind_stop = kwargs['inds'][0], kwargs['inds'][1]
        gals = gals[ind_start:ind_stop]

    n_gals = len(gals)
    size_deg = kwargs['size'] * 60. / 3600. #convert from arcminutes to degrees

    for i in range(n_gals):
        galname = gals['name'][i].replace(' ', '').upper()
        pgcname = gals['pgcname'][i]
        ra_ctr, dec_ctr = gals['ra_deg'][i], gals['dec_deg'][i]

        stamp_kwargs = {'ra_ctr': ra_ctr, 'dec_ctr': dec_ctr, 'size_deg': size_deg, 'name': galname, 
                        'pgcname': pgcname, 'model_bg': kwargs['model_bg'], 
                        'weight_ims': kwargs['weight_ims'], 'convert_mjysr': kwargs['convert_mjysr'], 
                        'imtype': kwargs['imtype'], 'wttype': kwargs['wttype'], 
                        'desired_pix_scale': kwargs['desired_pix_scale']}
        if wband == 'fuv':
            extract_stamp.galex(band='fuv', **stamp_kwargs)
        elif wband == 'nuv':
            extract_stamp.galex(band='nuv', **stamp_kwargs)
        else:
            extract_stamp.galex(band='fuv', **stamp_kwargs)
            extract_stamp.galex(band='nuv', **stamp_kwargs)



if __name__ == '__main__':
    args = get_args()
    main(**vars(args))
