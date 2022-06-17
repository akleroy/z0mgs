import astropy.io.fits
import os
import numpy as np

def gal_data(names=None, data=None, all=False, galdata_dir=None, tag=None):
    """ 
    Access the galbase 

    Parameters
    ----------
    names : str, list, optional
        A galaxy name or list of galaxy names in string format
    data : structured array, optional
        Structured array, such as read in from a FITS file, containing all of the galaxy data
    all : bool, optional
        Return all of the galaxy data (Default: False)
    galdata_dir : string, optional
        The path to the directory that contains the database
    tag : string, optional
        Select only galaxy data with a specific tag; e.g., 'SINGS'

    Returns
    -------
    data: structured array
        All of the galaxy data contained in the database that meets the desired tags
    """

    if not names and not all and not tag:
        print('Need a name to find a galaxy. Returning None.')
        return None

    # If galdata_dir is not specified, galdata_dir is contained within the working directory
    if not galdata_dir:
        galbase_dir, this_filename = os.path.split(__file__)
        galdata_dir = os.path.join(galbase_dir, "gal_data")

    # read in the data if it is not passed in
    if data is None:
        dbfile = os.path.join(galdata_dir, 'gal_base.fits')
        hdulist = astropy.io.fits.open(dbfile)
        data = hdulist[1].data
        hdulist.close()

    # all data are desried
    if all:
        return data

    # use only a specific survey
    if tag is not None:
        n_data = len(data)
        keep = np.ones(n_data)
        # survey_file = os.path.join(galdata_dir, 'survey_' + tag.lower() + '.txt')
        # gals = np.genfromtxt(survey_file, dtype='string')

        for i in range(n_data):
            this_tag = data['tags'][i].strip(';').split(';;')
            keep[i] = sum(np.in1d(tag, this_tag))

        if np.sum(keep) == 0:
            print('No targets found with that tag combination.')
            return None

        good_data = data[np.where(keep)]

        return good_data

    # build the alias
    aliases = {}
    fname = os.path.join(galdata_dir, "gal_base_alias.txt")
    f = open(fname)
    f.readline()
    f.readline()
    for line in f:
        s = [temp for temp in line.strip('\n').split(' ')]
        aliases[s[0].replace(' ', '').upper()] = s[1].replace(' ', '').upper()

    # name or names of galaxies -- put single galaxy name into list
    if type(names) == str:
        names = [names]
    keep = np.zeros(len(data), dtype=bool)

    # search for the specified galaxies
    for name in names:
        name_a = aliases[name.replace(' ', '').upper()]
        ind = data.field('NAME') == name_a

        if sum(ind) == 0:
            print('No match for ' + name)
        else:
            keep += ind

    return data[keep]