import numpy as np
import astropy.io.fits
import astropy.wcs
import montage_wrapper as montage
from matplotlib.path import Path
import os
import sys
import shutil
import glob
import time

import scipy.ndimage as sp

from pdb import set_trace


# directories that contain the input data
_TOP_DIR = '/data/tycho/0/leroy.42/allsky/'
_INDEX_DIR = os.path.join(_TOP_DIR, 'z0mgs/')
_WISE_DIR = os.path.join(_TOP_DIR, 'unwise', 'atlas')

# directories to do the work in
_WORK_DIR = '/data/tycho/0/leroy.42/allsky/galex/atlas/'
_MOSAIC_DIR = os.path.join(_WORK_DIR, 'cutouts')

# CALIBRATION FROM GALEX COUNTS TO ABMAG
FUV2AB = 18.82
NUV2AB = 20.08
UV2AB = {'fuv': FUV2AB, 'nuv': NUV2AB}

GALEX_PIX_AS = 1.5 ## galex pixel scale in arcseconds -- from documentation

class GalaxyHeader(object):
    def __init__(self, name, gal_dir, ra_ctr, dec_ctr, pix_len, pix_scale, factor=1):
        self.name = name
        self.gal_dir = gal_dir
        self.ra_ctr = ra_ctr
        self.dec_ctr = dec_ctr
        self.hdr, self.hdrfile = self._create_hdr_output(pix_len, pix_scale, factor=1)
        self.hdr_ext, self.hdrfile_ext = self._create_hdr_output(pix_len, pix_scale, factor=factor)

    def _create_hdr_obj(self, pix_len, pix_scale):
        """
        Create a FITS header

        Parameters
        ----------
        ra_ctr : float
            RA of center of galaxy
        dec_ctr : float
            Dec of center of galaxy
        pix_len : float
            Length of each axis (square, so the same for x and y)
        pix_scale : float
            Pixel scale in degrees

        Returns
        -------
        hdr : astropy Header() object
            Newly created header object
        """
        hdr = astropy.io.fits.Header()
        hdr['NAXIS'] = 2
        hdr['NAXIS1'] = pix_len
        hdr['NAXIS2'] = pix_len
        hdr['CTYPE1'] = 'RA---TAN'
        hdr['CRVAL1'] = float(self.ra_ctr)
        hdr['CRPIX1'] = (pix_len / 2.) * 1.
        hdr['CDELT1'] = -1.0 * pix_scale
        hdr['CTYPE2'] = 'DEC--TAN'
        hdr['CRVAL2'] = float(self.dec_ctr)
        hdr['CRPIX2'] = (pix_len / 2.) * 1.
        hdr['CDELT2'] = pix_scale
        hdr['EQUINOX'] = 2000
        return hdr


    def _create_hdr_output(self, size_degrees, pixel_scale, factor=1):
        """
        Create a header and write it to an ascii file for use in Montage

        Parameters
        ----------
        galname : str
            Name of the galaxy
        ra_ctr : float
            Central RA of galaxy
        dec_ctr : float
            Central Dec of galaxy
        size_degrees : float
            size of cutout, in degrees
        pixel_scale : float
            pixel scale of output in arcseconds per pixel
        factor : int, optional
            Number by which to multiply size_degrees to extend the size of the cutout for bg modeling. (Default: 1)

        Returns
        -------
        target_hdr : astropy.header object
            The output header object
        header_file : str
            Path to the ascii file containing the header information
        """
        pix_len = int(np.ceil(size_degrees * factor / pixel_scale))
        hdr = self._create_hdr_obj(pix_len, pixel_scale)
        ri_targ, di_targ = self._make_axes(hdr)
        sz_out = ri_targ.shape
        outim = ri_targ * np.nan

        prihdu = astropy.io.fits.PrimaryHDU(data=outim, header=hdr)
        target_hdr = prihdu.header

        suff = '_template.hdr'
        if factor != 1:
            suff = suff.replace('.hdr', '_ext.hdr')
        header_file = os.path.join(self.gal_dir, self.name + suff)
        self.write_headerfile(header_file, target_hdr)

        return target_hdr, header_file


    def _make_axes(self, hdr, quiet=False, novec=False, vonly=False, simple=False):
        """
        Create axes arrays for the new mosaiced image. This is a simple translation to Python of Adam's
        IDL routine of the same name.

        Parameters
        ----------
        hdr : FITS header object
            FITS header to hold astrometry of desired output image
        quiet : bool, optional
            NOT USED
        novec : bool
            Find RA and Dec for every point (Default: False)
        vonly : bool
            Return only velocity data (Default: False)
        simple : bool
            Do the simplest thing (Default: False)

        Returns
        -------
        rimg : array
            array for ouptut RA
        dimg : array
            array for output Dec
        """

        # PULL THE IMAGE/CUBE SIZES FROM THE HEADER
        naxis  = int(hdr['NAXIS'])
        naxis1 = int(hdr['NAXIS1'])
        naxis2 = int(hdr['NAXIS2'])
        if naxis > 2:
            naxis3 = hdr['NAXIS3']

        ## EXTRACT FITS ASTROMETRY STRUCTURE
        ww = astropy.wcs.WCS(hdr)

        #IF DATASET IS A CUBE THEN WE MAKE THE THIRD AXIS IN THE SIMPLEST WAY POSSIBLE (NO COMPLICATED ASTROMETRY WORRIES FOR FREQUENCY INFORMATION)
        if naxis > 3:
            #GRAB THE RELEVANT INFORMATION FROM THE ASTROMETRY HEADER
            cd = ww.wcs.cd
            crpix = ww.wcs.crpix
            cdelt = ww.wcs.crelt
            crval = ww.wcs.crval

        if naxis > 2:
        # MAKE THE VELOCITY AXIS (WILL BE M/S)
            v = np.arange(naxis3) * 1.0
            vdif = v - (hdr['CRPIX3']-1)
            vaxis = (vdif * hdr['CDELT3'] + hdr['CRVAL3'])

        # CUT OUT HERE IF WE ONLY WANT VELOCITY INFO
        if vonly:
            return vaxis

        #IF 'SIMPLE' IS CALLED THEN DO THE REALLY TRIVIAL THING:
        if simple:
            print('Using simple aproach to make axes.')
            print('BE SURE THIS IS WHAT YOU WANT! It probably is not.')
            raxis = np.arange(naxis1) * 1.0
            rdif = raxis - (hdr['CRPIX1'] - 1)
            raxis = (rdif * hdr['CDELT1'] + hdr['CRVAL1'])

            daxis = np.arange(naxis2) * 1.0
            ddif = daxis - (hdr['CRPIX1'] - 1)
            daxis = (ddif * hdr['CDELT1'] + hdr['CRVAL1'])

            rimg = raxis # (fltarr(naxis2) + 1.)
            dimg = (np.asarray(naxis1) + 1.) # daxis
            return rimg, dimg

        # OBNOXIOUS SFL/GLS THING
        glspos = ww.wcs.ctype[0].find('GLS')
        if glspos != -1:
            ctstr = ww.wcs.ctype[0]
            newtype = 'SFL'
            ctstr.replace('GLS', 'SFL')
            ww.wcs.ctype[0] = ctstr
            print('Replaced GLS with SFL; CTYPE1 now =' + ww.wcs.ctype[0])

        glspos = ww.wcs.ctype[1].find('GLS')
        if glspos != -1:
            ctstr = ww.wcs.ctype[1]
            newtype = 'SFL'
            ctstr.replace('GLS', 'SFL')
            ww.wcs.ctype[1] = ctstr
            print('Replaced GLS with SFL; CTYPE2 now = ' + ww.wcs.ctype[1])

        # CALL 'xy2ad' TO FIND THE RA AND DEC FOR EVERY POINT IN THE IMAGE
        if novec:
            rimg = np.zeros((naxis1, naxis2))
            dimg = np.zeros((naxis1, naxis2))
            for i in range(naxis1):
                j = np.asarray([0 for i in xrange(naxis2)])

                pixcrd = np.array([[zip(float(i), float(j))]], numpy.float_)
                ra, dec = ww.all_pix2world(pixcrd, 1)

                rimg[i, :] = ra
                dimg[i, :] = dec
        else:
            ximg = np.arange(naxis1) * 1.0
            yimg = np.arange(naxis1) * 1.0
            X, Y = np.meshgrid(ximg, yimg, indexing='xy')
            ss = X.shape
            xx, yy = X.flatten(), Y.flatten()

            pixcrd = np.array(zip(xx, yy), np.float_)
            img_new = ww.all_pix2world(pixcrd, 0)
            rimg_new, dimg_new = img_new[:,0], img_new[:,1]

            rimg = rimg_new.reshape(ss)
            dimg = dimg_new.reshape(ss)

        # GET AXES FROM THE IMAGES. USE THE CENTRAL COLUMN AND CENTRAL ROW
        raxis = np.squeeze(rimg[:, naxis2/2])
        daxis = np.squeeze(dimg[naxis1/2, :])

        return rimg, dimg


    def write_headerfile(self, header_file, header):
        """
        Write out the header for the output mosaiced image

        Parameters
        ----------
        header_file : str
            Path to file to which to write header
        header : array
            The header to which to write to ASCII file
        """
        f = open(header_file, 'w')
        for iii in range(len(header)):
            outline = str(header[iii:iii+1]).strip().rstrip('END').strip()+'\n'
            f.write(outline)
        f.close()


    def append2hdr(self, keyword=None, value=None, ext=False):
        """
        Append information to the header and write to ASCII file

        Parameters
        ----------
        headerfile : str
            The path to the ascii file containing the header information
        keyword : str, optional
            The keyword in the header that you want to create (Default: None)
        value : multiple, optional
            The value to apply to the keyword (Default: None)
        """
        if keyword is not None:
            if ext:
                self.hdr_ext[keyword] = value
                self.write_headerfile(self.hdrfile_ext, self.hdr_ext)
            else:
                self.hdr[keyword] = value
                self.write_headerfile(self.hdrfile, self.hdr)

def calc_tile_overlap(ra_ctr, dec_ctr, pad=0.0, min_ra=0., max_ra=180., min_dec=-90., max_dec=90.):
    """
    Find all tiles that fall within a given overlap (pad) of (ra_ctr, dec_ctr)

    Parameters
    ----------
    ra_ctr : float
        Central RA
    dec_ctr : float
        Central Dec
    pad : float, optional
        Size of region about center (Default: 0.0)
    min_ra : float. optional
        Min RA of box to search in for overlaps (Default: 0.)
    max_ra : float, optional
        Max RA of box to search in (Default 180.)
    min_dec : float, optional
        Min Dec of box to search in (Default: -90.)
    max_dec : float, optional
        Max Dec of box to search in (Default: 90.)

    Returns
    -------
    overlap : bool array
        Bool arrat indicatinng which tiles in the index file fall within the given region
    """
    overlap = ((min_dec - pad) < dec_ctr) & ((max_dec + pad) > dec_ctr)

    #TRAP HIGH LATITUDE CASE AND (I GUESS) TOSS BACK ALL TILES. DO BETTER LATER
    mean_dec = (min_dec + max_dec) * 0.5
    if np.abs(dec_ctr) + pad > 88.0:
        return overlap

    ra_pad = pad / np.cos(np.radians(mean_dec))

    # MERIDIAN CASES
    merid = np.where(max_ra < min_ra)
    overlap[merid] = overlap[merid] & ( ((min_ra-ra_pad) < ra_ctr) | ((max_ra+ra_pad) > ra_ctr) )[merid]

    # BORING CASE
    normal = np.where(max_ra > min_ra)
    overlap[normal] = overlap[normal] & ((((min_ra-ra_pad) < ra_ctr) & ((max_ra+ra_pad) > ra_ctr)))[normal]

    return overlap

def galex(band='fuv', ra_ctr=None, dec_ctr=None, size_deg=None, index=None, name=None, pgcname=None, model_bg=False, weight_ims=False, convert_mjysr=False, desired_pix_scale=GALEX_PIX_AS, imtype='int', wttype='rrhr'):
    """
    Create cutouts of a galaxy in a single GALEX band.

    Parameters
    ----------
    band : str
        GALEX band to use
    ra_ctr : float
        Central RA of galaxy
    dec_ctr : float
        Central Dec of galaxy
    size_deg : float
        Desired side length of each cutout, in degrees
    index : array, optional
        Structured array containing the galbase information. The default is to read it in inside this code. (Default: None)
    name : str, optional
        Name of the galaxy for which to generate a cutout
    model_bg : bool, optional
        Model the background of the mosaiced image (Default: False)
    weight_ims : bool, optional
         weight the input images with the weights images
    convert_mjysr : bool, optional
        convert input images from counts/sec to MJy/sr
    desired_pix_scale : float, optional
        Desired pixel scale of output image. Default is currently set to GALEX pixel scale (Default: 1.5)
    imtype : str, optional
        input image type to use from galex (Default: int)
    wttype : str, optional
        input weights image type to use from galex (Default: rrhr)
    """
    ttype = 'galex'
    data_dir = os.path.join(_TOP_DIR, ttype, 'sorted_tiles')
    problem_file = os.path.join(_WORK_DIR, 'problem_galaxies_{}.txt'.format(band))# 'problem_galaxies_' + band + '.txt')
    bg_reg_file = os.path.join(_WORK_DIR, 'galex_reprojected_bg.reg')
    numbers_file = os.path.join(_WORK_DIR, 'gal_reproj_info_{}.txt'.format(band))# 'gal_reproj_info_' + band + '.dat')

    galaxy_mosaic_file = os.path.join(_MOSAIC_DIR, '_'.join([name, band]).upper() + '.FITS')

    if not os.path.exists(galaxy_mosaic_file):
        start_time = time.time()
        print name, band.upper()

        # READ THE INDEX FILE (IF NOT PASSED IN)
        if index is None:
            indexfile = os.path.join(_INDEX_DIR, 'galex_index_file.fits')
            ext = 1
            index, hdr = astropy.io.fits.getdata(indexfile, ext, header=True)

        # CALCULATE TILE OVERLAP
        tile_overlaps = calc_tile_overlap(ra_ctr, dec_ctr, pad=size_deg,
                                          min_ra=index['MIN_RA'],
                                          max_ra=index['MAX_RA'],
                                          min_dec=index['MIN_DEC'],
                                          max_dec=index['MAX_DEC'])

        # FIND OVERLAPPING TILES WITH RIGHT BAND
        #  index file set up such that index['fuv'] = 1 where fuv and
        #                              index['nuv'] = 1 where nuv
        ind = np.where((index[band]) & tile_overlaps)

        # MAKE SURE THERE ARE OVERLAPPING TILES
        ct_overlap = len(ind[0])
        if ct_overlap == 0:
            with open(problem_file, 'a') as myfile:
                myfile.write(name + ': ' + 'No overlapping tiles\n')
            return

        pix_scale = desired_pix_scale / 3600.  # 1.5 arbitrary: how should I set it?

        try:
            # CREATE NEW TEMP DIRECTORY TO STORE TEMPORARY FILES
            gal_dir = os.path.join(_WORK_DIR, '_'.join([name, band]).upper())
            os.makedirs(gal_dir)

            # MAKE HEADER AND EXTENDED HEADER AND WRITE TO FILE
            gal_hdr = GalaxyHeader(name, gal_dir, ra_ctr, dec_ctr, size_deg, pix_scale, factor=3)

            # GATHER THE INPUT FILES
            input_dir = os.path.join(gal_dir, 'input')
            os.makedirs(input_dir)
            nfiles = get_input(index, ind, data_dir, input_dir, hdr=gal_hdr)
            im_dir, wt_dir = input_dir, input_dir

            # WRITE TABLE OF INPUT IMAGE INFORMATION
            input_table = os.path.join(im_dir, 'input.tbl')
            montage.mImgtbl(im_dir, input_table, corners=True)
           
            #set_trace()
            if convert_mjysr:
                converted_dir = os.path.join(gal_dir, 'converted')
                if not os.path.exists(converted_dir):
                    os.makedirs(converted_dir)
                convert_to_flux_input(im_dir, converted_dir, band, desired_pix_scale, imtype=imtype)
                im_dir = converted_dir


            # MASK IMAGES
            masked_dir = os.path.join(gal_dir, 'masked')
            im_masked_dir = os.path.join(masked_dir, imtype)
            wt_masked_dir = os.path.join(masked_dir, wttype)
            for outdir in [masked_dir, im_masked_dir, wt_masked_dir]:
                os.makedirs(outdir)

            mask_images(im_dir, wt_dir, im_masked_dir, wt_masked_dir, imtype=imtype, wttype=wttype)
            im_dir = im_masked_dir
            wt_dir = wt_masked_dir


            # REPROJECT IMAGES WITH EXTENDED HEADER
            reprojected_dir = os.path.join(gal_dir, 'reprojected')
            reproj_im_dir = os.path.join(reprojected_dir, imtype)
            reproj_wt_dir = os.path.join(reprojected_dir, wttype)
            for outdir in [reprojected_dir, reproj_im_dir, reproj_wt_dir]:
                os.makedirs(outdir)

            reproject_images(gal_hdr.hdrfile_ext, im_dir, reproj_im_dir, imtype)
            reproject_images(gal_hdr.hdrfile_ext, wt_dir, reproj_wt_dir, wttype)
            im_dir = reproj_im_dir
            wt_dir = reproj_wt_dir


            # MODEL THE BACKGROUND IN THE IMAGE FILES WITH THE EXTENDED HEADER
            if model_bg:
                bg_model_dir = os.path.join(gal_dir, 'background_model')
                diff_dir = os.path.join(bg_model_dir, 'differences')
                corr_dir = os.path.join(bg_model_dir, 'corrected')
                for outdir in [bg_model_dir, diff_dir, corr_dir]:
                    os.makedirs(outdir)
                bg_model(im_dir, bg_model_dir, diff_dir, corr_dir, gal_hdr.hdrfile_ext, im_type=imtype, level_only=False)
                im_dir = os.path.join(corr_dir, 'int')


            # WEIGHT IMAGES
            if weight_ims:
                weight_dir = os.path.join(gal_dir, 'weighted')
                im_weight_dir = os.path.join(weight_dir, imtype)
                wt_weight_dir = os.path.join(weight_dir, wttype)
                for outdir in [weight_dir, im_weight_dir, wt_weight_dir]:
                    os.makedirs(outdir)

                weight_images(im_dir, wt_dir, weight_dir, im_weight_dir, wt_weight_dir, imtype=imtype, wttype=wttype)
                im_dir = im_weight_dir
                wt_dir = wt_weight_dir


            # CREATE THE METADATA TABLES NEEDED FOR COADDITION
            weight_table = create_table(wt_dir, dir_type='weights')
            weighted_table = create_table(im_dir, dir_type='int')


            # COADD THE REPROJECTED, WEIGHTED IMAGES AND THE WEIGHT IMAGES WITH THE REGULAR HEADER FILE
            penultimate_dir = os.path.join(gal_dir, 'large_mosaic')
            final_dir = os.path.join(gal_dir, 'mosaic')
            for outdir in [penultimate_dir, final_dir]:
                os.makedirs(outdir)
            
            coadd(gal_hdr.hdrfile, penultimate_dir, im_dir, output='int', add_type='mean')
            coadd(gal_hdr.hdrfile, penultimate_dir, wt_dir, output='weights', add_type='mean')


            # DIVIDE OUT THE WEIGHTS AND CONVERT TO MJY/SR
            imagefile, wtfile = finish_weight(penultimate_dir)


            #if convert_mjysr:
            #    convert_to_flux_final(imagefile, band, desired_pix_scale)

            
            # SUBTRACT OUT THE BACKGROUND
            rm_overall_bg = False
            if rm_overall_bg:
                remove_background(final_dir, imagefile, bg_reg_file)
            else:
                outfile = os.path.join(final_dir, 'final_mosaic.fits')
                shutil.copy(imagefile, outfile)

            # copy weights mosaic to final directory
            shutil.copy(wtfile, os.path.join(final_dir, 'weights_mosaic.fits'))


            # COPY MOSAIC FILES TO CUTOUTS DIRECTORY
            mosaic_file = os.path.join(final_dir, 'final_mosaic.fits')
            weight_file = os.path.join(final_dir, 'weights_mosaic.fits')

            newsuffs = ['.FITS', '_weight.FITS']
            oldfiles = [mosaic_file, weight_file]
            newfiles = ['_'.join([name, band]).upper() + s for s in newsuffs]

            for files in zip(oldfiles, newfiles):
                shutil.copy(files[0], os.path.join(_MOSAIC_DIR, files[1]))


            # REMOVE TEMP GALAXY DIRECTORY AND EXTRA FILES
            shutil.rmtree(gal_dir, ignore_errors=True)


            # NOTE TIME TO FINISH
            stop_time = time.time()
            total_time = (stop_time - start_time) / 60.


            # WRITE OUT THE NUMBER OF TILES THAT OVERLAP THE GIVEN GALAXY
            out_arr = [name, band.upper(), nfiles, np.around(total_time, 2)]
            with open(numbers_file, 'a') as nfile:
                nfile.write('{0: >10}'.format(out_arr[0]))
                nfile.write('{0: >6}'.format(out_arr[1]))
                nfile.write('{0: >6}'.format(out_arr[2]))
                nfile.write('{0: >6}'.format(out_arr[3]) + '\n')
                #nfile.write(name + ': ' + str(len(infiles)) + '\n')

        # SOMETHING WENT WRONG -- WRITE ERROR TO FILE
        except Exception as inst:
            me = sys.exc_info()[0]
            with open(problem_file, 'a') as myfile:
                myfile.write(name + ': ' + str(me) + ': '+str(inst)+'\n')
            shutil.rmtree(gal_dir, ignore_errors=True)

    return


def get_input(index, ind, data_dir, input_dir, hdr=None):
    """
    Gather the input files for creating mosaics and copy them into a temporary directory

    Parameters
    ----------
    index : np.array
        structured array from galbase FITS file
    ind : np.array(np.dtype(int))
        An array of indices into the index locating the correct files
    data_dir : str
        Path to location of raw data downloaded from GALEX (or other) server
    input_dir : str
        Path to newly created temporary directory for storing temp files used in mosaicing

    Returns
    -------
    len(input_files) : int
        The number of files that will go into the mosaic.
    """
    infiles = index[ind[0]]['fname']
    wtfiles = index[ind[0]]['rrhrfile']
    flgfiles = index[ind[0]]['flagfile']
    infiles = [os.path.join(data_dir, f) for f in infiles]
    wtfiles = [os.path.join(data_dir, f) for f in wtfiles]
    flgfiles = [os.path.join(data_dir, f) for f in flgfiles]

    for i, infile in enumerate(infiles):
        basename = os.path.basename(infile)
        new_in_file = os.path.join(input_dir, basename)
        os.symlink(infile, new_in_file)
        if hdr is not None:
            keyw = 'INFILE{}'.format(str(i+1).zfill(2))
            hdr.append2hdr(keyword=keyw, value=basename, ext=False)

    for wtfile in wtfiles:
        basename = os.path.basename(wtfile)
        new_wt_file = os.path.join(input_dir, basename)
        os.symlink(wtfile, new_wt_file)

    for flgfile in flgfiles:
        basename = os.path.basename(flgfile)
        new_flg_file = os.path.join(input_dir, basename)
        os.symlink(flgfile, new_flg_file)

    return len(infiles)
        

def mask_images(im_dir, wt_dir, im_masked_dir, wt_masked_dir, imtype='int', wttype='rrhr'):
    """
    Mask pixels in the input images

    Parameters
    ----------
    im_dir : str
        Path to directory containing the images
    wt_dir : str
        Path to directory containing the weights
    im_masked_dir : str
        Path to temp directory for this galaxy in which to store masked image files
    wt_masked_dir : str
        Path to temp directory for this galaxy in which to store masked weight files
    """
    #int_suff, rrhr_suff = '*_mjysr.fits', '*-rrhr.fits'
    int_suff, rrhr_suff = '*-{}.fits'.format(imtype), '*-{}.fits'.format(wttype)
    int_images = sorted(glob.glob(os.path.join(im_dir, int_suff)))
    rrhr_images = sorted(glob.glob(os.path.join(wt_dir, rrhr_suff)))

    for i in range(len(int_images)):
        image_infile = int_images[i]
        wt_infile = rrhr_images[i]

        image_outfile = os.path.join(im_masked_dir, os.path.basename(image_infile))
        wt_outfile = os.path.join(wt_masked_dir, os.path.basename(wt_infile))

        mask_galex(image_infile, wt_infile, out_intfile=image_outfile, out_wtfile=wt_outfile)


def mask_galex(intfile, wtfile, chip_rad=1400, chip_x0=1920, chip_y0=1920, out_intfile=None, out_wtfile=None):
    """
    The actual masking routine. Selects pixels that are close to the edges of the chips
    or that have bad values, and masks them.

    Parameters
    ----------
    intfile : str
        input image file
    wtfile : str
        input weight file
    chip_rad : int
        Radius of the GALEX chip to use. The actual radius of data is ~1500 pixels. There are known edge effects. (Default: 1400)
    chip_x0 : int
        Center of GALEX chip on the x-axis (Default: 1920)
    chip_y0 : int
        Center of GALEX chip on the y-axis (Default: 1920)
    out_intfile : str, optional
        Path to output, masked image file. If not included, it will default to replacing the input file name as
        '.fits' --> '_masked.fits' (Default: None)
    out_wtfile : str, optional
        Path to output, masked weight file. If not included, it will default to replacing the input file name as
        '.fits' --> '_masked.fits' (Default: None)
    """
    if out_intfile is None:
        out_intfile = intfile.replace('.fits', '_masked.fits')
    if out_wtfile is None:
        out_wtfile = wtfile.replace('.fits', '_masked.fits')

    if not os.path.exists(out_intfile):
        # read in the data
        data, hdr = astropy.io.fits.getdata(intfile, header=True)
        wt, whdr = astropy.io.fits.getdata(wtfile, header=True)

        # determine distance of each pixel from the center
        x = np.arange(data.shape[1]).reshape(1, -1) + 1
        y = np.arange(data.shape[0]).reshape(-1, 1) + 1
        r = np.sqrt((x - chip_x0)**2 + (y - chip_y0)**2)

        # make pixel selections for masking
        i = (r > chip_rad)
        j = (wt == -1.1e30)

        # mask pixels that meet eithr of the above criterion
        data = np.where(i | j, np.nan, data)  #0
        wt = np.where(i | j, np.nan, wt) #1e-20

        # write new data to file
        astropy.io.fits.writeto(out_intfile, data, hdr)
        astropy.io.fits.writeto(out_wtfile, wt, whdr)


def reproject_images(template_header, input_dir, reproj_dir, imtype, whole=True, exact=True, corners=True, img_list=None):
    """
    Reproject input images to a new WCS as given by a template header

    Parameters
    ----------
    template_header : ascii file
        ASCII file containing the WCS to which you want to reproject. This is what Montage requires.
    input_dir : str
        Path to directory containing input data
    reproj_imtype_dir : 
        Path to new directory for storing reprojected data
    imtype : str
        The type of image you are reprojecting; one of [int, rrhr]
    whole : bool, optional
        Montage argument: Force reprojection of whole images, even if they exceed the area of the FITS 
        header template (Default: True)
    exact : bool, optional
        Montage argument: Flag indicating output image should exactly match the FITS header template, 
        and not crop off blank pixels (Default: True)
    corners : bool, optional
        Montage argument: Adds 8 columns for the RA and Dec of the image corners to the output metadata table 
        (Default: True)
    img_list : list of strs, optional 
        Montage argument: only process files with names specified in table img_list, ignoring any other files
        in the directory. (Default: None)
    """

    # get image metadata from input images
    input_table = os.path.join(input_dir, imtype + '_input.tbl')
    montage.mImgtbl(input_dir, input_table, corners=corners, img_list=img_list)

    # Create reprojection directory, reproject, and get image metadata
    stats_table = os.path.join(reproj_dir, imtype+'_mProjExec_stats.log')
    montage.mProjExec(input_table, template_header, reproj_dir, stats_table, raw_dir=input_dir, 
                      whole=whole, exact=exact)
    reprojected_table = os.path.join(reproj_dir, imtype + '_reprojected.tbl')
    montage.mImgtbl(reproj_dir, reprojected_table, corners=corners)


def bg_model(reprojected_dir, bg_model_dir, diff_dir, corr_dir, template_header, im_type='int', level_only=True):
    """
    Model the background for the mosaiced image

    Parameters
    ----------
    reprojected_dir : str
        Path to temp directory containing reprojected images 
    bg_model_dir : str
        Path to directory inside gal_dir to hold the background modeling information
    diff_dir : str
        Path to directory inside bg_model_dir to hold the difference images
    corr_dir : str
        Path to directory inside bg_model_dir to hold the background corrected images
    template_header : ascii file
        Path to file containing the WCS to which we want to reproject our images
    level_only : bool, optional
        Montage argument: Adjust background levels only, don't try to fit the slope (Default: True)
    """
    # FIND OVERLAPS
    diff_dir = os.path.join(diff_dir, im_type)
    os.makedirs(diff_dir)
    reprojected_table = os.path.join(reprojected_dir, im_type + '_reprojected.tbl')
    diffs_table = os.path.join(diff_dir, 'differences.tbl')
    montage.mOverlaps(reprojected_table, diffs_table)

    # CALCULATE DIFFERENCES BETWEEN OVERLAPPING IMAGES
    montage.mDiffExec(diffs_table, template_header, diff_dir,
                      proj_dir=reprojected_dir)

    # BEST-FIT PLANE COEFFICIENTS
    fits_table = os.path.join(diff_dir, 'fits.tbl')
    montage.mFitExec(diffs_table, fits_table, diff_dir)

    # CALCULATE CORRECTIONS
    corr_dir = os.path.join(corr_dir, im_type)
    os.makedirs(corr_dir)
    corrections_table = os.path.join(corr_dir, 'corrections.tbl')
    montage.mBgModel(reprojected_table, fits_table, corrections_table,
                     level_only=level_only)

    # APPLY CORRECTIONS
    montage.mBgExec(reprojected_table, corrections_table, corr_dir,
                    proj_dir=reprojected_dir)


def weight_images(im_dir, wt_dir, weight_dir, im_weight_dir, wt_weight_dir, imtype='-int', wttype='-rrhr'):
    """
    Weight the input images by a set of weights images

    Parameters
    ----------
    im_dir : str
        Path to directory containing the images
    wt_dir : str
        Path to directory containing the weights
    weight_dir : str
        Path to directory for the newly weighted images
    im_weight_dir : str
        Path to subdirectory containing the weighted images
    wt_weight_dir : str
        Path to subdirectory containgn the weights images (same as before, they haven't changed)    
    """
    #im_suff, wt_suff = '*_mjysr.fits', '*-rrhr.fits'
    im_suff, wt_suff = '*-{}.fits'.format(imtype), '*-{}.fits'.format(wttype)
    imfiles = sorted(glob.glob(os.path.join(im_dir, im_suff)))
    wtfiles = sorted(glob.glob(os.path.join(wt_dir, wt_suff)))    

    # weight each image
    for i in range(len(imfiles)):
        # read in the data
        imfile = imfiles[i]
        wtfile = os.path.join(os.path.dirname(wtfiles[i]), os.path.basename(imfile).replace('-int', '-rrhr'))
        im, hdr = astropy.io.fits.getdata(imfile, header=True)
        rrhr, rrhrhdr = astropy.io.fits.getdata(wtfile, header=True)

        # weight the data by the exposure time
        wt = rrhr
        newim = im * wt

        # write data to new files and copy the *_area.fits files created by Montage to have the same naming convention
        newfile = os.path.join(im_weight_dir, os.path.basename(imfile))
        astropy.io.fits.writeto(newfile, newim, hdr)
        old_area_file = imfile.replace('.fits', '_area.fits')
        if os.path.exists(old_area_file):
            new_area_file = newfile.replace('.fits', '_area.fits')
            shutil.copy(old_area_file, new_area_file)

        weightfile = os.path.join(wt_weight_dir, os.path.basename(wtfile))
        astropy.io.fits.writeto(weightfile, wt, rrhrhdr)
        old_area_file = wtfile.replace('.fits', '_area.fits')
        if os.path.exists(old_area_file):
            new_area_file = weightfile.replace('.fits', '_area.fits')
            shutil.copy(old_area_file, new_area_file)


def create_table(in_dir, dir_type=None):
    """
    Create a metadata table using Montage for all the files in a given directory

    Parameters
    ----------
    in_dir : str
        Path to directory containing the files
    dir_type : str, optional
        type of file you are creating a table for, e.g., 'int, rrhr, wt' (Default: None)

    Returns
    -------
    reprojected_table 
        Path to the table containing the metadata
    """
    if dir_type is None:
        reprojected_table = os.path.join(in_dir, 'reprojected.tbl')
    else:
        reprojected_table = os.path.join(in_dir, dir_type + '_reprojected.tbl')
    montage.mImgtbl(in_dir, reprojected_table, corners=True)
    return reprojected_table


def coadd(template_header, output_dir, input_dir, output=None, add_type=None):
    """
    Coadd input images to create mosaic.

    Parameters
    ----------
    template_header : ascii file
        File containing new WCS
    output_dir : str
        Path to directory contianing the output image
    input_dir : str
        path to directory containing the input images
    output : str, optional
        Type of mosaic you're making: e.g., int, wt, count (Default: None)
    add_type : str, optional
        Montage argument -- type of coadding to perform (Default: None -- defaults to Montage's default)
    """
    img_dir = input_dir
    # output is either 'weights' or 'int'
    if output is None:
        reprojected_table = os.path.join(img_dir, 'reprojected.tbl')
        out_image = os.path.join(output_dir, 'mosaic.fits')
    else:
        reprojected_table = os.path.join(img_dir, output + '_reprojected.tbl')
        out_image = os.path.join(output_dir, output + '_mosaic.fits')
    montage.mAdd(reprojected_table, template_header, out_image, img_dir=img_dir, exact=True, type=add_type)


def finish_weight(output_dir, convert_mjysr=True, band='fuv', gal_hdr=None, pix_as=None):
    """
    Divide out the weights from the final image to get back to flux density units

    Parameters
    ----------
    output_dir : str
        Path to directory containing the output image

    Returns
    -------
    newfile : str
        Path to new, mosaiced file
    """
    image_file = os.path.join(output_dir, 'int_mosaic.fits')
    wt_file = os.path.join(output_dir, 'weights_mosaic.fits')
    
    im, hdr = astropy.io.fits.getdata(image_file, header=True)
    wt = astropy.io.fits.getdata(wt_file)
    newim = im / wt

    newfile = os.path.join(output_dir, 'image_mosaic.fits')
    astropy.io.fits.writeto(newfile, newim, hdr)

    return newfile, wt_file


def counts2jy_galex(counts, cal, pix_as):
    """
    Convert GALEX counts/s to MJy/sr

    Parameters
    ----------
    counts : float
        Array containing counts data to be converted
    cal : float
        Calibration value from counts to AB mag for desired band (FUV or NUV)
    pix_as : float
        Pixel scale in arcseconds

    Returns
    -------
    val : float
        Count rate converted to MJy/sr
    """
    # first convert to abmag
    abmag = -2.5 * np.log10(counts) + cal

    # then convert to Jy
    f_nu = 10**(abmag/-2.5) * 3631.

    # then to MJy
    f_nu *= 1e-6

    # then to MJy/sr
    pix_rad = np.radians(pix_as / 3600.) # pixel scale coverted from arcsec to radians
    pix_sr = pix_rad ** 2. # pixel scale converted from radians to steradians
    val = f_nu / pix_sr

    return val


def convert_to_flux_input(indir, outdir, band, pix_as, imtype='int'):
    infiles = sorted(glob.glob(os.path.join(indir, '*-{}.fits'.format(imtype))))
    for infile in infiles:
        data, hdr = astropy.io.fits.getdata(infile, header=True)
        newdata = counts2jy_galex(data, UV2AB[band.lower()], pix_as)
        hdr['BUNIT'] = 'MJY/SR'
        outfile = os.path.join(outdir, os.path.basename(infile))
        astropy.io.fits.writeto(outfile, newdata, hdr)


def convert_to_flux_final(mosaicfile, band, pix_as):
    data, hdr = astropy.io.fits.getdata(mosaicfile, header=True)
    newim = counts2jy_galex(data, UV2AB[band.lower()], pix_as)

    # APPEND UNIT INFORMATION TO NEW HEADER AND WRITE OUT HEADER FILE
    hdr['BUNIT'] = 'MJY/SR' #gal_hdr.append2hdr(keyword='BUNIT', value='MJY/SR', ext=False)

    newfile = os.path.join(os.path.dirname(mosaicfile), 'image_mosaic_mjysr.fits')
    astropy.io.fits.writeto(newfile, newim, hdr)



# ------------------ #
## UNUSED functions ##
# ------------------ #
def remove_background(final_dir, imfile, bgfile):
    """
    Remove a background from the mosaiced image

    Parameters
    ----------
    final_dir : str
        Path to directory that will contain the final image
    imfile : str
        Path to file from which you want to remove a background
    bgfile : str 
        Path to file from which to calculate a background
    """
    # read in the data
    data, hdr = astropy.io.fits.getdata(imfile, header=True)
    box_inds = read_bg_regfile(bgfile)

    # calculate a background from the background file
    allvals = []
    sample_means = []
    for box in box_inds:
        rectangle = zip(box[0::2], box[1::2])
        sample = get_bg_sample(data, hdr, rectangle)
        for s in sample:
            allvals.append(s)
        sample_mean = np.nanmean(sample)
        sample_means.append(sample_mean)
    this_mean = np.around(np.nanmean(sample_means), 8)

    # subtract the background from the data
    final_data = data - this_mean

    # write the subtracted background file to the header for future use
    hdr['BG'] = this_mean
    hdr['comment'] = 'Background has been subtracted.'

    # write out the background-subtracted file
    outfile = os.path.join(final_dir, 'final_mosaic.fits')
    astropy.io.fits.writeto(outfile, final_data, hdr)


def read_bg_regfile(regfile):
    """
    Read the background region file (ds9 .reg file) and return the regions in which to calculate a background

    Parameters
    ----------
    regfile : str
        Path to file containing the regions in which to calculate a background

    Returns
    -------
    box_list : list
        List of data points that define the region in which to calculate a background
    """
    f = open(regfile, 'r')
    boxes = f.readlines()
    f.close()
    box_list = []
    for b in boxes:
        this_box = []
        box = b.strip('polygon()\n').split(',')
        [this_box.append(int(np.around(float(bb), 0))) for bb in box]
        box_list.append(this_box)
    return box_list


def get_bg_sample(data, hdr, box):
    """
    Determine the background sample in each reigon

    Parameters
    ----------
    data : float
        Input data
    hdr : FITS WCS
        header of input file containing WCS information
    box : np.2darray
        Numpy 2D array containing x,y coordinates of the box in which to determine the background

    Returns
    -------
    sample : float
        The data that lie within the region in which you want to calculate the background
    """
    wcs = astropy.wcs.WCS(hdr, naxis=2)
    x, y = np.arange(data.shape[0]), np.arange(data.shape[1])
    X, Y = np.meshgrid(x, y, indexing='ij')
    xx, yy = X.flatten(), Y.flatten()
    pixels = np.array(zip(yy, xx))
    box_coords = box
    sel = Path(box_coords).contains_points(pixels)
    sample = data.flatten()[sel]
    return sample


