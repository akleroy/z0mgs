#!/usr/bin/env python3

import os, glob, sys
import numpy as np
import warnings

import urllib.request
import urllib.error

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy import units as u
from astropy.stats import sigma_clipped_stats
import astropy.constants as const

# Astroquery for IRSA access
from astroquery.ipac.irsa import Irsa

# Reproject for image alignment
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to query data from IRSA
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

def check_spherex_collections():
    """Query IRSA's collection list to see which collections hold spherex
    data. Use this to inform the image search.
    """
    
    collections = Irsa.list_collections(servicetype='SIA', filter='spherex')

    print(f"Found {len(collections)} SPHEREX collections:")
    for collection in collections['collection']:
        print(f"  - {collection}")

    return()


def search_spherex_images(
        target: str = None,
        coordinates: SkyCoord = None,
        radius: u.Quantity = 10*u.arcmin,
        collection: str = 'spherex_qr2',
        verbose: bool = True) -> Table:
    
    """
    From Adam Ginsburg (pared down a lot)

    Search for SPHEREX images in IRSA archive.

    Parameters
    ----------
        
    target : str, optional
    Target name (e.g., 'M31', 'NGC 1234')
    
    coordinates : SkyCoord, optional
    Target coordinates

    radius : Quantity
    Search radius

    collection : str (default 'spherex_qr2')
    Collection holding the images, will need to be managed    

    Returns
    -------
    Table
    
    Table of available SPHEREX images
    """
    
    if verbose:
        print(f"Searching for SPHEREX images...")

    # Make sure we have target coordinates
    if coordinates is None and target is not None:
        coordinates = SkyCoord.from_name(target)
    elif coordinates is None:
        raise ValueError("Either target name or coordinates must be provided")

    if verbose:
        print(f"Target coordinates: {coordinates.to_string('hmsdms')}")
        print(f"Search radius: {radius}")

    # Search for images in the specified collection
    try:
        images = Irsa.query_sia(pos=(coordinates, radius), collection=collection)
    except Exception as e:
        if verbose:
            print(f"SIA query failed: {e}")

    if verbose:
        print(f"Found {len(images)} images")
        if len(images) > 0:
            print("Available columns (first 10): ", images.colnames[:10])  # Show first 10 columns

    return(images)

def download_images(
        images: Table,
        max_images: int = None,
        outdir: str = '',
        incremental: bool = True,
        verbose: bool = True) -> list[str]:
    """Download SPHEREX images from a table produced by the IRSA query above.

    Modified from version by Adam Ginsburg.
    
    Parameters
    ----------
    
    images : Table
             Table of images from search_spherex_images
        
    max_images : int, optional
                 Maximum number of images to download

    outdir : str, optional
             Location where the download directs

    incremental : bool, default True
                  Check if file is present before download

    verbose : bool, default True
              Turn on extra print statements
    
    Returns
    -------
    
    list[str] : List of downloaded file paths

    """

    # If an empty table, return
    if len(images) == 0:
        if verbose:
            print("No images to download")
        return([])

    # If we're capped in files to download impose that
    if max_images is not None:
        images = images[:max_images]

    # Initialize record
    downloaded_files = []

    # Loop and download
    if verbose:
        print(f"Downloading {len(images)} images...")
        if len(images) > 10:
            print("  (This may take a while for large datasets...)")

    # Loop over the list of images
    for ii, this_row in enumerate(images):

        # Determine the access URL        
        try:
            this_url = None
            # ... try the different potential column names
            for url_col in ['access_url', 'cloud_access', 'download_url', 'url']:
                if url_col in this_row.colnames and this_row[url_col]:
                    this_url = str(this_row[url_col])
                    break

            # ... catch the error case
            if this_url is None:
                if verbose:
                    print(f"  Skipping image {ii+1}: No access URL found")
                continue

            # Generate filename

            # AKL: Note that you cannot just use the obsid because two
            # images with different filters are obtained with the same
            # obsid - parse the full file name.

            # Observation ID + detector should be unique, but
            # processing time in the ID stamp is also potentially of
            # interest.
           
            #obs_id = this_row.get('obs_id', f'spherex_{ii:04d}')
            
            obs_fname = this_url.split('/')[-1]            
            this_fname = outdir + obs_fname

            # Check if the filename is present already and either delete
            # the file or proceed with a notification.
            if os.path.isfile(this_fname):
                if incremental:
                    if verbose:
                        print(f"  Image {ii+1}/{len(images)}: {this_fname} already exists.")
                        print(f"  Incremental is set to TRUE.")
                        downloaded_files.append(str(this_fname))
                    continue
                else:
                    if verbose:
                        print(f"  Image {ii+1}/{len(images)}: {this_fname} already exists.")
                        print(f"  Incremental is set to FALSE. Will proceed and overwrite.")
                    os.system('rm -rf '+this_fname)

            # Download the file
            if verbose:
                print(f"  Downloading image {ii+1}/{len(images)}: {this_fname}")

            # Download using urllib first, then open with astropy
            try:                
                # Download to a temporary location first
                temp_fname = this_fname + '.tmp'
                urllib.request.urlretrieve(this_url, temp_fname)
            
                # Verify it's a valid FITS file by opening it
                with fits.open(temp_fname) as this_hdu:                
                    # Save as final file
                    this_hdu.writeto(this_fname, overwrite=True)

                # Remove temporary file
                os.system('rm -rf '+temp_fname)

            except urllib.error.URLError as e:
                if verbose:
                    print(f"    URL error: {e}")
                continue            
            except Exception as fits_error:            
                # Clean up temp file if it exists
                temp_fname = this_fname + '.tmp'
                if os.path.isfile(temp_fname):
                    os.system('rm -rf '+temp_fname)
                    # Will crash the program?
                raise fits_error

            # If we made it here we're successful
            downloaded_files.append(str(this_fname))
        
        except Exception as e:
            if verbose:
                print(f"  Failed to download image {ii+1}: {e}")
            continue

    if verbose:
        print(f"Successfully downloaded {len(downloaded_files)} images")

    return(downloaded_files)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Routines to build a cube
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

