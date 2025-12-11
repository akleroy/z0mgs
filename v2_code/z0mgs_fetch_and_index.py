# Wrapper function to handle fetching and indexing of underlying data
# from z0mgs work. The core routines are in utils_fetch_and_index.

# Imports
from utils_fetch_and_index import *

# TBD could build this out to a fancy command line argument, but these
# need to be run only a limited number of times so this is probably
# fine.

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# What task?
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ... FETCH data from the web
# ... CHECK whether expected data are present, fetch missing
# ... FLIST compile the list of files to be indexed
# ... INDEX the coordinates and details of each file

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# What survey?
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ... UNWISE (has three subsurveys)

# ... ... CUSTOM_UNWISE calls for individual galaxies from 2022
# ... ... (these are important for WISE3 and 4 no-median cases)

# ... ... ALLWISE initial cool whole-sky survey with WISE 3 and 4

# ... ... NEOWISE through release 9 with whole sky WISE 1 and 2

# ... SDSS DR12 imaging

# ... GALEX whole sky tiles

# ... GAIA DR3 whole sky catalog

do_fetch = False
do_check = False
do_flist = False
do_index = False

do_unwise = False
do_sdss = False
do_galex = False
do_gaia = False

do_custom_unwise = False
do_allwise = False
do_neowise = False

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Commands
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Fetching these surveys tends to be a large job, so most of these
# calls will just print a short WGET or RSYNC command to be run by the
# user.

if do_fetch:

    if do_unwise:
        if do_unwise:
            unwise_fetch(dry_run=True, incremental=False, version='neo9')
        if do_allwise:
            unwise_fetch(dry_run=True, incremental=False, version='allwise')

    if do_galex:
        print("Not implemented, a GALEX WGET script is checked in.")
        galex_fetch(dry_run=True, incremental=False)

    if do_sdss:
        sdss_fetch(dry_run=True)

    if do_gaia:
        print("Not implemented, GAIA DR3 obtained via GLOBUS transfer manually.")
        gaia_fetch(dry_run=True)

# Check that we successfully pulled down all of the relevant files.
        
if do_check:

    # UNWISE provides file lists, we will check those and wget anything missing.
    
    if do_unwise:
        if do_neowise:
            unwise_fetch(dry_run=False, incremental=True, version='neo9')
        if do_allwise:
            unwise_fetch(dry_run=False, incremental=True, version='allwise')

    # GALEX to be implemented
            
    # Gaia GLOBUS transfer is assumed to be correct

    # SDSS rsync is assumed to be correct

# These obtain and write lists of files in a table for each
# survey. That table is the basis for the indexing in the next step.
    
if do_flist:
    
    if do_unwise:
        if do_custom_unwise:
            test = compile_list_of_images(survey='unwise_custom')
        if do_allwise:
            test = compile_list_of_images(survey='unwise_allwise')
        if do_neowise:
            test = compile_list_of_images(survey='unwise_neowise')

    if do_galex:
        test = compile_list_of_images(survey='galex')

    if do_sdss:
        test = compile_list_of_images(survey='sdss')
        
    if do_gaia:
        test = compile_list_of_images(survey='gaia')

# These read the files and record their astrometry, filter, and some
# other metadata. They take O(1 hour - 1 day) to run for GALEX and
# UNWISE. For SDSS this takes a week (even indexing only one band) and
# so the call is run with an incremental tag that writes results every
# 1k or so galaxies and can be resumed.
        
if do_index:
    
    if do_unwise:
        if do_custom_unwise:
            test = index_image_list(survey='unwise_custom')
        if do_allwise:
            test = index_image_list(survey='unwise_allwise')
        if do_neowise:
            test = index_image_list(survey='unwise_neowise')

    if do_galex:
        test = index_image_list(survey='galex')

    if do_sdss:
        #test = index_image_list(survey='sdss')
        test = index_image_list(survey='sdss', incremental=True)

    if do_gaia:        
        test = index_image_list(survey='gaia')

