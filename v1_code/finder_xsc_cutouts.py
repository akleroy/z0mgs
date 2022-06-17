#! /usr/bin/env python

# This creates URLs for 2MASS XSC Image Cutouts covering a point.
# It takes a .csv table containing ra,dec. The output is curl calls
# in a file named get_xsc.csh.  When executed, that will dump the fits files in the 
# current directory.  Need write permission, of course.
#
#   ./finder_xsc_cutouts.py table_radec.csv
#   ./get_xsc.csh

import os
import sys
import re
import subprocess
import errno
import argparse
import csv

def path_safe(path):
        original_umask = os.umask(0)
        try:
            os.makedirs(path, mode=0o775)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        os.umask(original_umask)
        
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("filecsv", help="CSV table with ra and dec (deg, J2000)")
    args = parser.parse_args()

    filecsv = args.filecsv

# Check that table is valid 

    with open(filecsv) as csvfile:
        reader = csv.DictReader(csvfile)
        i = 0
        for line in reader:
            ra = line['ra']
            try:
                float(ra)
            except ValueError:
                print('')
                print "ra not a float at line ", i, "with value = ", ra
                exit()
            if not (float(ra) >= 0.0 and float(ra) <= 360.0):
                print('')
                print "Bad or missing ra line ", i, "with value = ", ra
                exit()
            dec = line['dec']
            try:
                float(dec)
            except ValueError:
                print('')
                print "dec not a float at line ", i, "with value = ", dec
                exit()
            if not (float(dec) >= -90.0 and float(dec) <= 90.0):
                print('')
                print "Bad or missing dec line ", i, "with value = ", dec
                exit()
            i += 1

# 2MASS XSC Image cutouts

    curlname = os.getcwd() + '/get_xsc.csh'
    fh1 = open(curlname,'w')

    htmlname = os.getcwd() + '/temp.html'
    pattern = re.compile(r"/workspace.*?fits")
    with open(filecsv) as csvfile:
        reader = csv.DictReader(csvfile)
        for line in reader:
            ra = line['ra']
            dec = line['dec']
            irsa_call = '\"https://irsa.ipac.caltech.edu/cgi-bin/2MASS/PubGalPS/nph-galps?locstr=' + ra + '+' + dec + '&radius=3\"'
            os.system("curl -s -o %s %s" % (htmlname, irsa_call))
	    with open(htmlname) as htmlfile:
                for line in htmlfile:
                    if pattern.search(line) != None:
                         fh1.write("curl -O -s https://irsa.ipac.caltech.edu%s\n" % (pattern.findall(line)[0]))            
    fh1.close()
    os.remove(htmlname)

    print('')
    print('Done.')

if __name__ == '__main__':
    main()
