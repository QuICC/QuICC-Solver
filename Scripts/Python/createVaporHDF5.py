#!/usr/bin/env python

from __future__ import print_function

import sys, getopt

from math import *
import h5py
import numpy as np
import re

def main(argv):
    inputfile = ''
    outputfile = ''
    snapshots = 1 
    try:
        opts, args = getopt.getopt(argv,"hi:o:n:")
    except getopt.GetoptError:
        print('Single file: createXDMF.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Single file: createXDMF.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i"):
            inputfile = arg
        elif opt in ("-o"):
            outputfile = arg
        elif opt in ("-n"):
            snapshots = int(arg)
    # Extract file information
    h5_infile = h5py.File(inputfile, 'r')
    scheme = h5_infile['/'].attrs['type']
    g1D = 'x'
    g2D = 'y'
    g3D = 'z'
    n1D = h5_infile['mesh']['grid_'+g1D].size
    n2D = h5_infile['mesh']['grid_'+g2D].size
    n3D = h5_infile['mesh']['grid_'+g3D].size
    time = h5_infile['run']['time'].value
    sId = int(re.findall('\d+', inputfile.split('.')[0])[0])
    basename = inputfile.split(re.findall('\d+', inputfile.split('.')[0])[0])[0]
    # Set default output file
    if outputfile == '' and snapshots == 1:
        outputfile = basename+str(sId).zfill(4)+'_vapor.hdf5' 
    print("Input file: ", inputfile)
    print("Input ID: ", sId)
    print("Input base: ", basename)
    print('Output file: ', outputfile)
    print("Scheme: ", scheme)
    print("1D grid size: ", n1D)
    print("2D grid size: ", n2D)
    print("3D grid size: ", n3D)
    print("Time: ", time)

    # Create output file
    h5_outfile = h5py.File(outputfile,'w')

    # Copy general information (not used by VAPOR)
    #h5_outfile.copy(h5_infile['mesh'], '/mesh')
    #h5_outfile.copy(h5_infile['physical'], '/physical')
    #h5_outfile.copy(h5_infile['run'], '/run')
    #h5_outfile.copy(h5_infile['truncation'], '/truncation')

    # Create mesh files
    for ext in [g1D, g2D, g3D]:
        gridfile = basename+str(sId).zfill(4)+'_vapor_' + ext + '.dat'
        np.savetxt(gridfile, np.sort(h5_infile['mesh']['grid_' + ext]))

    # Copy and flatten scalars
    for s in list(h5_infile):
        if s in list(h5_infile[s]):
            h5_outfile.copy(h5_infile[s][s], s)

    # Copy and flatten vectors
    for v in list(h5_infile):
        for ext in [g1D, g2D, g3D]:
            if v +  '_' + ext in list(h5_infile[v]):
                h5_outfile.copy(h5_infile[v][v + '_' + ext], v + '_' + ext)

    # Close hdf5 files
    h5_outfile.close()
    h5_infile.close()

if __name__ == "__main__":
   main(sys.argv[1:])
