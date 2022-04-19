#!/usr/bin/env python

from __future__ import print_function

import sys, getopt

import h5py
import numpy as np

def main(argv):
    inputfile = ''
    outputfile = ''

    try:
        opts, args = getopt.getopt(argv,"hi:o:")
    except getopt.GetoptError:
        print('convert_to_epm_physical.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('convert_to_epm_physical.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i"):
            inputfile = arg
        elif opt in ("-o"):
            outputfile = arg
    # Extract file information
    quicc_file = h5py.File(inputfile, 'r')
    epm_file = h5py.File(outputfile, 'w')

    # Create OuterCore group
    root = epm_file.create_group('OuterCore')

    # Create mesh
    group = root.create_group('grid')
    group.create_dataset('r axis', data = quicc_file['mesh']['grid_r'])
    group.create_dataset('t axis', data = quicc_file['mesh']['grid_theta'])
    group.create_dataset('p axis', data = quicc_file['mesh']['grid_phi'])

    # Create field datasets
    for g in quicc_file:
        if g not in ['mesh','run','physical','truncation']:
            print("Converting " + g)
            makeGroup = True
            for d in quicc_file[g]:
                ext = ''
                if d == g+'_r':
                    ext = 'r'
                elif d == g+'_theta':
                    ext = 't'
                elif d == g+'_phi':
                    ext = 'p'

                if ext != '':
                    if makeGroup:
                        group = root.create_group(g)
                        makeGroup = False
                else:
                    group = root
                group.create_dataset(g+ext, data = quicc_file[g][d])

    # Finished
    epm_file.close()
    quicc_file.close()

if __name__ == "__main__":
    main(sys.argv[1:])
