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
        print('convert_epm_physical.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('convert_epm_physical.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i"):
            inputfile = arg
        elif opt in ("-o"):
            outputfile = arg
    # Extract file information
    epm_file = h5py.File(inputfile, 'r')
    quicc_file = h5py.File(outputfile, 'w')

    # Create header
    quicc_file.attrs['header'] = np.string_('StateFile'.encode('ascii')) 
    quicc_file.attrs['type'] = np.string_('WLFl'.encode('ascii'))
    quicc_file.attrs['version'] = np.string_('1.0'.encode('ascii'))

    # Create run group
    group = quicc_file.create_group('run')
    group.create_dataset('time', (), 'f8', data = 0)
    group.create_dataset('timestep', (), 'f8', data = 0)

    # Create mesh
    group = quicc_file.create_group('mesh')
    group.create_dataset('grid_r', data = epm_file['FSOuter']['grid']['r axis'])
    group.create_dataset('grid_theta', data = epm_file['FSOuter']['grid']['t axis'])
    group.create_dataset('grid_phi', data = epm_file['FSOuter']['grid']['p axis'])

    # Create field datasets
    for g in epm_file['FSOuter']:
        if g != 'grid':
            group = quicc_file.create_group(g)
            if g == 'Codensity':
                group.create_dataset(g, data = epm_file['FSOuter'][g])
            else:
                for d in epm_file['FSOuter'][g]:
                    ext = ''
                    if d == g+'r':
                        ext = '_r'
                    elif d == g+'t':
                        ext = '_theta'
                    elif d == g+'p':
                        ext = '_phi'

                    group.create_dataset(g+ext, data = epm_file['FSOuter'][g][d])

    # Finished
    epm_file.close()
    quicc_file.close()

if __name__ == "__main__":
    main(sys.argv[1:])
