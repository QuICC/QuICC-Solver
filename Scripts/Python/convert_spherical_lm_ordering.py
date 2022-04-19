import sys, getopt

import numpy as np
import h5py

input_filename = ''
output_filename = ''
argv = sys.argv[1:]

# get command line arguments
try:
    opts, args = getopt.getopt(argv,"hi:o:")
except getopt.GetoptError:
    print('convert_l_to_m.py -i <inputfile> -o <outputfile>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('convert_l_to_m.py -i <inputfile> -o <outputfile>')
        sys.exit()
    elif opt in ("-i"):
        input_filename = arg
    elif opt in ("-o"):
        output_filename = arg

# load input state
in_file = h5py.File(input_filename, 'r')
# create output state
out_file = h5py.File(output_filename, 'w')

# get mode ordering of input
in_ordering = in_file.attrs['type'].decode('ascii')[-1]
if in_ordering == 'm':
    out_ordering = 'l'
    print("Input file has M ordering")
else:
    out_ordering = 'm'
    print("Input file has L ordering")

# copy file attributes
for key, value in in_file.attrs.items():
    if key == "type":
        out_file.attrs[key] = np.string_((value.decode('ascii')[:-1]+out_ordering).encode('ascii'))
    else:
        out_file.attrs[key] = np.string_(value)

# copy all groups to output file
for g in in_file:
    in_file.copy(g, out_file)

# get truncation
resN = out_file['truncation']['spectral']['dim1D'][()]
resL = out_file['truncation']['spectral']['dim2D'][()]
resM = out_file['truncation']['spectral']['dim3D'][()]

l_idx = -1
ls = []
ms = []
for l in range(0,resL+1):
    for m in range(0,min(l+1,resM+1)):
        l_idx += 1
        m_idx = -1
        for mm in range(0, m):
            for ll in range(mm, resL+1):
                m_idx += 1
        for ll in range(m, l+1):
            m_idx += 1

        # set input and output file row index
        if out_ordering == 'm':
            out_idx = m_idx
            in_idx = l_idx
        else:
            out_idx = l_idx
            in_idx = m_idx

        # copy input spectral coefficients to new ordering in output
        for g in in_file:
            dataset_group = '/'+g + '/'+ g
            # scalar dataset
            if dataset_group in in_file:
                out_file[dataset_group][out_idx,:] = in_file[dataset_group][in_idx,:]
            # toroidal/poloidal dataset
            if (dataset_group + '_tor') in in_file:
                for c in ['_tor', '_pol']:
                    out_file[dataset_group + c][out_idx,:] = in_file[dataset_group + c][in_idx,:]

in_file.close()
out_file.close()
