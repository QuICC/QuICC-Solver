import h5py
import numpy as np
import sys, os


if __name__=="__main__":

    try:
        argv = sys.argv
        filename_inout = argv[1]
        filename_in = argv[2]
        assert sys.argv.__len__() == 3
    except:
        print('Supposed usage: python add_noise.py filename_inout filename_in')
        sys.exit()

    # open files
    finout = h5py.File(filename_inout,'r+')
    fin = h5py.File(filename_in,'r')

    for fl in ['/velocity/velocity_tor','/velocity/velocity_pol']:
        field1 = fin[fl]
        # copy the fields from in to out
        field2 = finout[fl]
        field2[:,:] = field1[:,:]

    finout.close()
    fin.close()
