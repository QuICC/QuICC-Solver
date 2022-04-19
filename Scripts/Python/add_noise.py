import h5py
import numpy as np
import sys, os


if __name__=="__main__":

    try:
        argv = sys.argv
        filename = argv[1]
        err_magnitude = float(argv[2])
    except:
        print('Supposed usage: python add_noise.py filename noise_size')
        sys.exit()

    fin = h5py.File(filename,'r+')

    for fl in ['/velocity/velocity_tor','/velocity/velocity_pol']:

        field = fin[fl]

        field[:,:]+=np.random.random(field[:,:].shape)*err_magnitude
        
        #print(field[:,:])

    #fin['/velocity/velocity_tor'][:,:]
    fin.close()
        
