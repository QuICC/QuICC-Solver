import h5py
import numpy as np

import getopt, sys

args = sys.argv[1:]
options = "hi:"

fname = "filenotset.hdf5"

try:
    arguments, values = getopt.getopt(args, options)
    for currentArg, currentVal in arguments:
        if currentArg in ("-h"):
            print("Usage: process_profile -i profile.hdf5")
        elif currentArg in ("-i"):
            fname = currentVal
except getopt.error as err:
    print(str(err))

h5file = h5py.File(fname, 'r')

refImpl = 'Flat'
variants = [20]
impls = [refImpl, 'View']
ops = ['Set', 'Add', 'Sub']
ds = [
        'Cross',
        'Dot',
        'MassFluxAdvection',
        'SphericalBuoyancy',
        'SphericalCoriolisAnelastic',
        'SphericalCoriolis',
        'SphericalHeatAdvection',
        'SphericalLorentzAnelastic',
        'SphericalOhmicDissipationAnelastic',
        'SphericalPoincare',
        'SphericalPrecession',
        'SphericalSComponent',
        'SphericalSelfAdvectionAnelastic',
        'SphericalViscousDissipationAnelastic',
        'SphericalZComponent',
        'StreamAdvection',
        'StreamHeatAdvection',
        'VelocityAdvection',
        'VelocityHeatAdvection',
     ]


sp = '  '
for dbase in ds:
    print(f'{dbase}:')
    for op in ops:
        print(f'{sp}{op}:')
        for variant in variants:
            print(f'{sp}{sp}variant {variant}:')
            for impl in impls:
                dname = f'{dbase}Tests::{impl}{op}_{variant}'
                if dname in h5file['timings'].keys():
                    ts = h5file['timings'][dname]['time'][()]
                    min = np.min(ts)
                    avg = np.average(ts)
                    max = np.max(ts)
                    if impl == refImpl:
                        refAvg = avg
                        speedup = 0.0
                    else:
                        speedup = 100.0*(avg - refAvg)/refAvg
                    print(f'{sp}{sp}{sp}{impl}: {min:.2E}/{max:.2E}/{avg:.2E} ({np.ceil(speedup):+} %)')
                else:
                    print(f'{sp}{sp}{sp}{impl}: NOT IMPLEMENTED')

