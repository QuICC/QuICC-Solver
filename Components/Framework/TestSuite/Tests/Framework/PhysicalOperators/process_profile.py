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
variants = [30]
impls = [refImpl, 'View']
ops = ['Set', 'Add', 'Sub', 'CSet', 'CAdd', 'CSub']
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
                        refMin = min
                        refAvg = avg
                        refMax = max
                        speedInfo =' '
                    else:
                        speedup = 100.0*(min - refMin)/refMin
                        minSpeedInfo = f'({np.ceil(speedup):+} %)'
                        speedup = 100.0*(avg - refAvg)/refAvg
                        avgSpeedInfo = f'({np.ceil(speedup):+} %)'
                        speedup = 100.0*(max - refMax)/refMax
                        maxSpeedInfo = f'({np.ceil(speedup):+} %)'
                        speedInfo = f' | {minSpeedInfo}/{maxSpeedInfo}/{avgSpeedInfo}'
                    print(f'{sp}{sp}{sp}{impl}: {min:.2E}/{max:.2E}/{avg:.2E}{speedInfo}')
                else:
                    print(f'{sp}{sp}{sp}{impl}: NOT IMPLEMENTED')

