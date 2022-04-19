import sys, getopt

import numpy as np
import h5py

input_file = ''
output_file = ''
argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv,"hi:o:")
except getopt.GetoptError:
    print('ExportEPMDynamo.py -i <inputfile> -o <outputfile>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('ExportEPMDynamo.py -i <inputfile> -o <outputfile>')
        sys.exit()
    elif opt in ("-i"):
        input_file = arg
    elif opt in ("-o"):
        output_file = arg

out_params = ['E', 'q', 'Ra', 'Ro']

quicc_epm_params = dict({
    'ekman':'E',
    'prandtl':'Pr',
    'magnetic_prandtl':'Pm',
    'rayleigh':'Ra',
    'magnetic_ekman':'Ro',
    'roberts':'q',
    })

epm_quicc_params = dict({
    'E':'ekman',
    'Pr':'prandtl',
    'Pm':'magnetic_prandtl',
    'Ra':'rayleigh',
    'q':'roberts',
    'Ro':'magnetic_ekman'
    })

# load QuICC state
quicc_file = h5py.File(input_file, 'r')
# Get truncation
epm_trunc = dict({
    'N':quicc_file['truncation']['spectral']['dim1D'][()],
    'L':quicc_file['truncation']['spectral']['dim2D'][()],
    'M':quicc_file['truncation']['spectral']['dim3D'][()],
    'Mp':1
    })
# Get physical parameters
quicc_params = dict()
for k in quicc_epm_params:
    if k in quicc_file['physical']:
        quicc_params[k] = quicc_file['physical'][k][()]

epm_params = dict()
if 'E' in out_params:
    epm_params['E'] = quicc_params['ekman']

if 'q' in out_params:
    if 'magnetic_prandtl' in quicc_params:
        epm_params['q'] = quicc_params['magnetic_prandtl']/quicc_params['prandlt']
    else:
        epm_params['q'] = quicc_params['prandtl']

if 'Ra' in out_params:
    epm_params['Ra'] = quicc_params['rayleigh']

if 'Ro' in out_params:
    if 'magnetic_prandtl' in quicc_params:
        epm_params['Ro'] = quicc_params['ekman']/quicc_params['magnetic_prandtl']
    else:
        epm_params['Ro'] = quicc_params['ekman']

print('QuICC physical parameters:')
print(epm_params)
for k, v in epm_params.items():
    print('\t'+ k + ' = ' + str(v))

# Get run parameters
runT = quicc_file['run']['time'][()]
runDt = quicc_file['run']['timestep'][()]

# create QuICC state file
epm_file = h5py.File(output_file,'w')
epm_file.attrs.create('StateFile', data=0, dtype='i4')
epm_file.attrs.create('Version', data=1, dtype='f4')

# Create run group
group = epm_file.create_group('RunParameters')
group.create_dataset('Time', (), 'f8', data = runT)
group.create_dataset('Step', (), 'f8', data = runDt)

# Create physical group
group = epm_file.create_group('PhysicalParameters')
for k, v in epm_params.items():
    group.create_dataset(k, (), 'f8', data = v)

# Create truncation group
group = epm_file.create_group('Truncation')
for k, v in epm_trunc.items():
    group.create_dataset(k, (), 'i4', data = v)

def rescale(d, l, m):
    out = np.zeros(d.shape)
    lm = -1
    for k in range(0, l+1):
        norm0 = 1.0/ np.sqrt(4.0*np.pi/(2.0*k + 1.0))
        normm = 1.0/np.sqrt(8.0*np.pi/(2.0*k + 1.0))
        for j in range(0, k+1):
            lm += 1
            if j == 0:
                norm = norm0
            else:
                norm = normm
            out[lm,:,:] = norm*d[lm,:,:]
    return out

def add_background(d):
    print("### Adding thermal background state ###")
    # Background state is: T = 1/2 (1 - r^2)
    # In normalized spherical harmonic space: T = 1/2 (1 - r^2) (2 sqrt(pi))
    # n = 0, W_0^0 = 1.0, norm: pi/2
    # n = 1 W_1^0 = r^2 - 1/2, norm: pi/16
    # t0 = 1.110720734539592
    # t1 = -0.7853981633974483
    import scipy.integrate as integrate
    t0 = integrate.quad(lambda r:0.5*np.sqrt(1.0-r**2),0,1,epsabs=-1,epsrel=50*np.finfo(float).eps)[0]
    t0 = t0*8.0**0.5
    t1 = integrate.quad(lambda r:0.5*np.sqrt(1.0-r**2)*(r**2 - 1.0/2.0),0,1,epsabs=-1,epsrel=70*np.finfo(float).eps)[0]
    t1 = t1*8.0
    d[0,0,0] += t0
    d[0,1,0] += t1
    return d


# Create temperature field
if 'temperature' in quicc_file:
    group = epm_file.create_group('Codensity')
    quicc_data = quicc_file['temperature']['temperature'][:]
    ds = group.create_dataset('Codensity', shape=quicc_data.shape[0:2], dtype = ('<f8', (2,)))
    ds[:] = rescale(add_background(quicc_data), epm_trunc['L'], epm_trunc['M'])

## Create velocity field
if 'velocity' in quicc_file:
    group = epm_file.create_group('Velocity')
    quicc_data = quicc_file['velocity']['velocity_tor'][:]
    ds = group.create_dataset('VelocityTor', shape=quicc_data.shape[0:2], dtype = ('<f8', (2,)))
    ds[:] = rescale(quicc_data, epm_trunc['L'], epm_trunc['M'])
    quicc_data = quicc_file['velocity']['velocity_pol'][:]
    ds = group.create_dataset('VelocityPol', shape=quicc_data.shape[0:2], dtype = ('<f8', (2,)))
    ds[:] = rescale(quicc_data, epm_trunc['L'], epm_trunc['M'])

## Create magnetic field
if 'magnetic' in quicc_file:
    group = epm_file.create_group('Magnetic')
    quicc_data = quicc_file['magnetic']['magnetic_tor'][:]
    ds = group.create_dataset('MagneticTor', shape=quicc_data.shape[0:2], dtype = ('<f8', (2,)))
    ds[:] = rescale(quicc_data, epm_trunc['L'], epm_trunc['M'])
    quicc_data = quicc_file['magnetic']['magnetic_pol'][:]
    ds = group.create_dataset('MagneticPol', shape=quicc_data.shape[0:2], dtype = ('<f8', (2,)))
    ds[:] = rescale(quicc_data, epm_trunc['L'], epm_trunc['M'])

# Finished
epm_file.close()
quicc_file.close()
