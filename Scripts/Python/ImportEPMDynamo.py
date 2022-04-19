import sys, getopt

import numpy as np
import h5py

input_file = ''
output_file = ''
argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv,"hi:o:")
except getopt.GetoptError:
    print('ImportEPMDynamo.py -i <inputfile> -o <outputfile>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('ImportEPMDynamo.py -i <inputfile> -o <outputfile>')
        sys.exit()
    elif opt in ("-i"):
        input_file = arg
    elif opt in ("-o"):
        output_file = arg

# Missing data flags
mdv = -4242.4242
mda = np.array(mdv)

# load EPMDynamo state
epm_file = h5py.File(input_file, 'r')
# Get truncation
spec1D = epm_file['Truncation']['N'][()]
spec2D = epm_file['Truncation']['L'][()]
spec3D = epm_file['Truncation']['M'][()]
trans1D = 3*(spec1D+1 + spec2D//2 + 8)//2
trans2D = 3*(spec2D+1)//2
trans3D =  3*(spec3D+1)
phys1D = 3*(spec1D+1 + spec2D//2 + 1)//2
phys2D = 3*(spec2D+1)//2
phys3D =  3*(spec3D+1)
# Get physical parameters
physE = epm_file['PhysicalParameters'].get('E',default=mda)[()]
physRa = epm_file['PhysicalParameters'].get('Ra',default=mda)[()]
physRo = epm_file['PhysicalParameters'].get('Ro',default=mda)[()]
physQ = epm_file['PhysicalParameters'].get('q',default=mda)[()]
physPm = physE/physRo
physPr = physPm/physQ
print('QuICC physical parameters:')
print('\t'+'E = ' + str(physE))
print('\t'+'Ra = ' + str(physRa))
print('\t'+'Pr = ' + str(physPr))
print('\t'+'Pm = ' + str(physPm))
# Get run parameters
runT = epm_file['RunParameters']['Time'][()]
runDt = epm_file['RunParameters']['Step'][()]

# create QuICC state file
quicc_file = h5py.File(output_file,'w')
quicc_file.attrs['header'] = np.string_('StateFile'.encode('ascii')) 
quicc_file.attrs['type'] = np.string_('WLFl'.encode('ascii'))
quicc_file.attrs['version'] = np.string_('1.0'.encode('ascii'))

# Create run group
group = quicc_file.create_group('run')
group.create_dataset('time', (), 'f8', data = runT)
group.create_dataset('timestep', (), 'f8', data = runDt)

# Create physical group
group = quicc_file.create_group('physical')
group.create_dataset('bc_magnetic', (), 'f8', data = 0)
group.create_dataset('bc_temperature', (), 'f8', data = 0)
group.create_dataset('bc_velocity', (), 'f8', data = 0)
group.create_dataset('ekman', (), 'f8', data = physE)
group.create_dataset('magnetic_prandtl', (), 'f8', data = physPm)
group.create_dataset('prandtl', (), 'f8', data = physPr)
group.create_dataset('rayleigh', (), 'f8', data = physRa)

# Create truncation group
group = quicc_file.create_group('truncation')
subgroup = group.create_group('physical')
subgroup.create_dataset('dim1D', (), 'i4', data = phys1D)
subgroup.create_dataset('dim2D', (), 'i4', data = phys2D)
subgroup.create_dataset('dim3D', (), 'i4', data = phys3D)
subgroup = group.create_group('transform')
subgroup.create_dataset('dim1D', (), 'i4', data = trans1D)
subgroup.create_dataset('dim2D', (), 'i4', data = trans2D)
subgroup.create_dataset('dim3D', (), 'i4', data = trans3D)
subgroup = group.create_group('spectral')
subgroup.create_dataset('dim1D', (), 'i4', data = spec1D)
subgroup.create_dataset('dim2D', (), 'i4', data = spec2D)
subgroup.create_dataset('dim3D', (), 'i4', data = spec3D)

def rescale(d, l, m):
    out = np.zeros(d.shape)
    lm = -1
    for k in range(0, l+1):
        norm0 = np.sqrt(4.0*np.pi/(2.0*k + 1.0))
        normm = np.sqrt(8.0*np.pi/(2.0*k + 1.0))
        for j in range(0, k+1):
            lm += 1
            if j == 0:
                norm = norm0
            else:
                norm = normm
            out[lm,:,:] = norm*d[lm,:,:]
    return out

def remove_background(d):
    print("### Removing thermal background state ###")
    d[0,0,0] -= 1.110720734539592
    d[0,1,0] -= -0.7853981633974483
    return d


# Create temperature field
if 'Codensity' in epm_file:
    group = quicc_file.create_group('temperature')
    epm_data = epm_file['Codensity']['Codensity'][:]
    ds = group.create_dataset('temperature', shape=epm_data.shape[0:2], dtype = ('<f8', (2,)))
    ds[:] = remove_background(rescale(epm_data, spec2D, spec3D))

# Create velocity field
if 'Velocity' in epm_file:
    group = quicc_file.create_group('velocity')
    epm_data = epm_file['Velocity']['VelocityTor'][:]
    ds = group.create_dataset('velocity_tor', shape=epm_data.shape[0:2], dtype = ('<f8', (2,)))
    ds[:] = rescale(epm_data, spec2D, spec3D)
    epm_data = epm_file['Velocity']['VelocityPol'][:]
    ds = group.create_dataset('velocity_pol', shape=epm_data.shape[0:2], dtype = ('<f8', (2,)))
    ds[:] = rescale(epm_data, spec2D, spec3D)

# Create magnetic field
if 'Magnetic' in epm_file:
    group = quicc_file.create_group('magnetic')
    epm_data = epm_file['Magnetic']['MagneticTor'][:]
    ds = group.create_dataset('magnetic_tor', shape=epm_data.shape[0:2], dtype = ('<f8', (2,)))
    ds[:] = rescale(epm_data, spec2D, spec3D)
    epm_data = epm_file['Magnetic']['MagneticPol'][:]
    ds = group.create_dataset('magnetic_pol', shape=epm_data.shape[0:2], dtype = ('<f8', (2,)))
    ds[:] = rescale(epm_data, spec2D, spec3D)

# Finished
epm_file.close()
quicc_file.close()
