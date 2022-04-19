#!/usr/bin/env python

from __future__ import print_function

import sys, getopt
import numpy as np
import numpy.fft as fft
import tables
import xml.etree.ElementTree as ET

argv = sys.argv[1:]
inputfile = ''
outputfile = ''
sim_type = None
viz_fields = False
try:
    opts, args = getopt.getopt(argv,"hi:o:t:v")
except getopt.GetoptError:
    print('generate_singlemode_state.py -i <inputfile> -o <outputfile> -t <simulation type>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('generate_singlemode_state.py -i <inputfile> -o <outputfile> -t <simulation type>')
        sys.exit()
    elif opt in ("-i"):
        inputfile = arg
    elif opt in ("-o"):
        outputfile = arg
    elif opt in ("-t"):
        sim_type = arg
    elif opt in ("-v"):
        viz_fields = True

if sim_type not in ["FPlane3DQG", "RRBCPlane", "RRBCPlaneMean", "RRBCPlaneDMean"]:
    print("Unknown simulation type")
    sys.exit(2)

if inputfile == '':
    print("Input file was not set")
    sys.exit(2)

if outputfile == '':
    print("Output file was not set")
    sys.exit(2)

cfg_file = open('parameters.cfg', 'r')
header = cfg_file.readline()
config = ET.fromstring(header + '<root>'+cfg_file.read()+'</root>')

# Extract truncation information
for h in config.find('framework').find('truncation'):
    if h.tag == "dim1D":
        nnz = int(h.text)+1
    elif h.tag == "dim2D":
        nnx = int(h.text)+1
    elif h.tag == "dim3D":
        nny = int(h.text)+1
    elif h.tag == "box2D":
        boxX = int(h.text)
    elif h.tag == "box3D":
        boxY = int(h.text)
    elif h.tag == "kc2D":
        kx = float(h.text)
    elif h.tag == "kc3D":
        ky = float(h.text)

# Extract physical information
phys_params = dict()
for h in config.find('simulation').find('physical'):
    phys_params[h.tag] = float(h.text)
for h in config.find('simulation').find('boundary'):
    phys_params['bc_' + h.tag] = int(h.text)
if 'ekman' in phys_params.keys():
    phys_params['rossby'] = phys_params['ekman']**(1./3.)

# read data and get Chebyshev expansion
tmp = np.genfromtxt(inputfile, skip_header = 1)
nz  = tmp.shape[0]
grid = np.array([np.cos(i*np.pi/(nz - 1.0)) for i in range(0,nz)])
fields = ["mean_temperature", "dz_meantemperature", "velocityz", "streamfunction", "temperature"]
data = dict()
for i,field in enumerate(fields):
    data[field] = tmp[:,i]

for f in fields:
    phys_data = np.zeros(2*(nz - 1))
    phys_data[0:nz] = data[f][0:nz][::-1]
    if f == "mean_temperature":
        phys_data[0:nz] -= 1.0 - 0.5*(grid + 1.0)
    if f == "mean_temperature" and sim_type in ["RRBCPlane", "RRBCPlaneMean", "RRBCPlaneDMean"]:
        phys_data[0:nz] /= phys_params.get('rossby', 1.0)
    phys_data[nz:] = phys_data[1:nz-1][::-1]
    cheb = fft.rfft(phys_data).real/phys_data.shape[0]
    data[f] = np.zeros((nnz,))
    data[f][0:min(cheb.shape[0],nnz)] = cheb[0:min(cheb.shape[0],nnz)]
    if viz_fields:
        import matplotlib.pylab as pl
        pl.subplot(121)
        pl.plot(grid, phys_data[0:nz])
        pl.title(f)
        pl.subplot(122)
        pl.semilogy(np.abs(data[f][0:min(cheb.shape[0],nnz)]))
        pl.title(f)
        pl.show()

if sim_type in ["RRBCPlane", "RRBCPlaneMean", "RRBCPlaneDMean"]:
    data["velocity_tor"] = -data["streamfunction"]
    data["velocity_pol"] = (1.0/kx**2)*data["velocityz"]

# Write HDF5 header attributes
hdf5_file = tables.open_file(outputfile, mode = 'w')
hdf5_file.set_node_attr('/', 'header', 'StateFile'.encode('ascii')) 
hdf5_file.set_node_attr('/', 'type', 'TFF'.encode('ascii')) 
hdf5_file.set_node_attr('/', 'version', '1.0'.encode('ascii')) 

# Write HDF5 physical parameters
group = hdf5_file.create_group("/","physical")
for k,v in phys_params.items():
    hdf5_file.create_array(group, k, v)

# Write HDF5 run parameters
group = hdf5_file.create_group("/", "run")
hdf5_file.create_array(group, "time", 0.0)
hdf5_file.create_array(group, "timestep", 1e-4)

# Write HDF5 truncation information
group = hdf5_file.create_group("/","truncation")
group2 = hdf5_file.create_group(group, "physical")
hdf5_file.create_array(group2, "dim1D", 1)
hdf5_file.create_array(group2, "dim2D", 1)
hdf5_file.create_array(group2, "dim3D", 1)
group2 = hdf5_file.create_group(group, "spectral")
hdf5_file.create_array(group2, "dim1D", nnz-1)
hdf5_file.create_array(group2, "dim2D", 2*(nnx-1))
hdf5_file.create_array(group2, "dim3D", nny-1)
group2 = hdf5_file.create_group(group, "transform")
hdf5_file.create_array(group2, "dim1D", nnz-1)
hdf5_file.create_array(group2, "dim2D", 2*(nnx-1))
hdf5_file.create_array(group2, "dim3D", nny-1)

def writeMean(name, c = 1.0):
    group = hdf5_file.create_group("/", name)
    tmp = np.zeros((2*nnx-1,nny,nnz), dtype=np.complex128)
    tmp[0,0,0:nnz].real = c*data[name][0:nnz]
    hdf5_file.create_array(group, name, tmp)

def writeScalar(name, c = 1.0, mean = None, cmean = 1.0):
    group = hdf5_file.create_group("/", name)
    tmp = np.zeros((2*nnx-1,nny,nnz), dtype=np.complex128)
    if mean is not None:
        tmp[0,0,0:nnz].real = cmean*data[mean][0:nnz]
    tmp[boxX,0,0:nnz].real = c*0.5*data[name][0:nnz]
    #tmp[-boxX,0,0:nnz].real = c*0.5*data[name][0:nnz]
    tmp[0,boxY,0:nnz].real = c*0.5*data[name][0:nnz]
    hdf5_file.create_array(group, name, tmp)

def writeVector(name, comp, coeff = None):
    group = hdf5_file.create_group("/",name)
    if coeff is None:
        coeff = (1.0,)*len(comp)
    for c,b in zip(comp,coeff):
        tmp = np.zeros((2*nnx-1,nny,nnz), dtype=np.complex128)
        tmp[boxX,0,0:nnz].real = b*0.5*data[name + '_' + c][0:nnz]
        #tmp[-boxX,0,0:nnz].real = b*0.5*data[name + '_' + c][0:nnz]
        tmp[0,boxY,0:nnz].real = b*0.5*data[name + '_' + c][0:nnz]
        hdf5_file.create_array(group, name + '_' + c, tmp)

############################################
# Write HDF5 FPlane3DQG fields
if sim_type == "FPlane3DQG":
    writeScalar("temperature")
    writeScalar("streamfunction")
    writeMean("dz_meantemperature")
    writeScalar("velocityz")

############################################
# Write HDF5 RRBCPlane fields
if sim_type == "RRBCPlane":
    writeScalar("temperature", mean = "mean_temperature")
    writeVector("velocity", ["tor","pol"])

############################################
# Write HDF5 RRBCPlaneMean fields
if sim_type == "RRBCPlaneMean":
    writeScalar("temperature")
    writeMean("meantemperature")
    writeVector("velocity", ["tor","pol"])

############################################
# Write HDF5 RRBCPlaneDMean fields
if sim_type == "RRBCPlaneDMean":
    writeScalar("temperature")
    writeMean("dz_meantemperature")
    writeVector("velocity", ["tor","pol"])

hdf5_file.close()
