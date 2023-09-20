# This reader is used to read the HDF5 output of the QuICC simulation code
#
# https://www.pasc-ch.org/projects/2021-2024/aqua-d/
#
# initial version provided by Philippe Marti (ETH-Zurich) to run serially
# This version was extended by Jean M. Favre, CSCS, to run in parallel
# Tested with ParaView v5.11
# Thu  6 Jul 10:42:48 CEST 2023
# This version was extended by Philippe Marti (ETH-Zurich) to deal with 
# missing phi point

# This is the module to import. It provides VTKPythonAlgorithmBase, the base class
# for all python-based vtkAlgorithm subclasses in VTK and decorators used to
# 'register' the algorithm with ParaView along with information about UI.
from paraview.util.vtkAlgorithm import *
from vtk.numpy_interface import algorithms as algs
from vtkmodules.numpy_interface import dataset_adapter as dsa
import vtk
import numpy as np

#------------------------------------------------------------------------------
# A reader example.
#------------------------------------------------------------------------------
def createModifiedCallback(anobject):
    import weakref
    weakref_obj = weakref.ref(anobject)
    anobject = None
    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()
    return _markmodified

#
# Base for QuICC Reader
#
class PythonQuICCReaderBase(VTKPythonAlgorithmBase):
    """Base for the a reader that reads QuICC's physical space HDF5 file."""
    def __init__(self, nInputPorts, nOutputPorts, outputType):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=nInputPorts, nOutputPorts=nOutputPorts, outputType=outputType)
        self._reset()

        # Setup array selection
        from vtkmodules.vtkCommonCore import vtkDataArraySelection
        self._arrayselection = vtkDataArraySelection()
        self._arrayselection.AddObserver("ModifiedEvent", createModifiedCallback(self))
        self._expand_phi = False

    def _get_raw_data(self, requested_time=None):
        if self._current_file is not None:
            if requested_time is not None:
                self._current_file = self._files[requested_time]
            return self._current_file

        if len(self._filenames) == 0:
            # Note, exceptions are totally fine!
            raise RuntimeError("No filename(s) specified")

        import h5py
        for fname in self._filenames:
            f = h5py.File(fname, 'r')
            time = f['run']['time'][()]
            self._files[time] = f
            self._timesteps.append(time)
        self._current_file = self._files[self._timesteps[0]]

        # Get spatial scheme
        self._scheme = self._current_file['/'].attrs['type']

        # Spherical schemes
        if self._scheme in [b'SLFl', b'SLFm', b'WLFl', b'WLFm']:
            self._components = ['r', 'theta', 'phi']
            self._make_mesh_coordinates = self._spherical_mesh
            self._vcomp = dict({'x':self._spherical_x_vcomp, 'y':self._spherical_y_vcomp, 'z':self._spherical_z_vcomp})
        # Plane layer
        elif self._scheme in [b'TFF']:
            self._components = ['z', 'x', 'y']
            self._make_mesh_coordinates = self._planelayer_mesh
            self._vcomp = dict({'x':self._planelayer_x_vcomp, 'y':self._planelayer_y_vcomp, 'z':self._planelayer_z_vcomp})
        else:
            raise RuntimeError('Unknown spatial scheme')

        # Add scalars to array selection
        for grp,s in self._get_scalar_fields():
            self._arrayselection.AddArray(s)

        # Add vector components to array selection
        for grp,v in self._get_vector_components():
            self._arrayselection.AddArray(v, False)

        # Add vectors to array selection
        for v in self._get_vector_fields():
            self._arrayselection.AddArray(v)

        return self._get_raw_data(requested_time)

    def _reset(self):
        """Reset members"""
        self._filenames = list()
        self._files = dict()
        self._current_file = None
        self._timesteps = list()
        self._scheme = None
        self._components = None
        self._make_mesh_coordinates = None
        self._vcomp = None

    def _get_timesteps(self):
        self._get_raw_data()
        return self._timesteps

    def _get_update_time(self, outInfo):
        executive = self.GetExecutive()
        timesteps = self._get_timesteps()
        if timesteps is None or len(timesteps) == 0:
            return None
        elif outInfo.Has(executive.UPDATE_TIME_STEP()) and len(timesteps) > 0:
            utime = outInfo.Get(executive.UPDATE_TIME_STEP())
            dtime = timesteps[0]
            for atime in timesteps:
                if atime > utime:
                    return dtime
                else:
                    dtime = atime
            return dtime
        else:
            assert(len(timesteps) > 0)
            return timesteps[0]

    def _get_update_extent(self, outInfo):
        executive = self.GetExecutive()
        if outInfo.Has(executive.UPDATE_EXTENT()):
            exts = outInfo.Get(executive.UPDATE_EXTENT())
            return exts

    def _get_whole_extent(self, outInfo):
        executive = self.GetExecutive()
        if outInfo.Has(executive.WHOLE_EXTENT()):
            dims = outInfo.Get(executive.WHOLE_EXTENT())
            return dims

    def _get_scalar_fields(self):
        """Get names of scalar fields"""

        assert(self._current_file is not None)
        raw_data = self._current_file
        names = list()
        for s in list(raw_data):
            if s in list(raw_data[s]):
                names.append((s,s))
        return names

    def _get_vector_components(self):
        """Get names of vector fields"""

        assert(self._current_file is not None)
        raw_data = self._current_file
        names = list()
        for v in list(raw_data):
            if f'{v}_{self._components[0]}' in list(raw_data[v]) and len(raw_data[v]) == 3:
                for c in self._components:
                    names.append((v,f'{v}_{c}'))
        return names

    def _get_vector_fields(self):
        """Get names of vector fields"""

        assert(self._current_file is not None)
        raw_data = self._current_file
        names = list()
        for v in list(raw_data):
            if f'{v}_{self._components[0]}' in list(raw_data[v]) and len(raw_data[v]) == 3:
                names.append(v)
        return names

    def _add_mesh_points(self, pvX, pvY, pvZ, mesh):
        """Add point coordinates to mesh"""

        coordinates = algs.make_vector(pvX.ravel(),pvY.ravel(),pvZ.ravel())
        pts = vtk.vtkPoints()
        pts.SetData(dsa.numpyTovtkDataArray(coordinates,'Points'))
        mesh.SetPoints(pts)

    def _add_mesh_cells(self, exts, dims, mesh):
        """Add cells to mesh"""

        raise RuntimeError('Adding cells to mesh is not implemented')

    def _has_phi_boundary(self, exts, raw_data):
        """Check if extents contain phi boundary"""

        val = False
        # Check if phi is expanded and boundary is in expand
        if self._expand_phi and exts[1] == raw_data['mesh'][f'grid_phi'].size:
            val = True

        return val

    def _add_scalar_pointdata(self, mesh, raw_data, exts):
        """Add scalar field point data to VTK grid"""

        for grp,s in self._get_scalar_fields():
            if self._arrayselection.ArrayIsEnabled(s):
                if self._has_phi_boundary(exts, raw_data):
                    # this is the partitioned subset selection for parallel reading
                    # and append 2pi data
                    data = np.append(
                            raw_data[grp][s][()][exts[4]:exts[5]+1, exts[2]:exts[3]+1, exts[0]:exts[1]],
                            raw_data[grp][s][()][exts[4]:exts[5]+1, exts[2]:exts[3]+1, 0:1],
                            2)
                else:
                    # this is the partitioned subset selection for parallel reading
                    data = raw_data[grp][s][()][exts[4]:exts[5]+1, exts[2]:exts[3]+1, exts[0]:exts[1]+1]
                mesh.GetPointData().AddArray(dsa.numpyTovtkDataArray(data.ravel(),
                                             s))

    def _add_vector_pointdata(self, mesh, raw_data, mgrid, exts):
        """Add vector field (X,Y,Z) point data to VTK mesh"""

        for v in self._get_vector_fields():
            if self._arrayselection.ArrayIsEnabled(v):
                tmp = np.zeros((mgrid[0].size, 3))
                comps = list()
                for i,c in enumerate(self._components):
                    if self._has_phi_boundary(exts, raw_data):
                        data = np.append(
                                raw_data[v][f'{v}_{c}'][()][exts[4]:exts[5]+1, exts[2]:exts[3]+1, exts[0]:exts[1]],
                                raw_data[v][f'{v}_{c}'][()][exts[4]:exts[5]+1, exts[2]:exts[3]+1, 0:1],
                                2)
                    else:
                        data = raw_data[v][f'{v}_{c}'][()][exts[4]:exts[5]+1, exts[2]:exts[3]+1, exts[0]:exts[1]+1]
                    comps.append(data)
                for i,c in enumerate(list(['z','x','y'])):
                    tmp[:,i] = self._vcomp[c](comps, mgrid).ravel()
                mesh.GetPointData().AddArray(dsa.numpyTovtkDataArray(tmp, v))

    def _add_vcomponent_pointdata(self, mesh, raw_data, exts):
        """Add vector field components from data geometry point data to VTK grid"""

        for grp,v in self._get_vector_components():
            if self._arrayselection.ArrayIsEnabled(v):
                if self._has_phi_boundary(exts, raw_data):
                    data = np.append(
                            raw_data[grp][v][()][exts[4]:exts[5]+1, exts[2]:exts[3]+1, exts[0]:exts[1]],
                            raw_data[grp][v][()][exts[4]:exts[5]+1, exts[2]:exts[3]+1, 0:1],
                           2)
                else:
                    # this is the partitioned subset selection for parallel reading
                    data = raw_data[grp][v][()][exts[4]:exts[5]+1, exts[2]:exts[3]+1, exts[0]:exts[1]+1]
                mesh.GetPointData().AddArray(dsa.numpyTovtkDataArray(data.ravel(),
                                             v))

    def _get_array_selection(self):
        return self._arrayselection

    def _AddFileName(self, name):
        """Specify filename for the file to read."""
        if name not in self._filenames:
            self._filenames.append(name)
            self.Modified()

    def RemoveAllFileNames(self):
        self._reset()
        self.Modified()

    def _GetTimestepValues(self):
        return self._get_timesteps()

    # Array selection API is typical with readers in VTK
    # This is intended to allow ability for users to choose which arrays to
    # load. To expose that in ParaView, simply use the
    # smproperty.dataarrayselection().
    # This method **must** return a `vtkDataArraySelection` instance.
    def _GetDataArraySelection(self):
        return self._get_array_selection()

    def RequestInformation(self, request, inInfoVec, outInfoVec):
        executive = self.GetExecutive()
        outInfo = outInfoVec.GetInformationObject(0)
        outInfo.Remove(executive.TIME_STEPS())
        outInfo.Remove(executive.TIME_RANGE())

        timesteps = self._get_timesteps()
        if timesteps is not None:
            for t in timesteps:
                outInfo.Append(executive.TIME_STEPS(), t)
            outInfo.Append(executive.TIME_RANGE(), timesteps[0])
            outInfo.Append(executive.TIME_RANGE(), timesteps[-1])

        data_time = self._get_update_time(outInfoVec.GetInformationObject(0))
        raw_data = self._get_raw_data(data_time)
        size_r = raw_data['mesh']['grid_r'].size
        size_theta = raw_data['mesh']['grid_theta'].size
        size_phi = raw_data['mesh']['grid_phi'].size
        # Expand phi grid by one to close domain
        if self._expand_phi:
            size_phi += 1
        dims = [size_phi, size_theta, size_r]
        outInfo.Set(executive.WHOLE_EXTENT(), 0, dims[0]-1 , 0, dims[1]-1 , 0, dims[2]-1)
        outInfo.Set(self.CAN_PRODUCE_SUB_EXTENT(), 1)

        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        data_time = self._get_update_time(outInfoVec.GetInformationObject(0))
        exts = self._get_update_extent(outInfoVec.GetInformationObject(0))
        raw_data = self._get_raw_data(data_time)
        # Add 2pi to grid if extent contains last grid point
        if self._has_phi_boundary(exts, raw_data):
            phi_grid = np.concatenate((raw_data['mesh'][f'grid_phi'][()][exts[0]:exts[1]], np.array([2*np.pi])))
        else:
            phi_grid = raw_data['mesh'][f'grid_phi'][()][exts[0]:exts[1]+1]
        grid = [raw_data['mesh'][f'grid_r'][()][exts[4]:exts[5]+1],
                raw_data['mesh'][f'grid_theta'][()][exts[2]:exts[3]+1],
                phi_grid]
        mgrid = np.meshgrid(*grid, indexing="ij")
        X, Y, Z = self._make_mesh_coordinates(*mgrid)

        # Create mesh
        output = self._create_mesh(outInfoVec)

        # Create points of mesh
        self._add_mesh_points(Z, X, Y, output)

        # Create cells of mesh
        whole_exts = self._get_whole_extent(outInfoVec.GetInformationObject(0))
        self._add_mesh_cells(exts, (whole_exts[1]+1, whole_exts[3]+1, whole_exts[5]+1), output)

        # Scalar data
        self._add_scalar_pointdata(output, raw_data, exts)

        # Vector field data
        self._add_vector_pointdata(output, raw_data, mgrid, exts)

        # Vector field components data
        self._add_vcomponent_pointdata(output, raw_data, exts)

        if data_time is not None:
            output.GetInformation().Set(output.DATA_TIME_STEP(), data_time)
        return 1

    def _spherical_x_vcomp(self, v, grid):
        """Compute X vector component from spherical field and grid"""

        vR, vT, vP = v
        gR, gT, gP = grid
        cT, sT = (np.cos(gT), np.sin(gT))
        cP, sP = (np.cos(gP), np.sin(gP))
        return vR*sT*cP + vT*cT*cP - vP*sP

    def _spherical_y_vcomp(self, v, grid):
        """Compute Y vector component from spherical field and grid"""

        vR, vT, vP = v
        gR, gT, gP = grid
        cT, sT = (np.cos(gT), np.sin(gT))
        cP, sP = (np.cos(gP), np.sin(gP))
        return vR*sT*sP + vT*cT*sP + vP*cP

    def _spherical_z_vcomp(self, v, grid):
        """Compute Z vector component from spherical field and grid"""

        vR, vT, vP = v
        gR, gT, gP = grid
        cT, sT = (np.cos(gT), np.sin(gT))
        return vR*cT - vT*sT

    def _spherical_mesh(self, gR, gT, gP):
        """Generate structured mesh for spherical geometry"""
        X = gR * np.sin(gT) * np.cos(gP)
        Y = gR * np.sin(gT) * np.sin(gP)
        Z = gR * np.cos(gT)

        return (X,Y,Z)

    def _planelayer_x_vcomp(self, v, grid):
        """Compute X vector component from spherical field and grid"""

        return v[1]

    def _planelayer_y_vcomp(self, v, grid):
        """Compute Y vector component from spherical field and grid"""

        return v[2]

    def _planelayer_z_vcomp(self, v, grid):
        """Compute Z vector component from spherical field and grid"""

        return v[0]

    def _planelayer_mesh(self, gZ, gX, gY):
        """Generate structured mesh for plane layer"""

        return (gZ,gX,gY)

#
# QuICC Reader generating a structured grid
#
@smproxy.reader(name="PythonQuICCStructuredReader", label="QuICC HDF5 Reader (structured grid)",
                extensions="hdf5", file_description="QuICC physical space HDF5 files")
class PythonQuICCStructuredReader(PythonQuICCReaderBase):
    """A reader that reads QuICC's physical space HDF5 files and creates a structured grid."""
    def __init__(self):
        PythonQuICCReaderBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkStructuredGrid')
        self._expand_phi = True

    @smproperty.stringvector(name="FileNames",
                                label="File Names",
                                animateable="1",
                                clean_command="RemoveAllFileNames",
                                command="AddFileName",
                                repeat_command="1",
                                number_of_elements="0",
                                panel_visibility="never"
                                )
    @smdomain.filelist()
    @smhint.filechooser(extensions="hdf5", file_description="QuICC's physical space HDF5 files")
    def AddFileName(self, name):
        self._AddFileName(name)

    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return self._GetTimestepValues()

    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        return self._GetDataArraySelection()

    def _add_mesh_cells(self, exts, dims, mesh):
        """Add cells to mesh"""

        # Nothing to do for structured grid

    def _create_mesh(self, outInfoVec):
        """Create the VTK mesh"""

        executive = self.GetExecutive()
        outInfo = executive.GetOutputInformation(0)
        exts = [executive.UPDATE_EXTENT().Get(outInfo, i) for i in range(6)]
        mesh = vtk.vtkStructuredGrid.GetData(outInfoVec,0)
        mesh.SetExtent(exts)

        return mesh


#
# QuICC Reader generating an unstructured grid
#
@smproxy.reader(name="PythonQuICCUnstructuredReader", label="QuICC HDF5 Reader (unstructured grid)",
                extensions="hdf5", file_description="QuICC physical space HDF5 files")
class PythonQuICCUnstructuredReader(PythonQuICCReaderBase):
    """A reader that reads QuICC's physical space HDF5 files and creates an unstructured grid."""
    def __init__(self):
        PythonQuICCReaderBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkUnstructuredGrid')
        self._expand_phi = False

    @smproperty.stringvector(name="FileNames",
                                label="File Names",
                                animateable="1",
                                clean_command="RemoveAllFileNames",
                                command="AddFileName",
                                repeat_command="1",
                                number_of_elements="0",
                                panel_visibility="never"
                                )
    @smdomain.filelist()
    @smhint.filechooser(extensions="hdf5", file_description="QuICC's physical space HDF5 files")
    def AddFileName(self, name):
        self._AddFileName(name)

    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return self._GetTimestepValues()

    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        return self._GetDataArraySelection()

    def _add_mesh_cells(self, exts, dims, mesh):
        """Add cells to mesh"""

        hex = vtk.vtkHexahedron()
        cells = vtk.vtkCellArray()
        nF, nM, nS = (exts[1]-exts[0]+1, exts[3]-exts[2]+1, exts[5]-exts[4]+1)
        nF_ = nF - 1
        if nF == dims[0]:
            nF_ = dims[0]
        for k in range(0, nS-1):
            for j in range(0, nM-1):
                for i in range(0, nF_):
                    hex.GetPointIds().SetId(0,i + j*nF + k*(nF*nM))
                    hex.GetPointIds().SetId(1,(i + 1)%nF + j*nF + k*(nF*nM))
                    hex.GetPointIds().SetId(2,(i + 1)%nF + (j+1)*nF + k*(nF*nM))
                    hex.GetPointIds().SetId(3,i + (j+1)*nF + k*(nF*nM))
                    hex.GetPointIds().SetId(4,i + j*nF + (k+1)*(nF*nM))
                    hex.GetPointIds().SetId(5,(i + 1)%nF + j*nF + (k+1)*(nF*nM))
                    hex.GetPointIds().SetId(6,(i + 1)%nF + (j+1)*nF + (k+1)*(nF*nM))
                    hex.GetPointIds().SetId(7,i + (j+1)*nF + (k+1)*(nF*nM))
                    cells.InsertNextCell(hex)
        mesh.SetCells(vtk.VTK_HEXAHEDRON, cells)

    def _create_mesh(self, outInfoVec):
        """Create the VTK mesh"""

        mesh = vtk.vtkUnstructuredGrid.GetData(outInfoVec,0)

        return mesh


# Test
def test_PythonQuICCStructuredReader(fname):
    reader = PythonQuICCStructuredReader()
    reader.AddFileName(fname)
    reader.Update()

def test_PythonQuICCUnstructuredReader(fname):
    reader = PythonQuICCUnstructuredReader()
    reader.AddFileName(fname)
    reader.Update()

if __name__ == "__main__":
    test_PythonQuICCStructuredReader("/tmp/visState0000.hdf5")
    test_PythonQuICCUnstructuredReader("/tmp/visState0000.hdf5")

    from paraview.detail.pythonalgorithm import get_plugin_xmls
    from xml.dom.minidom import parseString
    for xml in get_plugin_xmls(globals()):
        dom = parseString(xml)
        print(dom.toprettyxml(" ","\n"))
