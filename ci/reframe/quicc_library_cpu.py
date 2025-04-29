from quicc.reframe.library import splitTestTransform,testStateFile
import reframe as rfm

#
# Jones-Worland
#

# Quadrature algorithm

@rfm.simple_test
class TransformWorlandTests_Poly_P_base_t_projector(splitTestTransform):
    region = ['applyOperators']
    splitting = 8*36

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformWorlandTests_Poly_P_viewCpu_t_projector(splitTestTransform):
    region = ['applyOperators']
    splitting = 8*36

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformWorlandTests_Poly_P_base_t_integrator(splitTestTransform):
    region = ['applyOperators']
    splitting = 8*36

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformWorlandTests_Poly_P_viewCpu_t_integrator(splitTestTransform):
    region = ['applyOperators']
    splitting = 8*36

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

# FFT algorithm

@rfm.simple_test
class TransformWorlandTests_Fft_P_base_t_projector(splitTestTransform):
    region = ['transform']
    splitting = 8*36

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformWorlandTests_Fft_P_base_t_integrator(splitTestTransform):
    region = ['transform']
    splitting = 8*36

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

#
# Associated Legendre
#

@rfm.simple_test
class TransformALegendreTests_Poly_P_base_t_projector(splitTestTransform):
    region = ['applyOperators']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

        # Test id = 108, splitting = 8*12
        self.add_reference(108, 8*12, 'icelake',   'applyOperatorsAvg', (0.00610543, -0.25, 0.1, 's'))
        self.add_reference(108, 8*12, 'skylake',   'applyOperatorsAvg', (0.0303996, -0.25, 0.1, 's'))
        self.add_reference(108, 8*12, 'haswell',   'applyOperatorsAvg', (0.0181676, -0.25, 0.1, 's'))
        # Test id = 108, splitting = 8*36
        self.add_reference(108, 8*36, 'icelake',   'applyOperatorsAvg', (0.002163, -0.25, 0.1, 's'))
        self.add_reference(108, 8*36, 'skylake',   'applyOperatorsAvg', (0.0105246, -0.25, 0.1, 's'))
        # Test id = 108, splitting = 256*12
        self.add_reference(108, 256*12, 'icelake',   'applyOperatorsAvg', (0.00022611, -0.25, 0.1, 's'))
        self.add_reference(108, 256*12, 'skylake',   'applyOperatorsAvg', (0.00122683, -0.25, 0.1, 's'))
        self.add_reference(108, 256*12, 'broadwell', 'applyOperatorsAvg', (0.000509902, -0.25, 0.1, 's'))
        # Test id = 108, splitting = 256*36
        self.add_reference(108, 256*36, 'icelake',   'applyOperatorsAvg', (9.69103e-05, -0.25, 0.1, 's'))
        self.add_reference(108, 256*36, 'skylake',   'applyOperatorsAvg', (0.000536682, -0.25, 0.1, 's'))

@rfm.simple_test
class TransformALegendreTests_Poly_P_viewCpu_t_projector(splitTestTransform):
    region = ['applyOperators']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

        # Test id = 108, splitting = 8*12
        self.add_reference(108, 8*12, 'icelake',   'applyOperatorsAvg', (0.00657657, -0.25, 0.1, 's'))
        self.add_reference(108, 8*12, 'haswell',   'applyOperatorsAvg', (0.0186619, -0.25, 0.1, 's'))

@rfm.simple_test
class TransformALegendreTests_Poly_LlDivS1Dp_base_t_projector(splitTestTransform):
    region = ['applyOperators']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

        # Test id = 108, splitting = 8*12
        self.add_reference(108, 8*12, 'icelake',   'applyOperatorsAvg', (0.00613202, -0.25, 0.1, 's'))
        self.add_reference(108, 8*12, 'haswell',   'applyOperatorsAvg', (0.0185139, -0.25, 0.1, 's'))

@rfm.simple_test
class TransformALegendreTests_Poly_LlDivS1Dp_viewCpu_t_projector(splitTestTransform):
    region = ['applyOperators']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

        # Test id = 108, splitting = 8*12
        self.add_reference(108, 8*12, 'icelake',   'applyOperatorsAvg', (0.00677547, -0.25, 0.1, 's'))
        self.add_reference(108, 8*12, 'haswell',   'applyOperatorsAvg', (0.0186614, -0.25, 0.1, 's'))


@rfm.simple_test
class TransformALegendreTests_Poly_P_base_t_integrator(splitTestTransform):
    region = ['applyOperators']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

        # Test id = 108, splitting = 8*12
        self.add_reference(108, 8*12, 'icelake',   'applyOperatorsAvg', (0.00625215, -0.25, 0.1, 's'))
        self.add_reference(108, 8*12, 'haswell',   'applyOperatorsAvg', (0.019000, -0.25, 0.1, 's'))

@rfm.simple_test
class TransformALegendreTests_Poly_P_viewCpu_t_integrator(splitTestTransform):
    region = ['applyOperators']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

        # Test id = 108, splitting = 8*12
        self.add_reference(108, 8*12, 'icelake',   'applyOperatorsAvg', (0.00680607, -0.25, 0.1, 's'))
        self.add_reference(108, 8*12, 'haswell',   'applyOperatorsAvg', (0.019000, -0.25, 0.1, 's'))

#
# Fourier
#

@rfm.simple_test
class TransformFourierTests_Mixed_P_base_t_projector(splitTestTransform):
    region = ['transform']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

        # Test id = 108, splitting = 8*12
        self.add_reference(108, 8*12, 'icelake',   'transformAvg', (0.00607775, -0.25, 0.05, 's'))
        self.add_reference(108, 8*12, 'skylake',   'transformAvg', (0.00369302, -0.25, 0.05, 's'))

@rfm.simple_test
class TransformFourierTests_Mixed_P_viewCpu_t_projector(splitTestTransform):
    region = ['transform','DOp::applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

        # Test id = 108, splitting = 8*12
        self.add_reference(108, 8*12, 'icelake',   'transformAvg', (0.00629394, -0.25, 0.05, 's'))
        self.add_reference(108, 8*12, 'skylake',   'transformAvg', (0.00392791, -0.25, 0.05, 's'))

@rfm.simple_test
class TransformFourierTests_Mixed_D1_base_t_projector(splitTestTransform):
    region = ['transform']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

        # Test id = 108, splitting = 8*12
        self.add_reference(108, 8*12, 'icelake',   'transformAvg', (0.0114277, -0.25, 0.05, 's'))
        self.add_reference(108, 8*12, 'skylake',   'transformAvg', (0.0114620, -0.25, 0.05, 's'))

@rfm.simple_test
class TransformFourierTests_Mixed_D1_viewCpu_t_projector(splitTestTransform):
    region = ['transform','DOp::applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

        # Test id = 108, splitting = 8*12
        self.add_reference(108, 8*12, 'icelake',   'transformAvg', (0.00684981, -0.25, 0.05, 's'))
        self.add_reference(108, 8*12, 'skylake',   'transformAvg', (0.00435947, -0.25, 0.05, 's'))

@rfm.simple_test
class TransformFourierTests_Mixed_P_base_t_integrator(splitTestTransform):
    region = ['transform']
    splitting = 8*36

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformFourierTests_Mixed_P_viewCpu_t_integrator(splitTestTransform):
    region = ['transform','DOp::applyImpl']
    splitting = 8*36

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformFourierTests_Mixed_D1_base_t_integrator(splitTestTransform):
    region = ['transform']
    splitting = 8*36

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformFourierTests_Mixed_D1_viewCpu_t_integrator(splitTestTransform):
    region = ['transform','DOp::applyImpl']
    splitting = 8*36

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("cpu.json", self.__class__.__name__)

#
# HDF5 IO tests
# These tests are currently disable because they do not show the expected speedup.
#

# @rfm.simple_test
# class testStateFileTests_WLFl_read(testStateFile):
#     test = 'prof_FrameworkStateFileTests_WLFl_read_'
#     region = ['StateFileReader::read']
#     is_serial_test = False
#     refs = {   'icelake': {
#                     'perfIoAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoScalarsAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoVectorsAvg': (0.00831783, -0.25, 0.05, 's'),
#                 },
#                 'skylake': {
#                     'perfIoAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoScalarsAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoVectorsAvg': (0.00831783, -0.25, 0.05, 's'),
#                 },
#                 'broadwell': {
#                     'perfIoAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoScalarsAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoVectorsAvg': (0.00831783, -0.25, 0.05, 's'),
#                 },
#             }

# @rfm.simple_test
# class testStateFileTests_WLFl_write(testStateFile):
#     test = 'prof_FrameworkStateFileTests_WLFl_write_'
#     region = 'StateFileWriter::write'
#     is_serial_test = False
#     refs = {   'icelake': {
#                     'perfIoAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoScalarsAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoVectorsAvg': (0.00831783, -0.25, 0.05, 's'),
#                 },
#                 'skylake': {
#                     'perfIoAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoScalarsAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoVectorsAvg': (0.00831783, -0.25, 0.05, 's'),
#                 },
#                 'broadwell': {
#                     'perfIoAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoScalarsAvg': (0.00831783, -0.25, 0.05, 's'),
#                     'perfIoVectorsAvg': (0.00831783, -0.25, 0.05, 's'),
#                 },
#             }
