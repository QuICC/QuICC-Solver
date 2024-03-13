from quicc.reframe.library import splitTestTransform
import reframe as rfm

#
# Jones-Worland
#

@rfm.simple_test
class TransformWorlandTests_Poly_P_kokkos_t_projector(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyUnitOperator']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformWorlandTests_Poly_P_kokkos_t_integrator(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyUnitOperator']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformWorlandTests_Poly_P_viewGpu_t_projector(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformWorlandTests_Poly_P_viewGpu_t_integrator(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformWorlandTests_Fft_P_base_t_projector(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformWorlandTests_Fft_P_base_t_integrator(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

#
# Associated Legendre with cpu <-> gpu copy
#

@rfm.simple_test
class TransformALegendreTests_Poly_P_kokkos_t_projector(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyUnitOperator']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformALegendreTests_Poly_P_viewGpu_t_projector(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformALegendreTests_Poly_LlDivS1Dp_kokkos_t_projector(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyUnitOperator']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformALegendreTests_Poly_LlDivS1Dp_viewGpu_t_projector(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformALegendreTests_Poly_P_kokkos_t_integrator(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyUnitOperator']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

@rfm.simple_test
class TransformALegendreTests_Poly_P_viewGpu_t_integrator(splitTestTransform):
    splitting = 8
    region = ['applyOperators', 'applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__)

#
# Fourier with cpu <-> gpu copy
#

@rfm.simple_test
class TransformFourierTests_Mixed_P_viewGpu_t_projector(splitTestTransform):
    splitting = 96
    region = ['transform','DOp::applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__, filter='Min')

@rfm.simple_test
class TransformFourierTests_Mixed_P_viewGpuVkFFT_t_projector(splitTestTransform):
    splitting = 96
    region = ['transform','DOp::applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__, filter='Min')

@rfm.simple_test
class TransformFourierTests_Mixed_D1_viewGpu_t_projector(splitTestTransform):
    splitting = 96
    region = ['transform','DOp::applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__, filter='Min')

@rfm.simple_test
class TransformFourierTests_Mixed_D1_viewGpuVkFFT_t_projector(splitTestTransform):
    splitting = 96
    region = ['transform','DOp::applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__, filter='Min')

@rfm.simple_test
class TransformFourierTests_Mixed_P_viewGpu_t_integrator(splitTestTransform):
    splitting = 96
    region = ['transform','DOp::applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__, filter='Min')

@rfm.simple_test
class TransformFourierTests_Mixed_P_viewGpuVkFFT_t_integrator(splitTestTransform):
    splitting = 96
    region = ['transform','DOp::applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__, filter='Min')

@rfm.simple_test
class TransformFourierTests_Mixed_D1_viewGpu_t_integrator(splitTestTransform):
    splitting = 96
    region = ['transform','DOp::applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__, filter='Min')

@rfm.simple_test
class TransformFourierTests_Mixed_D1_viewGpuVkFFT_t_integrator(splitTestTransform):
    splitting = 96
    region = ['transform','DOp::applyImpl']

    def init_references(self):
        """ Initiallize references
        """

        self.read_references("gpu.json", self.__class__.__name__, filter='Min')

