from quicc.reframe.library import testBase
import reframe as rfm

#
# Associated Legendre with cpu <-> gpu copy
#

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_split8_0_cuda(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<kokkos_t\>_projector_id108_ulp.*_split8_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.0454438, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_split8_0_viewGpu(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<viewGpu_t\>_projector_id108_ulp.*_split8_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.0396912, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_split256_0_cuda(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<kokkos_t\>_projector_id108_ulp.*_split256_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.00190451, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_split256_0_viewGpu(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<viewGpu_t\>_projector_id108_ulp.*_split256_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.00130127, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_LlDivS1Dp_projector_id108_split8_0_cuda(testBase):
    test = 'prof_TransformALegendreTests_Poly_LlDivS1Dp_\<kokkos_t\>_projector_id108_ulp.*_split8_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.0454104, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_LlDivS1Dp_projector_id108_split8_0_viewGpu(testBase):
    test = 'prof_TransformALegendreTests_Poly_LlDivS1Dp_\<viewGpu_t\>_projector_id108_ulp.*_split8_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.0396904, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_integrator_id108_split8_0_cuda(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<kokkos_t\>_integrator_id108_ulp.*_split8_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.0433737, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_integrator_id108_split8_0_viewGpu(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<viewGpu_t\>_integrator_id108_ulp.*_split8_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.041728, -0.25, 0.1, 's'),
                },
            }

#
# Associated Legendre without cpu <-> gpu copy
#

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_split8_0_noCopy_cuda(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<kokkos_t\>_projector_id108_ulp.*_split8_0'
    region = 'applyUnitOperator'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.00249139, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_split8_0_noCopy_viewGpu(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<viewGpu_t\>_projector_id108_ulp.*_split8_0'
    region = 'ImplOp::applyImpl'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (1.18773e-05, -0.25, 0.1, 's'),
                },
            }

#
# Fourier with cpu <-> gpu copy
#

@rfm.simple_test
class testFourierTests_Mixed_D1_projector_id108_split96_0_viewGpu(testBase):
    test = 'prof_TransformFourierTests_Mixed_D1_\<viewGpu_t\>_projector_id108_ulp.*_split96_0'
    region = 'transform'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.00196887, -0.25, 0.2, 's'),
                },
            }
