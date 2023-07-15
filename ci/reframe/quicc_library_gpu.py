from quicc.reframe.library import testBase
import reframe as rfm

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
class testALegendreTests_Poly_P_projector_id108_split256_0_cuda(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<kokkos_t\>_projector_id108_ulp.*_split256_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.00190451, -0.25, 0.1, 's'),
                },
            }


@rfm.simple_test
class testFourierTests_Mixed_D1_projector_id108_split96_0_viewGpu(testBase):
    test = 'prof_TransformFourierTests_Mixed_D1_\<viewGpu_t\>_projector_id108_ulp.*_split96_0'
    region = 'transform'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.00196887, -0.25, 0.2, 's'),
                },
            }
