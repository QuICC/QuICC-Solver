from quicc.reframe.library import testBase
import reframe as rfm

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split8_0_cuda(testBase):
    test = 'prof_TransformALegendreTests_Poly_PP_\<CudaIALegendreOperatorTypes\>_projector_id108_ulp20800_split8_0'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.0919701, -0.1, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split256_0_cuda(testBase):
    test = 'prof_TransformALegendreTests_Poly_PP_\<CudaIALegendreOperatorTypes\>_projector_id108_ulp20800_split256_0'
    steps = 500
    refs =  {   'p100': {
                    'applyOperatorsAvg': (0.00351006, -0.1, 0.1, 's'),
                },
            }



