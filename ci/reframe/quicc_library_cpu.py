from quicc.reframe.library import testBase
import reframe as rfm


@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split96_0(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_projector_id108_ulp20800_split96_0'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.0303996, -0.1, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.0100948, -0.1, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split288_0(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_projector_id108_ulp20800_split288_0'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.0105246, -0.1, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.00395184, -0.1, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split3072_0(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_projector_id108_ulp20800_split3072_0'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00122683, -0.1, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.000509902, -0.1, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split9216_0(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_projector_id108_ulp20800_split9216_0'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.000536682, -0.1, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.000232762, -0.1, 0.1, 's'),
                },
            }


