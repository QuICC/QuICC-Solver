from quicc.reframe.library import testBase
import reframe as rfm


@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split96_0(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<\>_projector_id108_ulp20800_split96_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.0303996, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.0100948, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split288_0(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<\>_projector_id108_ulp20800_split288_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.0105246, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.00395184, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split3072_0(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<\>_projector_id108_ulp20800_split3072_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00122683, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.000509902, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_ulp20800_split9216_0(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<\>_projector_id108_ulp20800_split9216_0'
    region = 'applyOperators'
    steps = 1000
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.000536682, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.000232762, -0.25, 0.1, 's'),
                },
            }


@rfm.simple_test
class testFourierTests_Mixed_P_projector_id108_split96_0_base(testBase):
    test = 'prof_TransformFourierTests_Mixed_P_\<base_t\>_projector_id108_ulp.*_split96_0'
    region = 'transform'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00667155, -0.25, 0.05, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.00480215, -0.25, 0.05, 's'),
                },
            }


@rfm.simple_test
class testFourierTests_Mixed_P_projector_id108_split96_0_view(testBase):
    test = 'prof_TransformFourierTests_Mixed_P_\<viewCpu_t\>_projector_id108_ulp.*_split96_0'
    region = 'transform'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00650388, -0.25, 0.05, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.00475000, -0.25, 0.05, 's'),
                },
            }


@rfm.simple_test
class testFourierTests_Mixed_D1_projector_id108_split96_0_base(testBase):
    test = 'prof_TransformFourierTests_Mixed_D1_\<base_t\>_projector_id108_ulp.*_split96_0'
    region = 'transform'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.0114620, -0.25, 0.05, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.0142788, -0.25, 0.05, 's'),
                },
            }


@rfm.simple_test
class testFourierTests_Mixed_D1_projector_id108_split96_0_view(testBase):
    test = 'prof_TransformFourierTests_Mixed_D1_\<viewCpu_t\>_projector_id108_ulp.*_split96_0'
    region = 'transform'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00831783, -0.25, 0.05, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.00545372, -0.25, 0.05, 's'),
                },
            }

