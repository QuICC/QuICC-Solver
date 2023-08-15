from quicc.reframe.library import testBase
import reframe as rfm

#
# Associated Legendre
#

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_split96_0_base(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<base_t\>_projector_id108_ulp.*_split96_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00610543, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.0100948, -0.25, 0.1, 's'),
                },
                'haswell': {
                    'applyOperatorsAvg': (0.0181676, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_split96_0_viewCpu(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<viewCpu_t\>_projector_id108_ulp.*_split96_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00657657, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.0100948, -0.25, 0.1, 's'),
                },
                'haswell': {
                    'applyOperatorsAvg': (0.0186619, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_split288_0_base(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<base_t\>_projector_id108_ulp.*_split288_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.002163, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.00395184, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_split3072_0_base(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<base_t\>_projector_id108_ulp.*_split3072_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00022611, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.000509902, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_projector_id108_split9216_0_base(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<base_t\>_projector_id108_ulp.*_split9216_0'
    region = 'applyOperators'
    steps = 1000
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (9.69103e-05, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.000232762, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_LlDivS1Dp_projector_id108_split96_0_base(testBase):
    test = 'prof_TransformALegendreTests_Poly_LlDivS1Dp_\<base_t\>_projector_id108_ulp.*_split96_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00613202, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.0100948, -0.25, 0.1, 's'),
                },
                'haswell': {
                    'applyOperatorsAvg': (0.0185139, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_LlDivS1Dp_projector_id108_split96_0_viewCpu(testBase):
    test = 'prof_TransformALegendreTests_Poly_LlDivS1Dp_\<viewCpu_t\>_projector_id108_ulp.*_split96_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00677547, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.0100948, -0.25, 0.1, 's'),
                },
                'haswell': {
                    'applyOperatorsAvg': (0.0186614, -0.25, 0.1, 's'),
                },
            }


@rfm.simple_test
class testALegendreTests_Poly_P_integrator_id108_split96_0_base(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<base_t\>_integrator_id108_ulp.*_split96_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00625215, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.011000, -0.25, 0.1, 's'),
                },
                'haswell': {
                    'applyOperatorsAvg': (0.019000, -0.25, 0.1, 's'),
                },
            }

@rfm.simple_test
class testALegendreTests_Poly_P_integrator_id108_split96_0_viewCpu(testBase):
    test = 'prof_TransformALegendreTests_Poly_P_\<viewCpu_t\>_integrator_id108_ulp.*_split96_0'
    region = 'applyOperators'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00680607, -0.25, 0.1, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.011600, -0.25, 0.1, 's'),
                },
                'haswell': {
                    'applyOperatorsAvg': (0.019000, -0.25, 0.1, 's'),
                },
            }

#
# Fourier
#

@rfm.simple_test
class testFourierTests_Mixed_P_projector_id108_split96_0_base(testBase):
    test = 'prof_TransformFourierTests_Mixed_P_\<base_t\>_projector_id108_ulp.*_split96_0'
    region = 'transform'
    steps = 500
    refs =  {   'icelake': {
                    'applyOperatorsAvg': (0.00607775, -0.25, 0.05, 's'),
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
                    'applyOperatorsAvg': (0.00629394, -0.25, 0.05, 's'),
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
                    'applyOperatorsAvg': (0.0114277, -0.25, 0.05, 's'),
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
                    'applyOperatorsAvg': (0.00684981, -0.25, 0.05, 's'),
                },
                'broadwell': {
                    'applyOperatorsAvg': (0.00545372, -0.25, 0.05, 's'),
                },
            }

