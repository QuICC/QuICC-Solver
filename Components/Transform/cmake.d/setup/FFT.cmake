###################################################
#------------- FFT IMPLEMENTATION ----------------#
###################################################

set(_FFT_BACKENDS "FFTW" "cuFFT")

quicc_create_option(NAME QUICC_FFT_MIXED
                    OPTS ${_FFT_BACKENDS}
                    LABEL "Mixed FFT backend"
                    ADVANCED)
quicc_add_definition(QUICC_FFT_MIXED)

quicc_create_option(NAME QUICC_FFT_COMPLEX
                    OPTS ${_FFT_BACKENDS}
                    LABEL "Complex FFT backend"
                    ADVANCED)
quicc_add_definition(QUICC_FFT_COMPLEX)

quicc_create_option(NAME QUICC_FFT_CHEBYSHEV
                    OPTS ${_FFT_BACKENDS}
                    LABEL "Chebyshev FFT backend"
                    ADVANCED)
quicc_add_definition(QUICC_FFT_CHEBYSHEV)

quicc_create_option(NAME QUICC_FFT_WORLAND
                    OPTS ${_FFT_BACKENDS}
                    LABEL "Worland FFT backend"
                    ADVANCED)
quicc_add_definition(QUICC_FFT_WORLAND)
