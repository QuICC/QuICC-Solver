###################################################
#-------------- WORLAND TRANSFORM ----------------#
###################################################

# Select backend for Worland transform
quicc_create_option(NAME QUICC_WORLAND_BACKEND
                    OPTS "Eigen" "BLAS" "Scalar"
                    LABEL "Worland backend"
                    ADVANCED)
quicc_add_definition(QUICC_WORLAND_BACKEND)

# Select implementations for Worland projectors
quicc_create_option(NAME QUICC_WORLAND_PROJIMPL
                    OPTS "Matrix" "OTF"
                    LABEL "Worland projector"
                    ADVANCED)
quicc_add_definition(QUICC_WORLAND_PROJIMPL)

# Select implementations for Worland integrators
quicc_create_option(NAME QUICC_WORLAND_INTGIMPL
                    OPTS "Matrix" "OTF"
                    LABEL "Worland integrator"
                    ADVANCED)
quicc_add_definition(QUICC_WORLAND_INTGIMPL)

# Select implementations for Worland reductors
quicc_create_option(NAME QUICC_WORLAND_REDUIMPL
                    OPTS "Matrix" "OTF"
                    LABEL "Worland reductor"
                    ADVANCED)
quicc_add_definition(QUICC_WORLAND_REDUIMPL)
