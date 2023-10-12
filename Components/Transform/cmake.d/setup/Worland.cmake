###################################################
#-------------- WORLAND TRANSFORM ----------------#
###################################################

# Select backend for Worland transform
quicc_create_option(NAME QUICC_WORLAND_BACKEND
                    OPTS "Eigen" "BLAS" "Scalar"
                    LABEL "Worland backend"
                    ADVANCED)
quicc_add_definition(QUICC_WORLAND_BACKEND)
