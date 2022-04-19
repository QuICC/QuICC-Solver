###################################################
#-------- ASSOCIATED LEGENDRE TRANSFORM ----------#
###################################################

# Select implementations for ALegendre projectors
quicc_create_option(NAME QUICC_ALEGENDRE_PROJIMPL
                    OPTS "Matrix" "OTF"
                    LABEL "ALegendre projector"
                    ADVANCED)
quicc_add_definition(QUICC_ALEGENDRE_PROJIMPL)


# Select implementations for ALegendre integrators
quicc_create_option(NAME QUICC_ALEGENDRE_INTGIMPL
                    OPTS "Matrix" "OTF"
                    LABEL "ALegendre integrator"
                    ADVANCED)
quicc_add_definition(QUICC_ALEGENDRE_INTGIMPL)
