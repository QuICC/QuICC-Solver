# Normalization of Worland polynomials
quicc_create_option(NAME QUICC_WORLAND_NORM
                    OPTS "Unity" "Natural"
                    LABEL "Worland normalization"
                    ADVANCED)
quicc_add_definition(QUICC_WORLAND_NORM)

# Type of Worland polynomials
quicc_create_option(NAME QUICC_WORLAND_TYPE
                    OPTS "Chebyshev" "Legendre" "CylEnergy" "SphEnergy"
                    LABEL "Worland type"
                    ADVANCED)
quicc_add_definition(QUICC_WORLAND_TYPE)
