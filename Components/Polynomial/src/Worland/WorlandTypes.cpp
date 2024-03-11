/**
 * @file WorlandTypes.cpp
 * @brief Source of the Worland types types
 */

// Project includes
//
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   // Chebyshev type
   const Internal::MHDFloat worland_chebyshev_t::ALPHA= -MHD_MP(0.5);
   const Internal::MHDFloat worland_chebyshev_t::DBETA= -MHD_MP(0.5);

   // Legendre type
   const Internal::MHDFloat worland_legendre_t::ALPHA= MHD_MP(0.0);
   const Internal::MHDFloat worland_legendre_t::DBETA= -MHD_MP(0.5);

   // CylEnergy type
   const Internal::MHDFloat worland_cylenergy_t::ALPHA= MHD_MP(0.0);
   const Internal::MHDFloat worland_cylenergy_t::DBETA= MHD_MP(0.0);

   // SphEnergy type
   const Internal::MHDFloat worland_sphenergy_t::ALPHA= MHD_MP(0.0);
   const Internal::MHDFloat worland_sphenergy_t::DBETA= MHD_MP(0.5);

#if defined QUICC_WORLAND_TYPE_CHEBYSHEV
   #define QUICC_WORLAND_DEFAULT_TYPE worland_chebyshev_t
#elif defined QUICC_WORLAND_TYPE_LEGENDRE
   #define QUICC_WORLAND_DEFAULT_TYPE worland_legendre_t
#elif defined QUICC_WORLAND_TYPE_CYLENERGY
   #define QUICC_WORLAND_DEFAULT_TYPE worland_cylenergy_t
#elif defined QUICC_WORLAND_TYPE_SPHENERGY
   #define QUICC_WORLAND_DEFAULT_TYPE worland_sphenergy_t
#else
   #error "QUICC_WORLAND_TYPE_? is not defined"
#endif //QUICC_WORLAND_TYPE_CHEBYSHEV

   // default type
   const Internal::MHDFloat worland_default_t::ALPHA= QUICC_WORLAND_DEFAULT_TYPE::ALPHA;
   const Internal::MHDFloat worland_default_t::DBETA= QUICC_WORLAND_DEFAULT_TYPE::DBETA;

#undef QUICC_WORLAND_DEFAULT_TYPE
}
}
}

