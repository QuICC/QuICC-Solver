/**
 * @file D2.cpp
 * @brief Source of the implementation of boundary value of second derivative
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/Boundary/D2.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Boundary {

   D2::D2(const BesselKind type, const int l)
      : ICondition(type, l)
   {
   }

   D2::ACoeff_t D2::compute(const int maxN)
   {
      ACoeff_t val = ACoeff_t::Ones(maxN+1);
      return val;
   }

} // Boundary
} // Bessel
} // SparseSM
} // QuICC
