/** 
 * @file Value.cpp
 * @brief Source of the implementation of boundary value
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "QuICC/SparseSM/Worland/IDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

   Value::Value(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IDiags(alpha, dBeta, l)
   {
   }

   Value::ACoeff_t Value::compute(const int maxN, const int k, const bool normalized)
   {
      ACoeff_t val = ACoeff_t::Zero(maxN+1);
      val(0) = MHD_MP(1.0);

      for(int i = 1; i <= maxN; i++)
      {
         auto d2i = static_cast<Scalar_t>(2*i);
         auto d2k = static_cast<Scalar_t>(2*k);
         val(i) = val(i-1)*(d2i - MHD_MP(1.0) + d2k)/d2i;
      }

      ACoeff_t n = (ACoeffI::LinSpaced(maxN+1, 0, maxN)).cast<Scalar_t>();

      if(normalized)
      {
         val *= this->invnorm(n);
      }

      return val;
   }
 
} // Boundary
} // Worland
} // Polynomial
} // QuICC
