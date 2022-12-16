/** 
 * @file InsulatingSphere.cpp
 * @brief Source of the implementation of boundary value
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Worland/Boundary/InsulatingSphere.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

   InsulatingSphere::InsulatingSphere(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IDiags(alpha, dBeta, l, 0), mBCk0(alpha, dBeta, l), mBCk1(alpha, dBeta, l+1)
   {
   }

   InsulatingSphere::ACoeff_t InsulatingSphere::compute(const int maxN)
   {
      auto ab1 = this->alpha() + this->beta(this->l()) + MHD_MP(1.0);

      ACoeff_t val = ACoeff_t::Zero(maxN+1);

      ACoeff_t n = (ACoeffI::LinSpaced(maxN, 1, maxN)).cast<Scalar_t>();
      if(maxN > 0)
      {
         auto bcVal = this->mBCk1.compute(maxN-1, 1, false);
         val.bottomRows(maxN) += MHD_MP(2.0)*(ab1 + n)*bcVal;
      }

      auto bcVal = this->mBCk0.compute(maxN, 0, false);
      val += (2.0*this->l() + 1.0)*bcVal;

      n = (ACoeffI::LinSpaced(maxN+1, 0, maxN)).cast<Scalar_t>();
      return this->invnorm(n)*val;
   }
 
} // Boundary
} // Worland
} // Polynomial
} // QuICC
