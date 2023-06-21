/** 
 * @file D2.cpp
 * @brief Source of the implementation of boundary value of second derivative
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Boundary/D2.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

   D2::D2(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : ICondition(alpha, dBeta, l, 0), mBCk0(alpha, dBeta, l), mBCk1(alpha, dBeta, l+1), mBCk2(alpha, dBeta, l+2)
   {
   }

   D2::ACoeff_t D2::compute(const int maxN)
   {
      auto ab1 = this->alpha() + this->beta(this->l()) + MHD_MP(1.0);
      auto ab2 = this->alpha() + this->beta(this->l()) + MHD_MP(2.0);

      ACoeff_t n = (ACoeffI::LinSpaced(maxN-1, 2, maxN)).cast<Scalar_t>();

      ACoeff_t val = ACoeff_t::Zero(maxN+1);
      if(maxN > 2)
      {
         auto bcVal = this->mBCk2.compute(maxN-2, 2, false);
         val.bottomRows(maxN-1) = MHD_MP(4.0)*(ab1 + n)*(ab2 + n)*bcVal;
      }

      n = (ACoeffI::LinSpaced(maxN, 1, maxN)).cast<Scalar_t>();
      if(this->l() > 0)
      {
         if(maxN > 0)
         {
            auto d2l = MHD_MP(2.0)*this->l();
            auto bcVal = this->mBCk1.compute(maxN-1, 1, false);
            val.bottomRows(maxN) += MHD_MP(2.0)*(ab1 + n)*(MHD_MP(1.0) + d2l)*bcVal;
         }
      }
      else
      {
         if(maxN > 0)
         {
            auto bcVal = this->mBCk1.compute(maxN-1, 1, false);
            val.bottomRows(maxN) += MHD_MP(2.0)*(ab1 + n)*bcVal;
         }
      }

      if(this->l() > 0)
      {
         auto bcVal = this->mBCk0.compute(maxN, 0,false);
         auto ll1 = this->l()*(this->l() - MHD_MP(1.0));
         val += ll1*bcVal;
      }

      n = (ACoeffI::LinSpaced(maxN+1, 0, maxN)).cast<Scalar_t>();
      return this->invnorm(n)*val;
   }
 
} // Boundary
} // Worland
} // SparseSM
} // QuICC
