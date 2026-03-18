/**
 * @file Int.cpp
 * @brief Source of the implementation of Integral across domain
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Int.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

Int::Int(const Scalar_t lower, const Scalar_t upper) :
    ICondition(lower, upper, Position::BOTTOM) // for this condition, there is no position
                                               // ICondition expects pos, so we just give something 
{}

Int::ACoeff_t Int::compute(const int maxN)
{
   ACoeff_t val = ACoeff_t::Zero(maxN + 1);

   auto cnst = this->c() * this->a();
   //auto cnst = this->a();
   for (int i = 0; i < val.size(); i++)
   {
      if (i==0)
      {
         val(i) = cnst;// * MHD_MP(2.0);
      }
      else if (i==1)
      {
         val(i) = MHD_MP(0.0);
      }
      else
      {
         const auto n = static_cast<Scalar_t>(i);

         // Compute (-1)^n + 1 based on parity
         // For even n: (-1)^n = 1, so (-1)^n + 1 = 2
         // For odd n:  (-1)^n = -1, so (-1)^n + 1 = 0
         Scalar_t onePowNPlusOne = (i % 2 == 0) ? MHD_MP(2.0) : MHD_MP(0.0);

         val(i) = cnst * onePowNPlusOne / (MHD_MP(1.0) - n * n);
      }
   }

   return val;
}

} // namespace Boundary
} // namespace LinearMap
} // namespace Chebyshev
} // namespace SparseSM
} // namespace QuICC
