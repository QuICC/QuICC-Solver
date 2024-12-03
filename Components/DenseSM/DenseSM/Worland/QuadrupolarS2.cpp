/**
 * @file QuadrupolarS2.cpp
 * @brief Source of the implementation of the quadrupolar S2 field
 */

// System includes
//

// Project includes
//
#include "DenseSM/Worland/QuadrupolarS2.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

int QuadrupolarS2::nN() const
{
   return 3;
}

Internal::Array QuadrupolarS2::evaluate(const Internal::Array& r, const int l) const
{
   using namespace Internal::Literals;

   Internal::Array val;
   if(l != 2)
   {
      val = 0*r;
   }
   else
   {
      Internal::MHDFloat c2 = (5.0_mp/14.0_mp)*Internal::Math::sqrt(3.0_mp/182.0_mp);
      val = c2*r.array().pow(2)*(157.0_mp - 296.0_mp*r.array().pow(2) + 143.0_mp*r.array().pow(4));
   }

   return val;
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
