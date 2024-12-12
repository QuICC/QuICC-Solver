/**
 * @file DipolarS1.cpp
 * @brief Source of the implementation of the dipolar S1 field
 */

// System includes
//

// Project includes
//
#include "DenseSM/Worland/DipolarS1.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

int DipolarS1::nN() const
{
   return 2;
}

Internal::Array DipolarS1::evaluate(const Internal::Array& r, const int l, const int m) const
{
   using namespace Internal::Literals;

   Internal::Array val;
   if(l != 1)
   {
      val = 0*r;
   }
   else
   {
      Internal::MHDFloat c1 = 0.5_mp*Internal::Math::sqrt(7.0_mp/46.0_mp);
      val = c1*r.array()*(5.0_mp - 3.0_mp*r.array().pow(2));
   }

   return val;
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
