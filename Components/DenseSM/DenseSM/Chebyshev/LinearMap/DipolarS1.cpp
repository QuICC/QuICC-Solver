/**
 * @file DipolarS1.cpp
 * @brief Source of the implementation of the dipolar S1 field
 */

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/DipolarS1.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

int DipolarS1::nN() const
{
   return 3;
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
      Internal::MHDFloat c1 = (10.0_mp/13.0_mp)*Internal::Math::sqrt(15.0_mp/732251.0_mp);
      val = c1*(49.0_mp - 566.0_mp*r.array() + 400.0_mp*r.array().pow(2));
   }

   return val;
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
