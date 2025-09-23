/**
 * @file QuadrupolarS2.cpp
 * @brief Source of the implementation of the quadrupolar S2 field
 */

// System includes
//

// Project includes
//
#include "TestSuite/DenseSM/Chebyshev/LinearMap/QuadrupolarS2.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

int QuadrupolarS2::nN() const
{
   return 3;
}

Internal::Array QuadrupolarS2::evaluate(const Internal::Array& r, const int l, const int m) const
{
   using namespace Internal::Literals;

   Internal::Array val;
   if(l != 2)
   {
      val = 0*r;
   }
   else
   {
      Internal::MHDFloat c2 = (1.0_mp/(8.0_mp*Internal::Math::sqrt(587.0_mp)));
      val = c2*(35.0_mp - 200.0_mp*r.array() + 139.0_mp*r.array().pow(2));
   }

   return val;
}

std::vector<int> QuadrupolarS2::ls() const
{
   std::vector<int> l = {2};

   return l;
}

}
}
}
}
} // namespace QuICC
