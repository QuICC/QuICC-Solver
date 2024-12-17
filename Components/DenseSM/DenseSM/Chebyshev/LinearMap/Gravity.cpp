/**
 * @file Gravity.cpp
 * @brief Source of the implementation of the gravity fields
 */

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/Gravity.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

int Gravity::nN() const
{
   return 32; //manually set to the resolution, but we can change that. It modifies the internal resolution to compute the polynomial representaiton of this field.
}

Internal::Array Gravity::evaluate(const Internal::Array& r, const int l, const int m) const
{
   using namespace Internal::Literals;

   Internal::Array val;
   if(l != 0)
   {
      val = 0*r;
   }
   else
   {
      Internal::MHDFloat c1 = 1.0_mp;
      val = c1/r.array()/r.array();
   }

   return val;
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
