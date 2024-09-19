/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// System includes
//

// Project includes
//
#include "QuICC/TestSuite/Polynomial/Bessel/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

namespace Bessel {

Polynomial::TestArgs& args()
{
   static Polynomial::TestArgs a;

   return a;
}

} // namespace Bessel
} // namespace Polynomial
} // namespace TestSuite
} // namespace QuICC
