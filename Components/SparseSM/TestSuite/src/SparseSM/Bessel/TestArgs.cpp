/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// System includes
//

// Project includes
//
#include "QuICC/TestSuite/SparseSM/Bessel/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

namespace Bessel {

SparseSM::TestArgs& args()
{
   static SparseSM::TestArgs a;

   return a;
}

} // namespace Bessel
} // namespace SparseSM
} // namespace TestSuite
} // namespace QuICC
