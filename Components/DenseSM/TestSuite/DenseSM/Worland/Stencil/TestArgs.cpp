/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// System includes
//

// Project includes
//
#include "TestSuite/DenseSM/Worland/Stencil/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Worland {

namespace Stencil {

DenseSM::TestArgs& args()
{
   static DenseSM::TestArgs a;

   return a;
}

} // namespace Stencil
} // namespace Worland
} // namespace DenseSM
} // namespace TestSuite
} // namespace QuICC
