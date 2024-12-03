/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// System includes
//

// Project includes
//
#include "TestSuite/DenseSM/Worland/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Worland {

   DenseSM::TestArgs& args()
   {
      static DenseSM::TestArgs a;

      return a;
   }

}
}
}
}
