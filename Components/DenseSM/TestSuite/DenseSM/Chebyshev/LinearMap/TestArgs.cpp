/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// System includes
//

// Project includes
//
#include "TestSuite/DenseSM/Chebyshev/LinearMap/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

   DenseSM::TestArgs& args()
   {
      static DenseSM::TestArgs a;

      return a;
   }

}
}
}
}
}
