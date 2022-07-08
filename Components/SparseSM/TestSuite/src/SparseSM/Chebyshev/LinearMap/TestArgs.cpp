/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// Configuration includes
//

// System includes
//

// Project includes
//
#include "QuICC/TestSuite/SparseSM/Worland/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   SparseSM::TestArgs& args()
   {
      static SparseSM::TestArgs a;

      return a;
   }

}
}
}
}
}
