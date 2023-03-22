/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// System includes
//

// Project includes
//
#include "QuICC/TestSuite/SparseSM/Worland/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

namespace Worland {

   SparseSM::TestArgs& args()
   {
      static SparseSM::TestArgs a;

      return a;
   }

}
}
}
}
