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
#include "QuICC/TestSuite/Polynomial/Worland/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

namespace Worland {

   Polynomial::TestArgs& args()
   {
      static Polynomial::TestArgs a;

      return a;
   }

}
}
}
}
