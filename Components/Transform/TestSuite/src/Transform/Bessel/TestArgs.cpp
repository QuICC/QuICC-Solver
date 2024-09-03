/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// System includes
//

// Project includes
//
#include "QuICC/TestSuite/Transform/Bessel/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Transform {

namespace Bessel {

   Transform::TestArgs& args()
   {
      static Transform::TestArgs a;

      return a;
   }

}
}
}
}
