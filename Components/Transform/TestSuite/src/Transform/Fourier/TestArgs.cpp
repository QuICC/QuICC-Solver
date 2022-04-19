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
#include "QuICC/TestSuite/Transform/Fourier/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Transform {

namespace Fourier {

   Transform::TestArgs& args()
   {
      static Transform::TestArgs a;

      return a;
   }

}
}
}
}
