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
#include "QuICC/TestSuite/Transform/Chebyshev/LinearMap/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Transform {

namespace Chebyshev {

namespace LinearMap {

   Transform::TestArgs& args()
   {
      static Transform::TestArgs a;

      return a;
   }

} // LinearMap
} // Chebyshev
} // Transform
} // TestSuite
} // QuICC
