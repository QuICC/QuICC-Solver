/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// Configuration includes
//

// System includes
//
#include <limits>

// Project includes
// 
#include "QuICC/TestSuite/Polynomial/Quadrature/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

namespace Quadrature {
   
   Quadrature::TestArgs& args()
   {
      static Quadrature::TestArgs a;

      return a;
   }

}
}
}
}
