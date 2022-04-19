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
#include "QuICC/TestSuite/Polynomial/Jacobi/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

namespace Jacobi {
   
   Polynomial::TestArgs& args()
   {
      static Polynomial::TestArgs a;

      return a;
   }

}
}
}
}
