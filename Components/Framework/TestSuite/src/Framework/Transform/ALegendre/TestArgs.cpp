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
#include "QuICC/TestSuite/Framework/Transform/ALegendre/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Transform {

namespace ALegendre {
   
   bool TestArgs::useDefault = true;
   
   bool TestArgs::keepData = false;

   int TestArgs::specN = 0;

   int TestArgs::physN = 0;

   std::vector<int> TestArgs::harmonicM = std::vector<int>();

}
}
}
}
}
