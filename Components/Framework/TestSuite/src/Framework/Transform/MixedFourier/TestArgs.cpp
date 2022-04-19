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
#include "QuICC/TestSuite/Framework/Transform/MixedFourier/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Transform {

namespace MixedFourier {
   
   bool TestArgs::useDefault = true;
   
   bool TestArgs::keepData = false;

   int TestArgs::specN = 0;

   int TestArgs::physN = 0;

   int TestArgs::blockSize = 0;

   std::vector<std::pair<int,int> > TestArgs::idPairs = std::vector<std::pair<int,int> >();

}
}
}
}
}
