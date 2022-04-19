/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// Configuration includes
//

// System includes
//
#include <vector>

// Project includes
// 
#include "QuICC/TestSuite/Framework/Transform/CartesianChebyshev/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Transform {

namespace CartesianChebyshev {
   
   bool TestArgs::useDefault = true;
   
   bool TestArgs::keepData = false;

   int TestArgs::specN = 0;

   int TestArgs::physN = 0;

   int TestArgs::blockSize = 0;

   std::vector<std::pair<int,int> > TestArgs::idPairs = std::vector<std::pair<int,int> >();

   double TestArgs::lower = 0.0;

   double TestArgs::upper = 0.0;

}
}
}
}
}
