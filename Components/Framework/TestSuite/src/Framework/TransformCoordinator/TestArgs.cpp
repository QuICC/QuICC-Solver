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
#include "QuICC/TestSuite/Framework/TransformCoordinator/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   TestArgs::TestArgs()
      : useDefault(true), dumpData(false), timeOnly(false), ulp(11), dim1D(0), dim2D(0), dim3D(0), algorithm(""), grouper("")
   {
   }

   TestArgs& args()
   {
      static TestArgs a;

      return a;
   }

}
}
}
}
