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
#include "QuICC/TestSuite/Framework/Communicators/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Communicators {

   TestArgs::TestArgs()
      : useDefault(true), dumpData(false), dim1D(0), dim2D(0), dim3D(0), algorithm(""), grouper("")
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
