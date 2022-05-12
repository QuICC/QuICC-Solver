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
#include "QuICC/TestSuite/Framework/LoadSplitter/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace LoadSplitter {

   TestArgs::TestArgs()
      : useDefault(true), dumpData(false), op("P"), db(0), np(0), dim1D(0), dim2D(0), dim3D(0), stage(1)
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
