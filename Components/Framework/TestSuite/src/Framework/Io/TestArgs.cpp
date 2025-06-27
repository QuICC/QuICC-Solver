/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// System includes
//

// Project includes
//
#include "QuICC/TestSuite/Framework/Io/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Io {

   TestArgs::TestArgs()
      : useDefault(true), dumpData(false), timeOnly(false), ulp(11), iter(1)
   {
   }

   TestArgs& args()
   {
      static TestArgs a;

      return a;
   }

} // Io
} // Framework
} // TestSuite
} // QuICC
