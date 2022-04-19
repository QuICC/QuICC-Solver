/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_TRANSFORM_ALEGENDRE_TESTARGS_HPP
#define QUICC_TESTSUITE_FRAMEWORK_TRANSFORM_ALEGENDRE_TESTARGS_HPP

// Configuration includes
//

// System includes
//
#include <vector>

// Project includes
//

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Transform {

namespace ALegendre {

   struct TestArgs
   {
      /// Use default test setup
      static bool useDefault;

      /// Write output data to file
      static bool keepData;

      /// Spectral truncation
      static int specN;

      /// Physical grid size
      static int physN;

      /// Harmonic order
      static std::vector<int> harmonicM;
   };

}
}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_TRANSFORM_ALEGENDRE_TESTARGS_HPP
