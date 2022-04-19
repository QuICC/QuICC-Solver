/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_TRANSFORMCONFIGURATORS_TRANSFORMTREE_TESTARGS_HPP
#define QUICC_TESTSUITE_FRAMEWORK_TRANSFORMCONFIGURATORS_TRANSFORMTREE_TESTARGS_HPP

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

namespace TransformConfigurators {

namespace TransformTree {

   struct TestArgs
   {
      /// Use default test setup
      static bool useDefault;

      /// Write output data to file
      static bool keepData;
   };

}
}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_TRANSFORMCONFIGURATORS_TRANSFORMTREE_TESTARGS_HPP
