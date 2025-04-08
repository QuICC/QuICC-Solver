/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_IO_TESTARGS_HPP
#define QUICC_TESTSUITE_FRAMEWORK_IO_TESTARGS_HPP

// System includes
//
#include <vector>
#include <string>

// Project includes
//

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Io {

   struct TestArgs
   {
      /// Use default test setup
      bool useDefault;

      /// Write output data to file
      bool dumpData;

      /// Only time execution, don't check data
      bool timeOnly;

      /// Max ulp
      unsigned int ulp;

      /// Number of ierations
      unsigned int iter;

      /// ID of the tests
      std::vector<int> params;

      /**
       * @brief Constructor
       */
      TestArgs();

      /**
       * @brief Destructor
       */
      ~TestArgs() = default;
   };

   TestArgs& args();

} // Variable
} // Io
} // TestSuite
} // QuICC

#endif //QUICC_TESTSUITE_FRAMEWORK_IO_TESTARGS_HPP
