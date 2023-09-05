/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_STATEFILE_TESTARGS_HPP
#define QUICC_TESTSUITE_FRAMEWORK_STATEFILE_TESTARGS_HPP

// System includes
//
#include <vector>
#include <string>

// Project includes
//

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace StateFile {

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

      /// Dimension 1D
      unsigned int dim1D;

      /// Dimension 1D
      unsigned int dim2D;

      /// Dimension 1D
      unsigned int dim3D;

      /// Number of ierations
      unsigned int iter;

      /// Comm splitting algorithm
      std::string algorithm;

      /// Comm grouping algorithm
      std::string grouper;

      /// ID of the tests
      std::vector<int> params;

      /// Imposed CPU factors
      std::vector<int> factors;

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
} // StateFile
} // TestSuite
} // QuICC

#endif //QUICC_TESTSUITE_FRAMEWORK_STATEFILE_TESTARGS_HPP
