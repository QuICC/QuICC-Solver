/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_COMMUNICATORS_TESTARGS_HPP
#define QUICC_TESTSUITE_FRAMEWORK_COMMUNICATORS_TESTARGS_HPP

// Configuration includes
//

// System includes
//
#include <vector>
#include <string>

// Project includes
//

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Communicators {

   struct TestArgs
   {
      /// Use default test setup
      bool useDefault;

      /// Write output data to file
      bool dumpData;

      /// Dimension 1D
      unsigned int dim1D;

      /// Dimension 1D
      unsigned int dim2D;

      /// Dimension 1D
      unsigned int dim3D;

      /// Comm splitting algorithm
      std::string algorithm;

      /// Comm grouping algorithm
      std::string grouper;

      /// ???
      std::vector<double> params;

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

}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_COMMUNICATORS_TESTARGS_HPP
