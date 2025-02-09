/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_TIMESTEP_TESTARGS_HPP
#define QUICC_TESTSUITE_FRAMEWORK_TIMESTEP_TESTARGS_HPP

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

namespace Timestep {

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

      /// Timestepping scheme
      std::string scheme;

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

}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_TIMESTEP_TESTARGS_HPP
