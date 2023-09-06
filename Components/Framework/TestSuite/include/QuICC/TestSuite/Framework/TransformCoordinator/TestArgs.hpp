/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_TCOORD_TESTARGS_HPP
#define QUICC_TESTSUITE_FRAMEWORK_TCOORD_TESTARGS_HPP

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

namespace TCoord {

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

}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_TCOORD_TESTARGS_HPP
