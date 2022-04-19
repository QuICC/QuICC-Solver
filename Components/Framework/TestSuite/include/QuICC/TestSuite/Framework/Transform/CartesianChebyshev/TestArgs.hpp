/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_TRANSFORM_CARTESIANCHEBYSHEV_TESTARGS_HPP
#define QUICC_TESTSUITE_FRAMEWORK_TRANSFORM_CARTESIANCHEBYSHEV_TESTARGS_HPP

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

namespace CartesianChebyshev {

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

      /// Number of similar transform
      static int blockSize;

      /// ID pairs of individual blocks of transform
      static std::vector<std::pair<int,int> > idPairs;

      /// Lower bound of grid
      static double lower;

      /// Upper bound of grid
      static double upper;
   };

}
}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_TRANSFORM_CARTESIANCHEBYSHEV_TESTARGS_HPP
