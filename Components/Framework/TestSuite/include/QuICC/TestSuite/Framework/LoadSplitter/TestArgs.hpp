/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_LOADSPLITTER_TESTARGS_HPP
#define QUICC_TESTSUITE_FRAMEWORK_LOADSPLITTER_TESTARGS_HPP

// System includes
//
#include <vector>
#include <string>

// Project includes
//

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace LoadSplitter {

   struct TestArgs
   {
      /// Use default test setup
      bool useDefault;

      /// Write output data to file
      bool dumpData;

      /// operator name
      std::string op;

      /// algorithm name
      std::string algorithm;

      /// Truncation scheme
      std::string truncation;

      /// ID of database file
      unsigned int db;

      /// Number of MPI ranks used for splitting
      unsigned int np;

      /// Dimension 1D
      unsigned int dim1D;

      /// Dimension 1D
      unsigned int dim2D;

      /// Dimension 1D
      unsigned int dim3D;

      /** 
       * @brief MPI stage to look at
       * 0: radial, 1: AL, 2:FFT, 3: spectral
       */
      unsigned int stage;

      /**
       * @brief Check rank decomposition individually
       */
      bool checkRanks;

      /// MPI ranks for which to compute splitting
      std::vector<double> params;

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

#endif //QUICC_TESTSUITE_FRAMEWORK_LOADSPLITTER_TESTARGS_HPP
