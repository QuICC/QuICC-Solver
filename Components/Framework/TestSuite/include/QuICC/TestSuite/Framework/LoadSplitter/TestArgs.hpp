/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_LOADSPLITTER_TESTARGS_HPP
#define QUICC_TESTSUITE_FRAMEWORK_LOADSPLITTER_TESTARGS_HPP

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

namespace LoadSplitter {

   struct TestArgs
   {
      /// Use default test setup
      bool useDefault;

      /// Write output data to file
      bool dumpData;

      /// operator name
      std::string op;

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

      /// MPI stage to look at
      unsigned int stage;

      /// MPI ranks for which to compute splitting
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

#endif //QUICC_TESTSUITE_FRAMEWORK_LOADSPLITTER_TESTARGS_HPP
