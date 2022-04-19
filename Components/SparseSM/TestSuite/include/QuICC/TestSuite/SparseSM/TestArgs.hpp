/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_SPARSESM_TESTARGS_HPP
#define QUICC_TESTSUITE_SPARSESM_TESTARGS_HPP

// Configuration includes
//

// System includes
//
#include <vector>
#include <string>

// Project includes
//
#include "QuICC/TestSuite/SparseSM/TestType.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

   struct TestArgs
   {
      /// Use default test setup
      bool useDefault;

      /// Write output data to file
      bool dumpData;

      /// Test type
      TestType type;

      /// Max ulp
      unsigned int ulp;

      /// Test parameters
      std::vector<double> params;

      /**
       * @brief Constructor
       */
      TestArgs();

      /**
       * @brief Destructor
       */
      ~TestArgs() = default;

      /**
       * @brief Translate string test type to enum
       */
      void setType(const std::string& type);
   };

}
}
}

#endif //QUICC_TESTSUITE_SPARSESM_TESTARGS_HPP
