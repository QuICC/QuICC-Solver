/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_POLYNOMIAL_TESTARGS_HPP
#define QUICC_TESTSUITE_POLYNOMIAL_TESTARGS_HPP

// Configuration includes
//

// System includes
//
#include <array>
#include <vector>
#include <string>

// Project includes
//
#include "QuICC/TestSuite/Polynomial/TestType.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

   struct TestArgs
   {
      /// Use default test setup
      bool useDefault;

      /// Write output data to file
      bool dumpData;

      /// Test type
      TestType type;

      /// Spectral truncation
      int specN;

      /// Physical grid size
      int physN;

      /// Max ulp
      unsigned int ulp;

      /// Test metadata IDs
      std::vector<int> ids;

      /// Harmonic degree
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

#endif //QUICC_TESTSUITE_POLYNOMIAL_TESTARGS_HPP
