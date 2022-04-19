/**
 * @file TestArgs.hpp
 * @brief Command line arguments for tests
 */

#ifndef QUICC_TESTSUITE_POLYNOMIAL_QUADRATURE_TESTARGS_HPP
#define QUICC_TESTSUITE_POLYNOMIAL_QUADRATURE_TESTARGS_HPP

// Configuration includes
//

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/TestSuite/Polynomial/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

namespace Quadrature {

   struct TestArgs: public Polynomial::TestArgs
   {
      /// List of grid sizes
      std::vector<std::uint32_t> gridN;
      /// List of alpha/beta pairs sizes
      std::vector<std::array<double, 2>> ab;
   };

   Quadrature::TestArgs& args();

}
}
}
}

#endif //QUICC_TESTSUITE_POLYNOMIAL_QUADRATURE_TESTARGS_HPP
