/**
 * @file TestType.hpp
 * @brief Enum of available test types
 */

#ifndef QUICC_TESTSUITE_POLYNOMIAL_TESTTYPE_HPP
#define QUICC_TESTSUITE_POLYNOMIAL_TESTTYPE_HPP

// Configuration includes
//

// System includes
//

// Project includes
//

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

   /**
    * @brief Test types
    */
   enum class TestType {
      QUADRATURE = 0,
      MATRIX,
      WEIGHTED_MATRIX,
      OTF_INNER,
      OTF_OUTER,
      OTF_REDUCE
   };

}
}
}

#endif //QUICC_TESTSUITE_POLYNOMIAL_TESTTYPE_HPP
