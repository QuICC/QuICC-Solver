/**
 * @file TestType.hpp
 * @brief Enum of available test types
 */

#ifndef QUICC_TESTSUITE_SPARSESM_TESTTYPE_HPP
#define QUICC_TESTSUITE_SPARSESM_TESTTYPE_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

   /**
    * @brief Test types
    */
   enum class TestType {
      SPARSE = 0,
      BANDED,
      BOUNDARY,
      STENCIL,
   };

}
}
}

#endif //QUICC_TESTSUITE_SPARSESM_TESTTYPE_HPP
