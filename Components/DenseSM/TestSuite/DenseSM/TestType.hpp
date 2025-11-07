/**
 * @file TestType.hpp
 * @brief Enum of available test types
 */

#ifndef QUICC_TESTSUITE_DENSESM_TESTTYPE_HPP
#define QUICC_TESTSUITE_DENSESM_TESTTYPE_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

   /**
    * @brief Test types
    */
   enum class TestType {
      DENSE = 0,
      SPARSE,
   };

}
}
}

#endif //QUICC_TESTSUITE_DENSESM_TESTTYPE_HPP
