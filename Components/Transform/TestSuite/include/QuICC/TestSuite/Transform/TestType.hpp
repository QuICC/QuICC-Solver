/**
 * @file TestType.hpp
 * @brief Enum of available test types
 */

#ifndef QUICC_TESTSUITE_TRANSFORM_TESTTYPE_HPP
#define QUICC_TESTSUITE_TRANSFORM_TESTTYPE_HPP

// Configuration includes
//

// System includes
//

// Project includes
//

namespace QuICC {

namespace TestSuite {

namespace Transform {

   /**
    * @brief Test types
    */
   enum class TestType {
      PROJECTOR = 0,
      INTEGRATOR,
      REDUCTOR,
      BFLOOP,
   };

}
}
}

#endif //QUICC_TESTSUITE_TRANSFORM_TESTTYPE_HPP
