/**
 * @file TestUtils.hpp
 * @brief Some utils
 */

#ifndef QUICC_TESTSUITE_SPARSESM_TESTUTILS_HPP
#define QUICC_TESTSUITE_SPARSESM_TESTUTILS_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

   /**
    * @Brief Compute ULP
    */
   std::pair<bool, MHDFloat> computeUlp(const ::QuICC::MHDFloat data, const ::QuICC::MHDFloat ref, const ::QuICC::MHDFloat refMod, const ::QuICC::MHDFloat maxUlp, const ::QuICC::MHDFloat epsilon);

   /**
    * @brief Check stencils by computing boundary values
    */
   void checkStencil(const Matrix& bc, const Matrix& stencil, const MHDFloat maxUlp);

}
}
}

#endif //QUICC_TESTSUITE_SPARSESM_TESTUTILS_HPP
