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
std::pair<bool, MHDFloat> computeUlp(const MHDFloat data, const MHDFloat ref,
   const MHDFloat refMod, const MHDFloat maxUlp, const MHDFloat epsilon);

/**
 * @brief Check stencils by computing boundary values
 */
void checkStencil(const Matrix& bc, const Matrix& stencil,
   const MHDFloat maxUlp);

} // namespace SparseSM
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_SPARSESM_TESTUTILS_HPP
