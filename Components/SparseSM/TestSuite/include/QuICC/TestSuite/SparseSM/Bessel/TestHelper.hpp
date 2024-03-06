/**
 * @file TestHelper.hpp
 * @brief Some helpers for Worland tests
 */

#ifndef QUICC_TESTSUITE_SPARSESM_WORLAND_TESTHELPER_HPP
#define QUICC_TESTSUITE_SPARSESM_WORLAND_TESTHELPER_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

namespace Worland {

/**
 * @Brief Compute ULP
 */
void setJacobiParameters(Internal::MHDFloat& alpha, Internal::MHDFloat& dBeta,
   const std::string& wtype);

} // namespace Worland
} // namespace SparseSM
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_SPARSESM_WORLAND_TESTHELPER_HPP
