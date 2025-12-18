/**
 * @file Tools.hpp
 * @brief Tools specific to Chebyshev polynomial implementation
 */

#ifndef QUICC_POLYNOMIAL_CHEBYSHEV_OPERATORS_HPP
#define QUICC_POLYNOMIAL_CHEBYSHEV_OPERATORS_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Chebyshev {

namespace Operators {

/**
 * @brief Integrate r^p Tn over r
 */
void integrateRpTn(Matrix& iop, const int p,
   const int size, const MHDFloat ro, const MHDFloat ri);

} // namespace Operators
} // namespace Chebyshev
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_CHEBYSHEV_OPERATORS_HPP
