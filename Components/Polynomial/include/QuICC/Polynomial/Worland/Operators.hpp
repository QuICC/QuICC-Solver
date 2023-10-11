/**
 * @file Tools.hpp
 * @brief Tools specific to Worland polynomial implementation
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_OPERATORS_HPP
#define QUICC_POLYNOMIAL_WORLAND_OPERATORS_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

namespace Operators {

/**
* @brief Integrate r^p Wnl over r
*/
void integrateRpWnl(Internal::Matrix& iop, const int l, const int p, const int size);

}
}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_OPERATORS_HPP
