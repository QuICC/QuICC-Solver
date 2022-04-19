/** 
 * @file Tools.hpp
 * @brief Tools specific to Worland polynomial implementation
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_OPERATORS_HPP
#define QUICC_POLYNOMIAL_WORLAND_OPERATORS_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Precision.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

namespace Operators {

/**
* @brief Integrate r^p Wnl over r
*/
void integrateRpWnl(internal::Matrix& iop, const int l, const int p, const int size);

}
}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_OPERATORS_HPP
