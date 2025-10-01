/**
 * @file Utils.hpp
 * @brief Utilities specific to Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_UTILS_HPP
#define QUICC_POLYNOMIAL_WORLAND_UTILS_HPP

// System includes
//
#include <string>

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

namespace Utils {

void selectJacobi(const std::string& t, Internal::MHDFloat& a, Internal::MHDFloat& db, Internal::Array& igrid, Internal::Array& iweights);

} // namespace Utils
} // namespace Worland
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_WORLAND_UTILS_HPP
