/**
 * @file Literals.hpp
 * @brief Definition of multiple precision literals
 */

#ifndef QUICC_TYPES_MP_LITERALS_HPP
#define QUICC_TYPES_MP_LITERALS_HPP

// System includes
//

// Project includes
//
#include "Types/MP/BasicTypes.hpp"

namespace QuICC {
namespace MP {
/// @brief namespace for literal operators
namespace Literals {

/// @brief literal operator allow for parsing MP constants as string
/// @param str literal
/// @return MHDFloat
inline const MHDFloat operator""_mp(const char* str)
{
   return MHDFloat(str);
}

} // namespace Literals
} // namespace MP
} // namespace QuICC

#endif // QUICC_TYPES_MP_LITERALS_HPP
