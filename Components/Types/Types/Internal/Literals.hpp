/**
 * @file Literals.hpp
 * @brief Definition of internal literals
 */

#ifndef QUICC_TYPES_INTERNAL_LITERALS_HPP
#define QUICC_TYPES_INTERNAL_LITERALS_HPP

// System includes
//

// Project includes
//
#ifdef QUICC_MULTPRECISION
#include "Types/MP/Literals.hpp"
#endif

namespace QuICC {
namespace Internal {
/// @brief namespace for literal operators
namespace Literals {

#ifdef QUICC_MULTPRECISION
/// @brief propagate MP literals
using namespace MP::Literals;
#else
/// @brief literal operator allow for parsing MP constants as string
/// when internal is not MP
/// @param str literal
/// @return MHDFloat
inline const QuICC::MHDFloat operator""_mp(const char* str)
{
   return std::stod(str);
}
#endif


} // namespace Literals
} // namespace Internal
} // namespace QuICC

#endif // QUICC_TYPES_INTERNAL_LITERALS_HPP
