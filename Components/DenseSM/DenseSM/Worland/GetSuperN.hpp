/**
 * @file GetSuperN.hpp
 * @brief Utility for Worland Builder
 */

#ifndef QUICC_DENSESM_WORLAND_GETSUPERN_HPP
#define QUICC_DENSESM_WORLAND_GETSUPERN_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"
#include "QuICC/SparseSM/Worland/I6.hpp"

namespace QuICC {
namespace DenseSM {
namespace Worland {

/// @brief Get number of super diagonals based on builder type
/// @tparam TINBuilder
/// @return number of super diagonals
template <class TINBuilder> constexpr std::uint32_t getSuperN()
{
#ifdef QUICC_TRANSFORM_WORLAND_TRUNCATE_QI
   return 0;
#else
   if constexpr (std::is_same_v<TINBuilder, ::QuICC::SparseSM::Worland::I2>)
   {
      // I2 has 3 super diagonals
      return 3;
   }
   else if constexpr (std::is_same_v<TINBuilder,
                         ::QuICC::SparseSM::Worland::I4>)
   {
      // I4 has 6 super diagonals
      return 6;
   }
   else if constexpr (std::is_same_v<TINBuilder,
                         ::QuICC::SparseSM::Worland::I6>)
   {
      // I6 has 9 super diagonals
      return 9;
   }
   else
   {
      throw std::logic_error("Unknow operator.");
   }
#endif // QUICC_TRANSFORM_WORLAND_TRUNCATE_QI
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // define QUICC_DENSESM_WORLAND_GETSUPERN_HPP
