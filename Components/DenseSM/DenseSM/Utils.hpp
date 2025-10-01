/**
 * @file Utils.hpp
 * @brief Some utils shared by different basis
 */

#ifndef QUICC_DENSESM_UTILS_HPP
#define QUICC_DENSESM_UTILS_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Utils {

   /**
    * @brief Compute Gaunt's integral
    */
   MHDFloat gaunt(const int lA, const int mA, const int lB, const int mB,
         const int lG, const int mG);

   /**
    * @brief Compute Elsasser's integral
    */
   MHDFloat elsasser(const int lA, const int mA, const int lB, const int mB,
         const int lG, const int mG);

} // namespace Utils
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_UTILS_HPP
