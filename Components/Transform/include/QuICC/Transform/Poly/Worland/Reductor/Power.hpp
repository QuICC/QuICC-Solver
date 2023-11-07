/**
 * @file Power.hpp
 * @brief Implementation of the Worland based power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWER_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWER_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/Power.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/Power.hpp"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWER_HPP
