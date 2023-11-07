/**
 * @file RadialPower.hpp
 * @brief Implementation of the Worland based power spectrum operator on radial grid
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWER_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWER_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/RadialPower.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/RadialPower.hpp"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWER_HPP
