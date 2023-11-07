/**
 * @file RadialPowerDivR1D1R1.hpp
 * @brief Implementation of the Worland based 1/R D R power spectrum operator on radial grid
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWERDIVR1D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWERDIVR1D1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/RadialPowerDivR1D1R1.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/RadialPowerDivR1D1R1.hpp"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWERDIVR1D1R1_HPP
