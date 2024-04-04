/**
 * @file RadialPowerDivR1.hpp
 * @brief Implementation of the Worland based 1/R power spectrum operator on radial grid
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWERDIVR1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWERDIVR1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/RadialPowerDivR1.hpp"
#include "QuICC/Transform/Wrappers/Worland/Reductor/RadialPowerDivR1viewCpu_t.hpp.inc"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/RadialPowerDivR1.hpp"
#endif
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Reductor/RadialPowerDivR1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_RADIALPOWERDIVR1_HPP
