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
#include "QuICC/Transform/Wrappers/Worland/Reductor/PowerviewCpu_t.hpp.inc"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/Power.hpp"
#endif
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Reductor/PowerviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWER_HPP
