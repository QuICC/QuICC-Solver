/**
 * @file PowerR2.hpp
 * @brief Implementation of the Worland based R^2 power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWERR2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWERR2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/PowerR2.hpp"
#include "QuICC/Transform/Wrappers/Worland/Reductor/PowerR2viewCpu_t.hpp.inc"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/PowerR2.hpp"
#endif
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Reductor/PowerR2viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWERR2_HPP
