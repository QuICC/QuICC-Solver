/**
 * @file P_Zero.hpp
 * @brief Implementation of the Worland based integrator, zero for l = 0
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_P_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_P_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/P_Zero.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/P_Zero.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Integrator/P_ZeroviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Integrator/P_ZeroviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_P_ZERO_HPP
