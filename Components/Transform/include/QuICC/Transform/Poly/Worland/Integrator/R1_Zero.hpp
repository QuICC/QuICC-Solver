/**
 * @file R1_Zero.hpp
 * @brief Implementation of the Worland based R1_Zero integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_R1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_R1_ZERO_HPP
// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/R1_Zero.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/R1_Zero.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Integrator/R1_ZeroviewCpu_t.hpp.inc"

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_R1_ZERO_HPP
