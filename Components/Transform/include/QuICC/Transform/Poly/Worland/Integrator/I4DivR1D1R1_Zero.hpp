/**
 * @file I4DivR1D1R1_Zero.hpp
 * @brief Implementation of the Worland based I4 1/R1 D R1 integrator but 0 mode is zeroed
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I4DIVR1D1R1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I4DIVR1D1R1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/I4DivR1D1R1_Zero.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/I4DivR1D1R1_Zero.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Integrator/I4DivR1D1R1_ZeroviewCpu_t.hpp.inc"

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I4DIVR1D1R1_ZERO_HPP
