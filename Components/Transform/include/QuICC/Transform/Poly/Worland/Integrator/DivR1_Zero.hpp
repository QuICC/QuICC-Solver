/**
 * @file DivR1_Zero.hpp
 * @brief Implementation of the Worland based 1/R1 integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_DIVR1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_DIVR1_ZERO_HPP
// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/DivR1_Zero.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/DivR1_Zero.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Integrator/DivR1_ZeroviewCpu_t.hpp.inc"

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_DIVR1_ZERO_HPP
