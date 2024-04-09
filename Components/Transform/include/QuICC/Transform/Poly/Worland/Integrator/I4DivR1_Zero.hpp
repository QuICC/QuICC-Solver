/**
 * @file I4DivR1_Zero.hpp
 * @brief Implementation of the Worland based I4 1/R1 integrator but 0 mode is zeroed
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I4DIVR1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I4DIVR1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/I4DivR1_Zero.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/I4DivR1_Zero.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Integrator/I4DivR1_ZeroviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Integrator/I4DivR1_ZeroviewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I4DIVR1_ZERO_HPP
