/**
 * @file DivR1D1R1_Zero.hpp
 * @brief Implementation of the Worland based 1/R1 D R1 integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_DIVR1D1R1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_DIVR1D1R1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/DivR1D1R1_Zero.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/DivR1D1R1_Zero.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Integrator/DivR1D1R1_ZeroviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Integrator/DivR1D1R1_ZeroviewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_DIVR1D1R1_ZERO_HPP
