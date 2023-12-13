/**
 * @file DivR1D1R1_Zero.hpp
 * @brief Implementation of the Worland based 1/R D R projector and zero for l = 0
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1D1R1_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1D1R1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/DivR1D1R1_Zero.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/DivR1D1R1_Zero.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/DivR1D1R1_ZeroviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Projector/DivR1D1R1_ZeroviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1D1R1_ZERO_HPP
