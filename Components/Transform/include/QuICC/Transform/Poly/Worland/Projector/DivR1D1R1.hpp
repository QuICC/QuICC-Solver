/**
 * @file DivR1D1R1.hpp
 * @brief Implementation of the Worland based 1/R D R projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1D1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/DivR1D1R1.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/DivR1D1R1.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/DivR1D1R1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Projector/DivR1D1R1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1D1R1_HPP
