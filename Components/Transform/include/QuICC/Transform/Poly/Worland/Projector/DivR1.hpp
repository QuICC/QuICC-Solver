/**
 * @file DivR1.hpp
 * @brief Implementation of the Worland based 1/R projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/DivR1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/DivR1.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/DivR1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Projector/DivR1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1_HPP
