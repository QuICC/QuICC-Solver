/**
 * @file D1R1.hpp
 * @brief Implementation of the Worland based D R projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/D1R1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/D1R1.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/D1R1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Projector/D1R1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1R1_HPP
