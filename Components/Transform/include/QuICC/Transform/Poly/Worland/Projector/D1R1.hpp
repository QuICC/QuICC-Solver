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
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/D1R1.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/D1R1viewCpu_t.hpp.inc"

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1R1_HPP
