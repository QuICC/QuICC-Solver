/**
 * @file D1_P.hpp
 * @brief Implementation of the Worland based D projector but 0 mode is P projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1_P_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1_P_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/D1_P.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/D1_P.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/D1_PviewCpu_t.hpp.inc"

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1_P_HPP
