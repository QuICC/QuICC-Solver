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
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/D1_P.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/D1_PviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Projector/D1_PviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1_P_HPP
