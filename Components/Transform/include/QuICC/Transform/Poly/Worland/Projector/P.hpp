/**
 * @file P.hpp
 * @brief Implementation of the Worland based P projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_P_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_P_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/P.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/P.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/PviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Projector/PviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_P_HPP
