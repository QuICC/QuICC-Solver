/**
 * @file D1CylLaplh.hpp
 * @brief Implementation of the Worland based D of cylindrical horizontal laplacian projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1CYLLAPLH_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1CYLLAPLH_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/D1CylLaplh.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/D1CylLaplh.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/D1CylLaplhviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Projector/D1CylLaplhviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1CYLLAPLH_HPP
