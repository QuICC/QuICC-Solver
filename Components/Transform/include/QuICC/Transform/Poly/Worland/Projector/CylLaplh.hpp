/**
 * @file CylLaplh.hpp
 * @brief Implementation of the Worland based cylindrical horizontal laplacian projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_CYLLAPLH_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_CYLLAPLH_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/CylLaplh.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/CylLaplh.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/CylLaplhviewCpu_t.hpp.inc"

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_CYLLAPLH_HPP
