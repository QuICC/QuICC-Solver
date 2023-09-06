/**
 * @file D1.hpp
 * @brief Implementation of the associated Legendre based D projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_D1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_D1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/D1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/D1.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Projector/D1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Projector/D1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_D1_HPP
