/**
 * @file D2.hpp
 * @brief Implementation of the associated Legendre based D2 projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_D2_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_D2_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/D2.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/D2.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Projector/D2viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Projector/D2viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_D2_HPP
