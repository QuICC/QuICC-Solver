/**
 * @file D1DivS1.hpp
 * @brief Implementation of the associated Legendre based D1(1/sin P) projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_D1DIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_D1DIVS1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/D1DivS1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/D1DivS1.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Projector/D1DivS1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Projector/D1DivS1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_D1DIVS1_HPP
