/**
 * @file LlD1.hpp
 * @brief Implementation of the associated Legendre based l(l+1) D projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLD1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLD1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/LlD1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/LlD1.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Projector/LlD1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Projector/LlD1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLD1_HPP
