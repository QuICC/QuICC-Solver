/**
 * @file Ll.hpp
 * @brief Implementation of the associated Legendre based l(l+1) P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LL_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LL_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/Ll.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/Ll.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Projector/LlviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Projector/LlviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LL_HPP
