/**
 * @file Ll.hpp
 * @brief Implementation of the associated Legendre based l(l+1) P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LL_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LL_HPP
// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/Ll.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/Ll.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/LlviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/LlviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LL_HPP
