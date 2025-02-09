/**
 * @file Ll2.hpp
 * @brief Implementation of the associated Legendre based l(l+1)^2 P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LL2_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LL2_HPP
// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/Ll2.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/Ll2.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/Ll2viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/Ll2viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LL2_HPP
