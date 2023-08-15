/**
 * @file D1.hpp
 * @brief Implementation of the associated Legendre based D integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_D1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_D1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/D1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/D1.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/D1viewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/D1viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_D1_HPP
