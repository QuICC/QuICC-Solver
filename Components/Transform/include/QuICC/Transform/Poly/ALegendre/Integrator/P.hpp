/**
 * @file P.hpp
 * @brief Implementation of the associated Legendre based P integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_P_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_P_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/P.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/P.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/PviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/PviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_P_HPP
