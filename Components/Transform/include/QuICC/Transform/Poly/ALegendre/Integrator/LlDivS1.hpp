/**
 * @file LlDivS1.hpp
 * @brief Implementation of the associated Legendre based l(l+1) 1/Sin integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LLDIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LLDIVS1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/LlDivS1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/LlDivS1.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/LlDivS1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/LlDivS1viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LLDIVS1_HPP
