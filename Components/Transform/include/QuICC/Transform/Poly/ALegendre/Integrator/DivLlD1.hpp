/**
 * @file DivLlD1.hpp
 * @brief Implementation of the associated Legendre based 1/l(l+1) D integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLLD1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLLD1_HPP
// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/DivLlD1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/DivLlD1.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/DivLlD1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/DivLlD1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLLD1_HPP
