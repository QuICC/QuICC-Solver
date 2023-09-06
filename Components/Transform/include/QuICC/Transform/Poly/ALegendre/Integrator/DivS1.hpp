/**
 * @file DivS1.hpp
 * @brief Implementation of the associated Legendre based 1/Sin integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVS1_HPP
// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/DivS1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/DivS1.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/DivS1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/DivS1viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVS1_HPP
