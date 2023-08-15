/**
 * @file DivLlDivS1.hpp
 * @brief Implementation of the associated Legendre based 1/l(l+1) 1/Sin integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLLDIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLLDIVS1_HPP
// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/DivLlDivS1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/DivLlDivS1.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/DivLlDivS1viewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/DivLlDivS1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVLLDIVS1_HPP
