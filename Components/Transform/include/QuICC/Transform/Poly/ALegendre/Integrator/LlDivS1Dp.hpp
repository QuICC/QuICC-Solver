/**
 * @file LlDivS1Dp.hpp
 * @brief Implementation of the associated Legendre based l(l+1)/Sin D_phi integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LLDIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LLDIVS1DP_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/LlDivS1Dp.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/LlDivS1Dp.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/LlDivS1DpviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/LlDivS1DpviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_LLDIVS1DP_HPP
