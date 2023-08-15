/**
 * @file DivS1Dp.hpp
 * @brief Implementation of the associated Legendre based 1/Sin D_phi integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVS1DP_HPP
// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Base/DivS1Dp.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/DivS1Dp.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/DivS1DpviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/DivS1DpviewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_DIVS1DP_HPP
