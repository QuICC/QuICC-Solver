/**
 * @file LlDivS1Dp.hpp
 * @brief Implementation of the associated Legendre based 1/sin l(l+1) P d_phi P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLDIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLDIVS1DP_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/LlDivS1Dp.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/LlDivS1Dp.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Projector/LlDivS1DpviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/ALegendre/Projector/LlDivS1DpviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLDIVS1DP_HPP
