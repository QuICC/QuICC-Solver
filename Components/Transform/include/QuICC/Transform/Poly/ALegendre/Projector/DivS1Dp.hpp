/**
 * @file DivS1Dp.hpp
 * @brief Implementation of the associated Legendre based 1/sin P d_phi projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_DIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_DIVS1DP_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/DivS1Dp.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/DivS1Dp.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Projector/DivS1DpviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/ALegendre/Projector/DivS1DpviewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_DIVS1DP_HPP
