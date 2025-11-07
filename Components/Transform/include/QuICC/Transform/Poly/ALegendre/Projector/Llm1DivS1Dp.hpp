/**
 * @file Llm1DivS1Dp.hpp
 * @brief Implementation of the associated Legendre based 1/sin [l(l+1)-1] P d_phi  projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLM1DIVS1DP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLM1DIVS1DP_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/Llm1DivS1Dp.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/Llm1DivS1Dp.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Projector/Llm1DivS1DpviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/ALegendre/Projector/Llm1DivS1DpviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLM1DIVS1DP_HPP
