/**
 * @file LlDivS1.hpp
 * @brief Implementation of the associated Legendre based 1/sin l(l+1) P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLDIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLDIVS1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/LlDivS1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/LlDivS1.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Projector/LlDivS1viewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/ALegendre/Projector/LlDivS1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_LLDIVS1_HPP
