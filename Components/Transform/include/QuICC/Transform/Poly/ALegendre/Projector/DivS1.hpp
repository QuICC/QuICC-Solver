/**
 * @file DivS1.hpp
 * @brief Implementation of the associated Legendre based 1/sin P projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_DIVS1_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_DIVS1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Base/DivS1.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/DivS1.hpp"
#endif
#include "QuICC/Transform/Wrappers/ALegendre/Projector/DivS1viewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/ALegendre/Projector/DivS1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_DIVS1_HPP
