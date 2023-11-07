/**
 * @file D1CylLaplh_D1DivR1D1R1.hpp
 * @brief Implementation of the Worland based D of cylindrical horizontal laplacian projector but 0 mode is D 1/R D R
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1CYLLAPLH_D1DIVR1D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1CYLLAPLH_D1DIVR1D1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/D1CylLaplh_D1DivR1D1R1.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/D1CylLaplh_D1DivR1D1R1.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/D1CylLaplh_D1DivR1D1R1viewCpu_t.hpp.inc"

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_D1CYLLAPLH_D1DIVR1D1R1_HPP
