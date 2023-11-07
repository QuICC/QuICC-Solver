/**
 * @file CylLaplh_DivR1D1R1.hpp
 * @brief Implementation of the Worland based cylindrical horizontal laplacian projector but 0 mode 1/R D R projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_CYLLAPLH_DIVR1D1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_CYLLAPLH_DIVR1D1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/CylLaplh_DivR1D1R1.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/CylLaplh_DivR1D1R1.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/CylLaplh_DivR1D1R1viewCpu_t.hpp.inc"

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_CYLLAPLH_DIVR1D1R1_HPP
