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
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/CylLaplh_DivR1D1R1.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/CylLaplh_DivR1D1R1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Projector/CylLaplh_DivR1D1R1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_CYLLAPLH_DIVR1D1R1_HPP
