/**
 * @file DivR1CylLaplh_Zero.hpp
 * @brief Implementation of the Worland based 1/R cylindrical horizontal laplacian projector but 0 mode is zeroed
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1CYLLAPLH_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1CYLLAPLH_ZERO_HPP
// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/DivR1CylLaplh_Zero.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/DivR1CylLaplh_Zero.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/DivR1CylLaplh_ZeroviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Projector/DivR1CylLaplh_ZeroviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_DIVR1CYLLAPLH_ZERO_HPP
