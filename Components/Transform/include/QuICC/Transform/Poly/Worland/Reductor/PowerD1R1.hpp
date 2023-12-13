/**
 * @file PowerD1R1.hpp
 * @brief Implementation of the Worland based D R power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWERD1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWERD1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/PowerD1R1.hpp"
#include "QuICC/Transform/Wrappers/Worland/Reductor/PowerD1R1viewCpu_t.hpp.inc"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/PowerD1R1.hpp"
#endif
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Reductor/PowerD1R1viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWERD1R1_HPP
