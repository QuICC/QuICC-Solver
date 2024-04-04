/**
 * @file EnergyD1R1.hpp
 * @brief Implementation of the Worland based D R energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYD1R1_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYD1R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/EnergyD1R1.hpp"
#include "QuICC/Transform/Wrappers/Worland/Reductor/EnergyD1R1viewCpu_t.hpp.inc"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/EnergyD1R1.hpp"
#endif
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Reductor/EnergyD1R1viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYD1R1_HPP
