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
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/EnergyD1R1.hpp"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYD1R1_HPP
