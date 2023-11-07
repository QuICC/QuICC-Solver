/**
 * @file EnergyR2.hpp
 * @brief Implementation of the Worland based R^2 energy operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYR2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYR2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/EnergyR2.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/EnergyR2.hpp"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_ENERGYR2_HPP
