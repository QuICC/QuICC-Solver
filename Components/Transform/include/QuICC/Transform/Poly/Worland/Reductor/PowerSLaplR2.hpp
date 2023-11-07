/**
 * @file PowerSLaplR2.hpp
 * @brief Implementation of the Worland based spherical Laplacian R^2 power spectrum operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWERSLAPLR2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWERSLAPLR2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/PowerSLaplR2.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/PowerSLaplR2.hpp"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWERSLAPLR2_HPP
