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
#include "QuICC/Transform/Wrappers/Worland/Reductor/PowerSLaplR2viewCpu_t.hpp.inc"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/PowerSLaplR2.hpp"
#endif
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Reductor/PowerSLaplR2viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_POWERSLAPLR2_HPP
