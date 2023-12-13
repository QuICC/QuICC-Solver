/**
 * @file I2.hpp
 * @brief Implementation of the Worland based I2 integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I2_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/I2.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/I2.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Integrator/I2viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Integrator/I2viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I2_HPP
