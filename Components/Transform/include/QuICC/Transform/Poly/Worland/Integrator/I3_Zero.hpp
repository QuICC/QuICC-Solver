/**
 * @file I3_Zero.hpp
 * @brief Implementation of the Worland based I3 integrator but 0 mode is zeroed
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I3_ZERO_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I3_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/I3_Zero.hpp"
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/I3_Zero.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Integrator/I3_ZeroviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Integrator/I3_ZeroviewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_I3_ZERO_HPP
