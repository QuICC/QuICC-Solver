/**
 * @file P.hpp
 * @brief Implementation of the Worland based integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_P_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_P_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/Base/P.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/P.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Integrator/PviewCpu_t.hpp.inc"

#endif // QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_P_HPP
