/**
 * @file SphLapl.hpp
 * @brief Implementation of the Worland based spherical laplacian projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_SPHLAPL_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_SPHLAPL_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/SphLapl.hpp"
#if defined(QUICC_USE_KOKKOS_CUDA) || defined(QUICC_USE_KOKKOS_HIP)
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/SphLapl.hpp"
#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/SphLaplviewCpu_t.hpp.inc"

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_SPHLAPL_HPP
