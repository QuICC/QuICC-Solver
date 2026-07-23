/**
 * @file rGrad_r_th.hpp
 * @brief Implementation of the Worland based projector of radial derivatives in Grad_r_th
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_RGRAD_R_TH_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_RGRAD_R_TH_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Base/rGrad_r_th.hpp"
//#ifdef QUICC_USE_KOKKOS
//#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/rGrad_r_th.hpp"
//#endif
#include "QuICC/Transform/Wrappers/Worland/Projector/rGrad_r_thviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "QuICC/Transform/Wrappers/Worland/Projector/rGrad_r_thviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_RGRAD_R_TH_HPP
