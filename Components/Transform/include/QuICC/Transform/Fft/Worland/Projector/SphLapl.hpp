/**
 * @file SphLapl.hpp
 * @brief Implementation of the Worland based spherial laplacian projector
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_SPHLAPL_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_SPHLAPL_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Projector/Base/SphLapl.hpp"

#ifdef QUICC_USE_PFSOLVE
#include "QuICC/Transform/Wrappers/Worland/Projector/SphLaplviewGpuParallalt_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_SPHLAPL_HPP