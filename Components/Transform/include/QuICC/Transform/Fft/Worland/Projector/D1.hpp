/**
 * @file D1.hpp
 * @brief Implementation of the Worland based D1 projector
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_D1_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_D1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Projector/Base/D1.hpp"

#ifdef QUICC_USE_PFSOLVE
#include "QuICC/Transform/Wrappers/Worland/Projector/D1viewGpuParallalt_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_D1_HPP