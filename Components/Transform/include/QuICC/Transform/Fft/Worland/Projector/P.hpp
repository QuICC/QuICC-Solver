/**
 * @file P.hpp
 * @brief Implementation of the Worland based P projector
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_P_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_P_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Projector/Base/P.hpp"

#ifdef QUICC_USE_PFSOLVE
#include "QuICC/Transform/Wrappers/Worland/Projector/PviewGpuParallalt_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_P_HPP
