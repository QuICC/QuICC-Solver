/**
 * @file D1.hpp
 * @brief Implementation of the Fourier based D projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D1_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D1Base.hpp"
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Projector/D1viewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Projector/D1viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D1_HPP
