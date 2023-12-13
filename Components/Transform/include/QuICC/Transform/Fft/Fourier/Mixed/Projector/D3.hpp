/**
 * @file D3.hpp
 * @brief Implementation of the Fourier based D projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D3_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D3_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D3Base.hpp"
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Projector/D3viewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Projector/D3viewGpu_t.hpp.inc"
#endif
#ifdef QUICC_USE_VKFFT
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Projector/D3viewGpuVkFFT_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_D3_HPP
