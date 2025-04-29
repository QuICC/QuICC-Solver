/**
 * @file D1.hpp
 * @brief Implementation of the Fourier based mixed D1 integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D1Base.hpp"
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Integrator/D1viewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Integrator/D1viewGpu_t.hpp.inc"
#endif
#ifdef QUICC_USE_VKFFT
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Integrator/D1viewGpuVkFFT_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_HPP
