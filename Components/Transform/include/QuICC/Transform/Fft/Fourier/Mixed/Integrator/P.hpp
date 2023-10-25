/**
 * @file P.hpp
 * @brief Implementation of the Fourier based mixed P integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_P_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_P_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/PBase.hpp"
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Integrator/PviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Integrator/PviewGpu_t.hpp.inc"
#endif
#ifdef QUICC_USE_VKFFT
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Integrator/PviewGpuVkFFT_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_P_HPP
