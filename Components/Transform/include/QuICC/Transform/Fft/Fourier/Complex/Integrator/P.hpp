/**
 * @file P.hpp
 * @brief Implementation of the Fourier based P Integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_P_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_P_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/PBase.hpp"
#include "QuICC/Transform/Wrappers/Fourier/Complex/Integrator/PviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Fourier/Complex/Integrator/PviewGpu_t.hpp.inc"
#endif
#ifdef QUICC_USE_VKFFT
#include "QuICC/Transform/Wrappers/Fourier/Complex/Integrator/PviewGpuVkFFT_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_P_HPP
