/**
 * @file D2.hpp
 * @brief Implementation of the Fourier based mixed D2 integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D2_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D2_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D2Base.hpp"
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Integrator/D2viewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Fourier/Mixed/Integrator/D2viewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D2_HPP