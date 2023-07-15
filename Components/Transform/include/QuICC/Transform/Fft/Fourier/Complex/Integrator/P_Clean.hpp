/**
 * @file P_Clean.hpp
 * @brief Implementation of the Fourier based P_Clean Integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_P_CLEAN_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_P_CLEAN_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/P_CleanBase.hpp"
#include "QuICC/Transform/Wrappers/Fourier/Complex/Integrator/P_CleanviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Fourier/Complex/Integrator/P_CleanviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_P_CLEAN_HPP
