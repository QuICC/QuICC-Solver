/**
 * @file InvLapl2D.hpp
 * @brief Implementation of the Fourier based InvLapl2D Integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_INVLAPL2D_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_INVLAPL2D_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/InvLapl2DBase.hpp"
#include "QuICC/Transform/Wrappers/Fourier/Complex/Integrator/InvLapl2DviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Fourier/Complex/Integrator/InvLapl2DviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_INVLAPL2D_HPP
