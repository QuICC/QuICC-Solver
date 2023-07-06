/**
 * @file Lapl2D.hpp
 * @brief Implementation of the Fourier based Lapl2D Integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_LAPL2D_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_LAPL2D_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/Lapl2DBase.hpp"
#include "QuICC/Transform/Wrappers/Fourier/Complex/Integrator/Lapl2DviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Fourier/Complex/Integrator/Lapl2DviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_INTEGRATOR_LAPL2D_HPP
