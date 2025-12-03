/**
 * @file I2DivR1D1R1_Zero.hpp
 * @brief Implementation of the Worland based I2DivR1D1R1 integrator and zero l = 0 mode
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2DIVR1D1R1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2DIVR1D1R1_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Integrator/Base/I2DivR1D1R1_Zero.hpp"

#ifdef QUICC_USE_PFSOLVE
#include "QuICC/Transform/Wrappers/Worland/Integrator/I2DivR1D1R1_ZeroviewGpuParallalt_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2DIVR1D1R1_ZERO_HPP
