/**
 * @file I2_Zero.hpp
 * @brief Implementation of the Worland based I2_Zero integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2_ZERO_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Integrator/Base/I2_Zero.hpp"

#ifdef QUICC_USE_PFSOLVE
#include "QuICC/Transform/Wrappers/Worland/Integrator/I2_ZeroviewGpuParallalt_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2_ZERO_HPP
