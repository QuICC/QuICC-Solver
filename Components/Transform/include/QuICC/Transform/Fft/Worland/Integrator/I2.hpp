/**
 * @file I2.hpp
 * @brief Implementation of the Worland based I2 integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Integrator/Base/I2.hpp"

#ifdef QUICC_USE_PFSOLVE
#include "QuICC/Transform/Wrappers/Worland/Integrator/I2viewGpuParallalt_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_I2_HPP
