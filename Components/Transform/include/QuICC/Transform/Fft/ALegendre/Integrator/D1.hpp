/**
 * @file DivR1.hpp
 * @brief Implementation of the Worland based DivR1 integrator
 */

#ifndef QUICC_TRANSFORM_FFT_ALEGENDRE_INTEGRATOR_D1_HPP
#define QUICC_TRANSFORM_FFT_ALEGENDRE_INTEGRATOR_D1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/ALegendre/Integrator/Base/D1.hpp"

#ifdef QUICC_USE_PFSOLVE
#include "QuICC/Transform/Wrappers/ALegendre/Integrator/D1viewGpuParallalt_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_DIVR1_HPP
