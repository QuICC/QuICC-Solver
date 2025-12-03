/**
 * @file R1.hpp
 * @brief Implementation of the Worland based R1 integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_R1_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_R1_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Integrator/Base/R1.hpp"

#ifdef QUICC_USE_PFSOLVE
#include "QuICC/Transform/Wrappers/Worland/Integrator/R1viewGpuParallalt_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_R1_HPP
