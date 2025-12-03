/**
 * @file R1_Zero.hpp
 * @brief Implementation of the Worland based R1_Zero integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_R1_ZERO_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_R1_ZERO_HPP
// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Integrator/Base/R1_Zero.hpp"

#ifdef QUICC_USE_PFSOLVE
#include "QuICC/Transform/Wrappers/Worland/Integrator/R1_ZeroviewGpuParallalt_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_R1_ZERO_HPP
