/**
 * @file D1_Neg.hpp
 * @brief Implementation of the Fourier based mixed D1_Neg integrator
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_NEG_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_NEG_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D1_NegBase.hpp"
#include "QuICC/Transform/Wrappers/Mixed/Integrator/D1_NegviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Mixed/Integrator/D1_NegviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_INTEGRATOR_D1_NEG_HPP