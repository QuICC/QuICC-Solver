/**
 * @file P.hpp
 * @brief Implementation of the Fourier based mixed P projector (i.e. from phyiscal to modal space)
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_P_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_P_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/PBase.hpp"
#include "QuICC/Transform/Wrappers/Mixed/Projector/PviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Mixed/Projector/PviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_PROJECTOR_P_HPP
