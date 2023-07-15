/**
 * @file P.hpp
 * @brief Implementation of the Fourier based P projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_P_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_P_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/PBase.hpp"
#include "QuICC/Transform/Wrappers/Fourier/Complex/Projector/PviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Fourier/Complex/Projector/PviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_P_HPP
