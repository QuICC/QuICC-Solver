/**
 * @file Mean.hpp
 * @brief Implementation of the Fourier based D projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_MEAN_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_MEAN_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/MeanBase.hpp"
#include "QuICC/Transform/Wrappers/Fourier/Complex/Projector/MeanviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Fourier/Complex/Projector/MeanviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_MEAN_HPP
