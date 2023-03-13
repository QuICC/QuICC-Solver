/**
 * @file Lapl2D.hpp
 * @brief Implementation of the Fourier based 2D laplacian projector
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_LAPL2D_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_LAPL2D_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/Lapl2DBase.hpp"
#include "QuICC/Transform/Wrappers/Complex/Projector/Lapl2DviewCpu_t.hpp.inc"
#ifdef QUICC_USE_CUFFT
#include "QuICC/Transform/Wrappers/Complex/Projector/Lapl2DviewGpu_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_PROJECTOR_LAPL2D_HPP
