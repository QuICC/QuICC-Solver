/**
 * @file P.hpp
 * @brief Implementation of the ALegendre based P projector
 */

#ifndef QUICC_TRANSFORM_FFT_ALEGENDRE_PROJECTOR_P_HPP
#define QUICC_TRANSFORM_FFT_ALEGENDRE_PROJECTOR_P_HPP

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/ALegendre/Projector/Base/P.hpp"

#ifdef QUICC_USE_PFSOLVE
#include "QuICC/Transform/Wrappers/ALegendre/Projector/PviewGpuParallalt_t.hpp.inc"
#endif

#endif // QUICC_TRANSFORM_FFT_ALegendre_PROJECTOR_P_HPP
