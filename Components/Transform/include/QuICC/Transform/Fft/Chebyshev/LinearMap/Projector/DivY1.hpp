/**
 * @file DivY1.hpp
 * @brief Implementation of the Chebyshev based DivY1 projector, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/Base/DivY1.hpp"
#include "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Projector/DivY1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
//#include "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Projector/DivY1viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1_HPP
