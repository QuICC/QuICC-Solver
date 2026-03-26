/**
 * @file DivY1D1.hpp
 * @brief Implementation of the Chebyshev based DivY1D1 projector, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1D1_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1D1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/Base/DivY1D1.hpp"
#include "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Projector/DivY1D1viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
//#include "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Projector/DivY1D1viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1D1_HPP
