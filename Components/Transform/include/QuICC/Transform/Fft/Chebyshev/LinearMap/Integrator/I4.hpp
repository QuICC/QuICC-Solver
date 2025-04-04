/**
 * @file I4.hpp
 * @brief Implementation of the Chebyshev based I4 integrator, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Base/I4.hpp"
#include "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Integrator/I4viewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
//#include "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Integrator/I4viewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I4_HPP
