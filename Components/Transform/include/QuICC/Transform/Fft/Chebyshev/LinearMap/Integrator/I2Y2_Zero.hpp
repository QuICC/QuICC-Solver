/**
 * @file I2Y2_Zero.hpp
 * @brief Implementation of the Chebyshev based I2Y2_Zero integrator, with linear map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I2Y2_ZERO_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I2Y2_ZERO_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Base/I2Y2_Zero.hpp"
#include "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Integrator/I2Y2_ZeroviewCpu_t.hpp.inc"
#ifdef QUICC_HAS_CUDA_BACKEND
//#include "QuICC/Transform/Wrappers/Chebyshev/LinearMap/Integrator/I2Y2_ZeroviewGpu_t.hpp.inc"
#endif


#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_INTEGRATOR_I2Y2_ZERO_HPP
